#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DMFB synthesis tool aligned to SimBioSys Biochip_V.py.

Why this version is conservative:
- Biochip_V checks every move against the CURRENT chip state.
- A move is illegal if its destination is within Chebyshev distance <= 1 of
  any other droplet currently present on the chip.
- Active mixers remain represented as two droplets at their endpoints until the
  verifier expires them, so routing near an active mixer can trigger D-7 even
  if only one explicit move appears on that line.

Design choices in this tool:
1. Reagent dispense nodes are reusable/on-demand.
2. Only one droplet move is emitted per clock tick.
3. Only one mix / waste progression is handled at a time.
4. Future reservations use counted footprints so overlapping protections are not
   accidentally removed.
5. Waste routing is checked before a mixer placement is committed.

This is the correct / verifier-friendly version. It is intentionally not the
most cycle-optimized version.
"""
from __future__ import print_function

import sys
from collections import defaultdict, deque

CYCLE_SEC = 6
HORIZON = 1000
PLAN_HOLD_HORIZON = 250
MAX_PATH = 120


def ch_dist(p1, p2):
    return max(abs(p1[0] - p2[0]), abs(p1[1] - p2[1]))


def reserve_radius(table, ts, r, c, radius):
    for dr in range(-radius, radius + 1):
        for dc in range(-radius, radius + 1):
            key = (ts, r + dr, c + dc)
            table[key] += 1


def reserve_horizon(table, r, c, start_t, radius=2, horizon=HORIZON):
    for dt in range(start_t, start_t + horizon):
        reserve_radius(table, dt, r, c, radius)


def unreserve_horizon(table, r, c, start_t, radius=2, horizon=HORIZON):
    removed = defaultdict(int)
    for dt in range(start_t, start_t + horizon):
        for dr in range(-radius, radius + 1):
            for dc in range(-radius, radius + 1):
                key = (dt, r + dr, c + dc)
                if table.get(key, 0) > 0:
                    table[key] -= 1
                    removed[key] += 1
                    if table[key] <= 0:
                        del table[key]
    return removed


def restore_reservations(table, removed):
    for key, cnt in removed.items():
        table[key] += cnt


def bfs_space_time(src, dst, start_t, rows, cols, rtab, overlay=None, max_path=MAX_PATH):
    def blocked(t, r, c):
        key = (t, r, c)
        return rtab.get(key, 0) > 0 or (overlay is not None and overlay.get(key, 0) > 0)

    if src == dst:
        return [src]
    visited = set([(start_t, src[0], src[1])])
    q = deque([(src, start_t, [src])])
    while q:
        (r, c), ct, path = q.popleft()
        if len(path) > max_path:
            continue
        nt = ct + 1
        # Manhattan-only: DMFB droplets can only step to one of the 4
        # orthogonally-adjacent cells per tick, never diagonally.
        for nr, nc in (
            (r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)
        ):
            if not (1 <= nr <= rows and 1 <= nc <= cols):
                continue
            if blocked(nt, nr, nc):
                continue
            nxt = (nr, nc)
            if nxt == dst:
                return path + [nxt]
            st = (nt, nr, nc)
            if st not in visited:
                visited.add(st)
                q.append((nxt, nt, path + [nxt]))
    return []


def reserve_path_wide(table, path, start_t, rows, cols):
    for si in range(len(path) - 1):
        r2, c2 = path[si + 1]
        ct = start_t + si
        reserve_radius(table, ct, r2, c2, 2)
        reserve_radius(table, ct + 1, r2, c2, 2)


def emit_path(path, start_t, cmds, rtab):
    ct = start_t
    for si in range(len(path) - 1):
        r1, c1 = path[si]
        r2, c2 = path[si + 1]
        cmds[ct].append("move({},{},{},{})".format(r1, c1, r2, c2))
        reserve_radius(rtab, ct, r2, c2, 2)
        reserve_radius(rtab, ct + 1, r2, c2, 2)
        ct += 1
    return ct


def overlay_reserve_source_hold(overlay, pos, start_t, horizon=PLAN_HOLD_HORIZON):
    reserve_horizon(overlay, pos[0], pos[1], start_t, radius=2, horizon=horizon)


def line_for_reagents(rmeta):
    return " ".join(
        "reagent({},{},{})".format(rr, cc, rn)
        for _, (rr, cc, rn) in rmeta.items()
    )


def parse_architecture(fp):
    """Architecture files only declare the physical PRESENCE and LOCATION
    of reservoir/waste/output sites on the chip - never which reagent (if
    any) a reservoir will hold. That assignment is made per-run by the
    synthesis tool from the assay being compiled (see map_arch), the same
    way a microprocessor assigns registers per-program rather than the
    hardware spec hardcoding them. So each keyword here takes only a
    row/col pair, e.g.:

        GRID 14 14
        RESERVOIR 1 2
        RESERVOIR 1 6
        WASTE 14 2
        OUTPUT 8 14

    Order in the file is preserved (as plain lists) since that order is
    what map_arch uses to assign sites to this assay's dispense/output
    nodes.
    """
    rows = cols = 15
    reservoirs = []
    waste_ports = []
    output_ports = []
    with open(fp) as f:
        for raw in f:
            line = raw.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            kw = parts[0].upper()
            if kw == 'GRID':
                rows, cols = int(parts[1]), int(parts[2])
            elif kw == 'RESERVOIR':
                reservoirs.append((int(parts[1]), int(parts[2])))
            elif kw == 'WASTE':
                waste_ports.append((int(parts[1]), int(parts[2])))
            elif kw == 'OUTPUT':
                output_ports.append((int(parts[1]), int(parts[2])))
    return rows, cols, reservoirs, waste_ports, output_ports


def parse_assay(fp):
    nodes = {}
    edges = []
    with open(fp) as f:
        for raw in f:
            line = raw.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            kw = parts[0].upper()
            if kw == 'NODE':
                nid = parts[1]
                nt = parts[2].lower()
                rg = parts[3] if nt == 'dispense' and len(parts) > 3 else None
                sec = int(parts[3]) if nt == 'mix' and len(parts) > 3 else None
                nodes[nid] = {'type': nt, 'reagent': rg, 'seconds': sec}
            elif kw == 'EDGE':
                edges.append((parts[1], parts[2]))
    return nodes, edges


def build_graph(nodes, edges):
    adj = defaultdict(list)
    inc = defaultdict(list)
    ideg = dict((n, 0) for n in nodes)
    for a, b in edges:
        adj[a].append(b)
        inc[b].append(a)
        ideg[b] = ideg.get(b, 0) + 1
        ideg.setdefault(a, 0)
    return adj, ideg, inc


def topo_sort(nodes, adj, ideg):
    lid = dict(ideg)
    q = deque(sorted(n for n in nodes if lid.get(n, 0) == 0))
    order = []
    while q:
        n = q.popleft()
        order.append(n)
        for m in sorted(adj[n]):
            lid[m] -= 1
            if lid[m] == 0:
                q.append(m)
    return order, len(order) == len(nodes)


def map_arch(nodes, order, reservoirs, outputs):
    """Assign this assay's dispense/output nodes to physical sites.

    reservoirs/outputs are now plain location lists (no reagent names
    baked into the architecture file), so the assignment is made here,
    dynamically, per assay - round-robin over the declared sites in the
    order dispense/output nodes appear in topological order. This is the
    synthesis tool's job, not the architecture file's.
    """
    disp = [n for n in order if nodes[n]['type'] == 'dispense']
    out_ids = [n for n in order if nodes[n]['type'] == 'output']

    rc = {}
    if reservoirs:
        for i, n in enumerate(disp):
            rc[n] = reservoirs[i % len(reservoirs)]

    oc = {}
    if outputs:
        for i, n in enumerate(out_ids):
            oc[n] = outputs[i % len(outputs)]
    return rc, oc


def generate_mix_zones(rows, cols, no_park_zone):
    zones = []
    for c in range(2, cols + 1):
        for r in range(2, rows - 2):
            a = (r, c)
            b = (r + 3, c)
            if a in no_park_zone or b in no_park_zone:
                continue
            zones.append((a, b))
    for r in range(2, rows + 1):
        for c in range(2, cols - 2):
            a = (r, c)
            b = (r, c + 3)
            if a in no_park_zone or b in no_park_zone:
                continue
            zones.append((a, b))
    return zones


def zone_sort_key(zone, pa, pb):
    ta, tb = zone
    direct = abs(ta[0] - pa[0]) + abs(ta[1] - pa[1]) + abs(tb[0] - pb[0]) + abs(tb[1] - pb[1])
    swap = abs(tb[0] - pa[0]) + abs(tb[1] - pa[1]) + abs(ta[0] - pb[0]) + abs(ta[1] - pb[1])
    return min(direct, swap)


def zone_key_id(zone):
    # Canonical, orientation-independent identity for a physical mix
    # zone. Two mixes pinned to the exact same electrodes can never run
    # concurrently no matter how compaction reorders time, so we track
    # how often each site has been assigned (zone_use_count) and steer
    # later mixes to unused real estate instead of always taking the
    # nearest site.
    return frozenset(zone)


def choose_ready_mixes(pending_mixes, inc, nodes, parked_products, topo_index):
    ready = []
    for mid in pending_mixes:
        ok = True
        for src in inc[mid]:
            if nodes[src]['type'] == 'dispense':
                continue
            if src not in parked_products:
                ok = False
                break
        if ok:
            ready.append(mid)
    ready.sort(key=lambda x: topo_index[x])
    return ready


def source_descriptor(src, nodes, parked_products, reservoirs):
    if nodes[src]['type'] == 'dispense':
        return {'kind': 'dispense', 'name': src, 'pos': reservoirs[src]}
    return {'kind': 'parked', 'name': src, 'pos': parked_products[src]}


def reserve_zone_window(rtab, ta, tb, start_t, end_t, radius=2):
    for ts in range(start_t, end_t + 1):
        reserve_radius(rtab, ts, ta[0], ta[1], radius)
        reserve_radius(rtab, ts, tb[0], tb[1], radius)


def reserve_pos_window(table, pos, start_t, end_t, radius=2):
    for ts in range(start_t, end_t + 1):
        reserve_radius(table, ts, pos[0], pos[1], radius)


def try_plan_source(srcd, dst, ready_t, rows, cols, rtab, overlay=None):
    src = srcd['pos']
    if srcd['kind'] == 'dispense':
        dispense_t = ready_t
        start_move_t = ready_t + 1
        if rtab.get((dispense_t, src[0], src[1]), 0) > 0 or (
            overlay is not None and overlay.get((dispense_t, src[0], src[1]), 0) > 0
        ):
            return None
    else:
        dispense_t = None
        start_move_t = ready_t

    path = bfs_space_time(src, dst, start_move_t, rows, cols, rtab, overlay=overlay)
    if not path:
        return None

    arrival_t = start_move_t + len(path) - 1
    return {
        'srcd': srcd,
        'dispense_t': dispense_t,
        'start_t': start_move_t,
        'path': path,
        'arrival_t': arrival_t,
        'dst': dst,
    }


def commit_source_plan(plan, cmds, rtab, t_ref):
    srcd = plan['srcd']
    if srcd['kind'] == 'dispense':
        rr, cc = srcd['pos']
        cmds[plan['dispense_t']].append("dispense({},{})".format(rr, cc))
        reserve_horizon(rtab, rr, cc, plan['dispense_t'], radius=2, horizon=2)
    end_t = emit_path(plan['path'], plan['start_t'], cmds, rtab)
    return max(t_ref, end_t)


def simulate_candidate(mid, srcA, srcB, ta, tb, mix_seconds, t0,
                       rows, cols, rtab, waste_port):
    assignments = [((srcA, ta), (srcB, tb)), ((srcA, tb), (srcB, ta))]
    for (sa, da), (sb, db) in assignments:
        overlay = defaultdict(int)

        lifted_a = None
        if sa['kind'] == 'parked':
            lifted_a = unreserve_horizon(rtab, sa['pos'][0], sa['pos'][1], t0, radius=2)

        planA = try_plan_source(sa, da, t0, rows, cols, rtab, overlay=overlay)
        if not planA:
            if lifted_a is not None:
                restore_reservations(rtab, lifted_a)
            continue

        reserve_path_wide(overlay, planA['path'], planA['start_t'], rows, cols)
        overlay_reserve_source_hold(overlay, da, planA['arrival_t'])
        if planA['dispense_t'] is not None:
            reserve_radius(overlay, planA['dispense_t'], sa['pos'][0], sa['pos'][1], 2)

        lifted_b = None
        if sb['kind'] == 'parked':
            lifted_b = unreserve_horizon(rtab, sb['pos'][0], sb['pos'][1], planA['arrival_t'], radius=2)

        planB = try_plan_source(sb, db, planA['arrival_t'], rows, cols, rtab, overlay=overlay)
        if not planB:
            if lifted_a is not None:
                restore_reservations(rtab, lifted_a)
            if lifted_b is not None:
                restore_reservations(rtab, lifted_b)
            continue

        reserve_path_wide(overlay, planB['path'], planB['start_t'], rows, cols)
        overlay_reserve_source_hold(overlay, db, planB['arrival_t'])
        if planB['dispense_t'] is not None:
            reserve_radius(overlay, planB['dispense_t'], sb['pos'][0], sb['pos'][1], 2)

        st_mix = planB['arrival_t']
        fin_mix = st_mix + mix_seconds + 1

        blocked = False
        for ts in range(st_mix, fin_mix + 1):
            if rtab.get((ts, da[0], da[1]), 0) > 0 or rtab.get((ts, db[0], db[1]), 0) > 0:
                blocked = True
                break
        if blocked:
            if lifted_a is not None:
                restore_reservations(rtab, lifted_a)
            if lifted_b is not None:
                restore_reservations(rtab, lifted_b)
            continue

        waste_start = fin_mix
        waste_overlay = defaultdict(int)
        reserve_path_wide(waste_overlay, planA['path'], planA['start_t'], rows, cols)
        reserve_path_wide(waste_overlay, planB['path'], planB['start_t'], rows, cols)
        if planA['dispense_t'] is not None:
            reserve_radius(waste_overlay, planA['dispense_t'], sa['pos'][0], sa['pos'][1], 2)
        if planB['dispense_t'] is not None:
            reserve_radius(waste_overlay, planB['dispense_t'], sb['pos'][0], sb['pos'][1], 2)
        reserve_pos_window(waste_overlay, da, planA['arrival_t'], st_mix, radius=2)
        reserve_zone_window(waste_overlay, da, db, st_mix, fin_mix, radius=2)
        reserve_horizon(waste_overlay, da[0], da[1], fin_mix, radius=2)
        for ts in range(fin_mix, waste_start):
            reserve_radius(waste_overlay, ts, db[0], db[1], 2)

        waste_path = bfs_space_time(db, waste_port, waste_start, rows, cols, rtab, overlay=waste_overlay)
        if not waste_path:
            if lifted_a is not None:
                restore_reservations(rtab, lifted_a)
            if lifted_b is not None:
                restore_reservations(rtab, lifted_b)
            continue

        return {
            'mid': mid,
            'planA': planA,
            'planB': planB,
            'product_pos': da,
            'waste_pos': db,
            'st_mix': st_mix,
            'fin_mix': fin_mix,
            'waste_start': waste_start,
            'waste_path': waste_path,
        }
    return None


def commit_candidate(bundle, cmds, rtab, parked_products, nodes):
    planA = bundle['planA']
    planB = bundle['planB']
    da = bundle['product_pos']
    db = bundle['waste_pos']
    st_mix = bundle['st_mix']
    fin_mix = bundle['fin_mix']
    waste_start = bundle['waste_start']
    waste_path = bundle['waste_path']
    mid = bundle['mid']

    for plan in (planA, planB):
        srcd = plan['srcd']
        if srcd['kind'] == 'parked' and srcd['name'] in parked_products:
            del parked_products[srcd['name']]

    t = 1
    t = commit_source_plan(planA, cmds, rtab, t)
    reserve_pos_window(rtab, da, planA['arrival_t'], st_mix, radius=2)
    t = commit_source_plan(planB, cmds, rtab, t)

    mix_time = nodes[mid]['seconds'] or CYCLE_SEC
    cmds[st_mix].append(
        "mix_split({},{},{},{},{})".format(da[0], da[1], db[0], db[1], mix_time)
    )
    reserve_zone_window(rtab, da, db, st_mix, fin_mix, radius=2)

    reserve_horizon(rtab, da[0], da[1], fin_mix, radius=2)
    reserve_horizon(rtab, db[0], db[1], fin_mix, radius=2)

    unreserve_horizon(rtab, db[0], db[1], waste_start, radius=2)
    end_waste = emit_path(waste_path, waste_start, cmds, rtab)
    cmds[end_waste].append("waste({},{})".format(waste_path[-1][0], waste_path[-1][1]))

    parked_products[mid] = da
    return end_waste + 1


def synthesize_graph(arch_file, assay_file, output_file):
    rows, cols, reservoirs, waste_ports, output_ports = parse_architecture(arch_file)
    nodes, edges = parse_assay(assay_file)
    adj, ideg, inc = build_graph(nodes, edges)
    order, is_dag = topo_sort(nodes, adj, ideg)
    if not is_dag:
        print('[ERROR] Graph is not a DAG - cycles detected')
        return False

    rc, oc = map_arch(nodes, order, reservoirs, output_ports)
    rmeta = dict((n, (r, c, nodes[n].get('reagent') or n)) for n, (r, c) in rc.items())
    rline = line_for_reagents(rmeta)
    wh = ' '.join('waste_reservoir({},{})'.format(r, c) for (r, c) in waste_ports)
    oh = ' '.join('output_reservoir({},{})'.format(r, c) for (r, c) in output_ports)
    header_extra = ' '.join(x for x in [wh, oh] if x).strip()

    if not waste_ports:
        print('[ERROR] No waste port defined')
        return False
    waste_port = sorted(waste_ports)[0]

    no_park_zone = set()
    for pr, pc in list(reservoirs) + list(waste_ports) + list(output_ports):
        for dr in range(-2, 3):
            for dc in range(-2, 3):
                no_park_zone.add((pr + dr, pc + dc))

    all_zones = generate_mix_zones(rows, cols, no_park_zone)
    topo_index = dict((n, i) for i, n in enumerate(nodes.keys()))

    cmds = defaultdict(list)
    rtab = defaultdict(int)
    parked_products = {}
    pending_mixes = set([n for n in nodes if nodes[n]['type'] == 'mix'])
    zone_use_count = defaultdict(int)

    t = 1
    guard = 0
    while pending_mixes:
        guard += 1
        if guard > 10000:
            print('[ERROR] Scheduling guard triggered; unresolved mixes:', len(pending_mixes))
            return False

        ready = choose_ready_mixes(pending_mixes, inc, nodes, parked_products, topo_index)
        if not ready:
            print('[ERROR] No ready mix could be found, but mixes remain unresolved:', sorted(pending_mixes))
            return False

        scheduled = False
        for mid in ready:
            srcs = inc[mid]
            if len(srcs) != 2:
                print('[ERROR] Mix node {} does not have exactly 2 inputs'.format(mid))
                return False
            sa = source_descriptor(srcs[0], nodes, parked_products, rc)
            sb = source_descriptor(srcs[1], nodes, parked_products, rc)
            pa, pb = sa['pos'], sb['pos']
            # Reused physical zones are the reason compact_from_safe.py
            # could never parallelize this schedule: two mixes pinned to
            # the identical (ta,tb) rectangle are the same electrodes and
            # can never run concurrently no matter how time is reordered.
            # Bias strongly toward zones this run hasn't used yet
            # (falling back to reuse only once the grid genuinely runs low
            # on unused real estate), so independent branches of the assay
            # land on distinct physical sites and can actually overlap in
            # time once compacted.
            zones = sorted(
                all_zones,
                key=lambda z: (zone_use_count[zone_key_id(z)], zone_sort_key(z, pa, pb))
            )[:300]

            bundle = None
            for ta, tb in zones:
                bundle = simulate_candidate(
                    mid, sa, sb, ta, tb, nodes[mid]['seconds'] or CYCLE_SEC,
                    t, rows, cols, rtab, waste_port
                )
                if bundle is not None:
                    break

            if bundle is None:
                continue

            t = commit_candidate(bundle, cmds, rtab, parked_products, nodes)
            zone_use_count[zone_key_id((bundle['product_pos'], bundle['waste_pos']))] += 1
            pending_mixes.remove(mid)
            scheduled = True
            break

        if not scheduled:
            print('[ERROR] No ready mix could be placed/routed safely at time {}.'.format(t))
            print('        Remaining ready mixes:', ready)
            return False

    for nid in order:
        if nodes[nid]['type'] != 'output':
            continue
        preds = inc[nid]
        if len(preds) != 1:
            print('[ERROR] Output node {} must have exactly one input'.format(nid))
            return False
        src = preds[0]
        if src not in parked_products:
            print('[ERROR] Output source {} for {} is not available'.format(src, nid))
            return False
        if nid not in oc:
            print('[ERROR] No output port mapped for {}'.format(nid))
            return False
        sp = parked_products[src]
        op = oc[nid]
        removed = unreserve_horizon(rtab, sp[0], sp[1], t, radius=2)
        opath = bfs_space_time(sp, op, t, rows, cols, rtab)
        if not opath:
            restore_reservations(rtab, removed)
            print('[ERROR] Failed to route output {} from {} to {}'.format(nid, sp, op))
            return False
        t = emit_path(opath, t, cmds, rtab)
        cmds[t].append('output({},{})'.format(opath[-1][0], opath[-1][1]))
        del parked_products[src]
        t += 1

    with open(output_file, 'w') as f:
        f.write('dimension {} {}\n'.format(rows, cols))
        f.write('accuracy 5\n')
        header = rline
        if header_extra:
            header = (header + ' ' + header_extra).strip()
        f.write(header + '\n\n')
        for ts in sorted(cmds.keys()):
            if cmds[ts]:
                f.write('{} {}\n'.format(ts, ' '.join(cmds[ts])))
        f.write('{} end\n'.format(t))

    print('[SUCCESS] Compiled cleanly to: {}'.format(output_file))
    return True


if __name__ == '__main__':
    arch = sys.argv[1] if len(sys.argv) > 1 else 'arch6.txt'
    assay = sys.argv[2] if len(sys.argv) > 2 else 'graph6.txt'
    out = sys.argv[3] if len(sys.argv) > 3 else 'graph6_output.dmfb'
    ok = synthesize_graph(arch, assay, out)
    sys.exit(0 if ok else 1)