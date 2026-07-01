"""
DMFB Synthesis Tool - Graph/DAG aware (Phase 6)
Compatible with new_driver.py + Biochip.py (Python 2.7 style syntax)

Changes from v58:
    1. Two separate input files instead of one:
         - architecture file: GRID / RESERVOIR / WASTE / OUTPUT
           (fixed physical port locations - hardware description)
         - assay file: NODE / EDGE (dispense/mix/output operations,
           mix nodes now carry a duration in seconds)
       assign_reservoir_coords() / assign_output_coords() (auto
       placement) are removed - physical coordinates now come straight
       from the architecture file, and a mapping step binds the
       assay's logical dispense/output node IDs onto those physical
       ports (so a physical port can be designated waste/output while
       mapping the assay, per the architecture file's declared role).
    2. Per-mix-operation timing: seconds parsed from the assay file are
       converted to cycles (1 cycle = 6 sec) and rounded UP to the
       nearest multiple of 6, replacing the old flat MIX_DUR constant.
    3. Waste and output are now real hardware ports:
         - header line gets waste_reservoir(r,c) / output_reservoir(r,c)
           entries (matches new_driver.py's header() parser)
         - a droplet reaching the output port emits an output(r,c)
           command (matches new_driver.py's gen(), calls delete_droplet)
         - the unused mix_split daughter droplet is BFS-routed to the
           nearest waste port and emits a waste(r,c) command, instead
           of being left as a permanent static routing obstacle
    Rest operations are intentionally NOT implemented in this version.
"""
import sys, os
from collections import defaultdict, deque

CYCLE_SEC = 6  # 1 cycle = 6 seconds


def seconds_to_cycles(seconds):
    """Convert a mix duration in seconds to cycles, rounded UP to the
    nearest multiple of CYCLE_SEC (6). E.g. 12 sec -> 2 cycles -> 6;
    30 sec -> 5 cycles -> 6; 42 sec -> 7 cycles -> 12."""
    raw_cycles = (seconds + CYCLE_SEC - 1) // CYCLE_SEC  # ceil division
    if raw_cycles <= 0:
        raw_cycles = 1
    return ((raw_cycles + 5) // 6) * 6


# ---------------------------------------------------------------------------
# Routing engine (unchanged from v58)
# ---------------------------------------------------------------------------

def bfs(src, dst, rows, cols, blocked):
    if dst in blocked:
        return []
    if src == dst:
        return [src]
    visited = {src}
    queue = deque([(src, [src])])
    while queue:
        (r, c), path = queue.popleft()
        for nr, nc in [(r-1, c), (r+1, c), (r, c-1), (r, c+1)]:
            nxt = (nr, nc)
            if nxt == dst:
                return path + [nxt]
            if (1 <= nr <= rows and 1 <= nc <= cols
                    and nxt not in visited and nxt not in blocked):
                visited.add(nxt)
                queue.append((nxt, path + [nxt]))
    return []


def fluidic_ok(all_positions):
    """Check no two droplets are within 1 cell of each other (Rules 1, 2, and 3)."""
    pos_list = list(all_positions)
    for i in range(len(pos_list)):
        for j in range(i + 1, len(pos_list)):
            r1, c1 = pos_list[i]
            r2, c2 = pos_list[j]
            if abs(r1 - r2) <= 1 and abs(c1 - c2) <= 1:
                return False
    return True


def post_mix_positions(r1, c1, r2, c2):
    """Predicts final coordinates of split droplets based on Biochip.py state 5."""
    return (r1, c1), (r2, c2)


# ---------------------------------------------------------------------------
# NEW: Architecture file parsing (hardware description)
# ---------------------------------------------------------------------------

def parse_architecture(filepath):
    """Parses GRID / RESERVOIR / WASTE / OUTPUT lines.
    Returns (rows, cols, reservoirs, waste_ports, output_ports) where:
        reservoirs   : id -> (row, col)
        waste_ports  : id -> (row, col)
        output_ports : id -> (row, col)
    """
    rows = cols = 15
    reservoirs = {}
    waste_ports = {}
    output_ports = {}
    with open(filepath) as f:
        for raw in f:
            line = raw.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            kw = parts[0].upper()
            if kw == 'GRID':
                rows, cols = int(parts[1]), int(parts[2])
            elif kw == 'RESERVOIR':
                port_id, r, c = parts[1], int(parts[2]), int(parts[3])
                reservoirs[port_id] = (r, c)
            elif kw == 'WASTE':
                port_id, r, c = parts[1], int(parts[2]), int(parts[3])
                waste_ports[port_id] = (r, c)
            elif kw == 'OUTPUT':
                port_id, r, c = parts[1], int(parts[2]), int(parts[3])
                output_ports[port_id] = (r, c)
    return rows, cols, reservoirs, waste_ports, output_ports


# ---------------------------------------------------------------------------
# NEW: Assay file parsing (logical operation graph)
# ---------------------------------------------------------------------------

def parse_assay(filepath):
    """Parses NODE / EDGE lines. mix nodes carry a trailing seconds value.
    Returns (nodes, edges) where nodes[id] = {'type', 'reagent', 'seconds'}."""
    nodes = {}
    edges = []
    with open(filepath) as f:
        for raw in f:
            line = raw.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            kw = parts[0].upper()
            if kw == 'NODE':
                node_id = parts[1]
                node_type = parts[2].lower()
                reagent = None
                seconds = None
                if node_type == 'dispense':
                    reagent = parts[3] if len(parts) > 3 else None
                elif node_type == 'mix':
                    seconds = int(parts[3]) if len(parts) > 3 else None
                nodes[node_id] = {'type': node_type, 'reagent': reagent, 'seconds': seconds}
            elif kw == 'EDGE':
                edges.append((parts[1], parts[2]))
    return nodes, edges


def build_graph(nodes, edges):
    adj = defaultdict(list)
    incoming = defaultdict(list)
    indeg = dict((n, 0) for n in nodes)
    for (a, b) in edges:
        adj[a].append(b)
        incoming[b].append(a)
        indeg[b] = indeg.get(b, 0) + 1
        indeg.setdefault(a, 0)
    return adj, indeg, incoming


def topo_sort(nodes, adj, indeg):
    """Kahn's algorithm. Returns (order, is_dag)."""
    local_indeg = dict(indeg)
    queue = deque(sorted(n for n in nodes if local_indeg.get(n, 0) == 0))
    order = []
    while queue:
        n = queue.popleft()
        order.append(n)
        for m in sorted(adj[n]):
            local_indeg[m] -= 1
            if local_indeg[m] == 0:
                queue.append(m)
    is_dag = (len(order) == len(nodes))
    return order, is_dag


def validate_graph(nodes, adj, incoming):
    """Structural sanity checks, independent of cycle detection."""
    errors = []
    for node_id, meta in nodes.items():
        ntype = meta['type']
        n_in = len(incoming[node_id])
        n_out = len(adj[node_id])
        if ntype == 'dispense':
            if meta['reagent'] is None:
                errors.append("Node {}: type 'dispense' has no reagent name".format(node_id))
            if n_in != 0:
                errors.append("Node {}: dispense node should have 0 incoming edges, has {}".format(node_id, n_in))
            if n_out == 0:
                errors.append("Node {}: dispensed droplet is never consumed".format(node_id))
        elif ntype == 'mix':
            if n_in != 2:
                errors.append("Node {}: mix node must have exactly 2 incoming edges, has {}".format(node_id, n_in))
            if n_out == 0:
                errors.append("Node {}: mix result is never consumed".format(node_id))
            if meta['seconds'] is None:
                errors.append("Node {}: mix node has no duration (seconds) specified".format(node_id))
        elif ntype == 'output':
            if n_in != 1:
                errors.append("Node {}: output node must have exactly 1 incoming edge, has {}".format(node_id, n_in))
            if n_out != 0:
                errors.append("Node {}: output node should be a sink, has {} outgoing edges".format(node_id, n_out))
        else:
            errors.append("Node {}: unknown node type '{}'".format(node_id, ntype))
    return errors


def report_dag(nodes, edges, order, is_dag, errors):
    print("=" * 60)
    print("DAG VALIDATION REPORT")
    print("=" * 60)
    print("Nodes: {}   Edges: {}".format(len(nodes), len(edges)))
    if is_dag:
        print("[OK] Graph is a valid DAG (no cycles detected).")
        print("Topological execution order: {}".format(" -> ".join(order)))
    else:
        unresolved = [n for n in nodes if n not in order]
        print("[FAIL] Graph contains a CYCLE. Could not fully order: {}".format(unresolved))
    if errors:
        print("Structural warnings/errors:")
        for e in errors:
            print("  - {}".format(e))
    else:
        print("[OK] No structural errors (dispense/mix/output edge counts all correct).")
    print("=" * 60)


# ---------------------------------------------------------------------------
# NEW: Mapping step - bind assay logical nodes onto architecture ports
# ---------------------------------------------------------------------------

def map_assay_to_architecture(nodes, order, reservoirs, output_ports):
    """Binds each assay 'dispense' node onto a physical RESERVOIR port
    (matched by node id) and each assay 'output' node onto a physical
    OUTPUT port (matched by node id, falling back to the next unused
    OUTPUT port if ids don't line up 1:1). This is the "designate
    reservoirs as waste/output while mapping the assay" step - the
    physical port identity/coords are fixed in the architecture file,
    the assay only decides which logical operation uses which port.
    Returns (res_coords, out_coords)."""
    dispense_ids = [n for n in order if nodes[n]['type'] == 'dispense']
    output_ids = [n for n in order if nodes[n]['type'] == 'output']

    res_coords = {}
    for node_id in dispense_ids:
        if node_id in reservoirs:
            res_coords[node_id] = reservoirs[node_id]
        else:
            print("  [ERROR] Assay dispense node '{}' has no matching RESERVOIR "
                  "in architecture file".format(node_id))

    out_coords = {}
    unused_ports = [p for p in output_ports if p not in output_ids]
    for node_id in output_ids:
        if node_id in output_ports:
            out_coords[node_id] = output_ports[node_id]
        elif unused_ports:
            port_id = unused_ports.pop(0)
            out_coords[node_id] = output_ports[port_id]
        else:
            print("  [ERROR] Assay output node '{}' has no matching OUTPUT "
                  "port in architecture file".format(node_id))

    return res_coords, out_coords


def graph_to_mix_ops(nodes, order, incoming):
    """Walk the DAG in topological order and emit (src_a, src_b, dest,
    cycles) tuples - cycles derived from the assay's per-node seconds."""
    mix_ops = []
    for node_id in order:
        if nodes[node_id]['type'] == 'mix':
            ins = incoming[node_id]
            if len(ins) != 2:
                continue
            src_a, src_b = ins[0], ins[1]
            seconds = nodes[node_id]['seconds'] or CYCLE_SEC
            cycles = seconds_to_cycles(seconds)
            mix_ops.append((src_a, src_b, node_id, cycles))
    return mix_ops


# ---------------------------------------------------------------------------
# Synthesis
# ---------------------------------------------------------------------------

def synthesize_graph(arch_file, assay_file, output_file):
    rows, cols, reservoirs, waste_ports, output_ports = parse_architecture(arch_file)
    nodes, edges = parse_assay(assay_file)
    adj, indeg, incoming = build_graph(nodes, edges)
    order, is_dag = topo_sort(nodes, adj, indeg)
    errors = validate_graph(nodes, adj, incoming)
    report_dag(nodes, edges, order, is_dag, errors)

    if not is_dag:
        print("\n[ABORT] Cannot synthesize: graph is not a DAG (cycle present).")
        return

    res_coords, out_coords = map_assay_to_architecture(nodes, order, reservoirs, output_ports)

    # Build reservoirs dict in the shape the rest of the pipeline expects:
    # res_id -> (row, col, reagent_name)
    reservoir_meta = {}
    for node_id, (r, c) in res_coords.items():
        reservoir_meta[node_id] = (r, c, nodes[node_id]['reagent'] or node_id)

    reagent_line = " ".join(
        "reagent({},{},{})".format(rr, rc, rn)
        for _, (rr, rc, rn) in reservoir_meta.items()
    )

    # Header also carries waste/output PORT declarations (paints the
    # cell + registers its type in Biochip.py, matches new_driver.py's
    # header() opcodes 'waste_reservoir' / 'output_reservoir')
    waste_header = " ".join(
        "waste_reservoir({},{})".format(r, c) for (r, c) in waste_ports.values()
    )
    output_header = " ".join(
        "output_reservoir({},{})".format(r, c) for (r, c) in output_ports.values()
    )
    header_extra = " ".join(x for x in [waste_header, output_header] if x)

    res_lookup = {}
    for res_id, (rr, rc, rn) in reservoir_meta.items():
        res_lookup[res_id] = (rr, rc)
        res_lookup[rn] = (rr, rc)

    mix_ops = graph_to_mix_ops(nodes, order, incoming)

    cmds = defaultdict(list)
    droplet_pos = defaultdict(list)
    waste_droplets = []  # unused mix-split daughters awaiting routing to a waste port
    t = 1

    for (src_a, src_b, dest, mix_cycles) in mix_ops:
        for src in (src_a, src_b):
            if not droplet_pos[src] and src in res_lookup:
                rrow, rcol = res_lookup[src]
                cmds[t].append("dispense({},{})".format(rrow, rcol))
                droplet_pos[src].append((rrow, rcol))

        if cmds[t]:
            t += 1

        if not droplet_pos[src_a] or not droplet_pos[src_b]:
            print("  [ERROR] Dependency missing for MIX {} + {} -> {}".format(src_a, src_b, dest))
            continue

        pos_a = droplet_pos[src_a].pop(0)
        pos_b = droplet_pos[src_b].pop(0)

        zone_pool = []
        for zc in range(3, cols, 2):
            for zr in range(3, rows - 3, 2):
                zone_pool.append((zr, zc))

        routed = False
        path_a, path_b = [], []
        for zi in range(len(zone_pool)):
            zr, zc = zone_pool[zi]
            tgt_a = (zr, zc)
            tgt_b = (zr + 3, zc)

            others = []
            for name, positions in droplet_pos.items():
                others.extend(positions)
            others.extend(waste_droplets)

            static_blocks = set()
            for (or_, oc) in others:
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        static_blocks.add((or_ + dr, oc + dc))

            path_a = bfs(pos_a, tgt_a, rows, cols, static_blocks | {pos_b})
            path_b = bfs(pos_b, tgt_b, rows, cols, static_blocks | set(path_a))

            if not path_a or not path_b:
                path_a = bfs(pos_a, tgt_b, rows, cols, static_blocks | {pos_b})
                path_b = bfs(pos_b, tgt_a, rows, cols, static_blocks | set(path_a))
                if path_a and path_b:
                    tgt_a, tgt_b = tgt_b, tgt_a

            if path_a and path_b:
                len_a, len_b = len(path_a), len(path_b)
                max_len = max(len_a, len_b)

                full_a = path_a + [tgt_a] * (max_len - len_a)
                full_b = path_b + [tgt_b] * (max_len - len_b)

                def close(p, q):
                    return abs(p[0] - q[0]) <= 1 and abs(p[1] - q[1]) <= 1

                valid_run = True
                for step_idx in range(max_len):
                    curr_a = full_a[step_idx]
                    curr_b = full_b[step_idx]
                    is_final_step = (step_idx == max_len - 1)

                    if close(curr_a, curr_b) and not is_final_step:
                        valid_run = False
                        break

                    if step_idx + 1 < max_len:
                        next_a = full_a[step_idx + 1]
                        next_b = full_b[step_idx + 1]
                        next_is_final = (step_idx + 1 == max_len - 1)
                        if not next_is_final:
                            if close(next_a, curr_b) or close(curr_a, next_b):
                                valid_run = False
                                break

                if valid_run:
                    routed = True
                    break

        if not routed:
            print("  [ERROR] Fluidic routing deadlock for MIX {} + {} -> {}".format(src_a, src_b, dest))
            continue

        steps = max(len(path_a), len(path_b)) - 1
        for step in range(steps):
            step_cmds = []
            if step + 1 < len(path_a):
                r1, c1 = path_a[step]
                r2, c2 = path_a[step + 1]
                step_cmds.append("move({},{},{},{})".format(r1, c1, r2, c2))
            if step + 1 < len(path_b):
                r1, c1 = path_b[step]
                r2, c2 = path_b[step + 1]
                step_cmds.append("move({},{},{},{})".format(r1, c1, r2, c2))
            cmds[t].extend(step_cmds)
            t += 1

        r1, c1 = path_a[-1] if path_a else pos_a
        r2, c2 = path_b[-1] if path_b else pos_b
        cmds[t].append("mix_split({},{},{},{},{})".format(r1, c1, r2, c2, mix_cycles))
        t += mix_cycles + 1

        post_d1, post_d2 = post_mix_positions(r1, c1, r2, c2)
        droplet_pos[dest].append(post_d1)
        # Second mix_split daughter is real but unused by the assay graph.
        # Route it to the nearest waste port instead of parking it as a
        # permanent obstacle.
        waste_droplets.append(post_d2)

    # --- Route leftover waste daughters into WASTE ports ----------------
    waste_port_list = list(waste_ports.values())
    if waste_port_list:
        still_pending = []
        for wd_pos in waste_droplets:
            others = []
            for name, positions in droplet_pos.items():
                others.extend(positions)
            others.extend([w for w in waste_droplets if w != wd_pos])
            static_blocks = set()
            for (or_, oc) in others:
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        static_blocks.add((or_ + dr, oc + dc))

            best_path = []
            for wport in waste_port_list:
                candidate = bfs(wd_pos, wport, rows, cols, static_blocks)
                if candidate and (not best_path or len(candidate) < len(best_path)):
                    best_path = candidate

            if not best_path:
                print("  [ERROR] Could not route waste droplet at {} to any waste port".format(wd_pos))
                still_pending.append(wd_pos)
                continue

            for step in range(len(best_path) - 1):
                r1, c1 = best_path[step]
                r2, c2 = best_path[step + 1]
                cmds[t].append("move({},{},{},{})".format(r1, c1, r2, c2))
                t += 1
            wr, wc = best_path[-1]
            cmds[t].append("waste({},{})".format(wr, wc))
            t += 1
        waste_droplets = still_pending
    elif waste_droplets:
        print("  [WARNING] {} waste droplet(s) produced but no WASTE port declared "
              "in architecture file - left on board as static obstacles".format(len(waste_droplets)))

    # --- Route final droplets into OUTPUT ports --------------------------
    for node_id in order:
        if nodes[node_id]['type'] != 'output':
            continue
        src_ids = incoming[node_id]
        if not src_ids:
            continue
        src_id = src_ids[0]
        if not droplet_pos[src_id]:
            print("  [ERROR] No droplet available to route into output {}".format(node_id))
            continue
        if node_id not in out_coords:
            continue
        start_pos = droplet_pos[src_id].pop(0)
        out_pos = out_coords[node_id]

        others = []
        for name, positions in droplet_pos.items():
            others.extend(positions)
        others.extend(waste_droplets)
        static_blocks = set()
        for (or_, oc) in others:
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    static_blocks.add((or_ + dr, oc + dc))

        out_path = bfs(start_pos, out_pos, rows, cols, static_blocks)
        if not out_path:
            print("  [ERROR] Could not route droplet from {} to output {}".format(src_id, node_id))
            continue
        for step in range(len(out_path) - 1):
            r1, c1 = out_path[step]
            r2, c2 = out_path[step + 1]
            cmds[t].append("move({},{},{},{})".format(r1, c1, r2, c2))
            t += 1
        # Droplet has arrived at the output port - emit output(r,c) so
        # new_driver.py's gen() removes it from the board (delete_droplet)
        or_, oc = out_path[-1]
        cmds[t].append("output({},{})".format(or_, oc))
        t += 1
        droplet_pos[node_id].append(out_path[-1])

    end_t = t + 1
    with open(output_file, 'w') as f:
        f.write("dimension {} {}\n".format(rows, cols))
        f.write("accuracy 5\n")
        line = reagent_line
        if header_extra:
            line = (line + " " + header_extra).strip()
        f.write(line + "\n\n")
        for ts in sorted(cmds.keys()):
            if cmds[ts]:
                f.write("{} ".format(ts) + " ".join(cmds[ts]) + "\n")
        f.write("{} end\n".format(end_t))

    print("\n[SUCCESS] Compiled graph-based protocol to {}".format(output_file))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python synthesis_tool.py architecture.txt assay.txt [output.dmfb]")
        sys.exit(0)
    arch = sys.argv[1]
    assay = sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else "output.dmfb"
    synthesize_graph(arch, assay, out)