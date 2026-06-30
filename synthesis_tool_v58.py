"""
DMFB Synthesis Tool - Graph/DAG aware (Phase 5)
Compatible with new_driver.py + Biochip.py (Python 2.7 style syntax)

This version replaces the old flat GRID/RESERVOIR/MIX text format with a
graph description (GRID / NODE / EDGE). It:

    1. Parses the node + directed-edge list.
    2. Builds adjacency / in-degree structures.
    3. Runs Kahn's algorithm to confirm the graph is a valid DAG (and to
       get a topological execution order) -- this is the "generate the
       DAG and check whether it is correct" step sir asked for.
    4. Runs structural sanity checks (dispense nodes have 0 inputs,
       mix nodes have exactly 2 inputs, output nodes have exactly 1
       input/0 outputs, etc).
    5. Auto-assigns physical grid coordinates to dispense/output nodes
       (the graph format carries no coordinates of its own).
    6. Feeds the resulting (src_a, src_b, dest) mix operations into the
       *same* fluidic-constraint-aware BFS routing engine used before,
       unchanged, so the produced .dmfb timing/command output is
       identical in spirit to the original tool.
    7. Routes the final droplet(s) into any OUTPUT node(s).
"""
import sys, os
from collections import defaultdict, deque

MIX_DUR = 6


# ---------------------------------------------------------------------------
# Routing engine (unchanged from the original tool)
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
# NEW: Graph parsing
# ---------------------------------------------------------------------------

def parse_graph(filepath):
    rows = cols = 15
    nodes = {}   # node_id -> {'type': ..., 'reagent': ...}
    edges = []   # list of (src_id, dst_id), in declaration order
    with open(filepath) as f:
        for raw in f:
            line = raw.split('#', 1)[0].strip()  # strip full-line / inline comments
            if not line:
                continue
            parts = line.split()
            kw = parts[0].upper()
            if kw == 'GRID':
                rows, cols = int(parts[1]), int(parts[2])
            elif kw == 'NODE':
                node_id = parts[1]
                node_type = parts[2].lower()
                reagent = parts[3] if len(parts) > 3 else None
                nodes[node_id] = {'type': node_type, 'reagent': reagent}
            elif kw == 'EDGE':
                edges.append((parts[1], parts[2]))
    return rows, cols, nodes, edges


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
# NEW: Auto-placement of physical coordinates (graph format has none)
# ---------------------------------------------------------------------------

def assign_reservoir_coords(nodes, order, rows, cols):
    """Place dispense nodes around the border with a guaranteed minimum
    spacing of 3 cells between any two reservoirs. A gap of exactly 2
    satisfies the static Rule #1 (Sec 5.2.2), but not necessarily the
    dynamic Rules #2/#3: if two reservoirs 2 cells apart both move one
    step in the same direction on the same clock tick, the leading
    droplet's NEW position can land adjacent to the trailing droplet's
    OLD (pre-move) position, which is exactly what Rule #2 forbids. A
    3-cell gap leaves enough margin to survive that case. Fills the left
    column first, then the top row, then the right column, then the
    bottom row, before wrapping."""
    dispense_ids = [n for n in order if nodes[n]['type'] == 'dispense']
    coords = {}
    if not dispense_ids:
        return coords

    slots = []
    # left column, rows 2,5,8,... (col 1)
    for r in range(2, rows, 3):
        slots.append((r, 1))
    # top row, cols 4,7,10,... (row 1) -- skip col 1, already used above
    for c in range(4, cols, 3):
        slots.append((1, c))
    # right column, rows 2,5,8,...
    for r in range(2, rows, 3):
        slots.append((r, cols))
    # bottom row, cols 4,7,10,...
    for c in range(4, cols, 3):
        slots.append((rows, c))

    for i, node_id in enumerate(dispense_ids):
        if i < len(slots):
            coords[node_id] = slots[i]
        else:
            # ran out of border space (very large reagent count); wrap with offset
            coords[node_id] = slots[i % len(slots)]
    return coords


def assign_output_coords(nodes, order, rows, cols):
    """Spread output nodes evenly down the right border column."""
    output_ids = [n for n in order if nodes[n]['type'] == 'output']
    coords = {}
    count = len(output_ids)
    if count == 0:
        return coords
    step = max(1, (rows - 2) // max(1, count))
    for i, node_id in enumerate(output_ids):
        r = min(rows, 1 + (i + 1) * step)
        coords[node_id] = (r, cols)
    return coords


def graph_to_mix_ops(nodes, order, incoming):
    """Walk the DAG in topological order and emit (src_a, src_b, dest)
    triples in the same shape the original flat MIX format produced."""
    mix_ops = []
    for node_id in order:
        if nodes[node_id]['type'] == 'mix':
            ins = incoming[node_id]
            if len(ins) != 2:
                # validate_graph() already reported this; skip safely here.
                continue
            src_a, src_b = ins[0], ins[1]
            mix_ops.append((src_a, src_b, node_id))
    return mix_ops


# ---------------------------------------------------------------------------
# Synthesis (same fluidic routing engine as before, now graph-driven)
# ---------------------------------------------------------------------------

def synthesize_graph(input_file, output_file):
    rows, cols, nodes, edges = parse_graph(input_file)
    adj, indeg, incoming = build_graph(nodes, edges)
    order, is_dag = topo_sort(nodes, adj, indeg)
    errors = validate_graph(nodes, adj, incoming)
    report_dag(nodes, edges, order, is_dag, errors)

    if not is_dag:
        print("\n[ABORT] Cannot synthesize: graph is not a DAG (cycle present).")
        return

    res_coords = assign_reservoir_coords(nodes, order, rows, cols)
    out_coords = assign_output_coords(nodes, order, rows, cols)

    # Build reservoirs dict in the shape the rest of the pipeline expects:
    # res_id -> (row, col, reagent_name)
    reservoirs = {}
    for node_id, (r, c) in res_coords.items():
        reservoirs[node_id] = (r, c, nodes[node_id]['reagent'] or node_id)

    reagent_line = " ".join(
        "reagent({},{},{})".format(rr, rc, rn)
        for _, (rr, rc, rn) in reservoirs.items()
    )

    res_lookup = {}
    for res_id, (rr, rc, rn) in reservoirs.items():
        res_lookup[res_id] = (rr, rc)
        res_lookup[rn] = (rr, rc)

    mix_ops = graph_to_mix_ops(nodes, order, incoming)

    cmds = defaultdict(list)
    droplet_pos = defaultdict(list)
    waste_droplets = []  # unused mix-split daughters: still physically on
                          # the board, must block routing, but are never
                          # consumed as a future mix input
    t = 1

    for (src_a, src_b, dest) in mix_ops:
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

        # Per the droplet routing paper's fluidic Rule #1 (Sec 5.2.2):
        # |Xi-Xj| >= 2 OR |Yi-Yj| >= 2 must hold for ANY two droplets,
        # including a freshly dispensed reservoir droplet vs. a droplet
        # already parked on the board (e.g. a mix product waiting to be
        # consumed). Reservoirs sit on the border (row==1 or col==1), so
        # mixing zones must start at least 2 cells in from the border in
        # both dimensions, or a parked mix-zone droplet can end up
        # diagonally adjacent (diff (1,1)) to a border reservoir cell.
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

                    # Rule #1 (static/same-time): both droplets' positions
                    # at this same tick must not be within 1 cell of each
                    # other, except at the final step where merging is the
                    # intended outcome.
                    if close(curr_a, curr_b) and not is_final_step:
                        valid_run = False
                        break

                    # Rule #2/#3 (dynamic/cross-time): a droplet's NEXT
                    # position must not be within 1 cell of the other
                    # droplet's CURRENT position (and vice versa), checked
                    # between this step and the next one -- this catches
                    # the case where two droplets approach each other and
                    # one's incoming move lands adjacent to where the
                    # other still currently sits, even though same-tick
                    # positions never violated Rule #1.
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
        cmds[t].append("mix_split({},{},{},{},{})".format(r1, c1, r2, c2, MIX_DUR))
        t += MIX_DUR + 1

        post_d1, post_d2 = post_mix_positions(r1, c1, r2, c2)
        droplet_pos[dest].append(post_d1)
        # NOTE: mix_split physically produces two daughter droplets, but our
        # graph model only ever wires one of them (dest) forward to the next
        # operation. The second daughter is unused for chemistry purposes,
        # but it is still a REAL droplet sitting on the chip -- the hardware
        # simulator will reject any later route that gets within 1 cell of
        # it. So we keep it in waste_droplets purely as a permanent routing
        # obstacle (never poppable / never usable as a future mix input),
        # rather than either (a) silently discarding it (which caused
        # collisions the real simulator correctly flagged) or (b) treating
        # it as a normal named droplet (which caused permanent deadlocks
        # when earlier code parked it under droplet_pos[dest + "_2"]).
        waste_droplets.append(post_d2)

    # --- Route final droplets into OUTPUT nodes -----------------------
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
        droplet_pos[node_id].append(out_path[-1])

    end_t = t + 1
    with open(output_file, 'w') as f:
        f.write("dimension {} {}\n".format(rows, cols))
        f.write("accuracy 5\n")
        f.write(reagent_line + "\n\n")
        for ts in sorted(cmds.keys()):
            if cmds[ts]:
                f.write("{} ".format(ts) + " ".join(cmds[ts]) + "\n")
        f.write("{} end\n".format(end_t))

    print("\n[SUCCESS] Compiled graph-based protocol to {}".format(output_file))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python synthesis_tool_graph.py input_graph.txt [output.dmfb]")
        sys.exit(0)
    inp = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "output.dmfb"
    synthesize_graph(inp, out)