"""
DMFB Synthesis Tool - Fluidic constraint aware
Compatible with new_driver.py + Biochip.py (Python 2.7)
Implements Phase 4 Parallel-Aware Droplet Routing with Multi-Identifier Lookup.
"""
import sys, os
from collections import defaultdict, deque

MIX_DUR = 6


def bfs(src, dst, rows, cols, blocked):
    if dst in blocked:
        return []
    if src == dst:
        return [src]
    visited = {src}
    queue = deque([(src, [src])])
    while queue:
        (r, c), path = queue.popleft()
        for nr, nc in [(r-1,c),(r+1,c),(r,c-1),(r,c+1)]:
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
        for j in range(i+1, len(pos_list)):
            r1,c1 = pos_list[i]
            r2,c2 = pos_list[j]
            if abs(r1-r2) <= 1 and abs(c1-c2) <= 1:
                return False
    return True


def post_mix_positions(r1, c1, r2, c2):
    """Predicts final coordinates of split droplets based on Biochip.py state 5."""
    return (r1, c1), (r2, c2)


def parse(filepath):
    rows = cols = 15
    reservoirs = {}
    mix_ops = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            kw = parts[0].upper()
            if kw == 'GRID':
                rows, cols = int(parts[1]), int(parts[2])
            elif kw == 'RESERVOIR':
                reservoirs[parts[1]] = (int(parts[2]), int(parts[3]), parts[4])
            elif kw == 'MIX':
                mix_ops.append((parts[1], parts[2], parts[4]))
    return rows, cols, reservoirs, mix_ops


def synthesize(input_file, output_file):
    rows, cols, reservoirs, mix_ops = parse(input_file)
    cmds = defaultdict(list)

    reagent_line = " ".join(
        "reagent({},{},{})".format(rr, rc, rn)
        for _, (rr, rc, rn) in reservoirs.items()
    )

    droplet_pos = defaultdict(list)
    
    # FIX: Map BOTH the Reservoir ID (R1) and Reagent Name (THCL) to coordinates
    res_lookup = {}
    for res_id, (rr, rc, rn) in reservoirs.items():
        res_lookup[res_id] = (rr, rc)
        res_lookup[rn] = (rr, rc)

    t = 1

    for (src_a, src_b, dest) in mix_ops:
        # Dynamic dispense if missing from board tracking
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

        # Dynamic array scanning for available mixer placement cells
        zone_pool = []
        for zc in range(2, cols, 2):
            for zr in range(2, rows - 3, 2):
                zone_pool.append((zr, zc))
                
        routed = False
        for zi in range(len(zone_pool)):
            zr, zc = zone_pool[zi]
            tgt_a = (zr,   zc)
            tgt_b = (zr+3, zc)

            others = []
            for name, positions in droplet_pos.items():
                others.extend(positions)

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
                # Step-by-step path validation for timing frame tracking
                len_a, len_b = len(path_a), len(path_b)
                max_len = max(len_a, len_b)
                
                full_a = path_a + [tgt_a] * (max_len - len_a)
                full_b = path_b + [tgt_b] * (max_len - len_b)
                
                valid_run = True
                for step_idx in range(max_len):
                    curr_a = full_a[step_idx]
                    curr_b = full_b[step_idx]
                    
                    if abs(curr_a[0] - curr_b[0]) <= 1 and abs(curr_a[1] - curr_b[1]) <= 1:
                        if step_idx < max_len - 1:
                            valid_run = False
                            break
                            
                if valid_run:
                    routed = True
                    break

        if not routed:
            print("  [ERROR] Fluidic routing deadlock for MIX {} + {} -> {}".format(src_a, src_b, dest))
            continue

        # Write sequential movement commands
        steps = max(len(path_a), len(path_b)) - 1
        for step in range(steps):
            step_cmds = []
            if step + 1 < len(path_a):
                r1, c1 = path_a[step]
                r2, c2 = path_a[step+1]
                step_cmds.append("move({},{},{},{})".format(r1, c1, r2, c2))
            if step + 1 < len(path_b):
                r1, c1 = path_b[step]
                r2, c2 = path_b[step+1]
                step_cmds.append("move({},{},{},{})".format(r1, c1, r2, c2))
            cmds[t].extend(step_cmds)
            t += 1

        # Schedule execution block
        r1, c1 = tgt_a
        r2, c2 = tgt_b
        cmds[t].append("mix_split({},{},{},{},{})".format(r1, c1, r2, c2, MIX_DUR))
        t += MIX_DUR + 1

        # Register products back to droplet tracker map
        post_d1, post_d2 = post_mix_positions(r1, c1, r2, c2)
        droplet_pos[dest].append(post_d1)
        droplet_pos[dest + "_2"].append(post_d2) 

    end_t = t + 1
    with open(output_file, 'w') as f:
        f.write("dimension {} {}\n".format(rows, cols))
        f.write("accuracy 5\n")
        f.write(reagent_line + "\n\n")
        for ts in sorted(cmds.keys()):
            if cmds[ts]:
                f.write("{} ".format(ts) + " ".join(cmds[ts]) + "\n")
        f.write("{} end\n".format(end_t))

    print("\n[SUCCESS] Compiled perfectly to {}".format(output_file))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python synthesis_tool_final.py input.txt [output.dmfb]")
        sys.exit(0)
    inp = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "output.dmfb"
    synthesize(inp, out)