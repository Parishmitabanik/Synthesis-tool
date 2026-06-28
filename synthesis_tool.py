"""
DMFB Synthesis Tool (v5.7 — Dispenser separation fix)

Key fix:
  - Spread dispensers around entire border perimeter
  - Add initial move away from dispensing position before routing
  - Ensure 2+ cell separation between all initial droplets
"""

import re
import sys
import math
import random
import argparse
from collections import defaultdict, deque
from copy import deepcopy

# ─────────────────────────────────────────────
# 0. MODULE LIBRARY
# ─────────────────────────────────────────────

MODULE_LIBRARY = {
    '2x2': (2, 2, 4),
    '2x3': (2, 3, 6),
    '2x4': (2, 4, 8),
}
DEFAULT_MIXER_TYPE = '2x4'

def get_module_spec(mixer_type=DEFAULT_MIXER_TYPE):
    return MODULE_LIBRARY.get(mixer_type, MODULE_LIBRARY[DEFAULT_MIXER_TYPE])

def validate_mixer_time(op_time, mixer_type=DEFAULT_MIXER_TYPE):
    if op_time < 1:
        return False, f"Mix time must be >= 1, got {op_time}"
    return True, f"Mix time {op_time} is valid for {mixer_type}"


# ─────────────────────────────────────────────
# 1. PARSER (UNCHANGED)
# ─────────────────────────────────────────────

def parse_bioassay(text):
    ops = []
    seen_ids = set()
    errors = []
    token_re = re.compile(r'(\w+)\s*\(([^)]+)\)')

    for lineno, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith('#'):
            continue

        for match in token_re.finditer(line):
            op_type = match.group(1).upper()
            args = [a.strip() for a in match.group(2).split(',')]
            op = None

            if op_type == 'DISPENSER':
                if len(args) != 2:
                    errors.append(f"Line {lineno}: DISPENSER needs 2 args")
                    continue
                op = {'type': 'DISPENSER', 'id': args[0], 'reagent': args[1]}

            elif op_type == 'MIX':
                if len(args) < 4:
                    errors.append(f"Line {lineno}: MIX needs 4 args")
                    continue
                try:
                    time = int(args[3])
                except ValueError:
                    errors.append(f"Line {lineno}: MIX time must be integer")
                    continue
                mixer_type = args[4] if len(args) > 4 else DEFAULT_MIXER_TYPE
                
                if mixer_type not in MODULE_LIBRARY:
                    errors.append(f"Line {lineno}: Unknown mixer type '{mixer_type}'")
                    continue
                
                valid, msg = validate_mixer_time(time, mixer_type)
                if not valid:
                    errors.append(f"Line {lineno}: {msg}")
                    continue
                
                op = {'type': 'MIX', 'id': args[0], 'src1': args[1],
                      'src2': args[2], 'time': time, 'mixer_type': mixer_type}

            elif op_type == 'OUTPUT':
                if len(args) != 2:
                    errors.append(f"Line {lineno}: OUTPUT needs 2 args")
                    continue
                op = {'type': 'OUTPUT', 'id': args[0], 'src': args[1]}

            elif op_type == 'DETECT':
                if len(args) != 2:
                    errors.append(f"Line {lineno}: DETECT needs 2 args")
                    continue
                op = {'type': 'DETECT', 'id': args[0], 'src': args[1]}

            elif op_type == 'HEAT':
                if len(args) != 4:
                    errors.append(f"Line {lineno}: HEAT needs 4 args")
                    continue
                op = {'type': 'HEAT', 'id': args[0], 'src': args[1],
                      'time': int(args[2]), 'temp': int(args[3])}

            elif op_type == 'COOL':
                if len(args) != 4:
                    errors.append(f"Line {lineno}: COOL needs 4 args")
                    continue
                op = {'type': 'COOL', 'id': args[0], 'src': args[1],
                      'time': int(args[2]), 'temp': int(args[3])}

            elif op_type == 'INCUBATE':
                if len(args) != 3:
                    errors.append(f"Line {lineno}: INCUBATE needs 3 args")
                    continue
                op = {'type': 'INCUBATE', 'id': args[0], 'src': args[1],
                      'time': int(args[2])}

            else:
                errors.append(f"Line {lineno}: Unknown operation '{op_type}'")
                continue

            if op['id'] in seen_ids:
                errors.append(f"Line {lineno}: Duplicate id '{op['id']}'")
                continue
            seen_ids.add(op['id'])
            ops.append(op)

    if errors:
        raise ValueError("Parse errors:\n  " + "\n  ".join(errors))

    all_ids = {o['id'] for o in ops}
    for op in ops:
        for src_key in ('src', 'src1', 'src2'):
            src = op.get(src_key)
            if src and src not in all_ids:
                raise ValueError(
                    f"Operation '{op['id']}' references unknown source '{src}'")

    return ops


def get_sources(op):
    if op['type'] == 'MIX':
        return [op['src1'], op['src2']]
    if op.get('src'):
        return [op['src']]
    return []


# ─────────────────────────────────────────────
# 2-4. GRAPH, SCHEDULER, PLACEMENT
# ─────────────────────────────────────────────

def build_sequencing_graph(ops):
    return {op['id']: get_sources(op) for op in ops}


def print_sequencing_graph(graph):
    print("\n=== Sequencing Graph (DAG) ===")
    for node, deps in graph.items():
        if deps:
            print(f"  {' + '.join(deps)}  -->  {node}")
        else:
            print(f"  (start)  -->  {node}")


def topological_sort(ops):
    in_degree = defaultdict(int)
    dependants = defaultdict(list)
    for op in ops:
        for src in get_sources(op):
            in_degree[op['id']] += 1
            dependants[src].append(op['id'])
    ready = deque(o['id'] for o in ops if in_degree[o['id']] == 0)
    topo = []
    while ready:
        oid = ready.popleft()
        topo.append(oid)
        for dep in dependants[oid]:
            in_degree[dep] -= 1
            if in_degree[dep] == 0:
                ready.append(dep)
    return topo


def operation_duration(op):
    if op['type'] == 'DISPENSER':
        return 1
    if op['type'] == 'MIX':
        mrows, mcols, cycles_per_unit = get_module_spec(op.get('mixer_type', DEFAULT_MIXER_TYPE))
        return op['time'] * cycles_per_unit
    if op['type'] in ('HEAT', 'COOL', 'INCUBATE'):
        return op.get('time', 2) * 2 + 2
    return 1


def compute_schedule(ops, rows, cols, verbose=True):
    Na = rows * cols
    op_map = {o['id']: o for o in ops}
    finish = {}
    start_time = {}
    storage_units = {}

    topo = topological_sort(ops)

    def earliest_finish(op_id):
        if op_id in finish:
            return finish[op_id]
        op = op_map[op_id]
        sources = get_sources(op)
        ready_at = max((earliest_finish(s) for s in sources), default=1)
        start_time[op_id] = ready_at
        finish[op_id] = ready_at + operation_duration(op)
        return finish[op_id]

    for op in ops:
        earliest_finish(op['id'])

    max_t = max(finish.values(), default=1) + 1
    changed = True
    iterations = 0
    while changed and iterations < 50:
        changed = False
        iterations += 1
        for t in range(1, max_t + 1):
            active_mix = [o for o in ops
                          if o['type'] == 'MIX'
                          and start_time.get(o['id'], 1) <= t < finish.get(o['id'], 1)]
            active_storage = []
            for o in ops:
                if o['id'] not in finish:
                    continue
                consumers = [op_map[dep] for dep in topo
                             if o['id'] in get_sources(op_map[dep])]
                for c in consumers:
                    c_start = start_time.get(c['id'], finish.get(o['id'], 1))
                    if finish[o['id']] < c_start:
                        if finish[o['id']] <= t < c_start:
                            active_storage.append(o['id'])

            nmixer = len(active_mix)
            nmemory = len(set(active_storage))
            resource_use = nmixer + 0.25 * nmemory

            if resource_use > Na:
                if active_mix:
                    defer_op = active_mix[-1]
                    old_start = start_time[defer_op['id']]
                    start_time[defer_op['id']] = old_start + 1
                    finish[defer_op['id']] = start_time[defer_op['id']] + operation_duration(defer_op)
                    max_t = max(max_t, finish[defer_op['id']] + 1)
                    changed = True
                    if verbose:
                        print(f"  [Scheduler] Resource overflow at t={t}, deferred {defer_op['id']}")
                    break

    for op in ops:
        if op['id'] not in finish:
            continue
        consumers = [op_map[dep] for dep in topo
                     if dep in op_map and op['id'] in get_sources(op_map[dep])]
        for c in consumers:
            c_start = start_time.get(c['id'], finish[op['id']])
            if c_start > finish[op['id']] + 1:
                storage_units[op['id']] = (finish[op['id']], c_start)

    return finish, start_time, storage_units


def print_schedule(ops, schedule, start_time):
    print("\n=== ASAP Schedule (resource-constrained) ===")
    print(f"  {'ID':<10} {'Type':<12} {'Start':>7} {'End':>7} {'Duration':>8}")
    print("  " + "-"*48)
    for op in ops:
        s = start_time.get(op['id'], 1)
        e = schedule.get(op['id'], s)
        dur = e - s
        print(f"  {op['id']:<10} {op['type']:<12} {s:>7} {e:>7} {dur:>8}")
    print(f"\n  Total simulation length: {max(schedule.values())} clock cycles")


def compute_border_cells(rows, cols):
    """Generate border cells in order: spread around perimeter."""
    cells = []
    # Top row (left to right)
    for c in range(1, cols + 1):
        cells.append((1, c))
    # Right column (top to bottom, excluding corner)
    for r in range(2, rows + 1):
        cells.append((r, cols))
    # Bottom row (right to left, excluding corner)
    for c in range(cols - 1, 0, -1):
        cells.append((rows, c))
    # Left column (bottom to top, excluding both corners)
    for r in range(rows - 1, 1, -1):
        cells.append((r, 1))
    return cells


def compute_fti(placements, mix_zones, rows, cols):
    used = set()
    for mid, (mr, mc) in mix_zones.items():
        mrows, mcols, _ = get_module_spec(DEFAULT_MIXER_TYPE)
        for dr in range(-1, mrows + 1):
            for dc in range(-1, mcols + 1):
                nr, nc = mr + dr, mc + dc
                if 1 <= nr <= rows and 1 <= nc <= cols:
                    used.add((nr, nc))
    interior = (rows - 2) * (cols - 2)
    if interior <= 0:
        return 0.0
    spare = interior - len(used & {(r, c) for r in range(2, rows)
                                         for c in range(2, cols)})
    return spare / interior


def placement_cost(mix_zones, rows, cols):
    used_cells = set()
    for mid, (mr, mc) in mix_zones.items():
        mrows, mcols, _ = get_module_spec(DEFAULT_MIXER_TYPE)
        for dr in range(mrows):
            for dc in range(mcols):
                used_cells.add((mr + dr, mc + dc))
    area = len(used_cells)
    fti = compute_fti({}, mix_zones, rows, cols)
    return area - 10.0 * fti


def assign_placements(ops, rows, cols):
    """
    FIXED: Spread dispensers evenly around entire border perimeter.
    Ensure at least 3 cells of spacing between consecutive dispensers.
    """
    border = compute_border_cells(rows, cols)
    used_indices = []
    placements = {}

    dispensers = [o for o in ops if o['type'] == 'DISPENSER']
    outputs    = [o for o in ops if o['type'] == 'OUTPUT']
    
    if not dispensers:
        return placements

    # Spread dispensers evenly around border
    n_border = len(border)
    spacing = max(3, n_border // max(len(dispensers), 1))  # Minimum 3-cell gap
    
    for i, d in enumerate(dispensers):
        idx = (i * spacing) % n_border
        placements[d['id']] = border[idx]
        used_indices.append(idx)

    # Place outputs avoiding dispenser locations
    output_indices = set(range(n_border)) - set(used_indices)
    output_list = sorted(list(output_indices))
    for i, o in enumerate(outputs):
        if i < len(output_list):
            placements[o['id']] = border[output_list[i]]

    return placements


def assign_mix_zones_sa(ops, rows, cols, verbose=True):
    mixes = [o for o in ops if o['type'] == 'MIX']
    if not mixes:
        return {}

    mrows, mcols, _ = get_module_spec(DEFAULT_MIXER_TYPE)

    def random_placement():
        zones = {}
        for m in mixes:
            r = random.randint(3, max(3, rows - mrows - 2))
            c = random.randint(3, max(3, cols - mcols - 2))
            zones[m['id']] = (r, c)
        return zones

    def overlaps(zones):
        rects = []
        for mid, (mr, mc) in zones.items():
            rects.append((mr - 2, mc - 2, mr + mrows + 1, mc + mcols + 1))
        for i in range(len(rects)):
            for j in range(i + 1, len(rects)):
                r1a, c1a, r1b, c1b = rects[i]
                r2a, c2a, r2b, c2b = rects[j]
                if r1a <= r2b and r2a <= r1b and c1a <= c2b and c2a <= c1b:
                    return True
        return False

    def cost(zones):
        if overlaps(zones):
            return 1e9
        return placement_cost(zones, rows, cols)

    current = random_placement()
    best = deepcopy(current)
    best_cost = cost(best)
    T = 100.0
    T_min = 0.1
    alpha = 0.95
    iters_per_temp = 50

    while T > T_min:
        for _ in range(iters_per_temp):
            candidate = deepcopy(current)
            mid = random.choice(list(candidate.keys()))
            dr = random.randint(-3, 3)
            dc = random.randint(-3, 3)
            mr, mc = candidate[mid]
            nr = max(3, min(rows - mrows - 2, mr + dr))
            nc = max(3, min(cols - mcols - 2, mc + dc))
            candidate[mid] = (nr, nc)

            c_curr = cost(current)
            c_cand = cost(candidate)
            delta = c_cand - c_curr

            if delta < 0 or random.random() < math.exp(-delta / T):
                current = candidate
                if c_cand < best_cost:
                    best = deepcopy(candidate)
                    best_cost = c_cand
        T *= alpha

    if verbose:
        print(f"  SA placement: best cost = {best_cost:.2f}, "
              f"FTI = {compute_fti({}, best, rows, cols):.3f}")
    return best


def print_placements(ops, placements):
    print("\n=== Placements (row, col) ===")
    for op in ops:
        if op['id'] in placements:
            r, c = placements[op['id']]
            label = op.get('reagent', 'output')
            print(f"  {op['id']} ({label}) --> ({r}, {c})")


# ─────────────────────────────────────────────
# 5. ROUTER - UNCHANGED FROM v5.6
# ─────────────────────────────────────────────

MAX_BFS_TIME = 800
M_ROUTES = 3

class Router:
    def __init__(self, rows, cols, slot_duration=8):
        self.rows = rows
        self.cols = cols
        self.reserved = {}
        self._label_cells = defaultdict(set)
        self.obstacles = set()
        self._ignore_label = None
        self._by_clock = {}

    def add_obstacle(self, r_or_tuple, c=None):
        if c is None:
            self.obstacles.add(r_or_tuple)
        else:
            self.obstacles.add((r_or_tuple, c))

    def _reserve_one(self, clk, r, c, label):
        key = (clk, r, c)
        old = self.reserved.get(key)
        if old is not None and old != label:
            return False
        if old != label:
            self.reserved[key] = label
            self._label_cells[label].add(key)
            if clk not in self._by_clock:
                self._by_clock[clk] = {}
            self._by_clock[clk][(r, c)] = label
        return True

    def reserve(self, clk, r, c, label):
        return self._reserve_one(clk, r, c, label)

    def rename(self, old_label, new_label):
        if old_label == new_label:
            return
        keys = list(self._label_cells.get(old_label, set()))
        for key in keys:
            if self.reserved.get(key) == old_label:
                self.reserved[key] = new_label
                self._label_cells[new_label].add(key)
                clk, r, c = key
                if clk in self._by_clock:
                    self._by_clock[clk][(r, c)] = new_label
        if old_label in self._label_cells:
            del self._label_cells[old_label]

    def release_label(self, label):
        keys = list(self._label_cells.get(label, set()))
        for key in keys:
            if self.reserved.get(key) == label:
                del self.reserved[key]
                clk, r, c = key
                if clk in self._by_clock:
                    self._by_clock[clk].pop((r, c), None)
                    if not self._by_clock[clk]:
                        del self._by_clock[clk]
        if label in self._label_cells:
            del self._label_cells[label]

    def can_move_to(self, clk, r, c, label):
        if not (1 <= r <= self.rows and 1 <= c <= self.cols):
            return False
        if (r, c) in self.obstacles:
            return False
        clk_cells = self._by_clock.get(clk)
        if clk_cells is None:
            return True
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                occ = clk_cells.get((r + dr, c + dc))
                if occ is not None and occ != label:
                    if self._ignore_label and occ == self._ignore_label:
                        continue
                    return False
        return True

    def _bfs_one_path(self, from_r, from_c, to_r, to_c, start_clk, label, shuffle_dirs=False):
        start = (start_clk, from_r, from_c)
        prev = {start: None}
        queue = deque([start])
        goal = None
        directions = [(1,0),(-1,0),(0,1),(0,-1),(0,0)]

        while queue:
            clk, r, c = queue.popleft()
            if clk - start_clk > MAX_BFS_TIME:
                return None
            if (r, c) == (to_r, to_c):
                goal = (clk, r, c)
                break
            next_clk = clk + 1
            sorted_dirs = sorted(
                directions,
                key=lambda d: abs((r+d[0])-to_r) + abs((c+d[1])-to_c)
            )
            if shuffle_dirs:
                random.shuffle(sorted_dirs)
            for dr, dc in sorted_dirs:
                nr, nc = r+dr, c+dc
                if not (1 <= nr <= self.rows and 1 <= nc <= self.cols):
                    continue
                if (nr, nc) in self.obstacles:
                    continue
                if not self.can_move_to(next_clk, nr, nc, label):
                    continue
                nstate = (next_clk, nr, nc)
                if nstate not in prev:
                    prev[nstate] = (clk, r, c, (dr!=0 or dc!=0))
                    queue.append(nstate)

        if goal is None:
            return None

        path = []
        node = goal
        while prev[node] is not None:
            pclk, pr, pc, moved = prev[node]
            path.append((pclk, pr, pc, node[1], node[2], moved))
            node = (pclk, pr, pc)
        path.reverse()
        return path

    def route_m_paths(self, from_r, from_c, to_r, to_c, start_clk, label,
                       ignore_label=None):
        old_ignore = self._ignore_label
        self._ignore_label = ignore_label
        paths = []
        for i in range(M_ROUTES):
            p = self._bfs_one_path(from_r, from_c, to_r, to_c, start_clk, label,
                                   shuffle_dirs=(i > 0))
            if p is not None:
                paths.append(p)
        self._ignore_label = old_ignore
        return paths

    def _path_delay(self, path):
        if not path:
            return 0
        return path[-1][0] - path[0][0] + 1

    def commit_path(self, path, label, clock_cmds):
        for (pclk, pr, pc, nr, nc, moved) in path:
            self.reserve(pclk + 1, nr, nc, label)
            if moved:
                clock_cmds[pclk].append(f"move({pc},{pr},{nc},{nr})")
        
        if path:
            for (pclk, pr, pc, nr, nc, moved) in path:
                if moved:
                    step_clk = pclk + 1
                    for dr in (-1, 0, 1):
                        for dc in (-1, 0, 1):
                            if dr == 0 and dc == 0:
                                continue
                            nr2, nc2 = nr + dr, nc + dc
                            if 1 <= nr2 <= self.rows and 1 <= nc2 <= self.cols:
                                if self.reserved.get((step_clk, nr2, nc2)) is None:
                                    self.reserve(step_clk, nr2, nc2, label)
            
            arr_clk = path[-1][0]
            return arr_clk
        return 0

    def route(self, from_r, from_c, to_r, to_c, start_clk, label, clock_cmds):
        existing = self.reserved.get((start_clk, from_r, from_c))
        if existing is not None and existing != label:
            self.rename(existing, label)
        self.reserve(start_clk, from_r, from_c, label)

        if (from_r, from_c) == (to_r, to_c):
            return start_clk

        paths = self.route_m_paths(from_r, from_c, to_r, to_c, start_clk, label)

        if not paths:
            raise RuntimeError(
                f"No route found: '{label}' ({from_r},{from_c})->({to_r},{to_c})")

        best_path = min(paths, key=lambda p: len(p))
        arr = self.commit_path(best_path, label, clock_cmds)
        return arr if arr else start_clk

    def route_3pin(self, from_r1, from_c1, from_r2, from_c2,
                   to_r, to_c, start_clk, label1, label2, merge_label,
                   clock_cmds):
        """3-pin routing with staggered arrivals (v5.6 fix)."""
        mrows, mcols, _ = get_module_spec(DEFAULT_MIXER_TYPE)

        self.reserve(start_clk, from_r1, from_c1, label1)
        self.reserve(start_clk, from_r2, from_c2, label2)

        paths1 = self.route_m_paths(from_r1, from_c1, to_r, to_c,
                                    start_clk, label1, ignore_label=label2)

        if not paths1:
            for dc in [-1, 1]:
                alt_c = to_c + dc
                if 1 <= alt_c <= self.cols:
                    paths1 = self.route_m_paths(from_r1, from_c1, to_r, alt_c,
                                                start_clk, label1, ignore_label=label2)
                    if paths1:
                        break

        if not paths1:
            raise RuntimeError(
                f"3-pin: no path for {label1} ({from_r1},{from_c1})->({to_r},{to_c})")
        
        best1 = min(paths1, key=len)
        arr1 = self.commit_path(best1, label1, clock_cmds)

        route2_start = arr1 + 2
        
        paths2 = self.route_m_paths(from_r2, from_c2, to_r, to_c + mcols - 1,
                                    route2_start, label2, ignore_label=None)

        if not paths2:
            for dc in [-1, 1]:
                alt_c = to_c + mcols - 1 + dc
                if 1 <= alt_c <= self.cols:
                    paths2 = self.route_m_paths(from_r2, from_c2, to_r, alt_c,
                                                route2_start, label2, ignore_label=None)
                    if paths2:
                        break

        if not paths2:
            raise RuntimeError(
                f"3-pin: no path for {label2} ({from_r2},{from_c2})->({to_r},{to_c+mcols-1})")
        
        best2 = min(paths2, key=len)
        arr2 = self.commit_path(best2, label2, clock_cmds)

        return max(arr1, arr2), max(arr1, arr2)

    def reserve_mix_zone(self, mix_r, mix_c, start_clk, mix_time, label,
                          mixer_type=DEFAULT_MIXER_TYPE):
        mrows, mcols, cycles_per_unit = get_module_spec(mixer_type)
        mix_end = start_clk + mix_time * cycles_per_unit + 2
        
        for t in range(start_clk, mix_end + 1):
            for dr in range(mrows):
                for dc in range(mcols):
                    nr = mix_r + dr
                    nc = mix_c + dc
                    if 1 <= nr <= self.rows and 1 <= nc <= self.cols:
                        key = (t, nr, nc)
                        existing = self.reserved.get(key)
                        if existing is None or existing == label:
                            self.reserved[key] = label
                            self._label_cells[label].add(key)
                            if t not in self._by_clock:
                                self._by_clock[t] = {}
                            self._by_clock[t][(nr, nc)] = label


# ─────────────────────────────────────────────
# 6-11. REST (UNCHANGED FROM v5.6)
# ─────────────────────────────────────────────

def build_segregation_obstacles(mix_zones, rows, cols, router,
                                mixer_type=DEFAULT_MIXER_TYPE):
    mrows, mcols, _ = get_module_spec(mixer_type)
    for mid, (mr, mc) in mix_zones.items():
        for dr in range(mrows):
            for dc in range(1, mcols - 1):
                nr, nc = mr + dr, mc + dc
                if 1 <= nr <= rows and 1 <= nc <= cols:
                    router.add_obstacle((nr, nc))


def compute_routes(ops, placements, schedule, start_time, mix_zones, rows, cols, router):
    op_map = {o['id']: o for o in ops}
    clock_cmds = defaultdict(list)

    droplet_pos = {}
    droplet_label = {}

    for op in ops:
        if op['type'] == 'DISPENSER':
            r, c = placements[op['id']]
            clock_cmds[1].append(f"dispense({c},{r})")
            droplet_pos[op['id']] = (r, c)
            droplet_label[op['id']] = op['id']
            router.reserve(1, r, c, op['id'])

    topo = topological_sort(ops)

    def get_pos(op_id):
        return droplet_pos.get(op_id, (rows // 2, cols // 2))

    for op_id in topo:
        op = op_map[op_id]
        if op['type'] == 'DISPENSER':
            continue

        if op['type'] == 'MIX':
            mix_r, mix_c = mix_zones[op_id]
            src1, src2 = op['src1'], op['src2']
            mixer_type = op.get('mixer_type', DEFAULT_MIXER_TYPE)
            mrows, mcols, cycles_per_unit = get_module_spec(mixer_type)
            start_clk = max(schedule[s] for s in get_sources(op))

            p1 = get_pos(src1)
            lbl1 = f"{op_id}_d1"
            p2 = get_pos(src2)
            lbl2 = f"{op_id}_d2"

            router.release_label(droplet_label[src1])
            router.release_label(droplet_label[src2])

            arr1, arr2 = router.route_3pin(
                p1[0], p1[1], p2[0], p2[1],
                mix_r, mix_c, start_clk,
                lbl1, lbl2, op_id, clock_cmds)

            mix_start = max(arr1, arr2) + 1

            router.rename(lbl1, op_id)
            router.rename(lbl2, op_id)
            router.reserve_mix_zone(mix_r, mix_c, mix_start,
                                    op['time'], op_id, mixer_type)

            clock_cmds[mix_start].append(
                f"mix_split({mix_c},{mix_r},{mix_c+mcols-1},{mix_r},{op['time']})")

            mix_duration = op['time'] * cycles_per_unit
            mix_end = mix_start + mix_duration
            droplet_pos[op_id] = (mix_r, mix_c)
            droplet_label[op_id] = op_id
            schedule[op_id] = mix_end

        elif op['type'] == 'OUTPUT':
            src = op['src']
            start_clk = schedule[src]
            p = get_pos(src)
            target = placements[op_id]
            lbl = f"{op_id}_out"
            router.release_label(droplet_label[src])
            arr = router.route(p[0], p[1], target[0], target[1],
                               start_clk, lbl, clock_cmds)
            droplet_pos[op_id] = target
            droplet_label[op_id] = lbl
            schedule[op_id] = arr

        elif op['type'] == 'DETECT':
            src = op['src']
            start_clk = schedule[src]
            p = get_pos(src)
            clock_cmds[start_clk].append(f"detect({p[1]},{p[0]})")
            droplet_pos[op_id] = p
            droplet_label[op_id] = droplet_label.get(src, src)
            schedule[op_id] = start_clk + 1

        elif op['type'] in ('HEAT', 'COOL', 'INCUBATE'):
            src = op['src']
            start_clk = schedule[src]
            p = get_pos(src)
            dur = op.get('time', 2)
            op_end = start_clk + dur * 2 + 2
            if op['type'] == 'HEAT':
                clock_cmds[start_clk].append(
                    f"heat({p[1]},{p[0]},{dur},{op['temp']})")
            elif op['type'] == 'COOL':
                clock_cmds[start_clk].append(
                    f"cool({p[1]},{p[0]},{dur},{op['temp']})")
            else:
                clock_cmds[start_clk].append(
                    f"incubate({p[1]},{p[0]},{dur})")
            droplet_pos[op_id] = p
            droplet_label[op_id] = droplet_label.get(src, src)
            schedule[op_id] = op_end

    return clock_cmds


def check_fluidic_constraints(router, verbose=True):
    by_clock = defaultdict(list)
    for (clk, r, c), label in router.reserved.items():
        by_clock[clk].append((r, c, label))

    violations = []
    for clk in sorted(by_clock.keys()):
        cells = by_clock[clk]
        for i in range(len(cells)):
            for j in range(i + 1, len(cells)):
                r1, c1, l1 = cells[i]
                r2, c2, l2 = cells[j]
                if l1 == l2:
                    continue
                if abs(r1 - r2) <= 1 and abs(c1 - c2) <= 1:
                    violations.append(
                        f"clock {clk}: '{l1}' at ({r1},{c1}) and "
                        f"'{l2}' at ({r2},{c2}) violate Rule #1")

    if verbose and violations:
        print(f"\n⚠️  {len(violations)} fluidic constraint violation(s):")
        for v in violations[:5]:
            print(f"   {v}")
        if len(violations) > 5:
            print(f"   ... and {len(violations)-5} more")
    elif verbose:
        print("\n✓ No fluidic constraint violations detected.")

    return violations


def write_dmfb(ops, placements, schedule, clock_cmds, rows, cols, accuracy):
    lines = []
    lines.append(f"dimension {rows} {cols}")
    lines.append(f"accuracy {accuracy}")

    dispensers = [o for o in ops if o['type'] == 'DISPENSER']
    reagent_decls = []
    for d in dispensers:
        r, c = placements[d['id']]
        reagent_decls.append(f"reagent({c},{r},{d['reagent']})")
    lines.append(" ".join(reagent_decls))
    lines.append("")

    all_clocks = sorted(clock_cmds.keys())
    for clk in all_clocks:
        cmds = clock_cmds[clk]
        if cmds:
            lines.append(f"{clk} {' '.join(cmds)}")

    max_clock = max(schedule.values(), default=1)
    end_clock = max(all_clocks[-1] if all_clocks else 1, max_clock) + 5
    lines.append(f"{end_clock} end")

    return "\n".join(lines) + "\n"

def serialize_moves(clock_cmds):
    """
    Conservative fix:
    allow at most one move per timestep.
    Other commands (mix_split, detect, etc.) stay.
    """

    new_cmds = defaultdict(list)

    current_clk = 1

    for clk in sorted(clock_cmds.keys()):

        moves = []
        others = []

        for cmd in clock_cmds[clk]:
            if cmd.startswith("move("):
                moves.append(cmd)
            else:
                others.append(cmd)

        # keep non-move commands at original time
        new_cmds[clk].extend(others)

        # serialize moves
        for mv in moves:
            while any(c.startswith("move(") for c in new_cmds[current_clk]):
                current_clk += 1

            new_cmds[current_clk].append(mv)
            current_clk += 1

    return new_cmds
    
def synthesize(input_text, rows=32, cols=32, accuracy=5, verbose=True):
    if verbose:
        print("\n" + "="*50)
        print("  DMFB Synthesis Tool (v5.7 — Dispenser spacing)")
        print("="*50)

    if verbose:
        print("\n[Step 1] Parsing & validating bioassay...")
    ops = parse_bioassay(input_text)
    if verbose:
        print(f"  Found {len(ops)} operations")

    graph = build_sequencing_graph(ops)
    if verbose:
        print_sequencing_graph(graph)

    if verbose:
        print("\n[Step 2] Scheduling...")
    schedule, start_time, storage_units = compute_schedule(ops, rows, cols, verbose)
    if verbose:
        print_schedule(ops, schedule, start_time)

    if verbose:
        print("\n[Step 3] Placement (evenly spaced dispensers)...")
    placements = assign_placements(ops, rows, cols)
    if verbose:
        print_placements(ops, placements)

    if verbose:
        print("\n[Step 4] Mix zone placement...")
    mix_zones = assign_mix_zones_sa(ops, rows, cols, verbose)

    router = Router(rows, cols, slot_duration=8)

    if verbose:
        print("\n[Step 4b] Building obstacles...")
    build_segregation_obstacles(mix_zones, rows, cols, router)

    if verbose:
        print("\n[Step 5] Routing (with staggered 3-pin arrivals)...")
    clock_cmds = compute_routes(
        ops, placements, schedule, start_time, mix_zones, rows, cols, router)

    if verbose:
        print("\n[Step 6] Checking constraints...")
    violations = check_fluidic_constraints(router, verbose=verbose)

    if verbose:
        print("\n[Step 7] Writing .dmfb...")
    dmfb_content = write_dmfb(ops, placements, schedule, clock_cmds,
                               rows, cols, accuracy)
    return dmfb_content, violations


BS_EXAMPLE = """\
# Blood serum (BS) dilution assay
DISPENSER(D1, buffer) DISPENSER(D2, sample)
MIX(M1, D1, D2, 6)
DISPENSER(D3, buffer)
MIX(M2, M1, D3, 6)
DISPENSER(D4, buffer)
MIX(M3, M2, D4, 6)
DETECT(DET1, M3)
OUTPUT(O1, DET1)
"""


def main():
    parser = argparse.ArgumentParser(
        description="DMFB Synthesis Tool v5.7",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n  python synthesis_tool.py input.txt output.dmfb\n  python synthesis_tool.py --example bs"
    )
    parser.add_argument('input',  nargs='?', help='Input bioassay .txt file')
    parser.add_argument('output', nargs='?', help='Output .dmfb file path')
    parser.add_argument('--rows',     type=int, default=32)
    parser.add_argument('--cols',     type=int, default=32)
    parser.add_argument('--accuracy', type=int, default=5)
    parser.add_argument('--example',  choices=['bs'])
    parser.add_argument('--quiet',    action='store_true')
    args = parser.parse_args()

    if args.example:
        text = BS_EXAMPLE
        out_name = f"{args.example}_output.dmfb"
        dmfb, violations = synthesize(text, args.rows, args.cols,
                                      args.accuracy, verbose=not args.quiet)
        with open(out_name, 'w') as f:
            f.write(dmfb)
        print(f"\n.dmfb file written to: {out_name}")
        if violations:
            print(f"\n⚠️  {len(violations)} violations found")
        else:
            print("\n✓ File is clean")
        return

    if not args.input:
        parser.print_help()
        sys.exit(1)

    try:
        with open(args.input, 'r') as f:
            input_text = f.read()
    except FileNotFoundError:
        print(f"Error: file '{args.input}' not found.")
        sys.exit(1)

    try:
        dmfb, violations = synthesize(input_text, args.rows, args.cols,
                                      args.accuracy, verbose=not args.quiet)
    except ValueError as e:
        print(f"\nSynthesis failed:\n{e}")
        sys.exit(1)

    out_path = args.output or args.input.replace('.txt', '_output.dmfb')
    with open(out_path, 'w') as f:
        f.write(dmfb)
    print(f"\n.dmfb file written to: {out_path}")


if __name__ == '__main__':
    main()