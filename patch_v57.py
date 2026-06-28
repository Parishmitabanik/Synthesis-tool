from pathlib import Path

p = Path("synthesis_tool.py")

txt = p.read_text()

# =========================================================
# 1. DEFAULT MIXER
# =========================================================

txt = txt.replace(
    "DEFAULT_MIXER_TYPE = '2x4'",
    "DEFAULT_MIXER_TYPE = '2x3'"
)

# =========================================================
# 2. DIMENSION ORDER
# =========================================================

txt = txt.replace(
    'lines.append(f"dimension {rows} {cols}")',
    'lines.append(f"dimension {cols} {rows}")'
)

# =========================================================
# 3. mix_split PHYSICAL CYCLES
# =========================================================

old = '''clock_cmds[mix_start].append(
                f"mix_split({mix_c},{mix_r},{mix_c+mcols-1},{mix_r},{op['time']})")'''

new = '''physical_cycles = op['time'] * cycles_per_unit

            clock_cmds[mix_start].append(
                f"mix_split("
                f"{mix_c},{mix_r},"
                f"{mix_c+mcols-1},{mix_r},"
                f"{physical_cycles})"
            )'''

txt = txt.replace(old, new)

# =========================================================
# 4. MIX SETTLING DELAY
# =========================================================

txt = txt.replace(
    "mix_start = max(arr1, arr2) + 1",
    "mix_start = max(arr1, arr2) + 4"
)

# =========================================================
# 5. ACTUALLY CALL serialize_moves()
# =========================================================

old = '''clock_cmds = compute_routes(
        ops, placements, schedule, start_time, mix_zones, rows, cols, router)'''

new = '''clock_cmds = compute_routes(
        ops, placements, schedule, start_time,
        mix_zones, rows, cols, router)

    clock_cmds = serialize_moves(clock_cmds)'''

txt = txt.replace(old, new)

# =========================================================
# 6. LAZY DISPENSING
# =========================================================

old = '''for op in ops:
        if op['type'] == 'DISPENSER':
            r, c = placements[op['id']]
            clock_cmds[1].append(f"dispense({c},{r})")
            droplet_pos[op['id']] = (r, c)
            droplet_label[op['id']] = op['id']
            router.reserve(1, r, c, op['id'])'''

new = '''for op in ops:

        if op['type'] != 'DISPENSER':
            continue

        consumers = [
            o for o in ops
            if op['id'] in get_sources(o)
        ]

        if consumers:
            first_use = min(
                start_time[c['id']]
                for c in consumers
            )
        else:
            first_use = 1

        t_disp = max(1, first_use - 40)

        r, c = placements[op['id']]

        clock_cmds[t_disp].append(
            f"dispense({c},{r})"
        )

        droplet_pos[op['id']] = (r, c)
        droplet_label[op['id']] = op['id']

        for t in range(t_disp, t_disp + 10):
            router.reserve(t, r, c, op['id'])'''

txt = txt.replace(old, new)

# =========================================================
# WRITE OUTPUT
# =========================================================

out = Path("synthesis_tool_v58.py")
out.write_text(txt)

print("✅ Wrote synthesis_tool_v58.py")