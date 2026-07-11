#!/usr/bin/env python3
from __future__ import annotations
import re, sys
from collections import defaultdict, deque
from verify_biochipv import BiochipV


def parse_dmfb(path):
    with open(path) as f:
        lines=[ln.rstrip('\n') for ln in f]
    return lines

class TokenRec:
    def __init__(self, tid, origin_type, origin_data):
        self.id=tid
        self.origin_type=origin_type     # 'dispense' or 'mix'
        self.origin_data=origin_data     # pos or parent mix id
        self.path=[]                     # list of positions incl source and subsequent dests
        self.terminal=None               # ('mix', mid) / ('waste',pos) / ('output',pos)

class MixRec:
    def __init__(self, mid, secs, inputs, out_pos):
        self.id=mid
        self.secs=secs
        self.inputs=inputs
        self.out_pos=out_pos
        self.outputs=[]
        self.children=[]
        self.rank=secs


def group_instructions(line):
    disp=[]; mix=[]; moves=[]; waste=[]; out=[]
    instr_line = re.compile(r'[a-z_]+\s*\([\s*\d+\s*,\s*]+\s*\w+\s*\)').findall(line)
    for instr in instr_line:
        obj = re.compile(r'[a-z_]+').match(instr)
        opcode = obj.group()
        operands = re.compile(r'(\d+|\w+)').findall(instr[obj.end()+1:])
        if opcode=='dispense' and len(operands)==2:
            disp.append((int(operands[0]), int(operands[1])))
        elif opcode=='mix_split' and len(operands)==5:
            mix.append(tuple(map(int, operands)))
        elif opcode=='move' and len(operands)==4:
            moves.append(tuple(map(int, operands)))
        elif opcode=='waste' and len(operands)==2:
            waste.append((int(operands[0]), int(operands[1])))
        elif opcode=='output' and len(operands)==2:
            out.append((int(operands[0]), int(operands[1])))
    return disp,mix,moves,waste,out


def init_sim_from_header(lines):
    dim = lines[0].split(); rows, cols = int(dim[1]), int(dim[2])
    acc = int(lines[1].split()[1])
    sim = BiochipV(rows, cols, acc)
    header = lines[2]
    instr_line = re.compile(r'[a-z_]+\s*\([\s*\d+\s*,\s*]+\s*\w+\s*\)').findall(header)
    for instr in instr_line:
        obj = re.compile(r'[a-z_]+').match(instr)
        opcode = obj.group()
        operands = re.compile(r'(\d+|\w+)').findall(instr[obj.end()+1:])
        if opcode == 'reagent':
            pts=[]
            for i in range(0, len(operands)-1, 2):
                pts.append((int(operands[i]), int(operands[i+1])))
            sim.verify_set_reagent_reservior(pts, operands[-1])
        elif opcode == 'waste_reservoir':
            sim.verify_set_resevior(int(operands[0]), int(operands[1]), 'waste_reservoir')
        elif opcode == 'output_reservoir':
            sim.verify_set_resevior(int(operands[0]), int(operands[1]), 'op_reservoir')
    return sim


def extract_skeleton(dmfb_path):
    lines = parse_dmfb(dmfb_path)
    sim = init_sim_from_header(lines)
    op_lines=[ln for ln in lines[4:] if ln.strip() and not ln.strip().endswith(' end')]

    tokens={}
    mixes={}
    pos_to_token={}
    active=[]
    next_tok=1
    next_mix=1

    def expire_until(t):
        nonlocal next_tok
        rem=[]
        for rec in active:
            if rec['avail_t'] <= t:
                for idx,pos in enumerate(rec['out_pos']):
                    tid=f"T{next_tok}"; next_tok+=1
                    tr=TokenRec(tid,'mix',rec['id'])
                    tr.path=[pos]
                    tokens[tid]=tr
                    rec['mix'].outputs.append(tid)
                    pos_to_token[pos]=tid
            else:
                rem.append(rec)
        return rem

    for line in op_lines:
        t=int(re.match(r'\d+', line).group())
        active = expire_until(t)
        sim.delete_expired_mixers(t)
        disp,mix,moves,waste,out=group_instructions(line)
        # mix first
        for r1,c1,r2,c2,secs in mix:
            ta=(r1,c1); tb=(r2,c2)
            ida=pos_to_token.pop(ta)
            idb=pos_to_token.pop(tb)
            mrec=MixRec(f"M{next_mix}", secs, [ida,idb], [ta,tb])
            next_mix += 1
            mixes[mrec.id]=mrec
            tokens[ida].terminal=('mix', mrec.id)
            tokens[idb].terminal=('mix', mrec.id)
            active.append({'id':mrec.id,'avail_t':t+secs+1,'out_pos':[ta,tb],'mix':mrec})
        # moves
        for a,b,c,d in moves:
            tid=pos_to_token.pop((a,b))
            pos_to_token[(c,d)]=tid
            tokens[tid].path.append((c,d))
        # dispenses
        for r,c in disp:
            tid=f"T{next_tok}"; next_tok+=1
            tr=TokenRec(tid,'dispense',(r,c))
            tr.path=[(r,c)]
            tokens[tid]=tr
            pos_to_token[(r,c)] = tid
        # waste/output
        for pos in waste:
            tid=pos_to_token.pop(pos)
            tokens[tid].terminal=('waste', pos)
        for pos in out:
            tid=pos_to_token.pop(pos)
            tokens[tid].terminal=('output', pos)
        sim.verify_line(line)
    # build child links and ranks
    for tid, tr in tokens.items():
        if tr.terminal and tr.terminal[0]=='mix':
            mixes[tr.terminal[1]].children.append(tid)
    # map token to produced by mix for dependencies
    parent_mix_of_token={}
    for mid,mr in mixes.items():
        for tid in mr.outputs:
            parent_mix_of_token[tid]=mid
    # build mix DAG children
    child_mix=defaultdict(list)
    for mid,mr in mixes.items():
        for otid in mr.outputs:
            term = tokens[otid].terminal
            if term and term[0]=='mix':
                child_mix[mid].append(term[1])
    # ranks
    for mid in list(mixes.keys())[::-1]:
        pass
    changed=True
    while changed:
        changed=False
        for mid,mr in mixes.items():
            nr = mr.secs + max([mixes[ch].rank for ch in child_mix[mid]] + [0])
            if nr != mr.rank:
                mr.rank=nr; changed=True
    return lines[:4], tokens, mixes, child_mix


def build_line(t, mix_ops, move_ops, disp_ops, waste_ops, out_ops):
    ops=[]
    ops += [f'mix_split({a},{b},{c},{d},{m})' for (a,b,c,d,m) in mix_ops]
    ops += [f'move({a},{b},{c},{d})' for (a,b,c,d) in move_ops]
    ops += [f'dispense({a},{b})' for (a,b) in disp_ops]
    ops += [f'waste({a},{b})' for (a,b) in waste_ops]
    ops += [f'output({a},{b})' for (a,b) in out_ops]
    return f"{t} " + ' '.join(ops)


def cheb(a, b):
    return max(abs(a[0]-b[0]), abs(a[1]-b[1]))


def pathfind_static(sim, src, dst):
    if src == dst:
        return [src]
    occ = set(occupied_positions(sim))
    occ.discard(src)
    if dst in occ:
        return []
    safe = set()
    occ_list = list(occ)
    for r in range(1, sim.row+1):
        for c in range(1, sim.col+1):
            if (r,c) in occ and (r,c) != dst:
                continue
            ok = True
            for p in occ_list:
                if cheb((r,c), p) <= 1:
                    ok = False
                    break
            if ok or (r,c) == dst:
                safe.add((r,c))
    q = deque([(src, [src])])
    seen = {src}
    while q:
        cur, path = q.popleft()
        # Manhattan-only: DMFB droplets can only step to one of the 4
        # orthogonally-adjacent cells per tick, never diagonally.
        for nr, nc in ((cur[0]-1,cur[1]),(cur[0]+1,cur[1]),(cur[0],cur[1]-1),(cur[0],cur[1]+1)):
            nxt = (nr, nc)
            if not (1 <= nr <= sim.row and 1 <= nc <= sim.col):
                continue
            if nxt != dst and nxt not in safe:
                continue
            if nxt in seen:
                continue
            if nxt == dst:
                return path + [nxt]
            seen.add(nxt)
            q.append((nxt, path + [nxt]))
    return []


def occupied_positions(sim):
    occ=[]
    for r in range(1, sim.row+1):
        for c in range(1, sim.col+1):
            if sim.biochip[r][c].id is not None:
                occ.append((r,c))
    # --- mixer physical-footprint protection -----------------------------
    # verify_biochipv.BiochipV only ever marks the two *endpoint* cells of
    # an active mixer as occupied (droplet_in_active_mixers / check_FC).
    # But the real Tkinter renderer (Biochip.py, __mixer_14/__mixer_41)
    # physically animates the droplet back and forth through EVERY cell
    # of the 1x4 / 4x1 zone between those endpoints for as long as the
    # mixer is active. If this router sends some unrelated droplet's move
    # through or next to that zone while the mixer is running, the move
    # is perfectly legal under BiochipV's simplified FC model (so
    # verify_dmfb happily accepts the compacted file) but collides with
    # the GUI's own oval bookkeeping - which is exactly the "droplets
    # stop moving" freeze reported at runtime. So we protect the whole
    # mixer rectangle (+1 cell halo, matching the FC rule everywhere
    # else) here, in the router's own occupancy view, even though the
    # verifier itself doesn't require it.
    for m in sim.mixer_table:
        if m.mtype == 14:
            r = m.d1r
            lo, hi = min(m.d1c, m.d2c), max(m.d1c, m.d2c)
        elif m.mtype == 41:
            c = m.d1c
            lo, hi = min(m.d1r, m.d2r), max(m.d1r, m.d2r)
        else:
            continue
        for pos in range(lo, hi + 1):
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if m.mtype == 14:
                        occ.append((r + dr, pos + dc))
                    else:
                        occ.append((pos + dr, c + dc))
    return occ


def compact(safe_dmfb, out_dmfb, dispense_frontier=1, max_onchip=8):
    header, token_defs, mix_defs, child_mix = extract_skeleton(safe_dmfb)
    sim = init_sim_from_header(header)

    token_state={}
    for tid,tr in token_defs.items():
        token_state[tid]={
            'def': tr,
            'onchip': False,
            'available': tr.origin_type=='dispense',
            'pos': None,
            'goal': tr.path[-1],
            'done': False,
            'consumed': False,
        }
    mix_state={}
    mix_order = sorted(mix_defs.keys(), key=lambda m: int(m[1:]))
    mix_order_idx = dict((m, i) for i, m in enumerate(mix_order))
    for mid,mr in mix_defs.items():
        mix_state[mid]={'def':mr,'started':False,'start_t':None,'avail_t':None,'done':False}

    # for spawn tokens, availability controlled by parent mix completion
    lines=[]
    t=1
    MAXT=3000
    # Tunable concurrency knobs. Higher values let more mixes'
    # inputs be dispensed/routed at once and more droplets occupy
    # the chip simultaneously, which is what actually lets
    # independent branches of the assay DAG overlap in time instead
    # of being serialized. synthesis_tool_v58.py now sweeps several
    # (dispense_frontier, max_onchip) combinations and keeps the
    # best one that still verifies cleanly, so these are just the
    # conservative defaults used when compact() is called directly.
    DISPENSE_FRONTIER = dispense_frontier
    MAX_ONCHIP = max_onchip
    # STALL_LIMIT bounds how many consecutive "nothing scheduled" loop
    # iterations we tolerate before concluding the greedy router has
    # deadlocked (e.g. two mix operations whose baseline coordinates
    # reuse the same physical site, or two droplets whose goal cells
    # are mutually FC-blocking). Without this, the loop silently grinds
    # all the way to MAXT and writes a truncated/broken schedule.
    STALL_LIMIT = 200
    stall_iters = 0
    deadlocked = False

    # LIVELOCK_LIMIT bounds a *different* failure mode than STALL_LIMIT.
    # STALL_LIMIT only fires when a loop iteration emits nothing at all.
    # But the greedy router can also get two (or more) tokens stuck
    # fighting over the same lane/chokepoint: each tick it re-runs
    # pathfind_static() from scratch with no memory of history, so a
    # token can be routed one step forward, then next tick routed one
    # step back because the "forward" cell is now blocked by the other
    # contending token - forever. That's real movement every tick (so
    # stall_iters keeps getting reset to 0), but zero net progress: no
    # token gets closer to its goal, no mixer starts/finishes, nothing
    # ever reaches waste/output. Left unchecked this grinds silently to
    # MAXT and writes a "successful" file that never performs its final
    # output/waste ops (this is exactly what happened on graph3: two
    # droplets shuttled back and forth for ~2900 ticks near the waste
    # port and the schedule was cut off before ever reaching output()).
    #
    # We detect this by tracking a monotonic progress signature -
    # (terminals reached, mixers finished, -distance remaining to goals)
    # - and bailing out if it hasn't strictly improved in LIVELOCK_LIMIT
    # consecutive productive iterations.
    LIVELOCK_LIMIT = 300
    best_progress = None
    no_progress_iters = 0

    def token_prio(tid):
        tr=token_defs[tid]
        term=tr.terminal
        if term[0]=='mix':
            return 1000 + mix_defs[term[1]].rank
        if term[0]=='output':
            return 900
        return 100

    def progress_signature():
        # Monotonic-when-things-are-actually-progressing signature.
        # done_terminals and mixes_done can only go up over the life of
        # the run. dist_remaining can go up or down tick to tick (e.g.
        # while one token detours around another), which is fine - what
        # matters is whether the *combined* signature ever beats its
        # previous best. Pure oscillation (a token bounces between two
        # cells) returns dist_remaining to the same value every 2 ticks
        # and never sets a new best, so no_progress_iters climbs
        # steadily instead of being reset every other tick.
        done_terminals = sum(1 for st in token_state.values() if st['done'])
        mixes_done = sum(1 for ms in mix_state.values() if ms['done'])
        dist_remaining = sum(
            cheb(st['pos'], st['goal'])
            for st in token_state.values()
            if st['onchip'] and not st['consumed'] and not st['done'] and st['pos'] is not None
        )
        return (done_terminals, mixes_done, -dist_remaining)

    def expire_to(cur_t):
        for mid,ms in mix_state.items():
            if ms['started'] and not ms['done'] and ms['avail_t'] <= cur_t:
                ms['done']=True
                for otid in ms['def'].outputs:
                    st=token_state[otid]
                    st['available']=True
                    st['onchip']=True
                    st['pos']=st['def'].path[0]

    while t < MAXT:
        expire_to(t)
        sim.delete_expired_mixers(t)

        mix_ops=[]; move_ops=[]; disp_ops=[]; waste_ops=[]; out_ops=[]

        # start mixes whose input tokens are present at endpoints
        ready=[]
        for mid,ms in mix_state.items():
            if ms['started']:
                continue
            inp=ms['def'].inputs
            if all(token_state[it]['onchip'] and not token_state[it]['consumed'] and token_state[it]['pos'] == token_state[it]['goal'] for it in inp):
                ready.append(mid)
        ready.sort(key=lambda m: mix_order_idx[m])
        for mid in ready:
            mr=mix_defs[mid]
            a=token_state[mr.inputs[0]]['pos']; b=token_state[mr.inputs[1]]['pos']
            mix_ops.append((a[0],a[1],b[0],b[1],mr.secs))
            mix_state[mid]['started']=True
            mix_state[mid]['start_t']=t
            mix_state[mid]['avail_t']=t+mr.secs+1
            for it in mr.inputs:
                token_state[it]['consumed']=True
                token_state[it]['onchip']=False

        occ = set(occupied_positions(sim))

        # candidate moves on fixed paths using exact verifier-derived parallel conditions
        candidates=[]
        for tid,st in token_state.items():
            if not st['onchip'] or st['consumed'] or st['done']:
                continue
            if st['pos'] != st['goal']:
                path = pathfind_static(sim, st['pos'], st['goal'])
                if len(path) >= 2:
                    src=st['pos']; dst=path[1]
                    candidates.append((token_prio(tid), tid, src, dst))
        candidates.sort(reverse=True)
        chosen_dsts=[]
        for _, tid, src, dst in candidates:
            if dst in occ:
                continue
            ok = True
            for q in occ:
                if q == src:
                    continue
                if cheb(dst, q) <= 1:
                    ok = False
                    break
            if not ok:
                continue
            for q in chosen_dsts:
                if cheb(dst, q) <= 1:
                    ok = False
                    break
            if not ok:
                continue
            move_ops.append((src[0], src[1], dst[0], dst[1]))
            chosen_dsts.append(dst)

        # dispenses after moves; source must be empty and FC-safe after moves commit
        moved_srcs = set((a,b) for a,b,_,_ in move_ops)
        occ_after_moves = set([p for p in occ if p not in moved_srcs]) | set(chosen_dsts)
        onchip_count = sum(1 for st in token_state.values() if st['onchip'])
        frontier = []
        for mid, ms in mix_state.items():
            if ms['started']:
                continue
            ok = True
            for it in ms['def'].inputs:
                tr_in = token_defs[it]
                st_in = token_state[it]
                if tr_in.origin_type == 'mix' and not (st_in['available'] or st_in['onchip'] or st_in['consumed']):
                    ok = False
                    break
            if ok:
                frontier.append(mid)
        frontier.sort(key=lambda m: mix_order_idx[m])
        frontier = set(frontier[:DISPENSE_FRONTIER])

        disp_cands=[]
        if onchip_count < MAX_ONCHIP:
            for tid,st in token_state.items():
                tr=st['def']
                if tr.origin_type!='dispense' or st['onchip'] or st['consumed'] or st['done'] or not st['available']:
                    continue
                term_mix = tr.terminal[1]
                if term_mix not in frontier:
                    continue
                src=tr.path[0]
                disp_cands.append((token_prio(tid), tid, src))
        disp_cands.sort(reverse=True)
        chosen_disp=[]
        used_src=set()
        for _, tid, src in disp_cands:
            if src in used_src or src in occ_after_moves:
                continue
            ok = True
            for q in occ_after_moves:
                if cheb(src, q) <= 1:
                    ok = False
                    break
            if not ok:
                continue
            for q in chosen_disp:
                if cheb(src, q) <= 1:
                    ok = False
                    break
            if not ok:
                continue
            disp_ops.append(src)
            chosen_disp.append(src)
            used_src.add(src)
            if onchip_count + len(chosen_disp) >= MAX_ONCHIP:
                break

        # predicted positions after moves/dispenses
        pred_pos={}
        for tid,st in token_state.items():
            if st['onchip']:
                pred_pos[tid]=st['pos']
        for a,b,c,d in move_ops:
            for tid,pos in list(pred_pos.items()):
                if pos==(a,b):
                    pred_pos[tid]=(c,d)
                    break
        for src in disp_ops:
            for tid,st in token_state.items():
                if st['def'].origin_type=='dispense' and (not st['onchip']) and st['def'].path[0]==src and st['available'] and not st['consumed'] and not st['done']:
                    pred_pos[tid]=src
                    break

        # waste/output terminals
        terminal_cands=[]
        for tid,st in token_state.items():
            tr=st['def']
            if tr.terminal is None or tid not in pred_pos:
                continue
            if st['pos'] == st['goal']:
                if tr.terminal[0]=='waste' and pred_pos[tid]==tr.terminal[1]:
                    terminal_cands.append((token_prio(tid), 'waste', tid, pred_pos[tid]))
                if tr.terminal[0]=='output' and pred_pos[tid]==tr.terminal[1]:
                    terminal_cands.append((token_prio(tid), 'output', tid, pred_pos[tid]))
        terminal_cands.sort(reverse=True)
        for _, kind, tid, pos in terminal_cands:
            if kind=='waste':
                waste_ops.append(pos)
            else:
                out_ops.append(pos)

        if mix_ops or move_ops or disp_ops or waste_ops or out_ops:
            line=build_line(t, mix_ops, move_ops, disp_ops, waste_ops, out_ops)
            sim.verify_line(line)
            lines.append(line)
            # updates
            for src in disp_ops:
                for tid,st in token_state.items():
                    tr=st['def']
                    if tr.origin_type=='dispense' and (not st['onchip']) and tr.path[0]==src and st['available'] and not st['consumed'] and not st['done']:
                        st['onchip']=True; st['pos']=src
                        break
            moved_sources={(a,b):(c,d) for a,b,c,d in move_ops}
            for tid,st in token_state.items():
                if st['onchip'] and st['pos'] in moved_sources:
                    old=st['pos']
                    st['pos']=moved_sources[old]
            for pos in waste_ops:
                for tid,st in token_state.items():
                    if st['onchip'] and st['pos']==pos and st['def'].terminal[0]=='waste':
                        st['onchip']=False; st['done']=True
                        break
            for pos in out_ops:
                for tid,st in token_state.items():
                    if st['onchip'] and st['pos']==pos and st['def'].terminal[0]=='output':
                        st['onchip']=False; st['done']=True
                        break
            t += 1
            stall_iters = 0

            progress = progress_signature()
            if best_progress is None or progress > best_progress:
                best_progress = progress
                no_progress_iters = 0
            else:
                no_progress_iters += 1
                if no_progress_iters > LIVELOCK_LIMIT:
                    # Tokens are moving every tick (stall_iters never
                    # trips) but nothing is net progressing: a livelock,
                    # not a stall. Treat it the same as a real deadlock -
                    # bail rather than grind on and silently drop the
                    # remaining waste/output ops (see comment above
                    # LIVELOCK_LIMIT's definition).
                    deadlocked = True
                    break
        else:
            unfinished = [mid for mid,ms in mix_state.items() if not ms['done']]
            if not unfinished and all(st['done'] or st['consumed'] for st in token_state.values()):
                break
            next_ts=[ms['avail_t'] for ms in mix_state.values() if ms['started'] and not ms['done']]
            if next_ts:
                t = max(t+1, min(next_ts))
            else:
                t += 1
            stall_iters += 1
            if stall_iters > STALL_LIMIT:
                # Real deadlock (not just "waiting on a mixer to finish"):
                # STALL_LIMIT consecutive loop iterations produced no
                # mix/move/dispense/waste/output at all. Bail out instead
                # of silently grinding to MAXT and emitting a truncated
                # schedule that will hang the GUI simulator.
                deadlocked = True
                break

    if deadlocked:
        # Do not write a broken/truncated schedule. Signal failure so the
        # caller can fall back to the verified-safe baseline schedule.
        return False

    with open(out_dmfb,'w') as f:
        for ln in header:
            f.write(ln+'\n')
        for ln in lines:
            f.write(ln+'\n')
        f.write(f'{t} end\n')
    return True

if __name__=='__main__':
    ok = compact(sys.argv[1], sys.argv[2])
    sys.exit(0 if ok else 1)