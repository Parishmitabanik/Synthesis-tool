#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel-safe DMFB synthesis tool.

Pipeline:
1. Generate a verifier-safe baseline schedule using dmfb_synthesis_tool_reduced_cycles.py
2. Compact that schedule with same-tick parallel moves using compact_from_safe.py
   under the exact Biochip_V semantics.

This version is intended for the user's SimBioSys/Biochip_V environment.
"""
from __future__ import print_function
import os
import re
import shutil
import sys

from dmfb_synthesis_tool_reduced_cycles import synthesize_graph as synthesize_safe
from compact_from_safe import compact as compact_schedule
from verify_biochipv import verify_dmfb
from shadow_gui import check_file as shadow_gui_check


def _passes_all_checks(path, label):
    """verify_dmfb only checks that the schedule is legal under
    BiochipV's simplified fluidic-constraint model. That model treats an
    active mixer as occupying just its two endpoint cells, but the real
    Tkinter renderer (Biochip.py) physically walks a droplet through the
    whole 1x4 / 4x1 zone between those endpoints for the entire mixing
    duration, and its move_droplet()/dispense_droplet() have no
    occupancy checks of their own - they just trust the schedule. A
    schedule that is "legal" per verify_dmfb can still park or route an
    unrelated droplet through that zone, which desyncs the two engines:
    the GUI either raises (missing droplet at a move source) or leaves
    an orphaned duplicate oval behind, and in both cases the droplet
    that should have kept moving just stops. That is the stuck-droplet
    freeze this project has hit before. shadow_gui.py replays the exact
    sequential bookkeeping Biochip.py performs and catches this before
    the file is ever handed to the real simulator.
    """
    try:
        verify_dmfb(path)
    except AssertionError as e:
        print('[FAIL] {} failed verify_dmfb: {}'.format(label, e))
        return False
    if not shadow_gui_check(path):
        print('[FAIL] {} would desync the Biochip.py GUI (see SHADOW-GUI FAIL above).'.format(label))
        return False
    return True


_OP_RE = re.compile(r'[a-z_]+\s*\([\s*\d+\s*,\s*]+\s*\w+\s*\)')


def _count_terminal_ops(path):
    """Count waste(...) and output(...) instructions in a .dmfb file.

    This is the completeness check that _looks_complete's timestamp-gap
    heuristic cannot provide. A compactor can grind out instructions
    right up until its tick budget runs out - e.g. two droplets stuck
    livelocked, shuttling back and forth every tick near a chokepoint -
    so the gap between the last emitted line and 'end' stays small even
    though the assay never actually finished. Every individual emitted
    move can be perfectly legal (passes verify_dmfb and shadow_gui) while
    the schedule as a whole silently drops its final waste/output events.
    Comparing terminal-op counts against the verified baseline catches
    that: a parallel-compacted schedule may reorder/interleave freely,
    but it must still perform the exact same set of terminal events as
    the baseline it was compacted from.
    """
    waste = 0
    output = 0
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.endswith(' end'):
                continue
            for instr in _OP_RE.findall(line):
                opcode = re.match(r'[a-z_]+', instr).group()
                if opcode == 'waste':
                    waste += 1
                elif opcode == 'output':
                    output += 1
    return waste, output


def _looks_complete(path):
    """Sanity check beyond verify_dmfb: the schedule must not silently
    stop early. verify_dmfb only checks that each executed instruction is
    legal - it happily accepts a file that stops after 10% of the assay
    and then jumps straight to 'N end'. That's exactly the bug that made
    graph6/graph1/graph2a hang the GUI at a low clock value: the compactor
    deadlocked, gave up, and stamped '3000 end' after the last real op.
    We detect that shape here: a huge gap between the last real
    instruction's timestamp and the terminal 'end' timestamp means the
    file was abandoned, not finished.
    """
    last_op_t = 0
    end_t = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.endswith(' end'):
                end_t = int(s.split()[0])
                continue
            tok = s.split()
            if tok and tok[0].isdigit():
                last_op_t = int(tok[0])
    if end_t is None:
        return False
    # A real schedule's 'end' marker sits right after the last operation
    # (verify_dmfb requires end_t >= last op time; in practice it's within
    # a handful of ticks). A gap of hundreds/thousands of ticks means the
    # compactor bailed out and padded the rest with nothing.
    return (end_t - last_op_t) < 50


def synthesize_graph(arch_file, assay_file, output_file):
    stage1 = output_file + '.stage1_safe.dmfb'
    ok = synthesize_safe(arch_file, assay_file, stage1)
    if not ok:
        return False

    # verify the baseline itself before trusting it as a fallback target
    if not _passes_all_checks(stage1, 'Baseline safe schedule'):
        return False

    compacted_ok = False
    try:
        compacted_ok = compact_schedule(stage1, output_file)
    except Exception as e:
        print('[WARN] Compaction raised an error ({}); falling back to baseline.'.format(e))
        compacted_ok = False

    if compacted_ok:
        # Belt-and-suspenders: re-verify the compacted file end-to-end,
        # make sure it isn't a truncated/deadlocked schedule padded out to
        # a huge 'end' timestamp (see _looks_complete for why that
        # matters), AND make sure it won't desync the real GUI engine.
        if not _passes_all_checks(output_file, 'Compacted schedule'):
            compacted_ok = False
        elif not _looks_complete(output_file):
            print('[WARN] Compacted schedule looks truncated (deadlocked mid-assay); falling back to baseline.')
            compacted_ok = False
        else:
            baseline_terms = _count_terminal_ops(stage1)
            compacted_terms = _count_terminal_ops(output_file)
            if compacted_terms != baseline_terms:
                print('[WARN] Compacted schedule is missing terminal operations '
                      '(baseline waste/output={} vs compacted={}); the compactor '
                      'likely livelocked and gave up before finishing the assay. '
                      'Falling back to baseline.'.format(baseline_terms, compacted_terms))
                compacted_ok = False

    if not compacted_ok:
        shutil.copyfile(stage1, output_file)
        print('[SUCCESS] Verified safe (sequential) schedule compiled to: {}'.format(output_file))
        print('[INFO] Parallel compaction was not possible for this assay/architecture combination,')
        print('[INFO] so the verified-safe baseline schedule was used as-is instead of a broken one.')
        return True

    print('[SUCCESS] Parallel-safe compiled to: {}'.format(output_file))
    print('[INFO] Baseline safe schedule kept at: {}'.format(stage1))
    return True


if __name__ == '__main__':
    arch = sys.argv[1] if len(sys.argv) > 1 else 'arch6.txt'
    assay = sys.argv[2] if len(sys.argv) > 2 else 'graph6.txt'
    out = sys.argv[3] if len(sys.argv) > 3 else 'graph6_parallel_safe.dmfb'
    ok = synthesize_graph(arch, assay, out)
    sys.exit(0 if ok else 1)