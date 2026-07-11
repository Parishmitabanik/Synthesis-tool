#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Greedy nearest-available-site placement
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
import tempfile

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


def _end_tick(path):
    """Return the tick number on the 'N end' line of a .dmfb file."""
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.endswith(' end'):
                return int(s.split()[0])
    return None


def synthesize_graph(arch_file, assay_file, output_file):
    stage1 = output_file + '.stage1_safe.dmfb'
    ok = synthesize_safe(arch_file, assay_file, stage1)
    if not ok:
        return False

    # verify the baseline itself before trusting it as a fallback target
    if not _passes_all_checks(stage1, 'Baseline safe schedule'):
        return False

    # compact_from_safe.compact() is a greedy, history-blind router: how
    # much parallelism it actually finds is very sensitive to two knobs -
    # DISPENSE_FRONTIER (how many not-yet-started mixes are allowed to
    # pull their input reagents onto the chip at once) and MAX_ONCHIP (how
    # many droplets are allowed on the chip concurrently). Too low and the
    # router serializes independent branches of the assay for no reason;
    # too high and it overcrowds the grid, droplets fight over the same
    # lanes, and the router stalls/livelocks (caught by _passes_all_checks
    # / _looks_complete / the terminal-op count below) or simply produces
    # a *worse* schedule than a more modest setting. There's no closed-form
    # best setting - it depends on the architecture's grid size and the
    # assay's DAG shape - so instead of hand-picking one value we try a
    # short list of candidates and keep whichever fully-verified result
    # finishes in the fewest cycles. (1, 8) is the original conservative
    # setting and is always included, so this can never do worse than the
    # previous behaviour.
    CANDIDATE_CONFIGS = [
        (1, 8), (2, 8), (2, 9), (2, 10), (3, 8), (3, 10),
        (4, 8), (4, 10), (5, 8), (5, 10), (6, 8), (6, 10),
    ]

    baseline_terms = _count_terminal_ops(stage1)
    best_path = None
    best_end_t = None
    best_config = None

    tmpdir = tempfile.mkdtemp(prefix='dmfb_compact_')
    try:
        for dispense_frontier, max_onchip in CANDIDATE_CONFIGS:
            cand_file = os.path.join(tmpdir, 'cand_{}_{}.dmfb'.format(dispense_frontier, max_onchip))
            try:
                ok = compact_schedule(stage1, cand_file, dispense_frontier=dispense_frontier, max_onchip=max_onchip)
            except Exception as e:
                print('[WARN] Compaction (frontier={}, onchip={}) raised an error ({}); skipping.'.format(
                    dispense_frontier, max_onchip, e))
                ok = False

            if not ok:
                continue
            if not _passes_all_checks(cand_file, 'Compacted schedule (frontier={}, onchip={})'.format(
                    dispense_frontier, max_onchip)):
                continue
            if not _looks_complete(cand_file):
                print('[WARN] Compacted schedule (frontier={}, onchip={}) looks truncated; skipping.'.format(
                    dispense_frontier, max_onchip))
                continue
            compacted_terms = _count_terminal_ops(cand_file)
            if compacted_terms != baseline_terms:
                print('[WARN] Compacted schedule (frontier={}, onchip={}) is missing terminal operations '
                      '(baseline waste/output={} vs compacted={}); skipping.'.format(
                          dispense_frontier, max_onchip, baseline_terms, compacted_terms))
                continue

            end_t = _end_tick(cand_file)
            print('[INFO] Candidate (frontier={}, onchip={}) verified OK, finishes at tick {}.'.format(
                dispense_frontier, max_onchip, end_t))
            if best_end_t is None or end_t < best_end_t:
                best_end_t = end_t
                best_config = (dispense_frontier, max_onchip)
                # Copy out immediately since cand_file itself lives in the
                # temp dir we're about to delete.
                best_path = output_file + '.best_candidate.dmfb'
                shutil.copyfile(cand_file, best_path)

        if best_path is None:
            shutil.copyfile(stage1, output_file)
            print('[SUCCESS] Verified safe (sequential) schedule compiled to: {}'.format(output_file))
            print('[INFO] Parallel compaction was not possible for this assay/architecture combination,')
            print('[INFO] so the verified-safe baseline schedule was used as-is instead of a broken one.')
            return True

        shutil.move(best_path, output_file)
        print('[SUCCESS] Parallel-safe compiled to: {} (frontier={}, onchip={}, finishes at tick {})'.format(
            output_file, best_config[0], best_config[1], best_end_t))
        print('[INFO] Baseline safe schedule kept at: {}'.format(stage1))
        return True
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == '__main__':
    arch = sys.argv[1] if len(sys.argv) > 1 else 'arch6.txt'
    assay = sys.argv[2] if len(sys.argv) > 2 else 'graph6.txt'
    out = sys.argv[3] if len(sys.argv) > 3 else 'graph6_parallel_safe.dmfb'
    ok = synthesize_graph(arch, assay, out)
    sys.exit(0 if ok else 1)