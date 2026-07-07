#!/usr/bin/env python3
"""
Shadow GUI: reproduces the *exact* sequential bookkeeping that Biochip.py
(the real Tkinter renderer) performs, so we can catch the "GUI desyncs
from the verifier and droplets get stuck / crash" class of bug without
needing Tkinter or Python 2.

Key fact this encodes, straight from Biochip.py:

  dispense_droplet(x,y):   creates ONE oval at (x,y). No occupancy check.
  move_droplet(a,b,c,d):   looks up the oval CURRENTLY at (a,b) via
                            electrode.find_all() and takes its FIRST id.
                            If there is no oval there -> IndexError -> the
                            worker thread dies -> the whole GUI freezes
                            (this is the "stuck droplet" symptom).
                            It does NOT delete/replace an oval that may
                            already exist at the destination (c,d); it
                            just draws a new one, so a stray leftover oval
                            silently piles up if the destination was not
                            actually empty.
  instantiate_mixer / mixer expiry: doesn't move ovals itself here; we
                            only need dispense/move/waste/output for this
                            check since those are the only ops that touch
                            oval placement directly in Biochip.py.

new_driver.py's gen() executes each regex-matched instruction in a line
ONE AT A TIME, strictly in the order they appear in the text of that
line (NOT re-grouped into mix/move/dispense/waste/output order the way
verify_biochipv.group_instructions + verify_line does internally). So if
our synthesis tool ever emits ops on one line in an order that assumes
atomic/simultaneous execution, but Biochip.py executes them one at a
time in that same textual order, this script will catch the mismatch.
"""
from __future__ import annotations
import re


def parse_ops(line):
    """Return the list of (opcode, operands) IN TEXTUAL ORDER, exactly as
    new_driver.gen() would iterate them."""
    out = []
    for instr in re.compile(r'[a-z_]+\s*\([\s*\d+\s*,\s*]+\s*\w+\s*\)').findall(line):
        obj = re.compile(r'[a-z_]+').match(instr)
        opcode = obj.group()
        operands = re.compile(r'(\d+|\w+)').findall(instr[obj.end() + 1:])
        out.append((opcode, operands))
    return out


class ShadowGUIError(Exception):
    pass


class ShadowGUI:
    def __init__(self):
        self.oval_at = {}  # (r,c) -> True if an oval currently sits there

    def dispense(self, r, c, line):
        # Biochip.py never checks occupancy before drawing a new oval.
        if self.oval_at.get((r, c)):
            raise ShadowGUIError(
                "line {!r}: dispense({},{}) drawn on top of an existing "
                "un-cleared oval -> orphaned/duplicate droplet".format(line, r, c))
        self.oval_at[(r, c)] = True

    def move(self, r1, c1, r2, c2, line):
        if not self.oval_at.get((r1, c1)):
            raise ShadowGUIError(
                "line {!r}: move({},{},{},{}) has NO oval at source "
                "-> find_all() returns [] -> IndexError -> worker thread "
                "dies -> GUI freezes (this is the stuck-droplet bug)"
                .format(line, r1, c1, r2, c2))
        del self.oval_at[(r1, c1)]
        if self.oval_at.get((r2, c2)):
            raise ShadowGUIError(
                "line {!r}: move({},{},{},{}) lands on a cell that still "
                "has an un-cleared oval -> stray duplicate droplet left "
                "on the board".format(line, r1, c1, r2, c2))
        self.oval_at[(r2, c2)] = True

    def clear(self, r, c, line, kind):
        if not self.oval_at.get((r, c)):
            raise ShadowGUIError(
                "line {!r}: {}({},{}) has no droplet to remove"
                .format(line, kind, r, c))
        del self.oval_at[(r, c)]

    def run_line(self, line):
        for opcode, ops in parse_ops(line):
            nums = [int(x) for x in ops if x.isdigit()] if opcode != 'dispense' else None
            if opcode == 'dispense':
                r, c = int(ops[0]), int(ops[1])
                self.dispense(r, c, line)
            elif opcode == 'move':
                r1, c1, r2, c2 = map(int, ops)
                self.move(r1, c1, r2, c2, line)
            elif opcode == 'waste':
                r, c = int(ops[0]), int(ops[1])
                self.clear(r, c, line, 'waste')
            elif opcode == 'output':
                r, c = int(ops[0]), int(ops[1])
                self.clear(r, c, line, 'output')
            elif opcode == 'mix_split':
                # mix_split doesn't touch oval placement in Biochip.py at
                # the moment it's issued (only later, on mixer expiry).
                pass


def check_file(path):
    gui = ShadowGUI()
    with open(path) as f:
        lines = f.readlines()
    # skip header (dimension/accuracy/reagent line + blank)
    body = [ln for ln in lines[3:] if ln.strip()]
    for ln in body:
        if ln.strip().endswith(' end'):
            break
        try:
            gui.run_line(ln)
        except ShadowGUIError as e:
            print('[SHADOW-GUI FAIL]', e)
            return False
    print('[SHADOW-GUI OK] no desync detected across {} instruction lines'.format(len(body)))
    return True


if __name__ == '__main__':
    import sys
    ok = check_file(sys.argv[1])
    sys.exit(0 if ok else 1)