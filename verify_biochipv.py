#!/usr/bin/env python3
from __future__ import annotations
import re, sys
from collections import defaultdict

class Droplet:
    def __init__(self):
        self.id = None
        self.conc = []

class Mixer:
    def __init__(self, mixer_type, start, end, r1, c1, r2, c2):
        self.mtype = mixer_type
        self.s_time = start
        self.e_time = end
        self.d1r = r1
        self.d1c = c1
        self.d2r = r2
        self.d2c = c2

class BiochipV:
    def __init__(self, r, c, n):
        self.row = r
        self.col = c
        self.accuracy = n
        self.op_reservior = []
        self.ip_reservior = []
        self.waste_reservior = []
        self.c_time = 0
        self.biochip = [[Droplet() for _ in range(c+2)] for _ in range(r+2)]
        self.mixer_table = []
        self._ID = r + c
        self._RID = -1

    def _get_new_ID(self):
        self._ID += 1
        return self._ID
    def _get_new_RID(self):
        self._RID += 1
        return self._RID
    def _valid(self, x, y):
        return x in range(1, self.row+1) and y in range(1, self.col+1)
    def _FC_between_2(self, x1,y1,x2,y2):
        return abs(x1-x2) >= 2 or abs(y1-y2) >= 2

    def verify_set_reagent_reservior(self, reagent_dispenser_list, name):
        RID = self._get_new_RID()
        for (r,c) in reagent_dispenser_list:
            assert self._valid(r,c)
            assert (r in [1,self.row]) or (c in [1,self.col])
            for reservior in self.op_reservior + self.waste_reservior + self.ip_reservior:
                assert self._FC_between_2(reservior[0],reservior[1],r,c)
            self.ip_reservior.append((r,c,'reagent',name,RID))

    def verify_set_resevior(self,r,c,reservior_type):
        assert self._valid(r,c)
        assert (r in [1,self.row]) or (c in [1,self.col])
        for reservior in self.op_reservior + self.waste_reservior + self.ip_reservior:
            assert self._FC_between_2(reservior[0],reservior[1],r,c)
        if reservior_type == 'waste_reservoir':
            ID = self._get_new_ID()
            self.waste_reservior.append((r,c,'waste',ID))
        elif reservior_type == 'op_reservoir':
            ID = self._get_new_ID()
            self.op_reservior.append((r,c,'output',ID))
        else:
            raise AssertionError('bad reservoir type')

    def droplet_in_active_mixers(self,r,c):
        for mixer in self.mixer_table:
            if (mixer.d1r,mixer.d1c)==(r,c) or (mixer.d2r,mixer.d2c)==(r,c):
                return True
        return False

    def check_FC(self,r,c):
        assert self._valid(r,c)
        for x in range(1,self.row+1):
            for y in range(1,self.col+1):
                if self.biochip[x][y].id is None:
                    continue
                elif (x,y)!=(r,c):
                    if not self._FC_between_2(x,y,r,c):
                        return False
        return True

    def verify_dispense_insrt(self, dispense_instr, line):
        for instr in dispense_instr:
            flag=0
            for reagent_dispr in self.ip_reservior:
                if instr == (reagent_dispr[0],reagent_dispr[1]):
                    flag=1
                    RID=reagent_dispr[4]
                    break
            assert flag==1, ('disp invalid', line)
            (r,c)=instr
            assert self.biochip[r][c].id is None, ('disp occupied', line)
            reagent_conc=[0]*(self._RID+1)
            reagent_conc[RID]=2**self.accuracy
            self.biochip[r][c].id = RID
            self.biochip[r][c].conc = reagent_conc[:]
            assert self.check_FC(r,c), ('disp FC', line)

    def verify_mix_split_insrt(self, mix_split_instr, line):
        for instr in mix_split_instr:
            (d1r,d1c,d2r,d2c,mixing_time) = instr
            assert self._valid(d1r,d1c) and self._valid(d2r,d2c), ('mix invalid loc', line)
            assert not self.droplet_in_active_mixers(d1r,d1c) and not self.droplet_in_active_mixers(d2r,d2c), ('mix in active', line)
            assert self.biochip[d1r][d1c].id is not None and self.biochip[d2r][d2c].id is not None, ('mix absent', line)
            assert mixing_time % 6 == 0, ('mix time', line)
            if d1r == d2r and abs(d1c-d2c)==3:
                if d1c < d2c: m = Mixer(14,self.c_time,self.c_time+mixing_time,d1r,d1c,d2r,d2c)
                else: m = Mixer(14,self.c_time,self.c_time+mixing_time,d2r,d2c,d1r,d1c)
                self.mixer_table.append(m)
            elif d1c == d2c and abs(d1r-d2r)==3:
                if d1r < d2r: m = Mixer(41,self.c_time,self.c_time+mixing_time,d1r,d1c,d2r,d2c)
                else: m = Mixer(41,self.c_time,self.c_time+mixing_time,d2r,d2c,d1r,d1c)
                self.mixer_table.append(m)
            else:
                raise AssertionError(('mix cannot instantiate', line))

    def verify_move_instr(self, move_instr, line):
        for instr in move_instr:
            (d1r,d1c,d2r,d2c)=instr
            assert self._valid(d1r,d1c) and self._valid(d2r,d2c), ('move invalid loc', line)
            assert not self.droplet_in_active_mixers(d1r,d1c), ('move src in active', line)
            assert self.biochip[d1r][d1c].id is not None, ('move src absent', line)
            assert abs(d1r-d2r) <= 1 and abs(d1c-d2c) <= 1, ('move step invalid', line)
            assert self.biochip[d2r][d2c].id is None, ('move dst occupied', line)
            assert not self.droplet_in_active_mixers(d2r,d2c), ('move dst active', line)
            ID = self.biochip[d1r][d1c].id
            self.biochip[d1r][d1c].id = None
            self.biochip[d2r][d2c].id = ID
            assert self.check_FC(d2r,d2c), ('D-7', line)
            self.biochip[d2r][d2c].id = None
            self.biochip[d1r][d1c].id = ID
        for instr in move_instr:
            (d1r,d1c,d2r,d2c)=instr
            self.biochip[d2r][d2c].id = self.biochip[d1r][d1c].id
            self.biochip[d2r][d2c].conc = self.biochip[d1r][d1c].conc[:]
            self.biochip[d1r][d1c].id = None
            self.biochip[d1r][d1c].conc = []
            assert self.check_FC(d2r,d2c), ('SFC', line)

    def verify_waste_insrt(self, waste_instr, line):
        for instr in waste_instr:
            flag=0
            for waste_dispr in self.waste_reservior:
                if instr == (waste_dispr[0],waste_dispr[1]):
                    flag=1
            assert flag==1, ('waste invalid port', line)
            (r,c)=instr
            assert self.biochip[r][c].id is not None, ('waste no droplet', line)
            self.biochip[r][c].id = None
            self.biochip[r][c].conc = []

    def verify_op_insrt(self, op_instr, line):
        for instr in op_instr:
            flag=0
            for op_dispr in self.op_reservior:
                if instr == (op_dispr[0],op_dispr[1]):
                    flag=1
            assert flag==1, ('output invalid port', line)
            (r,c)=instr
            assert self.biochip[r][c].id is not None, ('output no droplet', line)
            self.biochip[r][c].id = None
            self.biochip[r][c].conc = []

    def group_instructions(self,line):
        dispense_instr=[]; mix_split_instr=[]; move_instr=[]; waste_instr=[]; op_instr=[]
        instr_line = re.compile(r'[a-z_]+\s*\([\s*\d+\s*,\s*]+\s*\w+\s*\)').findall(line)
        for instr in instr_line:
            obj = re.compile(r'[a-z_]+').match(instr)
            opcode = obj.group()
            operands = re.compile(r'(\d+|\w+)').findall(instr[obj.end()+1:])
            if opcode == 'dispense' and len(operands)==2:
                dispense_instr.append((int(operands[0]), int(operands[1])))
            elif opcode == 'mix_split' and len(operands)==5:
                mix_split_instr.append(tuple(map(int, operands)))
            elif opcode == 'move' and len(operands)==4:
                move_instr.append(tuple(map(int, operands)))
            elif opcode == 'waste' and len(operands)==2:
                waste_instr.append((int(operands[0]), int(operands[1])))
            elif opcode == 'output' and len(operands)==2:
                op_instr.append((int(operands[0]), int(operands[1])))
        return dispense_instr, mix_split_instr, move_instr, waste_instr, op_instr

    def conc_ratio_after_mixing(self, l, m):
        return [(l[i]+m[i])//2 for i in range(len(l))]

    def delete_expired_mixers(self, curr_time):
        active=[]
        for m in self.mixer_table:
            if curr_time > m.e_time:
                new_id = self._get_new_ID()
                new_conc = self.conc_ratio_after_mixing(self.biochip[m.d1r][m.d1c].conc, self.biochip[m.d2r][m.d2c].conc)
                self.biochip[m.d1r][m.d1c].id = new_id
                self.biochip[m.d1r][m.d1c].conc = new_conc[:]
                self.biochip[m.d2r][m.d2c].id = new_id
                self.biochip[m.d2r][m.d2c].conc = new_conc[:]
            else:
                active.append(m)
        self.mixer_table = active

    def verify_line(self, line):
        time = re.compile(r'\d+').match(line)
        assert time is not None, ('bad syntax', line)
        op_time = int(time.group())
        assert self.c_time <= op_time, ('timing violation', line)
        self.c_time = op_time
        self.delete_expired_mixers(self.c_time)
        dispense_instr, mix_split_instr, move_instr, waste_instr, op_instr = self.group_instructions(line)
        self.verify_mix_split_insrt(mix_split_instr, line)
        self.verify_move_instr(move_instr, line)
        self.verify_dispense_insrt(dispense_instr, line)
        self.verify_waste_insrt(waste_instr, line)
        self.verify_op_insrt(op_instr, line)
        self.c_time += 1


def verify_dmfb(path):
    with open(path) as f:
        line = f.readline()
        while line and not line.strip(): line = f.readline()
        tok=line.split(); assert tok[0]=='dimension' and len(tok)==3
        row,col = int(tok[1]), int(tok[2])
        line=f.readline();
        while line and not line.strip(): line=f.readline()
        tok=line.split(); assert tok[0]=='accuracy' and len(tok)==2
        acc=int(tok[1])
        B = BiochipV(row,col,acc)
        line=f.readline()
        while line and not line.strip(): line=f.readline()
        instr_line = re.compile(r'[a-z_]+\s*\([\s*\d+\s*,\s*]+\s*\w+\s*\)').findall(line)
        for instr in instr_line:
            obj = re.compile(r'[a-z_]+').match(instr)
            opcode = obj.group()
            operands = re.compile(r'(\d+|\w+)').findall(instr[obj.end()+1:])
            if opcode=='reagent' and len(operands)>=3 and len(operands)%2==1:
                lst=[]
                for i in range(0, len(operands)-1, 2):
                    lst.append((int(operands[i]), int(operands[i+1])))
                B.verify_set_reagent_reservior(lst, operands[-1])
            elif opcode=='waste_reservoir' and len(operands)==2:
                B.verify_set_resevior(int(operands[0]), int(operands[1]), 'waste_reservoir')
            elif opcode=='output_reservoir' and len(operands)==2:
                B.verify_set_resevior(int(operands[0]), int(operands[1]), 'op_reservoir')
        for line in f:
            if not line.strip():
                continue
            if line.strip().endswith(' end'):
                break
            B.verify_line(line)
    return True

if __name__ == '__main__':
    try:
        verify_dmfb(sys.argv[1])
        print('OK')
    except AssertionError as e:
        print('FAIL', e)
        raise