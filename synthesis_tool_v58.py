#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DMFB Synthesis Tool - Graph/DAG aware - SimBioSys compatible (rev4 - Parallelized)
- No dummy move(r,c,r,c)
- True timestamps
- mix_split time = seconds (multiple of 6) to satisfy Biochip_V
- Parallel execution using space-time reservation table
"""
from __future__ import print_function
import sys
from collections import defaultdict, deque

CYCLE_SEC = 6

def seconds_to_cycles(s):
    return (s + CYCLE_SEC -1)//CYCLE_SEC or 1

class Router:
    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self.res = {} # (t, r, c) -> uid

    def is_free(self, t, r, c, ignore_id1=None, ignore_id2=None):
        if not (1 <= r <= self.rows and 1 <= c <= self.cols): return False
        for dr in (-1,0,1):
            for dc in (-1,0,1):
                occ = self.res.get((t, r+dr, c+dc))
                if occ and occ not in (ignore_id1, ignore_id2):
                    return False
        return True

    def block_path(self, path, start_t, uid, hold_at_end=True, max_t=3000):
        t = start_t
        for i in range(len(path)):
            r, c = path[i]
            self.res[(t, r, c)] = uid
            if i < len(path)-1:
                t += 1
        if hold_at_end and path:
            r, c = path[-1]
            for tt in range(t, max_t):
                self.res[(tt, r, c)] = uid
        return t

    def free_from(self, start_t, r, c, uid, max_t=3000):
        for t in range(start_t, max_t):
            if self.res.get((t, r, c)) == uid:
                del self.res[(t, r, c)]

    def bfs_time(self, src, dst, start_t, uid1, uid2=None, no_park=None):
        no_park = no_park or set()
        if not self.is_free(start_t, dst[0], dst[1], uid1, uid2):
            pass
        q = deque([(src, start_t, [src])])
        vis = {(src, start_t)}
        while q:
            (r, c), t, path = q.popleft()
            if (r, c) == dst: return path, t
            if t > start_t + 150: continue
            nxt_t = t + 1
            moves = [(r,c), (r-1,c), (r+1,c), (r,c-1), (r,c+1)]
            moves.sort(key=lambda m: abs(m[0]-dst[0]) + abs(m[1]-dst[1]))
            for nr, nc in moves:
                if (nr, nc) in no_park and (nr, nc) != dst:
                    if (nr, nc) == (r, c): pass
                    else: continue
                if self.is_free(nxt_t, nr, nc, uid1, uid2):
                    if ((nr, nc), nxt_t) not in vis:
                        vis.add(((nr, nc), nxt_t))
                        q.append(((nr, nc), nxt_t, path + [(nr, nc)]))
        return [], start_t

def parse_arch(p):
    rows=cols=15; res={}; waste={}; out={} 
    with open(p) as f:
        for line in f:
            line=line.split('#',1)[0].strip()
            if not line: continue
            t=line.split(); k=t[0].upper()
            if k=='GRID': rows,cols=int(t[1]),int(t[2])
            elif k=='RESERVOIR': res[t[1]]=(int(t[2]),int(t[3]))
            elif k=='WASTE': waste[t[1]]=(int(t[2]),int(t[3]))
            elif k=='OUTPUT': out[t[1]]=(int(t[2]),int(t[3]))
    return rows,cols,res,waste,out

def parse_assay(p):
    nodes={}; edges=[]
    with open(p) as f:
        for line in f:
            line=line.split('#',1)[0].strip()
            if not line: continue
            t=line.split(); k=t[0].upper()
            if k=='NODE':
                nid=t[1]; ntype=t[2].lower(); reagent=None; sec=None
                if ntype=='dispense': reagent=t[3] if len(t)>3 else None
                elif ntype=='mix': sec=int(t[3]) if len(t)>3 else None
                nodes[nid]={'type':ntype,'reagent':reagent,'seconds':sec}
            elif k=='EDGE': edges.append((t[1],t[2]))
    return nodes,edges

def build_graph(nodes,edges):
    adj=defaultdict(list); inc=defaultdict(list); indeg={n:0 for n in nodes}
    for a,b in edges:
        adj[a].append(b); inc[b].append(a); indeg[b]=indeg.get(b,0)+1; indeg.setdefault(a,0)
    return adj,indeg,inc

def topo(nodes,adj,indeg):
    d=dict(indeg); q=deque(sorted([n for n in nodes if d.get(n,0)==0])); order=[]
    while q:
        n=q.popleft(); order.append(n)
        for m in sorted(adj[n]):
            d[m]-=1
            if d[m]==0: q.append(m)
    return order, len(order)==len(nodes)

def main(arch,assay,out_file):
    rows,cols,res,waste,outp=parse_arch(arch)
    nodes,edges=parse_assay(assay)
    adj,indeg,inc=build_graph(nodes,edges)
    order,ok=topo(nodes,adj,indeg)
    if not ok: print("cycle!"); return
    
    res_coords={n:res[n] for n in order if nodes[n]['type']=='dispense' and n in res}
    out_coords={}
    for n in order:
        if nodes[n]['type']!='output': continue
        if n in outp: out_coords[n]=outp[n]
        else:
            if outp: out_coords[n]=list(outp.values())[0]

    reagent_line=' '.join('reagent({},{},{})'.format(r,c,nodes[nid]['reagent'] or nid) for nid,(r,c) in res_coords.items())
    header_extra=' '.join([' '.join('waste_reservoir({},{})'.format(r,c) for r,c in waste.values()),
                          ' '.join('output_reservoir({},{})'.format(r,c) for r,c in outp.values())]).strip()

    mix_ops=[]
    for nid in order:
        if nodes[nid]['type']=='mix':
            ins=inc[nid]
            if len(ins)==2:
                sec=nodes[nid]['seconds'] or 6
                mix_ops.append((ins[0],ins[1],nid,seconds_to_cycles(sec),sec))

    cmds=defaultdict(list); droplet_pos={}; t_ready=defaultdict(lambda: 1)
    res_lookup={nid:coord for nid,coord in res_coords.items()}
    for nid,(r,c) in res_coords.items():
        rn=nodes[nid]['reagent']
        if rn: res_lookup[rn]=(r,c)

    no_park=set()
    for pr,pc in list(res.values())+list(waste.values())+list(outp.values()):
        for dr in (-1,0,1):
            for dc in (-1,0,1):
                no_park.add((pr+dr,pc+dc))

    router = Router(rows, cols)

    for src_a,src_b,dest,mix_cycles,mix_sec in mix_ops:
        to_disp=[s for s in (src_a,src_b) if s not in droplet_pos and s in res_lookup]
        if to_disp:
            for s in to_disp:
                rr,rc=res_lookup[s]
                cmds[t_ready[s]].append('dispense({},{})'.format(rr,rc))
                droplet_pos[s] = (rr,rc)
                router.block_path([(rr,rc)], t_ready[s], s, hold_at_end=True)
                t_ready[s] += 1
                
        if src_a not in droplet_pos or src_b not in droplet_pos:
            print('missing',src_a,src_b); continue
            
        pos_a=droplet_pos[src_a]; pos_b=droplet_pos[src_b]
        
        zones=[]
        for zc in range(2,cols):
            for zr in range(2,rows-3):
                a=(zr,zc); b=(zr+3,zc)
                if a in no_park or b in no_park: continue
                zones.append((zr,zc))
                
        start_t = max(t_ready[src_a], t_ready[src_b])
        routed = False
        
        def try_route(ta, tb, st):
            res_backup = dict(router.res)
            router.free_from(st, pos_a[0], pos_a[1], src_a)
            path_a, end_a = router.bfs_time(pos_a, ta, st, src_a)
            if not path_a:
                router.res = res_backup; return None
            router.block_path(path_a, st, src_a, hold_at_end=True)
            
            router.free_from(st, pos_b[0], pos_b[1], src_b)
            path_b, end_b = router.bfs_time(pos_b, tb, st, src_b)
            if not path_b:
                router.res = res_backup; return None
            router.block_path(path_b, st, src_b, hold_at_end=True)
            return path_a, end_a, path_b, end_b

        while not routed and start_t < 2500:
            valid_zones = []
            for zr, zc in zones:
                is_free = True
                for mt in range(start_t, start_t + 20 + mix_sec):
                    for mr in range(zr-1, zr+4):
                        for mc in range(zc-1, zc+2):
                            if router.res.get((mt, mr, mc)) not in (None, src_a, src_b):
                                is_free = False; break
                        if not is_free: break
                    if not is_free: break
                if is_free: valid_zones.append((zr, zc))
                
            valid_zones.sort(key=lambda z: abs(z[0]-pos_a[0])+abs(z[1]-pos_a[1])+abs(z[0]+3-pos_b[0])+abs(z[1]-pos_b[1]))
            
            for zr, zc in valid_zones:
                ta=(zr,zc); tb=(zr+3,zc)
                res_route = try_route(ta, tb, start_t)
                if not res_route:
                    res_route = try_route(tb, ta, start_t)
                    if res_route: ta, tb = tb, ta
                if res_route:
                    path_a, end_a, path_b, end_b = res_route
                    routed = True
                    break
            
            if not routed:
                start_t += 5
                
        if not routed:
            print('routing fail',src_a,src_b); continue

        mix_start_t = max(end_a, end_b)
        
        for i in range(len(path_a)-1):
            cmds[start_t + i].append('move({},{},{},{})'.format(path_a[i][0], path_a[i][1], path_a[i+1][0], path_a[i+1][1]))
        for i in range(len(path_b)-1):
            cmds[start_t + i].append('move({},{},{},{})'.format(path_b[i][0], path_b[i][1], path_b[i+1][0], path_b[i+1][1]))
            
        cmds[mix_start_t].append('mix_split({},{},{},{},{})'.format(ta[0],ta[1],tb[0],tb[1],mix_sec))
        mix_end_t = mix_start_t + 1 + mix_sec
        
        for mt in range(mix_start_t, mix_end_t):
            for mr in range(ta[0], tb[0]+1):
                router.res[(mt, mr, ta[1])] = dest
                
        router.free_from(mix_start_t, ta[0], ta[1], src_a)
        router.free_from(mix_start_t, tb[0], tb[1], src_b)
        
        router.block_path([(ta[0], ta[1])], mix_end_t, dest, hold_at_end=True)
        droplet_pos[dest] = (ta[0], ta[1])
        t_ready[dest] = mix_end_t
        
        post2=(tb[0], tb[1])
        waste_id = dest + "_waste"
        router.block_path([post2], mix_end_t, waste_id, hold_at_end=True)
        
        if waste:
            best_path = None
            best_end = float('inf')
            router.free_from(mix_end_t, post2[0], post2[1], waste_id)
            for wp in waste.values():
                p, end_p = router.bfs_time(post2, wp, mix_end_t, waste_id)
                if p and end_p < best_end:
                    best_path, best_end = p, end_p
            if best_path:
                router.block_path(best_path, mix_end_t, waste_id, hold_at_end=False)
                for i in range(len(best_path)-1):
                    cmds[mix_end_t + i].append('move({},{},{},{})'.format(best_path[i][0], best_path[i][1], best_path[i+1][0], best_path[i+1][1]))
                cmds[best_end].append('waste({},{})'.format(best_path[-1][0], best_path[-1][1]))
            else:
                router.block_path([post2], mix_end_t, waste_id, hold_at_end=True)
                
    for nid in order:
        if nodes[nid]['type']!='output': continue
        srcs=inc[nid]
        if not srcs: continue
        src=srcs[0]
        if src not in droplet_pos: continue
        
        start = droplet_pos[src]
        if nid not in out_coords: continue
        target = out_coords[nid]
        start_t = t_ready[src]
        
        router.free_from(start_t, start[0], start[1], src)
        path, end_t = router.bfs_time(start, target, start_t, src)
        if path:
            router.block_path(path, start_t, src, hold_at_end=False)
            for i in range(len(path)-1):
                cmds[start_t + i].append('move({},{},{},{})'.format(path[i][0], path[i][1], path[i+1][0], path[i+1][1]))
            cmds[end_t].append('output({},{})'.format(path[-1][0], path[-1][1]))

    with open(out_file,'w') as f:
        f.write('dimension {} {}\n'.format(rows,cols))
        f.write('accuracy 5\n')
        hdr=reagent_line + (' '+header_extra if header_extra else '')
        f.write(hdr.strip()+'\n\n')
        max_t = 0
        for ts in sorted(cmds):
            if cmds[ts]:
                f.write('{} {}\n'.format(ts,' '.join(cmds[ts])))
                max_t = max(max_t, ts)
        f.write('{} end\n'.format(max_t+1))
    print('Wrote',out_file,'last',max_t+1)

if __name__=='__main__':
    if len(sys.argv)<3:
        print('usage: python synthesis_tool_v58.py arch assay [out]')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv)>3 else 'output.dmfb')
