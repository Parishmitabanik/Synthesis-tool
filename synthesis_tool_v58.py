#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DMFB Synthesis Tool - Graph/DAG aware - SimBioSys compatible (rev3)
- No dummy move(r,c,r,c)
- True timestamps
- mix_split time = seconds (multiple of 6) to satisfy Biochip_V
- idle gap = mix_seconds
"""
from __future__ import print_function
import sys
from collections import defaultdict, deque

CYCLE_SEC = 6

def seconds_to_cycles(s):
    return (s + CYCLE_SEC -1)//CYCLE_SEC or 1

def bfs(src,dst,rows,cols,blocked):
    from collections import deque
    if dst in blocked: return []
    if src==dst: return [src]
    q=deque([(src,[src])]); vis={src}
    while q:
        (r,c),path=q.popleft()
        for nr,nc in [(r-1,c),(r+1,c),(r,c-1),(r,c+1)]:
            nxt=(nr,nc)
            if nxt==dst: return path+[nxt]
            if 1<=nr<=rows and 1<=nc<=cols and nxt not in vis and nxt not in blocked:
                vis.add(nxt); q.append((nxt,path+[nxt]))
    return []

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
    from collections import deque
    d=dict(indeg); q=deque(sorted([n for n in nodes if d.get(n,0)==0])); order=[]
    while q:
        n=q.popleft(); order.append(n)
        for m in sorted(adj[n]):
            d[m]-=1
            if d[m]==0: q.append(m)
    return order, len(order)==len(nodes)

def main(arch,assay,out):
    rows,cols,res,waste,outp=parse_arch(arch)
    nodes,edges=parse_assay(assay)
    adj,indeg,inc=build_graph(nodes,edges)
    order,ok=topo(nodes,adj,indeg)
    if not ok: print("cycle!"); return
    # map reservoirs
    res_coords={n:res[n] for n in order if nodes[n]['type']=='dispense' and n in res}
    out_coords={}
    for n in order:
        if nodes[n]['type']!='output': continue
        if n in outp: out_coords[n]=outp[n]
        else:
            # take first free
            free=[p for p in outp if p not in out_coords.values()]
            # Actually outp keys are port IDs, simplify: use first output port
            if outp: out_coords[n]=list(outp.values())[0]
    # header
    reagent_line=' '.join('reagent({},{},{})'.format(r,c,nodes[nid]['reagent'] or nid) for nid,(r,c) in res_coords.items())
    header_extra=' '.join([' '.join('waste_reservoir({},{})'.format(r,c) for r,c in waste.values()),
                          ' '.join('output_reservoir({},{})'.format(r,c) for r,c in outp.values())]).strip()
    # mix ops
    mix_ops=[]
    for nid in order:
        if nodes[nid]['type']=='mix':
            ins=inc[nid]
            if len(ins)==2:
                sec=nodes[nid]['seconds'] or 6
                mix_ops.append((ins[0],ins[1],nid,seconds_to_cycles(sec),sec))
    # scheduler
    cmds=defaultdict(list); droplet_pos=defaultdict(list); waste_drops=[]; t=1
    res_lookup={nid:coord for nid,coord in res_coords.items()}
    # add reagent name lookup
    for nid,(r,c) in res_coords.items():
        rn=nodes[nid]['reagent']
        if rn: res_lookup[rn]=(r,c)
    # no-park
    no_park=set()
    for pr,pc in list(res.values())+list(waste.values())+list(outp.values()):
        for dr in (-1,0,1):
            for dc in (-1,0,1):
                no_park.add((pr+dr,pc+dc))
    for src_a,src_b,dest,mix_cycles,mix_sec in mix_ops:
        # dispense
        to_disp=[s for s in (src_a,src_b) if not droplet_pos[s] and s in res_lookup]
        if to_disp:
            for s in to_disp:
                rr,rc=res_lookup[s]
                cmds[t].append('dispense({},{})'.format(rr,rc))
                droplet_pos[s].append((rr,rc))
            t+=1
        if not droplet_pos[src_a] or not droplet_pos[src_b]:
            print('missing',src_a,src_b); continue
        pos_a=droplet_pos[src_a].pop(0); pos_b=droplet_pos[src_b].pop(0)
        # find zone
        zones=[]
        for zc in range(2,cols):
            for zr in range(2,rows-3):
                a=(zr,zc); b=(zr+3,zc)
                if a in no_park or b in no_park: continue
                zones.append((zr,zc))
        zones.sort(key=lambda z: abs(z[0]-pos_a[0])+abs(z[1]-pos_a[1])+abs(z[0]+3-pos_b[0])+abs(z[1]-pos_b[1]))
        routed=False
        for zr,zc in zones:
            ta=(zr,zc); tb=(zr+3,zc)
            others=[] 
            for v in droplet_pos.values(): others.extend(v)
            others.extend(waste_drops)
            static=set()
            for or_,oc in others:
                if (or_,oc) in res.values(): continue
                for dr in (-1,0,1):
                    for dc in (-1,0,1):
                        static.add((or_+dr,oc+dc))
            def buf(cells):
                s=set()
                for rr,cc in cells:
                    for dr in (-1,0,1):
                        for dc in (-1,0,1):
                            s.add((rr+dr,cc+dc))
                return s
            pa=bfs(pos_a,ta,rows,cols,static|buf([pos_b]))
            pb=bfs(pos_b,tb,rows,cols,static|buf(pa))
            if not pa or not pb:
                pa=bfs(pos_a,tb,rows,cols,static|buf([pos_b]))
                pb=bfs(pos_b,ta,rows,cols,static|buf(pa))
                if pa and pb: ta,tb=tb,ta
            if pa and pb:
                # fluidic check
                la,lb=len(pa),len(pb); ml=max(la,lb)
                fa=pa+[ta]*(ml-la); fb=pb+[tb]*(ml-lb)
                def close(p,q): return abs(p[0]-q[0])<=1 and abs(p[1]-q[1])<=1
                ok=True
                for i in range(ml):
                    if close(fa[i],fb[i]): ok=False; break
                    if i+1<ml and (close(fa[i+1],fb[i]) or close(fa[i],fb[i+1])): ok=False; break
                if ok:
                    routed=True; path_a=pa; path_b=pb; final_a=ta; final_b=tb; break
        if not routed: print('routing fail',src_a,src_b); droplet_pos[src_a].insert(0,pos_a); droplet_pos[src_b].insert(0,pos_b); continue
        # emit moves
        steps=max(len(path_a),len(path_b))-1
        for step in range(steps):
            cmds_step=[]
            if step+1 < len(path_a):
                r1,c1=path_a[step]; r2,c2=path_a[step+1]; cmds_step.append('move({},{},{},{})'.format(r1,c1,r2,c2))
            if step+1 < len(path_b):
                r1,c1=path_b[step]; r2,c2=path_b[step+1]; cmds_step.append('move({},{},{},{})'.format(r1,c1,r2,c2))
            if cmds_step:
                cmds[t].extend(cmds_step); t+=1
        r1,c1=final_a; r2,c2=final_b
        mix_time_arg=mix_sec  # seconds, multiple of 6 for Biochip_V
        cmds[t].append('mix_split({},{},{},{},{})'.format(r1,c1,r2,c2,mix_time_arg))
        t+=1
        t+=mix_time_arg  # idle
        droplet_pos[dest].append((r1,c1))
        # waste second droplet
        post2=(r2,c2)
        # try route waste
        if waste:
            # find nearest waste
            others=[]
            for v in droplet_pos.values(): others.extend(v)
            others.extend(waste_drops)
            static=set()
            for or_,oc in others:
                if (or_,oc) in res.values(): continue
                for dr in (-1,0,1):
                    for dc in (-1,0,1):
                        static.add((or_+dr,oc+dc))
            best=None
            for wp in waste.values():
                cand=bfs(post2,wp,rows,cols,static)
                if cand and (best is None or len(cand)<len(best)): best=cand
            if best:
                for i in range(len(best)-1):
                    a=best[i]; b=best[i+1]
                    cmds[t].append('move({},{},{},{})'.format(a[0],a[1],b[0],b[1])); t+=1
                cmds[t].append('waste({},{})'.format(best[-1][0],best[-1][1])); t+=1
            else:
                waste_drops.append(post2)
        else:
            waste_drops.append(post2)
    # output routing
    for nid in order:
        if nodes[nid]['type']!='output': continue
        srcs=inc[nid]
        if not srcs: continue
        src=srcs[0]
        if not droplet_pos[src]: continue
        start=droplet_pos[src].pop(0)
        if nid not in out_coords: continue
        target=out_coords[nid]
        # bfs
        others=[]
        for v in droplet_pos.values(): others.extend(v)
        others.extend(waste_drops)
        static=set()
        for or_,oc in others:
            if (or_,oc) in res.values(): continue
            for dr in (-1,0,1):
                for dc in (-1,0,1):
                    static.add((or_+dr,oc+dc))
        path=bfs(start,target,rows,cols,static)
        if path:
            for i in range(len(path)-1):
                a=path[i]; b=path[i+1]
                cmds[t].append('move({},{},{},{})'.format(a[0],a[1],b[0],b[1])); t+=1
            cmds[t].append('output({},{})'.format(path[-1][0],path[-1][1])); t+=1
    # write
    with open(out,'w') as f:
        f.write('dimension {} {}\n'.format(rows,cols))
        f.write('accuracy 5\n')
        hdr=reagent_line + (' '+header_extra if header_extra else '')
        f.write(hdr.strip()+'\n\n')
        for ts in sorted(cmds):
            if cmds[ts]:
                f.write('{} {}\n'.format(ts,' '.join(cmds[ts])))
        f.write('{} end\n'.format(t+1))
    print('Wrote',out,'last',t+1)

if __name__=='__main__':
    if len(sys.argv)<3:
        print('usage: python synthesis_tool_v3.py arch assay [out]')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv)>3 else 'output.dmfb')
