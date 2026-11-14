"""Verify the compact exact rational H6 sixth-jet certificate.

Default: full CF fourth-jet checks, global certificate weights, and selected
exact sixth-order hook/H2/raising/closure checks. --full expands every raw
coefficient and verifies the complete quintic equation. Optional python-flint
speeds up full verification; --fraction uses only standard Fraction arithmetic.
Finite Taylor data does not prove existence of an actual closed field.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import json
import time

import CF_h6_fourth_jet as v
import CM_h6_sixth_codec as codec
import CM_h6_sixth_geometry as g
import CM_h6_sixth_h2 as h2

HERE=Path(__file__).resolve().parent;tick=time.time()

def memory_profile():
    try:
        import ctypes
        class Counters(ctypes.Structure):
            _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong)]+[(name,ctypes.c_size_t) for name in ['peak_working_set','working_set','peak_paged','paged','peak_nonpaged','nonpaged','pagefile','peak_pagefile']]
        c=Counters();c.cb=ctypes.sizeof(c)
        okay=ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.c_void_p(-1),ctypes.byref(c),c.cb)
        return {'peak_working_set_bytes':int(c.peak_working_set)} if okay else {}
    except (AttributeError,OSError):return {}

def selected(certificate):
    aux=certificate.aux;labels=[aux['closure_labels'][i] for i in aux['default_closure_indices']];highest=[tuple(m) for m in aux['homogeneous_monomials']];needed=set(highest)
    for mon,a,m in labels:
        for i in v.BO.axes(m):needed.add(g.union(tuple(mon),(i,)))
    for mon in highest:
        for action in v.raising:
            for i,row in enumerate(action):
                for j,c in row:
                    if j in mon:needed.add(g.union(g.subtract(mon,(j,)),(i,)))
    poly=certificate.expand(Fraction,needed);v.constraints(poly,6);h2.verify(poly)
    for mon in highest:
        for action in v.raising:
            out=v.BO.root_action(poly.get(mon,{}),action)
            for i,row in enumerate(action):
                for j,c in row:
                    if j not in mon:continue
                    source=g.union(g.subtract(mon,(j,)),(i,));factor=source.count(i)*c
                    for key,z in poly.get(source,{}).items():v.BO.add(out,key,factor*z)
            assert not out,('selected raising',mon,next(iter(out.items())))
    rhs=g.RHS()
    for mon,a,m in labels:
        mon=tuple(mon);left=Fraction(0)
        for i in v.BO.axes(m):
            source=g.union(mon,(i,));left+=source.count(i)*v.BO.wedge_sign(i,m)*poly.get(source,{}).get((a,m^(1<<i)),0)
        assert left==rhs.row((mon,a,m)),('selected closure',mon,a,m)
    result=dict(selected_exact_checks_passed=True,degree6_monomials_reconstructed=len(poly),degree6_coefficients=sum(map(len,poly.values())),full_certificate_weight_check=True,raising_monomials=len(highest),quintic_closure_coordinates=len(labels),seconds=time.time()-tick,scope='Selected verification only; run --full for the complete sixth-jet certificate.')
    print(json.dumps(result,indent=2),flush=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--fraction', action='store_true')
    parser.add_argument('--flint-dir', type=Path)
    args = parser.parse_args()
    v.main();certificate=codec.Certificate(HERE,v)
    if args.full:
        import CM_h6_sixth_full as full
        full.main(certificate,args.fraction,args.flint_dir)
    else:selected(certificate)
    print(json.dumps(dict(total_seconds=time.time()-tick,**memory_profile())),flush=True)


if __name__ == '__main__':
    main()
