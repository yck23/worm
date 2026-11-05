#!/usr/bin/env python3
"""Scalar response regressions: periodic exact Hamiltonian and histogram identity."""
from pathlib import Path
import itertools
import json
import os
import subprocess
import tempfile
import h5py
import numpy as np


def exact_chi(potential, hopping, beta):
    basis=np.array([n for n in itertools.product(range(2),repeat=9) if sum(n)==2])
    lookup={tuple(n):i for i,n in enumerate(basis)}
    bonds=[(x+3*y,(x+1)%3+3*y) for y in range(3) for x in range(3)]
    bonds += [(x+3*y,x+3*((y+1)%3)) for y in range(3) for x in range(3)]
    H=np.diag(basis@potential)
    O=np.eye(len(basis))*2
    for col,n in enumerate(basis):
        for source in range(9):
            for dest in range(9):
                if source!=dest and n[source] and not n[dest]:
                    updated=n.copy();updated[source]-=1;updated[dest]+=1
                    O[lookup[tuple(updated)],col]+=1
        for a,b in bonds:
            for source,dest in [(a,b),(b,a)]:
                if n[source] and not n[dest]:
                    updated=n.copy();updated[source]-=1;updated[dest]+=1
                    H[lookup[tuple(updated)],col]-=hopping
    energy,vectors=np.linalg.eigh(H)
    weight=np.exp(-beta*(energy-energy.min()));weight/=weight.sum()
    rho=(vectors*weight)@vectors.T
    return np.sum(rho*O.T)/9


def main():
    root=Path(__file__).resolve().parents[1]
    work=root/'build/response-validation'
    work.mkdir(exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='check-',dir=work))
    programmable=root/'dist/windows-square-programmable/qmc_worm_canonical.exe'
    uniform=root/'dist/windows-square/qmc_worm_canonical.exe'
    rows=[]
    for name,potential,hopping,exe in [
        ('zero_hopping',np.linspace(-.6,.8,9),0.,programmable),
        ('nonuniform',np.linspace(-.6,.8,9),.45,programmable),
        ('uniform',np.zeros(9),.45,uniform)]:
        case=work/name;case.mkdir()
        np.savetxt(case/'potential.dat',potential,fmt='%.17g')
        params=dict(model='"BoseHubbard"',Lx=3,Ly=3,Lz=1,pbcx=1,pbcy=1,pbcz=0,
            t_hop=hopping,U_on=1,V_nn=0,mu=0,beta=.7,nmax=1,canonical=2,canonical_window=.5,
            thermalization=10000,sweeps=400000,runtimelimit=120,seed=653128,
            Nmeasure=1,Nmeasure2=20,Ntest=1000000,Nsave=100000000,C_worm=2,E_off=1,
            onsite_potential_file='"potential.dat"' if exe==programmable else '""',
            outputfile='"result.h5"',checkpoint='"state.h5"')
        (case/'parameters.ini').write_text('\n'.join(f'{k} = {v}' for k,v in params.items())+'\n')
        env=os.environ.copy();env['PATH']=str(exe.parent)+os.pathsep+env['PATH'];env['TMPDIR']=case.as_posix()+'/'
        def launch(args,log):
            with (case/log).open('w') as stream:
                done=subprocess.run([str(exe),*args],cwd=case,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=150)
            assert done.returncode==0,(name,log)
        launch(['parameters.ini'],'run.log')
        with h5py.File(case/'result.h5') as data:
            result=data['simulation/results']
            chi=float(result['Equal_time_susceptibility/mean/value'][()])
            error=float(result['Equal_time_susceptibility/mean/error'][()])
            exact=exact_chi(potential,hopping,.7)
            if hopping==0:
                # ALPS can report NaN errors for a strictly constant sequence.
                np.testing.assert_allclose(chi,exact,rtol=0,atol=1e-12)
            else:
                assert np.isfinite(error) and abs(chi-exact)<7*error+.006,(name,chi,error,exact)
            total=float(result['Winding_number_squared_total/mean/value'][()])
            np.testing.assert_allclose(total,np.sum(result['Winding_number_squared/mean/value'][()]),atol=1e-12)
            assert float(result['Number_of_particles/mean/value'][()])==2
            if name=='uniform':
                dm=result['Density_Matrix/mean/value'][()]
                np.testing.assert_allclose(chi,dm.sum()-dm[0]+2/9,rtol=2e-12,atol=2e-12)
            rows.append(dict(case=name,chi_s=chi,se=error,exact=float(exact)))
        if name=='nonuniform':
            launch(['state.h5','--reset-statistics'],'restart.log')
            with h5py.File(case/'result.h5') as data:
                assert int(data['simulation/results/Equal_time_susceptibility/count'][()])==20000
        print(f'PASS {name}: chi={chi:.6f} +/- {error:.6f}; exact={exact:.6f}',flush=True)
    (work/'report.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('Evidence:',work)


if __name__=='__main__':
    main()
