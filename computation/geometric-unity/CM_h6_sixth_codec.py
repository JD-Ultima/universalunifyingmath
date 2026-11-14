"""Portable integer-stream decoder for the CM factored rational jet."""
from pathlib import Path
from fractions import Fraction
from itertools import combinations_with_replacement
import gzip,json,hashlib,time

class Reader:
    def __init__(self,path,magic):
        self.data=gzip.decompress(path.read_bytes());assert self.data[:len(magic)]==magic;self.pos=len(magic)
    def unsigned(self):
        n=shift=0
        while True:
            b=self.data[self.pos];self.pos+=1;n|=(b&127)<<shift
            if b<128:return n
            shift+=7
    def signed(self):
        n=self.unsigned();return n//2 if n%2==0 else -(n//2)-1
    def integer(self):
        tag=self.unsigned();size=tag//2;n=int.from_bytes(self.data[self.pos:self.pos+size],'little');self.pos+=size
        return -n if tag%2 else n
    def end(self):assert self.pos==len(self.data)

class Certificate:
    def __init__(self,folder,geometry):
        self.folder=Path(folder);self.v=geometry;self.mons=list(combinations_with_replacement(range(14),6));self.loc={m:i for i,m in enumerate(self.mons)}
        self.weights=[tuple(sum(geometry.weights[i][z] for i in m) for z in range(2)) for m in self.mons]
        self.aux=json.loads(gzip.decompress((self.folder/'CM_h6_sixth_aux.json.gz').read_bytes()))
        manifest=json.loads((self.folder/'CM_h6_sixth_manifest.json').read_text())
        for name,rec in manifest['files'].items():
            data=(self.folder/name).read_bytes();assert len(data)==rec['bytes'] and hashlib.sha256(data).hexdigest()==rec['sha256']
    def load(self,monomials=None,Q=Fraction):
        selected=None if monomials is None else {self.loc[m] for m in monomials};records=[];tags=[];needed=set()
        r=Reader(self.folder/'CM_h6_sixth_coefficients.bin.gz',b'CM6C\x01');den=[r.integer() for _ in range(r.unsigned())];assert all(d>0 for d in den)
        count=r.unsigned();assert count==self.aux['coefficients']
        for _ in range(count):
            mi,vi,n,di=r.unsigned(),r.unsigned(),r.integer(),r.unsigned();assert mi<len(self.mons) and di<len(den) and n
            tags.append((mi,vi))
            if selected is None or mi in selected:records.append((mi,vi,n,den[di]));needed.add(vi)
        r.end();r=Reader(self.folder/'CM_h6_sixth_vectors.bin.gz',b'CM6V\x01');nv=r.unsigned();assert nv==self.aux['basis_vectors'];vectors={};weights=[];entries=0;fw={}
        for vi in range(nv):
            ki=r.unsigned();wt=(r.signed(),r.signed());copy,j=r.unsigned(),r.unsigned();D=r.integer();count=r.unsigned();assert D>0;weights.append(wt);tensor={}
            for _ in range(count):
                a,m,n=r.unsigned(),r.unsigned(),r.integer();assert a<14 and m.bit_count()==6 and m<1<<14 and n;entries+=1
                if m not in fw:fw[m]=tuple(sum(self.v.weights[i][z] for i in self.v.BO.axes(m)) for z in range(2))
                assert tuple(self.v.weights[a][z]+fw[m][z] for z in range(2))==wt
                if vi in needed:tensor[a,m]=Q(n,D)
            if vi in needed:vectors[vi]=tensor
        r.end();assert entries==self.aux['basis_entries']
        for mi,vi in tags:assert vi<nv and tuple(self.weights[mi][z]+weights[vi][z] for z in range(2))==(5,5)
        grouped={}
        for mi,vi,n,d in records:grouped.setdefault(self.mons[mi],[]).append((vi,Q(n,d)))
        return grouped,vectors
    def expand(self,Q=Fraction,monomials=None):
        tick=time.time();coords,vectors=self.load(monomials,Q);out={}
        for n,(mon,terms) in enumerate(coords.items()):
            tensor={}
            for vi,c in terms:
                for key,z in vectors[vi].items():self.v.BO.add(tensor,key,c*z)
            if tensor:out[mon]=tensor
            if monomials is None and n%1000==0:print('CM exact expansion',n,'/',len(coords),'seconds',time.time()-tick,flush=True)
        return out
