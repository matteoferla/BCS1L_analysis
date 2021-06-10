## Code

These snippets were run on a Jupyter notebook running off a cluster.
The `pyrosetta_help` is my collection of helper functions
([pip](https://pypi.org/project/pyrosetta-help/), [GitHub](https://github.com/matteoferla/pyrosetta_help))
and `rdkit_to_params` is my opensource topology generator ([pip](https://pypi.org/project/rdkit-to-params/), [GitHub](https://github.com/matteoferla/rdkit_to_params))

### Init
```python
import pyrosetta
from pyrosetta_help import make_option_string, configure_logger, get_log_entries

# capture to log
logger = configure_logger()
# give CLI attributes in a civilised way
pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                mute='all',
                                                ignore_unrecognized_res=False,
                                                load_PDB_components=False,
                                                ignore_waters=False)
                               )
```

### ASG
```python
from rdkit_to_params import Params
# https://www.rcsb.org/ligand/AGS
p = Params.from_smiles_w_pdbfile('AGS.pdb',
                             'c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)([O-])OP(=O)([O-])OP(=S)([O-])[O-])O)O)N',
                              name='AGS')

p.dump('AGS.params')
```

### Loading
```python
## PDBs
# 6UKP 6UKS
from pyrosetta_help import download_pdb

print(download_pdb('6UKP'), download_pdb('6UKS'))

## ED
# EMDB 20808, PDB 6UKP (Apo mBcs1) and EMDB 20811, PDB 6UKS (ATPγS-bound mBcs1)
from pyrosetta_help import download_map

download_map('EMD-20808'), download_map('EMD-20811')
ed_maps = {'6UKP': 'EMD-20808', '6UKS': 'EMD-20811'}

## gunzip
import gzip
import shutil
for m in ed_maps.values():
    with gzip.open(f'{m}.map.gz', 'rb') as f_in:
        with open(f'{m}.map', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
```

For each of the two:
```python
from pyrosetta_help import pose_from_file, prep_ED

code = '6UKS' # or '6UKP'
pose = pose_from_file(f'{code.lower()}.pdb', ['AGS.params'])
ed = prep_ED(pose, f'{ed_maps[code]}.map')
print('How good is the density match:', ed.matchPose(pose))
```
The ED match for 6UKS is 71% uniminimised and 69% minimised.
For 6UKP is 66% and 61%.

## Non-symmetrical
It is a heptamer, but it is not symmetrical.
```python
ds = pyrosetta.rosetta.protocols.symmetry.DetectSymmetry(5, .01) 
ds.apply(pose) # ERROR: rmsd between subunits higher than subunit tolerance
```

## Minimisation

Three initial cartesian cycles of FastRelax

```python
scorefxnED = pyrosetta.create_score_function('ref2015_cart')
elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
scorefxnED.set_weight(elec_dens_fast, 30)

movemap = pyrosetta.MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)
movemap.set_jump(True)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, 3)
relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
relax.set_movemap(movemap)
relax.minimize_bond_angles(True)
relax.minimize_bond_lengths(True)
relax.cartesian(True)
relax.apply(pose)
```
Then a descending weighted dihedral FastRelax
```python
scorefxnED = pyrosetta.create_score_function('ref2015')
elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
scorefxnED.set_weight(elec_dens_fast, 30)

for w in (30, 20, 10):
    scorefxnED.set_weight(elec_dens_fast, w)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, 3)
    relax.apply(pose)
    print(f'complete {w}')
reg_scorefxn = pyrosetta.get_fa_scorefxn()
print('How good is the density match:', ed.matchPose(pose))
print(reg_scorefxn.get_name(), reg_scorefxn(pose))
pose.dump_pdb(f'{code.lower()}.r.pdb')
print(f'{code.lower()}.r.pdb')
```

### Threading
```python
from pyrosetta_help import  (get_alignment, 
                             write_grishin, 
                             thread, 
                             rangify, 
                             get_nonprotein_pose, 
                             oligomer_thread,correct_numbering
                             )

# >sp|Q9Y276|BCS1_HUMAN Mitochondrial chaperone BCS1 OS=Homo sapiens OX=9606 GN=BCS1L PE=1 SV=1
bcs1l_seq = '''MPLSDFILALKDNPYFGAGFGLVGVGTALALARKGVQLGLVAFRRHYMITLEVPARDRSY
AWLLSWLTRHSTRTQHLSVETSYLQHESGRISTKFEFVPSPGNHFIWYRGKWIRVERSRE
MQMIDLQTGTPWESVTFTALGTDRKVFFNILEEARELALQQEEGKTVMYTAVGSEWRPFG
YPRRRRPLNSVVLQQGLADRIVRDVQEFIDNPKWYTDRGIPYRRGYLLYGPPGCGKSSFI
TALAGELEHSICLLSLTDSSLSDDRLNHLLSVAPQQSLVLLEDVDAAFLSRDLAVENPVK
YQGLGRLTFSGLLNALDGVASTEARIVFMTTNHVDRLDPALIRPGRVDLKEYVGYCSHWQ
LTQMFQRFYPGQAPSLAENFAEHVLRATNQISPAQVQGYFMLYKNDPVGAIHNAESLRR'''.replace('\n', '')

template = pose
pose, threaders, unaltered_vector = oligomer_thread(pose, bcs1l_seq)
correct_numbering(pose)
pose.dump_pdb(f'{code.lower()}.thr.pdb')
```

### Residue deletion
...

### Final minimisation

```python
relax = pyrosetta.rosetta.protocols.relax.FastRelax(pyrosetta.create_score_function('ref2015_cart'), 3)
relax.minimize_bond_angles(True)
relax.minimize_bond_lengths(True)
relax.cartesian(True)
relax.apply(pose)

relax = pyrosetta.rosetta.protocols.relax.FastRelax(pyrosetta.get_fa_scorefxn(), 3)
relax.apply(pose)
```
### Missing TM helices
In hindsight I should have done the maths for the rototranslations to place 7 pillars at a given spots perpendicular to each other in PyMOL.
Instead I did it in a more crude way.

Made pillars:
```python
import pymol2
missing = 'MPLSDFILALKDNPYFGAGFGLVGVGTALALARKGVQLGLVAFRRHYMIT'

with pymol2.PyMOL() as pymol:
    pymol.cmd.fetch('6UKS')
    pymol.cmd.alter('*', 'segi=""')
    for chain in 'ABCDEFG':
        pymol.cmd.fab(missing, name=f'pillar{chain}', resi=1, chain=chain, ss=1)
        pymol.cmd.align(f'resi 50 and pillar{chain}', f'resi 50 and chain {chain} and 6UKS')
    pymol.cmd.delete('6UKS')
    pymol.cmd.save('pillars.pdb')
```
Loaded in PyRosetta
```python
from pyrosetta_help import pose_from_file
reference = pose_from_file('6uks.pdb', ['AGS.params'])
pose = pyrosetta.pose_from_pdb('pillars.pdb')
```
Added constrains to make a spiderweb of constraints for each residue to its equivalents in the other chains.
To prevent an oblique heptagonal prism, the 1:X-50:X-1:Y angle was forced to 90°.

```python
from itertools import permutations

def get_AtomID(chain:str, resi:int, atomname:str) -> pyrosetta.rosetta.core.id.AtomID:
    r = pose.pdb_info().pdb2pose(res=resi, chain=chain)
    assert r != 0, f'{resi}:{chain} is absent'
    residue = pose.residue(r)
    return pyrosetta.rosetta.core.id.AtomID(atomno_in=residue.atom_index(atomname), rsd_in=r)

HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
cons = []
for chain_A, chain_B in permutations('ABCDEFG', 2):
    A_xyz = reference.residue(reference.pdb_info().pdb2pose(res=51, chain=chain_A)).xyz("CA")
    B_xyz = reference.residue(reference.pdb_info().pdb2pose(res=51, chain=chain_B)).xyz("CA")
    distance = (A_xyz - B_xyz).norm()
    for a in range(1, 50):
        cons.append(AtomPairConstraint(get_AtomID(chain_A, a, 'CA'),
                                       get_AtomID(chain_B, a, 'CA'),
                                       HarmonicFunc(x0_in=distance, sd_in=3))
                   )
    A_xyz = pose.residue(pose.pdb_info().pdb2pose(res=1, chain=chain_A)).xyz("CA")
    B_xyz = pose.residue(pose.pdb_info().pdb2pose(res=50, chain=chain_A)).xyz("CA")
    height = (A_xyz - B_xyz).norm()
    diagonal = (distance**2 + height**2)**0.5
    cons.append(AtomPairConstraint(get_AtomID(chain_A, 1, 'CA'),
                                       get_AtomID(chain_B, 50, 'CA'),
                                       HarmonicFunc(x0_in=diagonal, sd_in=3))
                   )
    cons.append(AtomPairConstraint(get_AtomID(chain_A, 50, 'CA'),
                                       get_AtomID(chain_B, 1, 'CA'),
                                       HarmonicFunc(x0_in=diagonal, sd_in=3))
                   )

cl = pyrosetta.rosetta.utility.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t()
cl.extend(cons)
cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
cs.add_constraints(cl)
setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
setup.constraint_set(cs)
setup.apply(pose)
```
The "pillars" were added to the 6UKS model prior to threading.
Remodel was considered but this is far faster.

### Scoring
NB. variant not in this repo

Doublecheck chains, using `model.has_interface(None, 'A_BCDEFGHIJKLMNOPQRST')`
```python
## chec
pi = pose.pdb_info()
possible_chains = {pi.chain(res+1) for res in range(pose.total_residue())}
possible_chains
```
Score...
```python
import pandas as pd
from pyrosetta_help import MutantScorer, extend_scores

model = MutantScorer(pose=pose, modelname=code)
variants = ['👾👾👾', '👾👾👾', '👾👾👾']
scoresx = []
data = model.score_mutations(variants,
                            chain='A',
                            interfaces=(('chainA', 'A_BCDEFG'),),
                            preminimise=False,
                            distance=12,
                            cycles=5)
scoresx.append(pd.DataFrame(data))
scores = pd.concat(scoresx)
extend_scores(scores)
scores.to_csv('scores.csv')
scores
```
