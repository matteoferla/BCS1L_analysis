# BCS1 analysis
An analysis of BCS1 (no pathogenic variants stored here)

## Modelling

_BCS1L_ is the gene name and BCS1 is the protein name.
It is a homologue of Paraplegin.
There are no structures of human BCS1, but there are of mouse.

Murine BCS1 is 92% identical to human BCS1. This has been solved by crystallography and cryoEM in 

> Tang, W.K., Borgnia, M.J., Hsu, A.L., Esser, L., Fox, T., de Val, N., Xia, D. (2020) Nat Struct Mol Biol 27: 202-209
> [10.1038/s41594-020-0373-0](http://dx.doi.org/10.1038/s41594-020-0373-0)

As the following:

* ΔNmBcs1-AMP-PNP PDB:6U1Y (2.2 Å)
* FLmBcs1-ADP PDB:6UKO (4.4 Å)
* Apo mBcs1 EMDB-20808 PDB:6UKP
* ATPγS-bound mBcs1 EMDB-20811 PDB:6UKS

6U1Y (ADP) has waters and is nice, but is very ∆N.
6UKO has a low resolution, but a nice ligand and an N-terminal transmembrane helix.
6UKP has no lingand, but has an N-terminal transmembrane helix.
6UKS has ATP mimic, but lacks the N-terminal transmembrane helix.

The helix is resolved at L29 onwards. There should be a cleaved MTS signalling peptide.
The N-terminal sequence, `MPLSDFILALKDNPYFGAGFGLVGVGTA`+`LALARKGVQLGLVAFRRHYM` may contain a MTS,
but it is a non-canonical form as ELM cannot find it. 
The distance between L29 and M48 is 30 Å and is kinked by 30° inwards.
The IMM is thicker than the regular plasma membrane, which is 30Å.
But let's assume that M1-A28 is the MTS and L29-M48 is the transmembrane span and the N-terminal nitrogen is okay,
either unhappy in the membrane or in contact with the phosphate heads.

Despite the 92% homology, the premade SwissModel has a Qmean of –2, possibly due to subpar fitting.
It is a heptamer, so I-Tasser or Phyre2 would not yield a nice fit. As a result, the structures were energy minimised and threaded all within PyRosetta.
Namely,

1. 6UKP and 6UKS PDB and cryoEM map were downloaded
2. ATPγS was parameterised properly
3. density map constrained energy minimisation (`FastRelax`)
4. one-to-one threading (`ThreadingMover` part of RosettaCM)
5. removal of 1-28, 288-309 and 417-418
6. energy minimisation
7. variant scoring

The symmetry feature could not be used due to excessive difference between chains.
The missing density 288-309 and 417-418 were not modelled due to time concerns.
Likewise in hindsight, the model 6UKO is actually good and I should have used that.
See [code notes](code.md) for more.

