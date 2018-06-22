# Literature summary

## Demographic inference w/ Kingman coalescent
### SFS-based (independent sites)
- Coventry.etal-2010 simulated coalescent trees to fit N(t) to SFS
- fastsimcoal Excoffier.etal-2013 jSFS mult subpops coalescent simulations to calculate likelihood
- fastNeutrino Bhaskar.etal-2015 (SFS, coalescent, based on Polanski and Kimmel)

### Linkage info
- PSMC Li.Durbin-2011 (patterns of heterozygocity in pairs of seqs, markov approx to ARG, pairwise coalescent HMM)
- Palamara.etal-2012 length distrib of IBD fragments for demo inference
- Harris.Nielsen-2013 haplotype sharing distribution (IBS) for demo inference
- diCal Sheehan.etal-2013 (piecewise-constant demo, multiple seqs, coalescent HMM)
- MSMC Schiffels.Durbin-2014 (coalescent hmm with multiple samples)
- SMC++ Terhorst.etal-2017 (pairwise coalescent HMM, emission is SFS conditioned on pairwise coal time of single diploid individual, spline regularization)

### Non-Kingman diffusion methods (check which allow selection)
- dadi Gutenkunst.etal-2009 (multi-population SFS, diffusion)
- Lukic.etal-2011 (SFS diffusion spectral methods, multi-population)
- Ragsdale.Gutenkunst-2017 two-locus statistics
- moments Jouganous.etal-2017 ODEs for moments of allele frequencies -> sfs

## Problems with neutral models
- Schrider.etal-2016 regular sweeps skew statistics and bias neutral estimators
- Beichman.etal-2017 models fit one aspect of the data, but can't always predict others accurately (similar to our approach)

### "A sizeable fraction of the genome is influenced by natural selection"
From Schrider: Hahn 2008, Sella et al 2009, Corbett-Detig et al 2015

## Theoretical background for SFS & coalescent:
- polanski et al. 2003
- polanski and kimmel 2003
- griffiths and tavare 1994

## Multiple mergers (definitions and causes)
### Reproductive process

### Large sample size

### Sweeps

### Interference selection

## Neutrality/multiple mergers tests and demography
SFS-based tests:
- Find a review of neutrality tests that use the SFS
-
Tests that compare against a putatively neutral background:
- bustamante et al 2001
- boyko et al 2008 (both from r&g)
## Two-site statistics
- From R&G 17:
    - Kimura '63
    - Hill and Robinson '66
    - Karlin and McGregor '68
    - Ohta and Kimura 69
    - Watterson 70

### Neutral
### Selection
### Multiple-mergers

### With recombination?
- Xie 2011 - diffusion two-locus linked neutral or equal-s selection
- Ferretti 2016 coalescent 2-locus sfs completely linked neutral locus
- kamm spence chan two-locus likelihoods under var pop size and rec rate estimation
- jenkins and song 2012 pade approximants exact two-locu samp distr.
- (same) 2010 asymptotic sampling coal. w recomb.
- (same) 2009 closed form two locus samp. dist.

## DPGP3 data set
Data description:
- ORIGINAL PAPER

Demographic models:
- Ragsdale and Gutenkunst (2017) neutral diffusion method, one- and two-site statistics, two- and three-epoch piecewise constant model
- Terhorst et al. (2017) SMC++, splines
- OTHERS?

Selection:
- Schiffels
- Garud
- OTHERS?

# Our contribution:
- SNP data
- No phasing needed
- No ancestral state needed
- No recombination map needed (although could be useful)
- Robust to population growth
- Easy to compute
- In principle, works with small sample size (maybe need to check this)
- Genome-wide background test for neutrality. (Extension: use locally)
