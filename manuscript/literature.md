# Literature summary

## Demographic inference w/ Kingman coalescent
### SFS-based (independent sites)
- Coventry.etal-2010 simulated coalescent trees to fit N(t) to SFS
- fastsimcoal Excoffier.etal-2013 jSFS mult subpops coalescent simulations to calculate likelihood
- fastNeutrino Bhaskar.etal-2015 (SFS, coalescent, based on Polanski and Kimmel)

### Use linkage info
- PSMC Li.Durbin-2011 (patterns of heterozygocity in pairs of seqs, markov approx to ARG, pairwise coalescent HMM)
- Palamara.etal-2012 length distrib of IBD fragments for demo inference
- Harris.Nielsen-2013 haplotype sharing distribution (IBS) for demo inference
- diCal Sheehan.etal-2013 (piecewise-constant demo, multiple seqs, coalescent HMM)
- MSMC Schiffels.Durbin-2014 (coalescent hmm with multiple samples)
- SMC++ Terhorst.etal-2017 (pairwise coalescent HMM, emission is SFS conditioned on pairwise coal time of single diploid individual, spline regularization)

### Non-coalescent diffusion methods (check which allow selection)
- dadi Gutenkunst.etal-2009 (multi-population SFS, diffusion)
- Lukic.etal-2011 (SFS diffusion spectral methods, multi-population)
- Ragsdale.Gutenkunst-2017 two-locus statistics
- moments Jouganous.etal-2017 ODEs for moments of allele frequencies -> sfs

## Problems with neutral models
- Ewing.Jensen-2015 background selection distorts demographic inference
- Schrider.etal-2016 selective sweeps skew statistics and bias neutral estimators
- Beichman.etal-2017 models fit one aspect of the data, but can't always predict others accurately (similar to our approach)
- Cvijovic.etal-2018 any level of bg sel strong enough to reduce the pairwise diversity will distort the SFS

### "A sizeable fraction of the genome is influenced by natural selection"
- Hahn-2008 Review of "neutral theory" and it's selective alternatives
- Sella.etal-2009 Review of evidence for pervasive selection in Drosophila
- Corbett-Detig.etal-2015 Pi/rho correlations in many species
- Elyashiv.etal-2016 linked selection reducing diversity in D. mel.

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
- Rafajlovic et al 2014

## Theory
### One-site SFS
- polanski et al. 2003
- polanski and kimmel 2003
- griffiths and tavare 1994
- griffiths and tavare 1998

### Two-site SFS, no recombination:
First two moments of SFS (can be used to derive completely linked 2-SFS):
- Fu-1995
- Zivkovic.Wiehe-2008 coalescent, second-order moments of segragating sites; extends Fu-1995 to time-varying population size
- Birkner.etal-2013 lambda coalescent

Related but not main citations:
- Eriksson.etal-2010 (coalescent, moments of total branch length with variable population size)
- Rafajlovic.etal-2014 (coalescent, first two SFS moments in piecewise constant, no recombination)
- Xie-2011 diffusion two-locus sfs perfectly linked neutral or equal-s selection
- Sargsyan-2014 coalescent non-recombining locus pairs of polymorphic sites "analytical formulas for the numbers of the topologies of genealogies with two mutation events" (relates to the 2-sfs but it's unclear to me exactly what was done)
- Ferretti.etal-2016 coalescent 2-sfs completely linked neutral locus, whole population. also conditional 1-SFS (conditioned on freq at other site)

### 2-SFS with recombination (?):


### CHECK THESE:
- kamm spence chan two-locus likelihoods under var pop size and rec rate estimation
- jenkins and song 2012 pade approximants exact two-locu samp distr.
- (same) 2010 asymptotic sampling coal. w recomb.
- (same) 2009 closed form two locus samp. dist.
- griffiths and tavare 2003
- Hobolth and Wiuf 2009

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
