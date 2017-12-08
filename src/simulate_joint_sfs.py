import sys
import numpy as np
import argparse

#FIXME: this is a hack to use my local version of msprime
# sys.path.insert(1, '/users/danielrice/msprime-lambda/')
sys.path.insert(1, '/home/dpr/mmc_genomics/src/msprime-lambda/')
import msprime
import jsfs

# Constant for the length of time for a bottleneck to last
BOTTLENECK_DURATION = 1e-6

# TODO: Check what happens at the lower level when alpha != 2 and you change the population size

parser = argparse.ArgumentParser()
# TODO: argument groups
# TODO: errors

# Whole-population parameters
parser.add_argument("--nSamples", type=int, default=2, help="number of samples to simulate DEFAULT=2")
# parser.add_argument("--nLoci", type=int, default=1, help="number of independent loci to simulate DEFAULT=1")
parser.add_argument("--nLoci", nargs='*', type=int, default=[1], help="number of independent loci to simulate DEFAULT=1")
parser.add_argument("-Tc", "--coalescent_time", nargs='*', type=float, default=[1], help="initial average coalescent time DEFAULT=1")

# FIXME: make this play nicely with the recombination rate
parser.add_argument("--unlinked", action='store_true', default=False,
                    help='Simulate unlinked trees, i.e. r->Infinity.')
parser.add_argument("-r", "--recombination_rate", type=float, default=0.0, help="per-locus recombination rate per coal unit")
parser.add_argument("-G", "--growth_rate", type=float, default=0.0, help="population growth rate DEFAULT=0.0") 
parser.add_argument("-a", "--alpha", type=float, default=2.0, help="Beta-coalescent parameter alpha. DEFAULT=2.0 (Kingman)")

# Population size changes
parser.add_argument("-T", "--time", type=float, nargs='+', help='Time of population-size change', default=[])
parser.add_argument("-S", "--size", type=float, nargs='+', help='(used with -T) size of population after size change', default=[])

# Bottleneck
parser.add_argument("-Tb", "--bottleneck_time", type=float, nargs='+', help='Time of population-size bottlenck', default=[])
parser.add_argument("-Sb", "--bottleneck_strength", type=float, nargs='+', help='(used with -Tb) strength of bottleneck', default=[])

# Demes and migration
# TODO: Make this more flexible to support assymetric migration, stepping-stone, etc.
parser.add_argument("-d", "--demes", type=int, help='Number of demes')
parser.add_argument("-M", "--migration_rate", type=float, default=1.0, help='Overall symmetric migration rate')
mig_type_group = parser.add_mutually_exclusive_group()
mig_type_group.add_argument("--migration_type", type=str, default='island', help='Type of migration model. DEFAULT=island')

# Mutations and errors
parser.add_argument("-U", "--mutation_rate", type=float, default=None, help='Per-locus mutation rate per coal time. DEFAULT=None (use trees instead of mutations to calculate SFS')
parser.add_argument("-e", "--error_rate", type=float, default=None, help='Per sample per site poisson error rate.')

# Errors and site-frequency spectrum output
parser.add_argument("--individual_loci", action='store_true', default=False,
                    help='Print individual loci. Default: print SFS moments only.')

args = parser.parse_args()

# Check whether inputs are good
if len(args.nLoci) != len(args.coalescent_time):
    sys.stderr.write('ERROR: Must specify the number of loci for each Tc.')
    sys.exit(1)
if args.recombination_rate < 0.0:
    sys.stderr.write('ERROR: Recombination rate must be non-negative.\n')
    sys.exit(1)
if args.growth_rate < 0.0:
    sys.stderr.write('ERROR: Growth rate must be non-negative.\n')
    sys.exit(1)
if args.alpha <= 1.0 or args.alpha > 2.0:
    sys.stderr.write('ERROR: Multiple-merger parameter alpha must be in (1.0, 2.0].\n')
    sys.exit(1)
if args.time and args.demes:
    sys.stderr.write('ERROR: Cannot specify step-changes and multiple demes.\n')
    sys.exit(1)
if args.bottleneck_time and args.demes:
    sys.stderr.write('ERROR: Cannot specify bottlenecks and multiple demes.\n')
    sys.exit(1)
if len(args.time) != len(args.size):
    sys.stderr.write('ERROR: Must specify the same number of population-size change times and sizes.\n')
    sys.exit(1)
if len(args.bottleneck_time) != len(args.bottleneck_strength):
    sys.stderr.write('ERROR: Must specify the same number of bottleneck times and strengths.\n')
    sys.exit(1)
if args.demes and args.nSamples % args.demes:
    sys.stderr.write("Error: must nSamples must be divisible by demes!\n")
    sys.exit(1)
if args.error_rate and not args.mutation_rate:
    sys.stderr.write("ERROR: must specify a mutation rate to have errors!\n")
    sys.exit(1)

# TODO: Make sure all the parameter outputs conform to the new standards
sys.stdout.write('#N_SAMPLES={}\n'.format(args.nSamples))
sys.stdout.write('#N_LOCI={}\n'.format(','.join([str(nL) for nL in args.nLoci])))
sys.stdout.write('#Tc={}\n'.format(','.join([str(Tc) for Tc in args.coalescent_time])))
if args.mutation_rate:
    sys.stdout.write('#SFS_FROM=MUTATIONS\n')
    sys.stdout.write('#MUTATION_RATE={}\n'.format(args.mutation_rate))
else:
    sys.stdout.write('#SFS_FROM=TREES\n')
if args.error_rate:
    sys.stdout.write('#ERROR_RATE={}\n'.format(args.error_rate))
sys.stdout.write('#RECOMBINATION_RATE={}\n'.format(args.recombination_rate))
sys.stdout.write('#GROWTH_RATE={}\n'.format(args.growth_rate))
sys.stdout.write('#ALPHA={}\n'.format(args.alpha))

# SET PARAMETERS

# Set defaults
migration_matrix = None 
population_configurations = None
demographic_events = []

# Population structure
if args.demes:
    # TODO: check migration rate normalizations
    if args.migration_type == 'island':
        # Normalize mutation rate to a per-generation per-deme basis
        m = args.migration_rate / (4.0 * (args.demes - 1.0))
        migration_matrix = np.ones((args.demes,args.demes)) - np.eye(args.demes)
    elif args.migration_type == 'steppingstone':
        # Normalize mutation rate to a per-generation per-deme basis
        m = args.migration_rate / (4.0 * 2.0)
        # Stepping-stone migraton
        migration_matrix = np.diag(np.ones(args.demes - 1), k=1)
        migration_matrix[0,-1] = 1
        migration_matrix += migration_matrix.T
    migration_matrix *= m
    migration_matrix = migration_matrix.tolist()
    samples_per_deme = args.nSamples / args.demes
    population_configurations = [msprime.PopulationConfiguration(
                                        sample_size=samples_per_deme,
                                        growth_rate=args.growth_rate,
                                        multiple_merger_para=args.alpha) 
                                    for d in range(args.demes)]
    sys.stdout.write('#DEMES={}\n#MIGRATION_RATE={}\n'.format(args.demes, args.migration_rate))
    sys.stdout.write('#MIGRATION_TYPE={}\n'.format(args.migration_type))

else:
    population_configurations = [msprime.PopulationConfiguration(sample_size=args.nSamples,
                                    growth_rate=args.growth_rate,
                                    multiple_merger_para=args.alpha)]
    # Population size changes only
    # FIXME: make sure this is compatable with multiple-mergers
    if args.time:
        for t, s in zip(args.time, args.size):
            demographic_events.append(msprime.PopulationParametersChange(t, initial_size=s/2.0))
        sys.stdout.write('#DEMOCHANGE_TIMES={}\n'.format(','.join([str(t) for t in args.time])))
        sys.stdout.write('#DEMOCHANGE_SIZES={}\n'.format(','.join([str(s) for s in args.size])))
    if args.bottleneck_time:
        # FIXME: don't hard-code the post-bottleneck size
        for tb, strength in zip(args.bottleneck_time, args.bottleneck_strength):
            b_size = BOTTLENECK_DURATION / strength
            demographic_events.append(msprime.PopulationParametersChange(tb, initial_size=b_size/2))
            demographic_events.append(msprime.PopulationParametersChange(tb+BOTTLENECK_DURATION, initial_size=0.5))
        sys.stdout.write('#BOTTLENECK_TIMES={}\n'.format(','.join([str(tb) for tb in args.bottleneck_time])))
        sys.stdout.write('#BOTTLENECK_STRENGTHS={}\n'.format(','.join([str(s) for s in args.bottleneck_strength])))
        sys.stdout.write('#BOTTLENECK_DURATION={}\n'.format(BOTTLENECK_DURATION))

# Print header
if args.individual_loci:
    header = 'N_SITES\t' + '\t'.join(['tau_{}'.format(i) for i in range(1, args.nSamples)]) + '\n'
    sys.stdout.write(header)

# RUN SIMULATION
mTAU  = np.zeros((args.nSamples - 1))
jTAU  = np.zeros((args.nSamples -1, args.nSamples - 1))
cumulative_loci = 0
for nLoci, Tc in zip(args.nLoci, args.coalescent_time):
    sys.stderr.write('{} loci with Tc={}\n'.format(nLoci, Tc))
    sys.stderr.write('Simulating genealogies...\n')
    if args.unlinked:
        simulations = msprime.simulate(
                recombination_rate=0,
                demographic_events=demographic_events,
                population_configurations=population_configurations,
                migration_matrix=migration_matrix,
                num_replicates=2*nLoci,
                Ne=Tc/2.0,
                mutation_rate=args.mutation_rate)
        for rep1 in simulations:
            rep2 = next(simulations)
            tau1, _ = jsfs.get_jsfs(rep1)
            tau2, _ = jsfs.get_jsfs(rep2)
            mTAU += tau1
            mTAU += tau2
            jTAU += tau1[None,:] * tau2[:,None]

    else:
        simulations = msprime.simulate(
                recombination_rate=args.recombination_rate,
                demographic_events=demographic_events,
                population_configurations=population_configurations,
                migration_matrix=migration_matrix,
                num_replicates=nLoci,
                Ne=Tc/2.0,
                mutation_rate=args.mutation_rate)

        # CALCULATE AND OUTPUT SFS
        sys.stderr.write('Calculating times...\n')
        for i_rep, rep in enumerate(simulations):
            try:
                if (i_rep+1) % (nLoci/10) == 0:
                    sys.stderr.write('Completed {} reps.\n'.format(i_rep+1))
            except ZeroDivisionError:
                pass

            tau1, tau2 = jsfs.get_jsfs(rep)
            mTAU += tau1
            mTAU += tau2
            jTAU += tau1[None,:] * tau2[:,None]

    cumulative_loci += nLoci

# if not args.individual_loci:
sys.stderr.write('Calculating average...\n')

mTAU /= 2.0*cumulative_loci
jTAU += jTAU.T
jTAU /= 2.0*cumulative_loci

sys.stdout.write(' '.join([str(x) for x in mTAU]) + '\n')
sys.stdout.write(' '.join([str(x) for x in jTAU[np.triu_indices(args.nSamples-1)]]) + '\n')
