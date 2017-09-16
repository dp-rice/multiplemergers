snakemake -j30 --cluster-config config/cluster.json --cluster "sbatch -J {cluster.name} -t {cluster.time} --mem-per-cpu {cluster.mem} -c {cluster.ncores} -o {cluster.out} -e {cluster.err} --partition={cluster.partition}" "$@"