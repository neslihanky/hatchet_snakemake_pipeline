rule final:
    input:
        expand("{patient}/bb/bulk.bb", patient=config['Patients'].keys()),
        expand("{patient}/bbc/bulk.bbc", patient=config['Patients'].keys()),
        expand("{patient}/bbc/bulk.seg", patient=config['Patients'].keys()),
        expand("{patient}/cluster_plots/bb.pdf", patient=config['Patients'].keys()),
#        expand("{patient}/results/best.bbc.ucn", patient=config['Patients'].keys()),


        
rule genome_binning:
    input:
        normal   = lambda wildcards: config["Patients"][wildcards.patient]["normal"],
    output:
        normalbin= "{patient}/bin/normal.bin",
        tumorbin = "{patient}/bin/bulk.bin",
        total_Rc = "{patient}/bin/total_read.counts",
    params:
        script=config['UTILS'],
        readqual = 20,
        binsize  = "50kb",
        g_build  = "hg19",
        tumors   = lambda wildcards: config["Patients"][wildcards.patient]["tumors"],
        allnames = lambda wildcards: config["Patients"][wildcards.patient]["allnames"]
    threads:
        22
    log:
        out = "{patient}/bin/binBAM.o.log",
        err = "{patient}/bin/binBAM.e.log"
    message: 
        "Binning genomes to 50kb windows"
    shell:
        "python {params.script}binBAM.py -N {input.normal} -T {params.tumors} -S {params.allnames} -t {output.total_Rc} -b {params.binsize} -g {params.g_build} -j {threads} -q {params.readqual} -O {output.normalbin} -o {output.tumorbin} 2> {log.err} 1> {log.out}" 

# genome_binning:
# bin size is selected as 50kb as suggested for WGS data. 200-250 kb is suggested for WES data. 
# -q is read quaility        
        
rule SNP_calling_and_BAF:
    input:
        ref      = config['REF'],
        normal   = lambda wildcards: config["Patients"][wildcards.patient]["normal"]
    output:
        normalbaf= "{patient}/baf/normal.baf",
        tumorbaf = "{patient}/baf/bulk.baf",
        SNPs= "{patient}/baf/selectedSNPs.csv"
    params:
        SAM      = config['SAM'],
        BCF      = config['BCF'],
        script   = config['UTILS'],
        min_cov  = 4,
        max_cov  = 200, 
        readqual = 20,
        basequal = 20,
        snpqual  = 20,
        tumors   = lambda wildcards: config["Patients"][wildcards.patient]["tumors"],
        allnames = lambda wildcards: config["Patients"][wildcards.patient]["allnames"]
    threads:
        22
    log:
        out = "{patient}/baf/deBAF.o.log",
        err = "{patient}/baf/deBAF.e.log"
    shell:
        """
        export PATH=$PATH:{params.SAM}
        export PATH=$PATH:{params.BCF}
        python {params.script}deBAF.py -N {input.normal} -T {params.tumors} -S {params.allnames} -r {input.ref} -j {threads} -l {output.SNPs} -q {params.readqual} -Q {params.basequal} -U {params.snpqual} -c {params.min_cov} -C {params.max_cov} -O {output.normalbaf} -o {output.tumorbaf} 2> {log.err} 1> {log.out}
        """
# SNP_calling_anf_BAF:
# -c and -C are min and max read depth thrseholds for a SNP (for max, twice higher than expected avg coverage is suggested).
# For WES, -e parameter can be used to provide path for BED file.

rule combine_bin_and_SNPs:
    input:
        normal_bin = "{patient}/bin/normal.bin",
        bulk_bin   = "{patient}/bin/bulk.bin",
        baffile    = "{patient}/baf/bulk.baf"
    output:
        "{patient}/bb/bulk.bb"
    params:
        script = config['UTILS'],
        seed   = 12,
        method = "MIRROR", 
        shift  = 0.08, 
    shell:
        "python {params.script}comBBo.py -c {input.normal_bin} -C {input.bulk_bin} -B {input.baffile} -d {params.shift} -m {params.method} -e {params.seed} > {output}"

#
# combine_bin_and_SNP
# -d is Maximum expected shift from 0.5 for BAF of diploid or tetraploid clusters. 

rule global_clustering:
    input: 
        "{patient}/bb/bulk.bb"
    output:
        seg = "{patient}/bbc/bulk.seg",
        bbc = "{patient}/bbc/bulk.bbc"
    params:
        bnpy   = config['BNPY'],
        script = config['UTILS'],
        shift  = 0.08,
        tlr_rdr= 0.15,
        tlr_baf= 0.03
    shell:
        "python {params.script}cluBB.py {input} -by {params.bnpy} -o {output.seg} -O {output.bbc} -e $RANDOM -tB {params.tlr_baf} -tR {params.tlr_rdr} -d {params.shift}"

# global_clustering
# -tR and -tB are Thresholds for clustering refinement. Used to merge clusters whose difference is no more than these values in all samples. In case of overclustering, can be increased. (In demo, they increase rdr tolerance to 0.5 from 0.15)
#Bootstrapping is provided and suggested for WES data (-u 20 -dR 0.002 -dB 0.002)

rule plot_clustering_results:
    input:
        "{patient}/bbc/bulk.bbc"
    output:
         dir = "{patient}/cluster_plots/",
         rdr = "{patient}/cluster_plots/readdepthratio.pdf",
         crd = "{patient}/cluster_plots/readdepthratio_clustered.pdf",
         baf = "{patient}/cluster_plots/ballelefrequency.pdf",
         bb  = "{patient}/cluster_plots/bb.pdf",
         cbb = "{patient}/cluster_plots/bb_clustered.pdf"
    params:
        script = config['UTILS']
    shell:
        """
        python {params.script}BBot.py -x {output.dir} -c RD --figsize 6,3 {input} 
        python   {params.script}BBot.py -x {output.dir} -c CRD --figsize 6,3 {input}
        python  {params.script}BBot.py  -x {output.dir} -c BAF --figsize 6,3 {input}
        python  {params.script}BBot.py -x {output.dir} -c BB --figsize 6,3 {input}
        python  {params.script}BBot.py -x {output.dir} -c CBB --figsize 6,3 {input}
        """ 
        
        
#rule Hatchet: CAN NOT BE SUBMIITED TO CLUSTER BECAUSE OF GUROBI LICENCE ISSUE!
#    input:
#        one = "{patient}/bbc/bulk.bbc",
#        two = "{patient}/bbc/bulk.seg"
#    params:
#        hatchet = config['HATCHET'], ##add hatchet.py to config file
#        solver  = config['SOLVER'], ##add solver to config
#        prefix  = "{patient}/bbc/bulk"
#    output:
#        dir = "{patient}/results/",
#        bucn= "{patient}/results/best.bbc.ucn",
#        sucn= "{patient}/results/best.seg.ucn"
#    threads:
#        22
#    shell:
#        "python {params.hatchet} {params.solver} -i {params.prefix} -n2,8 -p 400 -v 3 -u 0.03 -r $RANDOM -j {threads} -eD 6 -eT 12 -g 0.35 -l 0.5 -x {output.dir} > hatchet.log"
