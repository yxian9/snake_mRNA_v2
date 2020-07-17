shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())

gtf = config["MYGTF"]
STARINDEX = config['STARINDEX']
hisat = config["hisat"]
stringtie = config['stringtie']
TARGETS = []
NGmerge = config['NGmerge'] ## remove the adaptor by the overlap , if the adaptor contaimination happen, use the cutadaptor


## constructe the target if the inputs are fastqs
ALL_TRIMMED_FASTQ_1 = expand("01_trim_seq/{sample}_1.fastq", sample = SAMPLES)
ALL_TRIMMED_FASTQ_2 = expand("01_trim_seq/{sample}_2.fastq", sample = SAMPLES)
ALL_FASTQC  = expand("02_fqc/{sample}_1_fastqc.zip", sample = SAMPLES)
ALL_SAM = expand("03_bam/{sample}_Aligned.out.sam", sample = SAMPLES)
ALL_SORTED_BAM = expand("04_sortBam/{sample}.sorted.bam", sample = SAMPLES)
# ALL_stringtie_gtf = expand("05_stringtie/{sample}/{sample}.stringtie.gtf", sample = SAMPLES)
ALL_bw = expand("06_bigwig/{sample}.bw", sample = SAMPLES)
ALL_QC = ["07_multiQC/multiQC_log.html"]
# ball_grown = ['ballgown_gene_table.tsv']
ALL_feature_count = ("07_featurecount/featureCount.txt") ## using feature count to generate the table

# TARGETS.extend(ALL_TRIMMED_FASTQ_1) 
# TARGETS.extend(ALL_TRIMMED_FASTQ_2) 
# TARGETS.extend(ALL_BAM) ##append all list to 
# TARGETS.extend(ALL_SORTED_BAM)
TARGETS.extend(ALL_FASTQC) ## check later
# TARGETS.extend(ALL_QC)
TARGETS.extend(ALL_bw) ##
TARGETS.extend(ALL_feature_count)




localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

rule all:
	input: TARGETS


rule trim_fastqs: ## merge fastq
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		temp("01_trim_seq/{sample}_1.fastq" ), temp("01_trim_seq/{sample}_2.fastq")
	log: "00_log/{sample}_trim_adapter.log"
	params:
		jobname = "{sample}"
	threads : 8
	# group: "mygroup"
	message: "trim fastqs {input}: {threads} threads"
	shell:
		"""
		NGmerge  -a -y  -n {threads} -1 {input[0]} -2 {input[1]}  -o 01_trim_seq/{params.jobname} 2> {log} 
		"""

rule fastqc:
	input:  "01_trim_seq/{sample}_1.fastq" , "01_trim_seq/{sample}_2.fastq"
	output: "02_fqc/{sample}_1_fastqc.zip" 
	log:    "00_log/{sample}_fastqc"
	# group: "mygroup"
	params : jobname = "{sample}"
	message: "fastqc {input}: {threads}"
	shell:
	    """
	    module load fastqc
	    fastqc -o 02_fqc -f fastq --noextract {input}  2> {log}
	    """



rule hisat_mapping:
	input: 
		"01_trim_seq/{sample}_1.fastq", 
		"01_trim_seq/{sample}_2.fastq"
	output: temp("03_bam/{sample}_Aligned.out.sam")
	log: "00_log/{sample}_hisat_align"
	params: 
		jobname = "{sample}"
	threads: 12
	# group: "mygroup"
	message: "aligning {input} using hisat: {threads} threads"
	shell:
		"""
		{hisat} -p {threads} \
		--dta \
		-x {STARINDEX} \
		-1 {input[0]} \
		-2 {input[1]} \
		-S {output} \
		&> {log}
		"""
		#--rna-strandness R ## for stand strandness
		#with the genome_tran no need for the slicesite
		#--known-splicesite-infile {splicesite_index} \



		
rule featureCount_fq:
    input: ALL_SAM
    output: "07_featurecount/featureCount.txt"
    log: "00log/featureCount.log"
    threads: 12
    message: "feature-count {input} : {threads} threads"
    shell:
        """
        # -p for paried-end, counting fragments rather reads
        featureCounts -T {threads} -p  -Q 10 -t exon -g gene_id --extraAttributes gene_name -a {gtf} -o {output} 03_bam/*.sam 2> {log}
        """

	
rule sortBam:
	input: "03_bam/{sample}_Aligned.out.sam"
	output: "04_sortBam/{sample}.sorted.bam"
	log: "00_log/{sample}_sortbam.log"
	params:
		jobname = "{sample}"
	threads: 12
	# group: "mygroup"
	message: "sorting {input} : {threads} threads"
	shell:
		"""
		module load samtools
		samtools sort  -@ {threads} -T {output}.tmp -o {output} {input} 2> {log}
		"""

rule index_bam:
    input:  "04_sortBam/{sample}.sorted.bam"
    output: "04_sortBam/{sample}.sorted.bam.bai"
    log:    "00_log/{sample}.index_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "index_bam {input}: {threads} threads"
    shell:
        """
        module load samtools
        samtools index {input} 2> {log}
        """


rule make_bigwigs: ## included if need the coverage depth 
	input : "04_sortBam/{sample}.sorted.bam", "04_sortBam/{sample}.sorted.bam.bai"
	# input : expand("05_sortBam/{sample}.sorted.bam", sample = SAMPLES)
	# output: "07_bigwig/{sample}_forward.bw", "07_bigwig/{sample}_reverse.bw"
	output: "06_bigwig/{sample}.bw"
	log: "00_log/{sample}.makebw"
	threads: 10
	params: jobname = "{sample}"
	message: "making bigwig for {input} : {threads} threads"
	shell:
		"""
	# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
	bamCoverage -b {input[0]}  --binSize 100  --skipNonCoveredRegions --normalizeUsing RPKM -p {threads}  -o {output[0]} 2> {log}
		"""

## for the strand specifc coverage
# --effectiveGenomeSize 2864785220 for GRCh37
	# 		# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
	# bamCoverage -b {input[0]}  --skipNonCoveredRegions --normalizeUsing RPKM --samFlagExclude 16  -p {threads}  -o {output[0]} 2> {log}
	# bamCoverage -b {input[0]}  --skipNonCoveredRegions --normalizeUsing RPKM --samFlagInclude 16  -p {threads}  -o {output[1]} 2>> {log}


# rule multiQC:
#     input :
# 	ALL_FASTQC, ALL_SORTED_BAM
#     output: "07_multiQC/multiQC_log.html"
#     log: "00_log/multiqc.log"
#     message: "multiqc for all logs"
#     shell:
#         """
#         multiqc 02_fqc 00_log -o 07_multiQC -d -f -v -n multiQC_log 2> {log}
#         """



