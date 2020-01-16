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
# splicesite_index = config['splicesite_index']

## constructe the target if the inputs are fastqs
ALL_TRIMMED_FASTQ_1 = expand("01_trim_seq/{sample}_1.fastq", sample = SAMPLES)
ALL_TRIMMED_FASTQ_2 = expand("01_trim_seq/{sample}_2.fastq", sample = SAMPLES)
ALL_FASTQC  = expand("02_fqc/{sample}_1_fastqc.zip", sample = SAMPLES)
ALL_BAM = expand("03_bam/{sample}_Aligned.out.sam", sample = SAMPLES)
ALL_SORTED_BAM = expand("04_sortBam/{sample}.sorted.bam", sample = SAMPLES)
ALL_stringtie_gtf = expand("05_stringtie/{sample}/{sample}.stringtie.gtf", sample = SAMPLES)
ALL_bw = expand("06_bigwig/{sample}.bw", sample = SAMPLES)
ALL_QC = ["07_multiQC/multiQC_log.html"]
ball_grown = ['ballgown_gene_table.tsv']
TARGETS.extend(ALL_TRIMMED_FASTQ_1) 
TARGETS.extend(ALL_TRIMMED_FASTQ_2) 
TARGETS.extend(ALL_BAM) ##append all list to 
TARGETS.extend(ALL_SORTED_BAM)
TARGETS.extend(ALL_stringtie_gtf)
TARGETS.extend(ALL_FASTQC) ## check later
TARGETS.extend(ALL_QC)
TARGETS.extend(ball_grown)
TARGETS.extend(ALL_bw)


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
		("01_trim_seq/{sample}_1.fastq" ), ("01_trim_seq/{sample}_2.fastq")
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
	threads: 10
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

rule stringtie_FPKM_caculation:
	input: "04_sortBam/{sample}.sorted.bam"
	output: "05_stringtie/{sample}/{sample}.stringtie.gtf", directory("05_stringtie/{sample}")
	log: "00_log/{sample}_stringtie.log"
	params:
		jobname = "{sample}"
	threads: 10
	group: "mygroup"
	message: "stringtie {input} : {threads} threads"
	shell:
		"""
		# -p for paried-end, counting fragments rather reads
		{stringtie} -e -B -p {threads} -G {gtf} -o {output[0]} {input}
		"""
## add the bowgrow later

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
	bamCoverage -b {input[0]}  --binSize 100 --effectiveGenomeSize 2864785220 --skipNonCoveredRegions --normalizeUsing RPKM -p {threads}  -o {output[0]} 2> {log}
		"""

## for the strand specifc coverage
# --effectiveGenomeSize 2864785220 for GRCh37
	# 		# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
	# bamCoverage -b {input[0]}  --skipNonCoveredRegions --normalizeUsing RPKM --samFlagExclude 16  -p {threads}  -o {output[0]} 2> {log}
	# bamCoverage -b {input[0]}  --skipNonCoveredRegions --normalizeUsing RPKM --samFlagInclude 16  -p {threads}  -o {output[1]} 2>> {log}


rule multiQC:
    input :
        expand("00_log/{sample}_hisat_align", sample = SAMPLES),
        # expand("04aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES),
        expand("02_fqc/{sample}_1_fastqc.zip", sample = SAMPLES)
    output: "07_multiQC/multiQC_log.html"
    log: "00_log/multiqc.log"
    message: "multiqc for all logs"
    shell:
        """
        multiqc 02_fqc 00_log -o 07_multiQC -d -f -v -n multiQC_log 2> {log}
        """

rule ballgown:
    input: 
    	expand("05_stringtie/{sample}", sample = SAMPLES)
    output: 'ballgown_gene_table.tsv'
    shell: 
    	"""
        Rscript ballgown.R {input} {output}
    	""" 


