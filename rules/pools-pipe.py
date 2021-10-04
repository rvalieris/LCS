
# paired end fastq input
rule bwa_paired:
	input: f1='data/fastq/{x}_1.fastq.gz', f2='data/fastq/{x}_2.fastq.gz'
	output: 'outputs/pool_map/{x}.bam'
	threads: 8
	params: rgline='''"@RG\\tID:{x}\\tSM:{x}\\tLB:{x}"'''
	resources: mem_gb=3
	shell: 'bwa mem -t {threads} -R {params.rgline} {REF} {input.f1} {input.f2} | samtools view -F256 -F2048 -F4 -q60 -Sb - | samtools sort - -o {output} && samtools index {output}'

# cram input
rule cram_add_rg_tag:
	input: 'data/cram/{x}.cram'
	output: 'outputs/pool_map/{x}.bam'
	params: tag='{x}'
	shell: 'picard -Xmx4G AddOrReplaceReadGroups I={input} O={output} RGID={params.tag} RGLB={params.tag} RGSM={params.tag} RGPL=illumina RGPU=unit1 VALIDATION_STRINGENCY=SILENT && samtools index {output}'

# bam input
rule bam_add_rg_tag:
	input: 'data/bams/{x}.bam'
	output: 'outputs/pool_map/{x}.bam'
	params: tag='{x}'
	shell: 'picard -Xmx4G AddOrReplaceReadGroups I={input} O={output} RGID={params.tag} RGLB={params.tag} RGSM={params.tag} RGPL=illumina RGPU=unit1 VALIDATION_STRINGENCY=SILENT && samtools index {output}'

# single end fastq input
rule bwa_single:
	input: 'data/fastq/{x}.fastq.gz'
	output: 'outputs/pool_map/{x}.bam'
	threads: 8
	params: rgline='''"@RG\\tID:{x}\\tSM:{x}\\tLB:{x}"'''
	resources: mem_gb=3
	shell: 'bwa mem -t {threads} -R {params.rgline} {REF} {input} | samtools view -F256 -F2048 -F4 -q60 -Sb - | samtools sort - -o {output} && samtools index {output}'

rule pool_mutect:
	input: bam='outputs/pool_map/{x}.bam', sites=get_target_sites()
	output: 'outputs/pool_mutect/{x}.vcf.gz'
	resources: mem_gb=8
	log: 'outputs/pool_mutect/{x}.log'
	shell: '''gatk --java-options '-Xmx{resources.mem_gb}G' Mutect2 -R {REF} -O {output} -I {input.bam} --alleles {input.sites} --max-reads-per-alignment-start 500 -mnp-dist 0 >{log} 2>&1'''

rule pool_vcf_to_table:
	input: sites=get_target_sites(), vcfs=lambda wc: expand('outputs/pool_mutect/{x}.vcf.gz',x=get_tags(wc.z))
	output: tbl='outputs/variants_table/pool_samples_{z}.tsv'
	params: unused_path='outputs/pool_mutect_unused/'
	shell: 'mkdir -p {params.unused_path} && python scripts/pools-variants_to_table.py {output} {input.sites} {params.unused_path} {input.vcfs}'

rule decompose:
	input: markers=get_target_markers_table(), pool=get_target_pools_table()
	output: stats='outputs/decompose/{}.status'.format(config['dataset']), out='outputs/decompose/{}.out'.format(config['dataset'])
	threads: 8
	resources: mem_gb=16
	shell: 'python scripts/vg-decompose.py -m {input.markers} -p {input.pool} -o {output.out} -s {output.stats} --threads {threads}'
