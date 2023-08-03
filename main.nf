/*
========================================================================================
Variant-Calling Nextflow Workflow
========================================================================================
Github   : FHisam27, ene220, NathanielT99
Contact  : fatima.hisam@utdallas.edu ene220000@utdallas.edu, NathanielT99
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Pipeline Input parameters

params.outdir = 'results'
// (Done for the first read) TODO Find the urls for these files https://github.com/sateeshperi/nextflow_varcal/tree/master/data
params.genome = "/scratch/applied-genomics/nextflow_varcal/data/ref_genome/ecoli_rel606.fasta"
params.reads = "/scratch/applied-genomics/nextflow_varcal/data/trimmed_fastq/SRR2584863_{1,2}.trim.fastq.gz"

println """\
        V A R I A N T-C A L L I N G - N F   P I P E L I N E
        ===================================
        genome       : ${params.genome}
        reads        : ${params.reads}
        outdir       : ${params.outdir}
        """
        .stripIndent()

/*
========================================================================================
Create Channels
========================================================================================
*/

ref_ch = Channel.fromPath( params.genome, checkIfExists: true )
reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )

/*
========================================================================================
MAIN Workflow
========================================================================================
*/

workflow {

    FASTQC( reads_ch )
    BWA_INDEX( ref_ch )
    BWA_ALIGN( BWA_INDEX.out.bwa_index.combine(reads_ch) ) // https://www.nextflow.io/docs/latest/process.html#understand-how-multiple-input-channels-work
    SAMTOOLS_SORT( BWA_ALIGN.out.aligned_bam )
    SAMTOOLS_INDEX( SAMTOOLS_SORT.out.sorted_bam )
    BCFTOOLS_MPILEUP( SAMTOOLS_SORT.out.sorted_bam )
    BCFTOOLS_CALL( BCFTOOLS_MPILEUP.out.bcf )
    VCFUTILS( BCFTOOLS_CALL.out.vcf )

}

/*
========================================================================================
Processes
========================================================================================
*/

/*
 * Align reads to reference genome & create BAM file.
 */
process FASTQC {
    tag{"${reads}"}
    label 'process_low'
    conda 'fastqc'

    publishDir("${params.outdir}/fastqc_trim", mode: 'copy')

    input:
    tuple val( sample_id ), path( reads )

    output:
    path( "*_fastqc*" ), emit: fastqc_out

    script:
    """
    fastqc ${reads}
    """
}

/*
 * Index the reference genome for use by bwa and samtools.
 */
process BWA_INDEX {
    tag{"${genome}"}
    label 'process_low'
    conda 'bwa samtools'

    publishDir("${params.outdir}/bwa_index", mode: 'copy')

    input:
    path genome

    output:
    tuple path( genome ), path( "*" ), emit: bwa_index

    script:
    """
    bwa index ${genome}
    samtools faidx ${genome}
    """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BWA_ALIGN {
    tag{"${sample_id}"}
    label 'process_medium'
    conda 'bwa samtools'

    publishDir("${params.outdir}/bwa_align", mode: 'copy')

    input:
    tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.bam" ), emit: aligned_bam

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem \$INDEX ${reads} > ${sample_id}.aligned.sam
    samtools view -S -b ${sample_id}.aligned.sam > ${sample_id}.aligned.bam
    """
}

/*
 * Convert the format of the alignment to sorted BAM.
 */
process SAMTOOLS_SORT {
    tag{"${sample_id}"}
    label 'process_low'
    conda 'samtools'

    publishDir("${params.outdir}/bam_align", mode: 'copy')

    input:
    tuple val( sample_id ), path( bam )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam" ), emit: sorted_bam

    script:
    """
    samtools sort -o "${sample_id}.aligned.sorted.bam" ${bam}
    """
}

/*
 * Index the BAM file for visualization purpose
 */
process SAMTOOLS_INDEX {
    tag{"${sample_id}"}
    label 'process_low'
    conda 'samtools'

    publishDir("${params.outdir}/bam_align", mode: 'copy')

    input:
    tuple val( sample_id ), path( bam )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam.bai" ), emit: sorted_bam_index

    script:
    """
    samtools index ${bam}
    """
}

/*
 * Calculate the read coverage of positions in the genome.
*/
process BCFTOOLS_MPILEUP {
    tag{"${sample_id}"}
    label 'process_high'
    conda "bioconda::bcftools=1.17"

    publishDir("${params.outdir}/bcftools_mpileup", mode: 'copy')

    input:
    tuple val( sample_id ), path( bam )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam.bcf" ), emit: bcf

    script:
    """
    bcftools mpileup -O b -o ${sample_id}.aligned.sorted.bam.bcf -f ${params.genome} ${bam}
    """
}

/*
 * Detect the single nucleotide variants (SNVs).
*/
process BCFTOOLS_CALL {
    tag{"${sample_id}"}
    label 'process_high'
    conda "bioconda::bcftools=1.17"

    publishDir("${params.outdir}/bcftools_call", mode: 'copy')

    input:
    tuple val( sample_id ), path( bcf )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam.vcf" ), emit: vcf

    script:
    """
    bcftools call -vmO v -o ${sample_id}.aligned.sorted.bam.vcf ${bcf}
    """
}

//* Filter and report the SNVs in VCF (variant calling format).

process VCFUTILS {
    tag{"VCFUTILS ${sample_id}"}
    label 'process_high'
    conda 'vcftools'

    publishDir("${params.outdir}/vcfutils", mode: 'copy')

    input:
    tuple val( sample_id ), path( vcf )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam.vcf" ), emit: vcf

    script:
    """
    vcfutils.pl varFilter -Q 20 -d 10 -D 1000 ${vcf} > ${sample_id}.aligned.sorted.bam.vcf
    """
    
}


/*
========================================================================================
Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

    println ( workflow.success ? """
            Pipeline execution summary
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            """ : """
            Failed: ${workflow.errorReport}
            exit status : ${workflow.exitStatus}
            """
    )
}

/*
========================================================================================
THE END
========================================================================================
*/
