/*
========================================================================================
Variant-Calling Nextflow Workflow
========================================================================================
<<<<<<< HEAD
<<<<<<< HEAD
   Github   : // ene220
   Contact  : // ene220000@utdallas.edu
=======
   Github   : // FHisam27
   Contact  : // fatima.hisam@utdallas.edu
>>>>>>> e3e98307d88b511a237342b44df62e84382fb189
=======

   Github   : // FHisam27, NathanielT99
   Contact  : // fatima.hisam@utdallas.edu, NathanielT99

>>>>>>> e51e668d66dfc30620d7ee4010033e0dce802f11
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Pipeline Input parameters

params.outdir = 'results'
<<<<<<< HEAD
// (Done for the first read) TODO Find the urls for these files https://github.com/sateeshperi/nextflow_varcal/tree/master/data
params.genome = "/scratch/applied-genomics/nextflow_varcal/data/ref_genome/ecoli_rel606.fasta"
params.reads = "/scratch/applied-genomics/nextflow_varcal/data/trimmed_fastq/SRR2384863_{1,2}.trim.fastq.gz"
=======
// Find the urls for these files https://github.com/sateeshperi/nextflow_varcal/tree/master/data
params.genome = "https://github.com/sateeshperi/nextflow_varcal/raw/master/data/ref_genome/ecoli_rel606.fasta"
params.reads = "https://github.com/sateeshperi/nextflow_varcal/raw/master/data/trimmed_fastq/SRR2584863_1.trim.fastq.gz"

/* 
Something like the below code will replace params.reads once we confirm working for single file
params.reads = /scratch/applied-genomics/nextflow_varcal/data/data/trimmed_fastq/SRR*2584863{1,2}.trim.fastq.gz
*/

>>>>>>> e3e98307d88b511a237342b44df62e84382fb189

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
    // TODO Enter the rest of the processes for variant calling based on the bash script below

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
    tag{"FASTQC ${reads}"}
    label 'process_low'
    conda 'varcal'

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
    tag{"BWA_INDEX ${genome}"}
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
    tag{"BWA_ALIGN ${sample_id}"}
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
    tag{"SAMTOOLS_SORT ${sample_id}"}
    label 'process_low'
<<<<<<< HEAD
    conda 'samtools'
=======
    conda 'bioconda::samtools'
>>>>>>> e51e668d66dfc30620d7ee4010033e0dce802f11

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
<<<<<<< HEAD
=======
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
>>>>>>> e3e98307d88b511a237342b44df62e84382fb189
}

/*
 * Calculate the read coverage of positions in the genome.
 */
process BCFTOOLS_MPILEUP {
    tag{"BCFTOOLS_MPILEUP ${sample_id}"}
    label 'process_high'
    conda 'bcftools'

    publishDir("${params.outdir}/bcftools_mpileup", mode: 'copy')

    input:
    tuple val( sample_id ), path( bam )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam.bcf" ), emit: bcf

    script:
    """
    bcftools mpileup -O b -o ${sample_id}.aligned.sorted.bam.bcf ${bam}
    """
}

/*
 * Detect the single nucleotide variants (SNVs).
 */
process BCFTOOLS_CALL {
    tag{"BCFTOOLS_CALL ${sample_id}"}
    label 'process_high'
    conda 'bcftools'





}

/*
 * Filter and report the SNVs in VCF (variant calling format).
 */
process VCFUTILS {
    // TODO
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
