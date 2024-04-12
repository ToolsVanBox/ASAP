//
// BAM GERMLINE SHORT VARIANT DISCOVERY CONTROLFREEC

// Include local modules
include { CONTROLFREEC_MAKEKARYOTYPE } from '../../../modules/local/controlfreec/makekaryotype/main'
include { CONTROLFREEC_MAKEBAFPLOT} from '../../../modules/local/controlfreec/makebafplot/main'

// Include nf-core modules
include { SAMTOOLS_MPILEUP } from '../../../modules/nf-core/samtools/mpileup/main'
include { CONTROLFREEC_FREEC as FREEC_GERMLINE } from '../../../modules/nf-core/controlfreec/freec/main'
include { CONTROLFREEC_ASSESSSIGNIFICANCE } from '../../../modules/nf-core/controlfreec/assesssignificance/main'
include { CONTROLFREEC_MAKEGRAPH } from '../../../modules/nf-core/controlfreec/makegraph/main'
include { CONTROLFREEC_FREEC2BED } from '../../../modules/nf-core/controlfreec/freec2bed/main'
include { CONTROLFREEC_FREEC2CIRCOS } from '../../../modules/nf-core/controlfreec/freec2circos/main'

workflow BAM_GERMLINE_COPY_NUMBER_DISCOVERY_CONTROLFREEC {
  take:
    ch_bam_bai  // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    dbsnp // path(dbsnp)
    dbsnp_tbi // path(dbsnp_tbi)
    
  main:
    ch_versions = Channel.empty()

    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    
    def chr_files = file( params.genomes[params.genome].freec_chr_files, checkIfExists: true )
    def mappability = file( params.genomes[params.genome].freec_mappability, checkIfExists: true )
    def len = file( params.genomes[params.genome].freec_len, checkIfExists: true )
    def ploidy = params.bam_germline_copy_number_discovery_controlfreec.ploidy ?: 2

    ch_bams_mpileup = ch_bam_bai.map{ meta, bam, bai -> [ meta, bam, [] ] }
     
    SAMTOOLS_MPILEUP( ch_bams_mpileup, fasta )
    ch_versions = ch_versions.mix( SAMTOOLS_MPILEUP.out.versions )

    controlfreec_input = SAMTOOLS_MPILEUP.out.mpileup.map{ meta, pileup -> [ meta, [], pileup, [], [], [], [] ] }
    FREEC_GERMLINE(controlfreec_input, fasta, len, [], dbsnp, dbsnp_tbi, chr_files, mappability, [], [])
    ch_versions = ch_versions.mix( FREEC_GERMLINE.out.versions )

    controlfreec_assesssignificance_input = FREEC_GERMLINE.out.CNV.join( FREEC_GERMLINE.out.ratio )
    CONTROLFREEC_ASSESSSIGNIFICANCE( controlfreec_assesssignificance_input )
    ch_versions = ch_versions.mix( CONTROLFREEC_ASSESSSIGNIFICANCE.out.versions )
    
    controlfreec_makegraph_input = FREEC_GERMLINE.out.ratio.join( FREEC_GERMLINE.out.BAF ).map{ meta, ratio, baf -> [ meta, ratio, baf, ploidy ] }
    CONTROLFREEC_MAKEGRAPH( controlfreec_makegraph_input )
    ch_versions = ch_versions.mix( CONTROLFREEC_MAKEGRAPH.out.versions )

    CONTROLFREEC_FREEC2BED( FREEC_GERMLINE.out.ratio )
    ch_versions = ch_versions.mix( CONTROLFREEC_FREEC2BED.out.versions )

    CONTROLFREEC_FREEC2CIRCOS( FREEC_GERMLINE.out.ratio )
    ch_versions = ch_versions.mix( CONTROLFREEC_FREEC2CIRCOS.out.versions )

    CONTROLFREEC_MAKEKARYOTYPE( FREEC_GERMLINE.out.ratio )
    ch_versions = ch_versions.mix( CONTROLFREEC_MAKEKARYOTYPE.out.versions )

    CONTROLFREEC_MAKEBAFPLOT( FREEC_GERMLINE.out.BAF )
    ch_versions = ch_versions.mix( CONTROLFREEC_MAKEBAFPLOT.out.versions )
    
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


