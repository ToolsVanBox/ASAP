//
// BAM SOMATIC SHORT VARIANT DISCOVERY CONTROLFREEC

// Include local modules
include { CONTROLFREEC_MAKEKARYOTYPE } from '../../../modules/local/controlfreec/makekaryotype/main'
include { CONTROLFREEC_MAKEBAFPLOT} from '../../../modules/local/controlfreec/makebafplot/main'

// Include nf-core modules
include { SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_NORMAL } from '../../../modules/nf-core/samtools/mpileup/main'
include { SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_TUMOR } from '../../../modules/nf-core/samtools/mpileup/main'
include { CONTROLFREEC_FREEC as FREEC_SOMATIC } from '../../../modules/nf-core/controlfreec/freec/main'
include { CONTROLFREEC_ASSESSSIGNIFICANCE } from '../../../modules/nf-core/controlfreec/assesssignificance/main'
include { CONTROLFREEC_MAKEGRAPH } from '../../../modules/nf-core/controlfreec/makegraph/main'
include { CONTROLFREEC_FREEC2BED } from '../../../modules/nf-core/controlfreec/freec2bed/main'
include { CONTROLFREEC_FREEC2CIRCOS } from '../../../modules/nf-core/controlfreec/freec2circos/main'

workflow BAM_SOMATIC_COPY_NUMBER_DISCOVERY_CONTROLFREEC {
  take:
    ch_bam_bai_normal  // channel: [ meta, path(normal_bam), path(normal_bai) ]
    ch_bam_bai_tumor // channel: [ meta,  path(tumor_bam), path(tumor_bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    dbsnp // path(dbsnp)
    dbsnp_tbi // path(dbsnp_tbi)
    
  main:
    ch_versions = Channel.empty()

    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    
    def chr_files = file( params.genomes[params.genome].freec_chr_files, checkIfExists: true )
    def mappability = file( params.genomes[params.genome].freec_mappability, checkIfExists: true )
    def len = file( params.genomes[params.genome].freec_len, checkIfExists: true )
    def ploidy = params.bam_somatic_copy_number_discovery_controlfreec.ploidy ?: 2

    ch_bams_normal_mpileup = ch_bam_bai_normal.map{ meta, bam, bai -> [ meta, bam, [] ] }
    ch_bams_tumor_mpileup = ch_bam_bai_tumor.map{ meta, bam, bai -> [ meta, bam, [] ] }
     
    SAMTOOLS_MPILEUP_NORMAL( ch_bams_normal_mpileup, fasta )
    SAMTOOLS_MPILEUP_TUMOR( ch_bams_tumor_mpileup, fasta )

    ch_versions = ch_versions.mix( SAMTOOLS_MPILEUP_NORMAL.out.versions )
    ch_versions = ch_versions.mix( SAMTOOLS_MPILEUP_TUMOR.out.versions )
    
    ch_freec_somatic = SAMTOOLS_MPILEUP_NORMAL.out.mpileup
      .combine( SAMTOOLS_MPILEUP_TUMOR.out.mpileup )
      .map{ normal_meta, normal_pileup, tumor_meta, tumor_pileup ->
       meta = [:]
       meta.id = "${normal_meta.id}_${tumor_meta.id}"
       meta.normal_id = normal_meta.id
       meta.tumor_id = tumor_meta.id
       meta.run_id = "${normal_meta.run_id}"
       [ meta, normal_pileup, tumor_pileup, [], [], [], [] ]
      }

    FREEC_SOMATIC(ch_freec_somatic, fasta, len, [], dbsnp, dbsnp_tbi, chr_files, mappability, [], [])
    ch_versions = ch_versions.mix( FREEC_SOMATIC.out.versions )

    cnv_files = FREEC_SOMATIC.out.CNV
      .multiMap{ meta, cnv ->
          def meta_clone_tumor = meta.clone()
          def meta_clone_normal = meta.clone()
          meta_clone_tumor.id = meta.id + "_tumor" //updating meta id so that the p.value file is named differently 
          meta_clone_normal.id = meta.id + "_normal" //updating meta id so that the p.value file is named differently 

          def tumor_file = cnv instanceof List ? cnv.find { it.toString().endsWith("gz_CNVs") } : cnv //only find if its a list, else it returns only the filename without the path
          def normal_file = cnv instanceof List ? cnv.find { it.toString().endsWith("gz_normal_CNVs") } : null //only find if its a list, else it returns only the filename without the path

          normal: normal_file == null ? [] : [meta_clone_normal,normal_file] //only fill channel if file was found, else leave it empty
          tumor: tumor_file == null ? [] : [meta_clone_tumor,tumor_file] //only fill channel if file was found, else leave it empty
      }
    
    ratio_files = FREEC_SOMATIC.out.ratio
      .multiMap{ meta, ratio ->
          def meta_clone_tumor = meta.clone()
          def meta_clone_normal = meta.clone()
          meta_clone_tumor.id = meta.id + "_tumor" //updating meta id so that the p.value file is named differently 
          meta_clone_normal.id = meta.id + "_normal" //updating meta id so that the p.value file is named differently 

          def tumor_file = ratio instanceof List ? ratio.find { it.toString().endsWith("gz_ratio.txt") } : ratio //same here as cnv
          def normal_file = ratio instanceof List ? ratio.find { it.toString().endsWith("gz_normal_ratio.txt") } : null //same here as cnv
          
          normal: normal_file == null ? [] : [meta_clone_normal,normal_file] //same here as ratio
          tumor: tumor_file == null ? [] : [meta_clone_tumor,tumor_file] //same here as ratio
      }
    
    baf_files = FREEC_SOMATIC.out.BAF
      .multiMap{ meta, baf ->
          def meta_clone_tumor = meta.clone()
          def meta_clone_normal = meta.clone()
          meta_clone_tumor.id = meta.id + "_tumor" //updating meta id so that the p.value file is named differently 
          meta_clone_normal.id = meta.id + "_normal" //updating meta id so that the p.value file is named differently 

          def tumor_file = baf instanceof List ? baf.find { it.getName().toString().startsWith("${meta.tumor_id}") } : baf //same here as cnv
          def normal_file = baf instanceof List ? baf.find { it.getName().toString().startsWith("${meta.normal_id}") } : null //same here as cnv
          
          normal: normal_file == null ? [] : [meta_clone_normal,normal_file] //same here as ratio
          tumor: tumor_file == null ? [] : [meta_clone_tumor,tumor_file] //same here as ratio
      }

    //Join the pairs
    normal_files = cnv_files.normal.join(ratio_files.normal, failOnDuplicate: true, failOnMismatch: true)
    tumor_files = cnv_files.tumor.join(ratio_files.tumor, failOnDuplicate: true, failOnMismatch: true)
    
    //Mix all the pairs into input channel
    controlfreec_assesssignificance_input = tumor_files.mix(normal_files)

    CONTROLFREEC_ASSESSSIGNIFICANCE( controlfreec_assesssignificance_input )
    ch_versions = ch_versions.mix( CONTROLFREEC_ASSESSSIGNIFICANCE.out.versions )

    normal_makegraph_files = ratio_files.normal.join( baf_files.normal, failOnDuplicate: true, failOnMismatch: true)
    tumor_makegraph_files = ratio_files.tumor.join( baf_files.tumor, failOnDuplicate: true, failOnMismatch: true)

    controlfreec_makegraph_input = tumor_makegraph_files.mix( normal_makegraph_files ).map{ meta, ratio, baf -> [ meta, ratio, baf, ploidy ]}
    CONTROLFREEC_MAKEGRAPH( controlfreec_makegraph_input )
    ch_versions = ch_versions.mix( CONTROLFREEC_MAKEGRAPH.out.versions )

    CONTROLFREEC_FREEC2BED( FREEC_SOMATIC.out.ratio )
    ch_versions = ch_versions.mix( CONTROLFREEC_FREEC2BED.out.versions )

    CONTROLFREEC_FREEC2CIRCOS( FREEC_SOMATIC.out.ratio )
    ch_versions = ch_versions.mix( CONTROLFREEC_FREEC2CIRCOS.out.versions )

    CONTROLFREEC_MAKEKARYOTYPE( FREEC_SOMATIC.out.ratio)
    ch_versions = ch_versions.mix( CONTROLFREEC_MAKEKARYOTYPE.out.versions )

    makebafplot_input = baf_files.tumor.mix( baf_files.normal )
    CONTROLFREEC_MAKEBAFPLOT( makebafplot_input )
    ch_versions = ch_versions.mix( CONTROLFREEC_MAKEBAFPLOT.out.versions )
    
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


