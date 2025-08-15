//
// HLA type calling from bam files 
//

// Include local subworkflows
include { BAM_HLA_TYPE_CALLING_LILAC } from '../../../subworkflows/local/bam_hla_type_calling_lilac/main'      



workflow BAM_HLA_TYPE_CALLING {
    take: 
        ch_normal_bam       // channel: [mandatory] [ meta, bam, bai ]
        ch_tumor_bam        // channel: [mandatory] [ meta, bam, bai ]
        ch_tumor_rna_bam    // channel: [mandatory] [ meta, bam, bai ]
        ch_cnv_dir          // channel: [mandatory] [ meta, cnv_dir ]
        ch_somatic_vcf      // channel: [mandatory] [ meta, path(vcf) ]
        ch_fasta            // channel : [ meta, path(fasta) ]
        ch_fai              // channel: [ meta, path(fai) ]

    main: 
        ch_versions = Channel.empty()

        for ( tool in params.bam_hla_type_calling.tool ) {
            tool = tool.toLowerCase()      
            known_tool = false

            // HMFtools Lilac
            if ( tool == "lilac" ) {

                // Retrieve the resource variables  --> Check if this is still needed 
                def lilac_dir = file( params.genomes[params.genome].lilac_dir, checkIfExists: true )
                def lilac_slice = file( params.genomes[params.genome].lilac_slice, checkIfExists: true )
                def lilac_chr_version = params.genomes[params.genome].lilac_chr_version
                ch_input = Channel.empty()

                // Get a combined meta from your input bam files 
                ch_input = ch_normal_bam
                    .combine(ch_tumor_bam)
                    .map{ meta, normal_bam, normal_bai, meta2, tumor_bam, tumor_bai -> 
                        meta = meta - meta.subMap("id","sample_type")
                        meta = meta + [ id: meta.sample+"-"+meta2.sample ]
                        meta = meta - meta.subMap("sample")
                        [ meta ]
                    }

                // Run Lilac for HLA-type calling 
                BAM_HLA_TYPE_CALLING_LILAC(ch_input, ch_normal_bam, ch_tumor_bam, ch_tumor_rna_bam, ch_cnv_dir, ch_somatic_vcf,  ch_fasta, ch_fai, lilac_slice)
                ch_versions = ch_versions.mix( BAM_HLA_TYPE_CALLING_LILAC.out.versions )

                known_tool = true
            } 

            // Skip if tool is unknown
            if ( ! known_tool ) {
                println ("WARNING: Skip ${tool}, because it's not known as a hla type calling tool, this tool is not build in (yet).")
            }
        }

    emit: 
        versions  = ch_versions // channel: [ versions.yml ]
}


