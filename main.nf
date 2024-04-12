
nextflow.enable.dsl=2

include { ASAP } from './workflows/asap.nf'

workflow {
    myFile = file("${projectDir}/vanboxtellab_ascii_art.txt")
    allLines = myFile.readLines()
    for( line : allLines ) {
        println line
    }

    myFile = file("${projectDir}/asap_ascii_art.txt")
    allLines = myFile.readLines()
    for( line : allLines ) {
        println line
    }

    // ASAP()
}