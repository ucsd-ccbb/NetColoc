#!/usr/bin/env nextflow

// Include modules

/*
* NetColoc Parameters
*/

params {
    input_genes: Path
    network: Path
    batch: String
}

/*
* Pipeline
*/

workflow {

    main:
    // Import and process network
    network_ch = channel.fromPath()
    processNetwork(network_ch)
    
    // Import and process input genes
    input_gene_ch = channel.fromPath('all_input_gene_lists')

    // how to do this efficiently??
    getInputGenes(input_gene_ch, processNetwork.out.nodes)

    // perform network propagation
    networkPropagation(getInputGenes.out, processNetwork.out.nodes)

    // perform colocalization
    networkColocalization(networkPropagation.out.collect())

    // Get statistics
    colocStatistics(networkColocalization.out)


    publish:
    network_w_prime = processNetwork.out.network_w_prime
    nodes = processNetwork.out.nodes
    genes = getInputGenes.out

}


output {
    network_w_prime {
        path { "${params.batch}/inputs/" }
        mode 'copy'
    }
    nodes {
        path { "${params.batch}/inputs/" }
        mode 'copy'
    }
    genes {
        path { "${params.batch}/inputs/" }
        mode 'copy'
    }
}