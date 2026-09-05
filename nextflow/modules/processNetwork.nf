process processNetwork {

    input:
        path networkfile

    output:
        path "wprime_${networkfile}", emit: network_w_prime
        path "nodes_${networkfile}", emit: nodes

    script:
    """
    """
}