process SENTIEON_BWAINDEX {
    tag "$fasta"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    path fasta

    output:
    path "bwa_index"   , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    source sentieon_init.sh SENTIEON_LICENSE_BASE64

    mkdir bwa_index

    sentieon bwa index \\
        $args \\
        -p bwa_index/${fasta.baseName} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    """
    source sentieon_init.sh SENTIEON_LICENSE_BASE64

    mkdir bwa_index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
