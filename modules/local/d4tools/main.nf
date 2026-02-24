process D4TOOLS_CREATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/d4tools:0.3.11--h3ab6199_2':
        'biocontainers/d4tools:0.3.11--h3ab6199_2' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.d4"), emit: d4
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    d4tools create \\
        $args \\
        ${bam} \\
        ${prefix}.d4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        d4tools: \$(d4tools --version | head -n 1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '-Az'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.d4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        d4tools: "stub-version"
    END_VERSIONS
    """
}
