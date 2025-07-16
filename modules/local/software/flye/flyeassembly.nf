process FLYE_ASSEMBLE_MOD {
    
    label 'big_task'
    label 'flye'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9.5--py39hdf45acc_1':
        null }"
    
    input:
    tuple val(name), path(reads)
    
    output:
    tuple val(name), path("flye_${name}"), emit: flyeassembled
    tuple val(name), path("flye_${name}/${name}.fasta"), emit: fasta
    
    script:
    def arg_read_error = params.read_error
    def arg_genome = params.genome_size
    def arg_asm = params.asm_coverage
    
    """
    flye --nano-hq ${reads} \
    --read-error ${arg_read_error} \
    --genome-size ${arg_genome} \
    --asm-coverage ${arg_asm} \
    --threads ${task.cpus} \
    --out-dir flye_${name}

    mv flye_${name}/assembly.fasta flye_${name}/${name}.fasta

    cat <<-VERSIONS > versions.yml
    "${task.process}":
         flye: \$(flye --version | cut -d" " -2 2>&1)
    VERSIONS
    """
}