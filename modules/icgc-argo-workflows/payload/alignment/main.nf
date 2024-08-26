process PAYLOAD_ALIGNMENT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::multiqc=1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:  // input, make update as needed
      tuple val(meta), path(files_to_upload), path(metadata_analysis)
      path pipeline_yml
    //   val genome_build
    //   val genome_annotation

    output:  // output, make update as needed
      tuple val(meta), path("*.payload.json"), path("out/*"), emit: payload_files
      path "versions.yml", emit: versions

    script:
      // add and initialize variables here as needed
      def arg_pipeline_yml = pipeline_yml.name != 'NO_FILE' ? "-p $pipeline_yml" : ''
      """
      main.py \
        -f ${files_to_upload} \
        -a ${metadata_analysis} \
        -w "${workflow.manifest.name}" \
        -r ${workflow.runName} \
        -s "${workflow.sessionId}" \
        -v "${workflow.manifest.version}" \
        -c "${meta.read_groups_count}" \
        -b "${meta.genome_build}" \
        -n "${meta.genome_annotation}" \
        $arg_pipeline_yml

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$(python --version | sed 's/Python //g')
      END_VERSIONS
      """
  }
