process PREP_META {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biocontainers/pandas' :
        'docker.io/gfeng2023/python_pandas_requests:latest' }"

        // 'ghcr.io/icgc-argo/argo-data-submission.sanity-check:0.1.3' }"

  input:  // input, make update as needed
    path experiment_info_tsv
    val api_token
    val song_url
    val clinical_url
    val skip_duplicate_check

  output:  // output, make update as needed
    path "updated*tsv", emit: updated_experiment_info_tsv
    path "versions.yml", emit: versions

  script:
    // add and initialize variables here as needed
    def args_skip_duplicate_check = skip_duplicate_check==true ? "--force" : ""
    """
    main.py \
      -x ${experiment_info_tsv} \
      -t ${api_token} \
      -s ${song_url} \
      -c ${clinical_url} \
      ${args_skip_duplicate_check}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

