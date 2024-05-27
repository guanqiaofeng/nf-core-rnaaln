#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2021,  Ontario Institute for Cancer Research

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Linda Xiang
"""

import os
import sys
import argparse
import subprocess
import json
import re
import hashlib
import uuid
import tarfile
from datetime import date
import copy
from glob import glob
import yaml
import csv
import io
from math import log10, isnan

workflow_process_map = {
    'Pre Alignment QC': 'prealn',
    'DNA Alignment QC': 'aln'
}

tool_list = ['fastqc', 'cutadapt', 'CollectMultipleMetrics', 'CollectWgsMetrics', 'CollectHsMetrics', 'stats', 'mosdepth', 'CollectOxoGMetrics', 'contamination']

def calculate_size(file_path):
    return os.stat(file_path).st_size


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()

def get_files_info(file_to_upload, date_str, analysis_dict, process_indicator, multiqc={}):
    # sampleId = analysis_dict['samples'][0]['sampleId']
    file_info = {
        'fileSize': calculate_size(file_to_upload),
        'fileMd5sum': calculate_md5(file_to_upload),
        'fileAccess': 'controlled',
        'info': {
            'data_category': 'Quality Control Metrics',
            'data_subtypes': None,
            'files_in_tgz': []
        }
    }

    if re.match(r'.+?fastqc\.tgz$', file_to_upload):
        file_type = 'fastqc'
        file_info.update({'dataType': 'Sequencing QC'})
        file_info['info']['data_subtypes'] = ['Read Characteristics']
        file_info['info'].update({'analysis_tools': ['FastQC']})
        file_info['info'].update({'description': 'High level sequencing reads QC metrics generated by FastQC.'})

    elif re.match(r'.+?cutadapt\.tgz$', file_to_upload):
        file_type = 'cutadapt'
        file_info.update({'dataType': 'Sequencing QC'})
        file_info['info']['data_subtypes'] = ['Read Characteristics']
        file_info['info'].update({'analysis_tools': ['Cutadapt']})
        file_info['info'].update({'description': 'Perform cutadpat on sequencing reads and generate log file.'})

    elif re.match(r'.+?CollectMultipleMetrics\.tgz$', file_to_upload):
        file_type = 'picard_QualityYieldMetrics'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Read Characteristics']
        file_info['info'].update({'analysis_tools': ['Picard:CollectQualityYieldMetrics']})
        file_info['info'].update({'description': 'Picard tool to collect metrics about reads that pass quality thresholds and Illumina-specific filters.'})

    elif re.match(r'.+?CollectWgsMetrics\.tgz$', file_to_upload):
        file_type = 'picard_wgsmetrics'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Coverage']
        file_info['info'].update({'analysis_tools': ['Picard:CollectWgsMetrics']})
        file_info['info'].update({'description': 'Picard tool to collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.'})

    elif re.match(r'.+?CollectHsMetrics\.tgz$', file_to_upload):
        file_type = 'picard_hsmetrics'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Coverage']
        file_info['info'].update({'analysis_tools': ['Picard:CollectHsMetrics']})
        file_info['info'].update({'description': 'Picard tool to collect hybrid-selection (HS) metrics for targeted sequencing experiments.'})

    elif re.match(r'.+?CollectOxoGMetrics\.tgz$', file_to_upload):
        file_type = 'picard_OxoGMetrics'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Read Characteristics']
        file_info['info'].update({'analysis_tools': ['Picard:CollectOxoGMetrics']})
        file_info['info'].update({'description': 'Picard tool to collects metrics quantifying the error rate resulting from oxidative artifacts.'})

    elif re.match(r'.+?stats\.tgz$', file_to_upload):
        file_type = 'samtools_stats'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Library Quality', 'Read Characteristics']
        file_info['info'].update({'analysis_tools': ['Samtools:stats']})
        file_info['info'].update({'description': 'Samtools to collect comprehensive statistics from alignment file.'})

    elif re.match(r'.+?mosdepth\.tgz$', file_to_upload):
        file_type = 'mosdepth'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Coverage']
        file_info['info'].update({'analysis_tools': ['Mosdepth']})
        file_info['info'].update({'description': 'Mosdepth performs fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing.'})

    elif re.match(r'.+?contamination\.tgz$', file_to_upload):
        file_type = 'gatk_contamination'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Contamination']
        file_info['info'].update({'analysis_tools': ['GATK:CalculateContamination']})
        file_info['info'].update({'description': 'GATK4 tool to calculate the fraction of reads coming from cross-sample contamination.'})

    elif re.match(r'.+?RnaSeqMetrics.', file_to_upload): # to be more specific in the future to match other matching styles
        file_type = 'picard_RnaSeqMetrics'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Library Quality', 'Read Characteristics']
        file_info['info'].update({'analysis_tools': ['Picard:RnaSeqMetrics']})
        file_info['info'].update({'description': 'Picard tool to collect metrics describing the distribution of the bases within the transcripts.'})

    elif re.match(r'.+?hisat2.', file_to_upload): # to be more specific in the future to match other matching styles
        file_type = 'hisat2_summary'
        file_info.update({'dataType': 'Aligned Reads QC'})
        file_info['info']['data_subtypes'] = ['Library Quality', 'Read Characteristics']
        file_info['info'].update({'analysis_tools': ['Hisat2:summary']})
        file_info['info'].update({'description': 'Hisat2 alignment alignment summary file to collect metrics describing mapping rates.'})

    else:
        sys.exit('Error: unknown QC metrics file: %s' % file_to_upload)

    # retrieve qc metrics from multiqc_data
    metric_info = multiqc.get(file_type, [])
    metric_info_updated = []
    for metric_item in metric_info:
      metric_info_updated.append(metric_item)
    file_info['info'].update({'metrics': metric_info_updated})

    # file naming patterns:
    #   pattern:  <argo_study_id>.<argo_donor_id>.<argo_sample_id>.<experiment_strategy>.<date>.<process_indicator>.<file_type>.<file_ext>
    #   process_indicator: pre-alignment, alignment(aligner), post-alignment(caller)
    #   example: TEST-PR.DO250183.SA610229.rna-seq.20200319.star.genome_aln.cram
    new_fname = '.'.join([
        analysis_dict['studyId'],
        analysis_dict['samples'][0]['donor']['donorId'],
        analysis_dict['samples'][0]['sampleId'],
        analysis_dict['experiment']['experimental_strategy'].lower() if analysis_dict['experiment'].get('experimental_strategy') else analysis_dict['experiment']['library_strategy'],
        date_str,
        process_indicator,
        file_type,
        'tgz'
      ])

    file_info['fileName'] = new_fname
    file_info['fileType'] = new_fname.split('.')[-1].upper()

    with tarfile.open(file_to_upload, 'r') as tar:
      for member in tar.getmembers():
        file_info['info']['files_in_tgz'].append(member.name)

    new_dir = 'out'
    try:
      os.mkdir(new_dir)
    except FileExistsError:
      pass

    dst = os.path.join(os.getcwd(), new_dir, new_fname)
    os.symlink(os.path.abspath(file_to_upload), dst)

    return file_info

def get_basename(metadata):
    study_id = metadata['studyId']
    donor_id = metadata['samples'][0]['donor']['donorId']
    sample_id = metadata['samples'][0]['sampleId']

    if not sample_id or not donor_id or not study_id:
      sys.exit('Error: missing study/donor/sample ID in the provided metadata')

    return ".".join([study_id, donor_id, sample_id])

def get_sample_info(sample_list):
    samples = copy.deepcopy(sample_list)
    for sample in samples:
      for item in ['info', 'sampleId', 'specimenId', 'donorId', 'studyId']:
        sample.pop(item, None)
        sample['specimen'].pop(item, None)
        sample['donor'].pop(item, None)

    return samples

def prepare_tarball(sampleId, qc_files, tool_list):

    tgz_dir = 'tarball'
    try:
      os.mkdir(tgz_dir)
    except FileExistsError:
      pass

    files_to_tar = {}
    for tool in tool_list:
      if not tool in files_to_tar: files_to_tar[tool] = []
      for f in sorted(qc_files):
        if tool in f:
          files_to_tar[tool].append(f)

    for tool in tool_list:
      if not files_to_tar[tool]: continue
      tarfile_name = f"{tgz_dir}/{sampleId}.{tool}.tgz"
      with tarfile.open(tarfile_name, "w:gz", dereference=True) as tar:
        for f in files_to_tar[tool]:
          tar.add(f, arcname=os.path.basename(f))

def main():
    """
    Python implementation of tool: payload-gen-qc
    """

    parser = argparse.ArgumentParser(description='Tool: payload-gen-qc')
    parser.add_argument("-a", "--metatada-analysis", dest="metadata_analysis", required=True,
                        help="Input metadata analysis", type=str)
    parser.add_argument("-f", "--files_to_upload", dest="files_to_upload", type=str, required=True,
                        nargs="+", help="All files to upload")
    parser.add_argument("-w", "--wf-name", dest="wf_name", required=True, help="Workflow name")
    parser.add_argument("-s", "--wf-session", dest="wf_session", required=True, help="workflow session ID")
    parser.add_argument("-v", "--wf-version", dest="wf_version", required=True, help="Workflow version")
    parser.add_argument("-p", "--pipeline_yml", dest="pipeline_yml", required=False, help="Pipeline info in yaml")
    parser.add_argument("-m", "--multiqc", dest="multiqc", required=False, help="multiqc json file")

    args = parser.parse_args()

    with open(args.metadata_analysis, 'r') as f:
      analysis_dict = json.load(f)

    pipeline_info = {}
    if args.pipeline_yml:
      with open(args.pipeline_yml, 'r') as f:
        pipeline_info = yaml.safe_load(f)

    # get tool_specific & aggregated metrics from multiqc
    mqc_stats = {}
    if args.multiqc:
      with open(args.multiqc, 'r') as f:
        mqc_stats = json.load(f)

    payload = {
        'analysisType': {
            'name': 'qc_metrics'
        },
        'studyId': analysis_dict.get('studyId'),
        'info': {},
        'workflow': {
            'workflow_name': args.wf_name,
            'workflow_version': args.wf_version,
            'session_id': args.wf_session,
            'inputs': [
                {
                    'analysis_type': analysis_dict['analysisType']['name'],
                    'input_analysis_id': analysis_dict.get('analysisId')
                }
            ],
            'pipeline_info': pipeline_info,
            'metrics': mqc_stats.get('metrics', None)
        },
        'files': [],
        'experiment': analysis_dict.get('experiment'),
        'samples': get_sample_info(analysis_dict.get('samples'))
    }
    if analysis_dict.get('workflow'):
      if analysis_dict['workflow'].get('genome_build'):
         payload['workflow']['genome_build'] = analysis_dict['workflow'].get('genome_build')
      if analysis_dict['workflow'].get('genome_annotation'):
         payload['workflow']['genome_annotation'] = analysis_dict['workflow'].get('genome_annotation')

    # pass `info` dict from seq_experiment payload to new payload
    if 'info' in analysis_dict and isinstance(analysis_dict['info'], dict):
      payload['info'] = analysis_dict['info']
    else:
      payload.pop('info')

    if 'library_strategy' in payload['experiment']:
      experimental_strategy = payload['experiment'].pop('library_strategy')
      payload['experiment']['experimental_strategy'] = experimental_strategy

    new_dir = 'out'
    try:
        os.mkdir(new_dir)
    except FileExistsError:
        pass

    # generate date string
    date_str = date.today().strftime("%Y%m%d")

    # prepare tarball to include all QC files generated by one tool
    prepare_tarball(analysis_dict['samples'][0]['sampleId'], args.files_to_upload, tool_list)

    process_indicator = workflow_process_map.get(args.wf_name)
    for f in sorted(glob('tarball/*.tgz')):
      file_info = get_files_info(f, date_str, analysis_dict, process_indicator, mqc_stats)
      payload['files'].append(file_info)

    with open("%s.%s.payload.json" % (str(uuid.uuid4()), args.wf_name.replace(" ","_")), 'w') as f:
        f.write(json.dumps(payload, indent=2))



if __name__ == "__main__":
    main()

