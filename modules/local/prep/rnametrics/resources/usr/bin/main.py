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
    Guanqiao Feng
"""

import argparse
import json
import os

qc_metrics = [
   'CORRECT_STRAND_READS',
   'PCT_RIBOSOMAL_BASES',
   'PCT_CODING_BASES',
   'PCT_UTR_BASES',
   'PCT_INTRONIC_BASES',
   'PCT_INTERGENIC_BASES',
   'PCT_MRNA_BASES',
   'PCT_USABLE_BASES',
   'PCT_CORRECT_STRAND_READS',
   'MEDIAN_CV_COVERAGE',
   'MEDIAN_5PRIME_TO_3PRIME_BIAS',
   'paired_total',
   'unpaired_total',
   'overall_alignment_rate',
   'total_reads',
   'avg_input_read_length',
   'uniquely_mapped_percent',
   'avg_mapped_read_length',
   'num_splices',
   'num_annotated_splices',
   'mismatch_rate',
   'multimapped_percent'
  ]

def parse_metrics_file(metrics_file):
    with open(metrics_file) as f:
        lines = f.readlines()
    header = lines[0].strip().split('\t')
    values = lines[1].strip().split('\t')
    # Exclude 'Sample' from parsed metrics
    return {k: v for k, v in zip(header, values) if k != 'Sample'}

def aggregate_metrics(parsed_metrics):
    aggregated = {}
    for metrics in parsed_metrics:
        for key, value in metrics.items():
            if key not in aggregated:
                if key in qc_metrics:
                    try:
                        aggregated[key] = float(value)
                    except ValueError:
                        aggregated[key] = value
    return aggregated

def parse_metrics_file_by_type(file_path):
    if 'multiqc_hisat2' in file_path:
        return 'hisat2_summary', parse_metrics_file(file_path)
    elif 'multiqc_picard_RnaSeqMetrics' in file_path:
        return 'picard_RnaSeqMetrics', parse_metrics_file(file_path)
    elif 'multiqc_star' in file_path:
        return 'star_log', parse_metrics_file(file_path)
    else:
        raise ValueError(f"Unknown file type for {file_path}")

def main():
    parser = argparse.ArgumentParser(description="Process QC metrics")
    parser.add_argument("-s", "--sample_id", required=True, help="Sample ID")
    parser.add_argument("-m", "--metrics_dir", required=True, help="Directory containing the metrics files")
    parser.add_argument("-q", "--qc_files", required=True, help="Comma-separated list of QC files")
    args = parser.parse_args()

    qc_files = args.qc_files.split(",")
    parsed_metrics = {}
    detailed_metrics = {}

    for qc_file in qc_files:
        file_path = os.path.join(args.metrics_dir, qc_file)
        metrics_type, metrics = parse_metrics_file_by_type(file_path)
        parsed_metrics[metrics_type] = metrics
        detailed_metrics[metrics_type] = metrics

    aggregated_metrics = aggregate_metrics(parsed_metrics.values())

    output_data = {
        "sample_id": args.sample_id,
        "metrics": aggregated_metrics,
        **detailed_metrics
    }

    output_file = f"{args.sample_id}.argo_metrics.json"
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=4)

if __name__ == "__main__":
    main()
