import json

with open('test-job-4.payload.json') as jfile:
    data = json.load(jfile)

# Open a file to write the CSV output
with open('local_sequencing-4.csv', 'w') as file:
    # Write the header row
    file.write("analysis_type,study_id,patient,sex,status,sample,lane,fastq_1,fastq_2,library_name,platform_unit,platform,sequencing_center,sequencing_date,platform_model,single_end,read_group_count,experiment,analysis_json\n")

    # Iterate through each sample and its read groups
    for sample in data['samples']:
        for read_group in data['read_groups']:
            # Prepare the values for each column according to the specified headers
            row = [
                data['analysisType']['name'],  # analysis_type
                data['studyId'],  # study_id
                sample['donor']['donorId'],  # patient
                sample['donor']['gender'],  # sex
                sample['specimen']['tumourNormalDesignation'],  # status
                sample['sampleId'],  # sample
                read_group['submitter_read_group_id'],  # lane
                read_group['file_r1'],  # fastq_1
                read_group['file_r2'] if read_group['is_paired_end'] else '',  # fastq_2
                read_group['library_name'],  # library_name
                read_group['platform_unit'],  # platform_unit
                data['experiment']['platform'],  # platform
                data['experiment']['sequencing_center'],  # sequencing_center
                data['experiment']['sequencing_date'],  # sequencing_date
                data['experiment']['platform_model'],  # platform_model
                not read_group['is_paired_end'],  # single_end
                data['read_group_count'],  # read_group_count
                data['experiment']['submitter_sequencing_experiment_id'],  # experiment
                'path/to/your/analysis.json'  # analysis_json (placeholder path)
            ]
            # Write the row to the CSV file
            file.write(','.join(map(str, row)) + '\n')
