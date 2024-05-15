// Import your local module
include { NOVEL_SPLICE_MERGE } from './main.nf'

// Define a channel with test data
// txt_files_ch = Channel.from([ \
//     [analysis_type: 'sequencing_experiment', id: 'S01L1', study_id: 'TCRB-CA', patient: 'DO262498', sex: 'XY', status: 1, sample: 'SA622799', read_group: '@RG\tID:S01L1\tSM:SA622799\tLB:sim01\tPU:sim01\tPI:250\tCN:SIM\tPL:ILLUMINA\tPM:Polyester\tDT:2021-04-21\tDS:RNA-Seq|TCRB-CA|SP222784|DO262498|Primary tumour|Tumour', data_type: 'fastq', numLanes: 2, experiment: 'RNA-Seq', single_end: false, library_strandedness: 'UNSTRANDED', date: '2021-04-21'], \
//     "./TCRB-CA.DO262498.SA622799.S01L1.hisat2.novel_splicesites.txt"], \
//     [analysis_type: 'sequencing_experiment', id: 'S01L2', study_id: 'TCRB-CA', patient: 'DO262498', sex: 'XY', status: 1, sample: 'SA622799', read_group: '@RG\tID:S01L2\tSM:SA622799\tLB:sim01\tPU:sim01\tPI:250\tCN:SIM\tPL:ILLUMINA\tPM:Polyester\tDT:2021-04-21\tDS:RNA-Seq|TCRB-CA|SP222784|DO262498|Primary tumour|Tumour', data_type: 'fastq', numLanes: 2, experiment: 'RNA-Seq', single_end: false, library_strandedness: 'UNSTRANDED', date: '2021-04-21'], \
//     "./TCRB-CA.DO262498.SA622799.S01L2.hisat2.novel_splicesites.txt" \
// ])

txt_files_ch = Channel.from([ \
    [analysis_type: 'sequencing_experiment', id: 'S01L1', study_id: 'TCRB-CA', patient: 'DO262498', sex: 'XY', status: 1, sample: 'SA622799', read_group: '@RG\tID:S01L1\tSM:SA622799\tLB:sim01\tPU:sim01\tPI:250\tCN:SIM\tPL:ILLUMINA\tPM:Polyester\tDT:2021-04-21\tDS:RNA-Seq|TCRB-CA|SP222784|DO262498|Primary tumour|Tumour', data_type: 'fastq', numLanes: 2, experiment: 'RNA-Seq', single_end: false, library_strandedness: 'UNSTRANDED', date: '2021-04-21'], \
    file("./TCRB-CA.DO262498.SA622799.S01L1.hisat2.novel_splicesites.txt")], \
    [[analysis_type: 'sequencing_experiment', id: 'S01L2', study_id: 'TCRB-CA', patient: 'DO262498', sex: 'XY', status: 1, sample: 'SA622799', read_group: '@RG\tID:S01L2\tSM:SA622799\tLB:sim01\tPU:sim01\tPI:250\tCN:SIM\tPL:ILLUMINA\tPM:Polyester\tDT:2021-04-21\tDS:RNA-Seq|TCRB-CA|SP222784|DO262498|Primary tumour|Tumour', data_type: 'fastq', numLanes: 2, experiment: 'RNA-Seq', single_end: false, library_strandedness: 'UNSTRANDED', date: '2021-04-21'], \
    file("./TCRB-CA.DO262498.SA622799.S01L2.hisat2.novel_splicesites.txt") \
])

txt_files_ch.flatten().buffer( size: 2 )
    .map{
        meta,txt ->
        [
            [
            id:"${meta.study_id}.${meta.patient}.${meta.sample}",
            study_id:"${meta.study_id}",
            patient:"${meta.patient}",
            sex:"${meta.sex}",
            sample:"${meta.sample}",
            numLanes:"${meta.numLanes}",
            experiment:"${meta.experiment}",
            date:"${meta.date}",
            tool: "${meta.tool}"
            ],
            [
            read_group:"${meta.id}",
            data_type:"${meta.data_type}",
            size:"${meta.size}",
            ],
            txt
        ]
    }.groupTuple(by: 0).
    map{
        meta,info,txt ->
        [
            [
            id:"${meta.study_id}.${meta.patient}.${meta.sample}",
            study_id:"${meta.study_id}",
            patient:"${meta.patient}",
            sex:"${meta.sex}",
            sample:"${meta.sample}",
            numLanes:"${meta.numLanes}",
            experiment:"${meta.experiment}",
            date:"${meta.date}",
            read_group:"${info.read_group.collect()}",
            data_type:"${info.data_type.collect()}",
            size:"${info.size.collect()}",
            tool: "${meta.tool}"
            ],txt.collect()
        ]
    }.set{ch_txts}

// Run the sorting process wisth the test data
workflow {
    NOVEL_SPLICE_MERGE(ch_txts)
}
