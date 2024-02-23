#!/usr/bin/env cwl-runner
{
    "cwlVersion": "v1.0",
    "$graph": [
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--annot-r1",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AddtoBam.cwl/Annotation_R1"
                },
                {
                    "inputBinding": {
                        "prefix": "--cell-order"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/Cell_Order"
                },
                {
                    "inputBinding": {
                        "prefix": "--annot-mol-file"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/Molecular_Annotation"
                },
                {
                    "inputBinding": {
                        "prefix": "--r2-bam"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/R2_Bam"
                },
                {
                    "inputBinding": {
                        "prefix": "--seq-stats"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/Seq_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--tag-calls"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AddtoBam.cwl/Tag_Calls"
                },
                {
                    "inputBinding": {
                        "prefix": "--target-gene-mapping"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AddtoBam.cwl/Target_Gene_Mapping"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "Annotated_mapping_R2.BAM"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/Annotated_Bam"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AddtoBam.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_add_to_bam.py"
            ],
            "class": "CommandLineTool",
            "id": "#AddtoBam.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AlignR2.cwl/Extra_Seqs"
                },
                {
                    "inputBinding": {
                        "prefix": "--index"
                    },
                    "type": "File",
                    "id": "#AlignR2.cwl/Index"
                },
                {
                    "inputBinding": {
                        "prefix": "--r2-fastqs",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AlignR2.cwl/R2"
                },
                {
                    "inputBinding": {
                        "prefix": "--vdj-version"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#AlignR2.cwl/VDJ_Version"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "envDef": [
                        {
                            "envName": "CORES_ALLOCATED_PER_CWL_PROCESS",
                            "envValue": "$(String(runtime.cores))"
                        }
                    ],
                    "class": "EnvVarRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*zip"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AlignR2.cwl/Alignments"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AlignR2.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_align_R2.py"
            ],
            "id": "#AlignR2.cwl",
            "arguments": [
                {
                    "prefix": "--assay",
                    "valueFrom": "WTA"
                }
            ],
            "class": "CommandLineTool"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--umi-option"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateMolecules.cwl/AbSeq_UMI"
                },
                {
                    "inputBinding": {
                        "prefix": "--num-bc"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateMolecules.cwl/Barcode_Num"
                },
                {
                    "inputBinding": {
                        "prefix": "--use-dbec"
                    },
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#AnnotateMolecules.cwl/Use_DBEC"
                },
                {
                    "inputBinding": {
                        "prefix": "--valid-annot"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/Valids"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_GeneStatus.csv.*"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/Gene_Status_List"
                },
                {
                    "outputBinding": {
                        "glob": "*_Annotation_Molecule.csv.*"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/Mol_Annot_List"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AnnotateMolecules.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_annotate_molecules.py"
            ],
            "class": "CommandLineTool",
            "id": "#AnnotateMolecules.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--label-version"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateR1.cwl/Label_Version"
                },
                {
                    "inputBinding": {
                        "prefix": "--R1"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/R1"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_Annotation_R1.csv.gz"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/Annotation_R1"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AnnotateR1.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_annotate_R1.py"
            ],
            "class": "CommandLineTool",
            "id": "#AnnotateR1.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateR2.cwl/Extra_Seqs"
                },
                {
                    "inputBinding": {
                        "prefix": "--gtf"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateR2.cwl/GTF_Annotation"
                },
                {
                    "inputBinding": {
                        "prefix": "--R2-zip"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/R2_zip"
                },
                {
                    "inputBinding": {
                        "prefix": "--transcript-length"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateR2.cwl/Transcript_Length"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*Annotation_R2.csv.gz"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/Annot_R2"
                },
                {
                    "outputBinding": {
                        "glob": "*-annot.gtf"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateR2.cwl/GTF"
                },
                {
                    "outputBinding": {
                        "glob": "*mapping_R2.BAM"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/R2_Bam"
                },
                {
                    "outputBinding": {
                        "glob": "*_picard_quality_metrics.csv.gz"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/R2_Quality_Metrics"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AnnotateR2.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_annotate_R2.py"
            ],
            "id": "#AnnotateR2.cwl",
            "arguments": [
                {
                    "prefix": "--assay",
                    "valueFrom": "WTA"
                }
            ],
            "class": "CommandLineTool"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--umi-option"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateReads.cwl/AbSeq_UMI"
                },
                {
                    "inputBinding": {
                        "prefix": "--bam-input"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#AnnotateReads.cwl/Bam_Input"
                },
                {
                    "inputBinding": {
                        "prefix": "--extra-seqs",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/Extra_Seqs"
                },
                {
                    "inputBinding": {
                        "prefix": "--filtering-stats"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/Filter_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--label-version"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateReads.cwl/Label_Version"
                },
                {
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#AnnotateReads.cwl/Putative_Cell_Call"
                },
                {
                    "inputBinding": {
                        "prefix": "--annotR1",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/R1_Annotation"
                },
                {
                    "inputBinding": {
                        "prefix": "--annotR2",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/R2_Annotation"
                },
                {
                    "inputBinding": {
                        "prefix": "--r2-quality-metrics"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/R2_Quality_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--reference-panel-names"
                    },
                    "type": "File",
                    "id": "#AnnotateReads.cwl/Reference_Panel_Names"
                },
                {
                    "inputBinding": {
                        "prefix": "--sample-tags-version"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#AnnotateReads.cwl/Sample_Tags_Version"
                },
                {
                    "inputBinding": {
                        "prefix": "--subsample-tags"
                    },
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#AnnotateReads.cwl/Subsample_Tags"
                },
                {
                    "inputBinding": {
                        "prefix": "--target-gene-mapping"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/Target_Gene_Mapping"
                },
                {
                    "inputBinding": {
                        "prefix": "--vdj-version"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#AnnotateReads.cwl/VDJ_Version"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_Annotation_Read.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/Annotation_Read"
                },
                {
                    "outputBinding": {
                        "glob": "metadata.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).is_trueno)"
                    },
                    "type": "boolean",
                    "id": "#AnnotateReads.cwl/Is_Trueno"
                },
                {
                    "outputBinding": {
                        "glob": "metadata.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).sample)"
                    },
                    "type": "string",
                    "id": "#AnnotateReads.cwl/Sample_Name"
                },
                {
                    "outputBinding": {
                        "glob": "*_SeqMetrics.csv.gz"
                    },
                    "type": "File",
                    "id": "#AnnotateReads.cwl/Seq_Metrics"
                },
                {
                    "outputBinding": {
                        "glob": "*Sorted_Valid_Reads.csv.*"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#AnnotateReads.cwl/Valid_Reads"
                },
                {
                    "outputBinding": {
                        "glob": "metadata.json"
                    },
                    "type": "File",
                    "id": "#AnnotateReads.cwl/metadata"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#AnnotateReads.cwl/output"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_IG_Valid_Reads.fasta.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/validIgReads"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_TCR_Valid_Reads.fasta.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#AnnotateReads.cwl/validTcrReads"
                }
            ],
            "baseCommand": [
                "mist_annotate_reads.py"
            ],
            "class": "CommandLineTool",
            "id": "#AnnotateReads.cwl"
        },
        {
            "inputs": [
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#BundleLogs.cwl/log_files"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "id": "#BundleLogs.cwl/logs_dir"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */\n  function uuid() {\n    var uuid = \"\", i, random;\n    for (i = 0; i < 32; i++) {\n      random = Math.random() * 16 | 0;\n      if (i == 8 || i == 12 || i == 16 || i == 20) {\n        uuid += \"-\";\n      }\n      uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);\n    }\n    return uuid;\n  }\n  var listing = [];\n  for (var i = 0; i < inputs.log_files.length; i++) {\n    var log_file = inputs.log_files[i];\n    log_file.basename = uuid() + \"-\" + log_file.basename;\n    listing.push(log_file);\n  }\n  return ({\n    logs_dir: {\n      class: \"Directory\",\n      basename: \"Logs\",\n      listing: listing\n    }\n  });\n}",
            "id": "#BundleLogs.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--metrics-files",
                        "itemSeparator": ","
                    },
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#BundleMetrics.cwl/Metrics"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*.zip"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#BundleMetrics.cwl/Bundle_Metrics"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#BundleMetrics.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_bundle_metrics.py"
            ],
            "class": "CommandLineTool",
            "id": "#BundleMetrics.cwl",
            "hints": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ]
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 0
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Cell_Classifier.cwl/molsPerCellMatrix"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*cell_type_experimental.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Cell_Classifier.cwl/cellTypePredictions"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#Cell_Classifier.cwl/log"
                }
            ],
            "baseCommand": [
                "mist_cell_classifier.py"
            ],
            "class": "CommandLineTool",
            "id": "#Cell_Classifier.cwl"
        },
        {
            "inputs": [
                {
                    "doc": "The minimum size (megabytes) of a file that should get split into chunks of a size designated in NumRecordsPerSplit\n",
                    "inputBinding": {
                        "prefix": "--min-split-size"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CheckFastqs.cwl/MinChunkSize"
                },
                {
                    "inputBinding": {
                        "prefix": "--reads",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#CheckFastqs.cwl/Reads"
                },
                {
                    "inputBinding": {
                        "prefix": "--subsample"
                    },
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#CheckFastqs.cwl/Subsample"
                },
                {
                    "inputBinding": {
                        "prefix": "--subsample-seed"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CheckFastqs.cwl/UserInputSubsampleSeed"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "CheckFastqs does several quality control routines including: (1) ensuring that read pair file names are formatted correctly and contain a read pair mate; (2) disambiguating the \"Subsample Reads\" input and; (3) if not provided, generating a subsampling seed that the downstream instances can use.\n",
            "baseCommand": [
                "mist_check_fastqs.py"
            ],
            "id": "#CheckFastqs.cwl",
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "fastq_read_pairs.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).fastq_read_pairs)\n"
                    },
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/filename"
                                },
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/readFlag"
                                },
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/readPairId"
                                },
                                {
                                    "type": "string",
                                    "name": "#CheckFastqs.cwl/FastqReadPairs/library"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#CheckFastqs.cwl/FastqReadPairs"
                },
                {
                    "outputBinding": {
                        "glob": "files_to_skip_split_and_subsample.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).files_to_skip_split_and_subsample)\n"
                    },
                    "type": {
                        "items": "string",
                        "type": "array"
                    },
                    "id": "#CheckFastqs.cwl/FilesToSkipSplitAndSubsample"
                },
                {
                    "outputBinding": {
                        "glob": "subsampling_info.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).subsampling_seed)\n"
                    },
                    "type": "int",
                    "id": "#CheckFastqs.cwl/SubsampleSeed"
                },
                {
                    "outputBinding": {
                        "glob": "subsampling_info.json",
                        "loadContents": true,
                        "outputEval": "$(JSON.parse(self[0].contents).subsampling_ratio)\n"
                    },
                    "type": "float",
                    "id": "#CheckFastqs.cwl/SubsamplingRatio"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#CheckFastqs.cwl/log"
                }
            ],
            "class": "CommandLineTool"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--abseq-reference"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#CheckReference.cwl/AbSeq_Reference"
                },
                {
                    "inputBinding": {
                        "prefix": "--label-version"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CheckReference.cwl/Label_Version"
                },
                {
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#CheckReference.cwl/Putative_Cell_Call"
                },
                {
                    "inputBinding": {
                        "prefix": "--reference",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#CheckReference.cwl/Reference"
                },
                {
                    "inputBinding": {
                        "prefix": "--sample-tags-version"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#CheckReference.cwl/Sample_Tags_Version"
                },
                {
                    "inputBinding": {
                        "prefix": "--supplemental-reference"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#CheckReference.cwl/Supplemental_Reference"
                },
                {
                    "inputBinding": {
                        "prefix": "--vdj-version"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#CheckReference.cwl/VDJ_Version"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "combined_extra_seq.fasta"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Extra_Seqs"
                },
                {
                    "outputBinding": {
                        "glob": "full-gene-list.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Full_Genes"
                },
                {
                    "outputBinding": {
                        "glob": "*gtf",
                        "outputEval": "${\n    if (self.length == 1) { // modified GTF\n        return self;\n    } else if (self.length == 0){ // WTA without extra seqs or Targeted\n        for (var i = 0; i < inputs.Reference.length; i++) {\n            if (inputs.Reference[i].basename.toLowerCase().indexOf('gtf') !== -1) {\n                return inputs.Reference[i];\n            }\n        }\n        return null\n    }\n}\n"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/GTF"
                },
                {
                    "outputBinding": {
                        "glob": "*-annot.*",
                        "outputEval": "${\n    if (self.length == 1) { // Targeted\n        return self;\n    } else if (self.length == 0){ // WTA without extra seqs or targets\n        for (var i = 0; i < inputs.Reference.length; i++) {\n            if (inputs.Reference[i].basename.toLowerCase().indexOf('tar.gz') !== -1) {\n                return inputs.Reference[i];\n            }\n        }\n        return null\n    }\n}\n"
                    },
                    "type": "File",
                    "id": "#CheckReference.cwl/Index"
                },
                {
                    "outputBinding": {
                        "glob": "reference_panel_names.json"
                    },
                    "type": "File",
                    "id": "#CheckReference.cwl/Reference_Panel_Names"
                },
                {
                    "outputBinding": {
                        "glob": "target-gene.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Target_Gene_Mapping"
                },
                {
                    "outputBinding": {
                        "glob": "transcript_length.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#CheckReference.cwl/Transcript_Length"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#CheckReference.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_check_references.py"
            ],
            "class": "CommandLineTool",
            "id": "#CheckReference.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--cell-order"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/Cell_Order"
                },
                {
                    "inputBinding": {
                        "prefix": "--dense-data-table"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#DensetoSparse.cwl/Dense_Data_Table"
                },
                {
                    "inputBinding": {
                        "prefix": "--gene-list"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/Gene_List"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*.csv.gz"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/Data_Tables"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#DensetoSparse.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_dense_to_sparse.py"
            ],
            "class": "CommandLineTool",
            "id": "#DensetoSparse.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#DensetoSparseFile.cwl/GDT_cell_order"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "stdout": "cell_order.json",
            "outputs": [
                {
                    "type": "stdout",
                    "id": "#DensetoSparseFile.cwl/Cell_Order"
                }
            ],
            "baseCommand": "cat",
            "id": "#DensetoSparseFile.cwl",
            "class": "CommandLineTool"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--basic-algo-only"
                    },
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#GetDataTable.cwl/Basic_Algo_Only"
                },
                {
                    "inputBinding": {
                        "prefix": "--exact-cell-count"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#GetDataTable.cwl/Exact_Cell_Count"
                },
                {
                    "inputBinding": {
                        "prefix": "--full-gene-list"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Full_Genes"
                },
                {
                    "inputBinding": {
                        "prefix": "--gene-status",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Gene_Status_List"
                },
                {
                    "inputBinding": {
                        "prefix": "--mol-annot",
                        "itemSeparator": ","
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Molecule_Annotation_List"
                },
                {
                    "inputBinding": {
                        "prefix": "--putative-cell-call"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#GetDataTable.cwl/Putative_Cell_Call"
                },
                {
                    "inputBinding": {
                        "prefix": "--seq-stats"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Seq_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--tag-names",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#GetDataTable.cwl/Tag_Names"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "metrics-files.tar.gz"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Annot_Files"
                },
                {
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*.png"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#GetDataTable.cwl/Cell_Label_Filter"
                },
                {
                    "outputBinding": {
                        "glob": "cell_order.json"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Cell_Order"
                },
                {
                    "outputBinding": {
                        "glob": "*_Annotation_Molecule_corrected.csv.gz"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Corrected_Molecular_Annotation"
                },
                {
                    "outputBinding": {
                        "glob": "*PerCell_Dense.csv.gz"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Dense_Data_Tables"
                },
                {
                    "outputBinding": {
                        "glob": "*PerCell_Unfiltered_Dense.csv.gz"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#GetDataTable.cwl/Dense_Data_Tables_Unfiltered"
                },
                {
                    "outputBinding": {
                        "glob": "*_Expression_Data.st.gz"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Expression_Data"
                },
                {
                    "outputBinding": {
                        "glob": "*_Expression_Data_Unfiltered.st.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Expression_Data_Unfiltered"
                },
                {
                    "outputBinding": {
                        "glob": "gene_list.json"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Gene_List"
                },
                {
                    "outputBinding": {
                        "glob": "Annotations/*_Annotation_Molecule.csv.gz"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/Molecular_Annotation"
                },
                {
                    "outputBinding": {
                        "glob": "Cell_Label_Filtering/*_Putative_Cells_Origin.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Putative_Cells_Origin"
                },
                {
                    "outputBinding": {
                        "glob": "Annotations/*_Annotation_Molecule_Trueno.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Tag_Annotation"
                },
                {
                    "outputBinding": {
                        "glob": "Trueno/*_Calls.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/Tag_Calls"
                },
                {
                    "outputBinding": {
                        "glob": "Trueno/*"
                    },
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#GetDataTable.cwl/Trueno_out"
                },
                {
                    "outputBinding": {
                        "glob": "Annotations/*_UMI_Adjusted_CellLabel_Stats.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/UMI_Adjusted_CellLabel_Stats"
                },
                {
                    "outputBinding": {
                        "glob": "Annotations/*_UMI_Adjusted_Stats.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#GetDataTable.cwl/UMI_Adjusted_Stats"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#GetDataTable.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_get_datatables.py"
            ],
            "class": "CommandLineTool",
            "id": "#GetDataTable.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    },
                    "type": "File",
                    "id": "#IndexBAM.cwl/BamFile"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "samtools_index.log",
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*.bai"
                    },
                    "type": "File",
                    "id": "#IndexBAM.cwl/Index"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#IndexBAM.cwl/log"
                }
            ],
            "baseCommand": [
                "samtools",
                "index"
            ],
            "id": "#IndexBAM.cwl",
            "arguments": [
                {
                    "position": 2,
                    "valueFrom": "${\n    return inputs.BamFile.basename + \".bai\"\n}"
                }
            ],
            "class": "CommandLineTool"
        },
        {
            "inputs": [],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/AbSeq_UMI"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/Barcode_Num"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#InternalSettings.cwl/Extra_Seqs"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/Label_Version"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/MinChunkSize"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#InternalSettings.cwl/NumRecordsPerSplit"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#InternalSettings.cwl/Putative_Cell_Call"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Read_Filter_Off"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#InternalSettings.cwl/Seq_Run"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/Subsample_Tags"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Target_analysis"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#InternalSettings.cwl/Use_DBEC"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#InternalSettings.cwl/VDJ_VGene_Evalue"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  var internalInputs = [\n    '_Label_Version',\n    '_Read_Filter_Off',\n    '_Barcode_Num',\n    '_Seq_Run',\n    '_AbSeq_UMI',\n    '_Putative_Cell_Call',\n    '_Use_DBEC',\n    '_Extra_Seqs',\n    '_MinChunkSize',\n    '_NumRecordsPerSplit',\n    '_Target_analysis',\n    '_Subsample_Tags',\n    '_VDJ_VGene_Evalue',\n    '_VDJ_JGene_Evalue',\n  ];\n  var internalOutputs = {}\n  for (var i = 0; i < internalInputs.length; i++) {\n    var internalInput = internalInputs[i];\n    var internalOutput = internalInput.slice(1); // remove leading underscore\n    if (inputs.hasOwnProperty(internalInput)) {\n      internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output\n    } else {\n      internalOutputs[internalOutput] = null; // if input not specified, provide a null\n    }\n  }\n  return internalOutputs;\n}",
            "id": "#InternalSettings.cwl"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/AbSeq_Reference",
                    "label": "AbSeq Reference"
                },
                {
                    "doc": "Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set.",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#main/Basic_Algo_Only",
                    "label": "Disable Refined Putative Cell Calling"
                },
                {
                    "doc": "Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count",
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/Exact_Cell_Count",
                    "label": "Exact Cell Count"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#main/Reads",
                    "label": "Reads"
                },
                {
                    "type": "File",
                    "id": "#main/Reference_Genome",
                    "label": "Reference Genome"
                },
                {
                    "doc": "The sample multiplexing kit version.  This option should only be set for a multiplexed experiment.",
                    "type": [
                        "null",
                        {
                            "symbols": [
                                "#main/Sample_Tags_Version/Sample_Tags_Version/human",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/hs",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/mouse",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/mm",
                                "#main/Sample_Tags_Version/Sample_Tags_Version/custom"
                            ],
                            "type": "enum",
                            "name": "#main/Sample_Tags_Version/Sample_Tags_Version"
                        }
                    ],
                    "id": "#main/Sample_Tags_Version",
                    "label": "Sample Tags Version"
                },
                {
                    "doc": "Any number of reads >1 or a fraction between 0 < n < 1 to indicate the percentage of reads to subsample.\n",
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#main/Subsample",
                    "label": "Subsample Reads"
                },
                {
                    "doc": "For use when replicating a previous subsampling run only. Obtain the seed generated from the log file for the SplitFastQ node.\n",
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/Subsample_seed",
                    "label": "Subsample Seed"
                },
                {
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Supplemental_Reference",
                    "label": "Supplemental Reference"
                },
                {
                    "doc": "Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Do not use the special characters: &, (), [], {}, <>, ?, |\n",
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Tag_Names",
                    "label": "Tag Names"
                },
                {
                    "type": "File",
                    "id": "#main/Transcriptome_Annotation",
                    "label": "Transcriptome Annotation"
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "The BD Rhapsody\u2122 WTA Analysis Pipeline is used to create sequencing libraries from single cell transcriptomes without having to specify a targeted panel.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files, a reference genome file and a transcriptome annotation file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
            "label": "BD Rhapsody\u2122 WTA Analysis Pipeline",
            "steps": [
                {
                    "run": "#AddtoBam.cwl",
                    "scatter": [
                        "#main/AddtoBam/R2_Bam"
                    ],
                    "in": [
                        {
                            "source": "#main/AnnotateR1/Annotation_R1",
                            "id": "#main/AddtoBam/Annotation_R1"
                        },
                        {
                            "source": "#main/Dense_to_Sparse_File/Cell_Order",
                            "id": "#main/AddtoBam/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Corrected_Molecular_Annotation",
                            "id": "#main/AddtoBam/Molecular_Annotation"
                        },
                        {
                            "source": "#main/AnnotateR2/R2_Bam",
                            "id": "#main/AddtoBam/R2_Bam"
                        },
                        {
                            "source": "#main/AnnotateReads/Seq_Metrics",
                            "id": "#main/AddtoBam/Seq_Metrics"
                        },
                        {
                            "source": "#main/GetDataTable/Tag_Calls",
                            "id": "#main/AddtoBam/Tag_Calls"
                        },
                        {
                            "source": "#main/CheckReference/Target_Gene_Mapping",
                            "id": "#main/AddtoBam/Target_Gene_Mapping"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 16000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AddtoBam",
                    "out": [
                        "#main/AddtoBam/Annotated_Bam",
                        "#main/AddtoBam/output"
                    ]
                },
                {
                    "run": "#AlignR2.cwl",
                    "out": [
                        "#main/AlignR2/Alignments",
                        "#main/AlignR2/output"
                    ],
                    "requirements": [
                        {
                            "coresMin": 16,
                            "ramMin": 48000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AlignR2",
                    "in": [
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AlignR2/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/Index",
                            "id": "#main/AlignR2/Index"
                        },
                        {
                            "source": "#main/QualityFilter/R2",
                            "id": "#main/AlignR2/R2"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/AlignR2/VDJ_Version"
                        }
                    ]
                },
                {
                    "run": "#AnnotateMolecules.cwl",
                    "scatter": [
                        "#main/AnnotateMolecules/Valids"
                    ],
                    "in": [
                        {
                            "source": "#main/Internal_Settings/AbSeq_UMI",
                            "id": "#main/AnnotateMolecules/AbSeq_UMI"
                        },
                        {
                            "source": "#main/Internal_Settings/Barcode_Num",
                            "id": "#main/AnnotateMolecules/Barcode_Num"
                        },
                        {
                            "source": "#main/Internal_Settings/Use_DBEC",
                            "id": "#main/AnnotateMolecules/Use_DBEC"
                        },
                        {
                            "source": "#main/AnnotateReads/Valid_Reads",
                            "id": "#main/AnnotateMolecules/Valids"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 32000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AnnotateMolecules",
                    "out": [
                        "#main/AnnotateMolecules/Mol_Annot_List",
                        "#main/AnnotateMolecules/Gene_Status_List",
                        "#main/AnnotateMolecules/output"
                    ]
                },
                {
                    "id": "#main/AnnotateR1",
                    "out": [
                        "#main/AnnotateR1/Annotation_R1",
                        "#main/AnnotateR1/output"
                    ],
                    "run": "#AnnotateR1.cwl",
                    "scatter": [
                        "#main/AnnotateR1/R1"
                    ],
                    "in": [
                        {
                            "source": "#main/Internal_Settings/Label_Version",
                            "id": "#main/AnnotateR1/Label_Version"
                        },
                        {
                            "source": "#main/QualityFilter/R1",
                            "id": "#main/AnnotateR1/R1"
                        }
                    ]
                },
                {
                    "run": "#AnnotateR2.cwl",
                    "scatter": [
                        "#main/AnnotateR2/R2_zip"
                    ],
                    "in": [
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AnnotateR2/Extra_Seqs"
                        },
                        {
                            "source": "#main/CheckReference/GTF",
                            "id": "#main/AnnotateR2/GTF_Annotation"
                        },
                        {
                            "source": "#main/AlignR2/Alignments",
                            "id": "#main/AnnotateR2/R2_zip"
                        },
                        {
                            "source": "#main/CheckReference/Transcript_Length",
                            "id": "#main/AnnotateR2/Transcript_Length"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 10000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AnnotateR2",
                    "out": [
                        "#main/AnnotateR2/Annot_R2",
                        "#main/AnnotateR2/R2_Bam",
                        "#main/AnnotateR2/GTF",
                        "#main/AnnotateR2/output",
                        "#main/AnnotateR2/R2_Quality_Metrics"
                    ]
                },
                {
                    "run": "#AnnotateReads.cwl",
                    "out": [
                        "#main/AnnotateReads/Seq_Metrics",
                        "#main/AnnotateReads/Valid_Reads",
                        "#main/AnnotateReads/Annotation_Read",
                        "#main/AnnotateReads/Is_Trueno",
                        "#main/AnnotateReads/Sample_Name",
                        "#main/AnnotateReads/output",
                        "#main/AnnotateReads/validTcrReads",
                        "#main/AnnotateReads/validIgReads",
                        "#main/AnnotateReads/metadata"
                    ],
                    "requirements": [
                        {
                            "ramMin": 32000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/AnnotateReads",
                    "in": [
                        {
                            "source": "#main/Internal_Settings/AbSeq_UMI",
                            "id": "#main/AnnotateReads/AbSeq_UMI"
                        },
                        {
                            "source": "#main/CheckReference/Extra_Seqs",
                            "id": "#main/AnnotateReads/Extra_Seqs"
                        },
                        {
                            "source": "#main/Bundle_Filter_Metrics/Bundle_Metrics",
                            "id": "#main/AnnotateReads/Filter_Metrics"
                        },
                        {
                            "source": "#main/Internal_Settings/Label_Version",
                            "id": "#main/AnnotateReads/Label_Version"
                        },
                        {
                            "source": "#main/Internal_Settings/Putative_Cell_Call",
                            "id": "#main/AnnotateReads/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/AnnotateR1/Annotation_R1",
                            "id": "#main/AnnotateReads/R1_Annotation"
                        },
                        {
                            "source": "#main/AnnotateR2/Annot_R2",
                            "id": "#main/AnnotateReads/R2_Annotation"
                        },
                        {
                            "source": "#main/Bundle_R2_Quality_Metrics/Bundle_Metrics",
                            "id": "#main/AnnotateReads/R2_Quality_Metrics"
                        },
                        {
                            "source": "#main/CheckReference/Reference_Panel_Names",
                            "id": "#main/AnnotateReads/Reference_Panel_Names"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/AnnotateReads/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Internal_Settings/Subsample_Tags",
                            "id": "#main/AnnotateReads/Subsample_Tags"
                        },
                        {
                            "source": "#main/CheckReference/Target_Gene_Mapping",
                            "id": "#main/AnnotateReads/Target_Gene_Mapping"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/AnnotateReads/VDJ_Version"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/AnnotateVDJResults/vdjCellsDatatable",
                        "#main/AnnotateVDJResults/vdjCellsDatatableUnfiltered",
                        "#main/AnnotateVDJResults/vdjCellsDatatableCellCorrected",
                        "#main/AnnotateVDJResults/vdjCellsDatatableDBECCellCorrected",
                        "#main/AnnotateVDJResults/vdjCellChainDatatableUnfiltered",
                        "#main/AnnotateVDJResults/vdjValidReadsDatatable",
                        "#main/AnnotateVDJResults/vdjInvalidReadsDatatable",
                        "#main/AnnotateVDJResults/vdjMetricsJson",
                        "#main/AnnotateVDJResults/vdjMetricsCsv",
                        "#main/AnnotateVDJResults/vdjReadsAndMoleculesPerCellFigure",
                        "#main/AnnotateVDJResults/vdjReadsPerCellByChainTypeFigure"
                    ],
                    "run": "#VDJ_Annotate_Molecules.cwl",
                    "id": "#main/AnnotateVDJResults",
                    "in": [
                        {
                            "source": "#main/AnnotateReads/Sample_Name",
                            "id": "#main/AnnotateVDJResults/Sample_Name"
                        },
                        {
                            "source": "#main/CellClassifier/cellTypePredictions",
                            "id": "#main/AnnotateVDJResults/cellTypeMapping"
                        },
                        {
                            "valueFrom": "$([])",
                            "id": "#main/AnnotateVDJResults/chainsToIgnore"
                        },
                        {
                            "source": "#main/Internal_Settings/VDJ_JGene_Evalue",
                            "id": "#main/AnnotateVDJResults/evalueJgene"
                        },
                        {
                            "source": "#main/Internal_Settings/VDJ_VGene_Evalue",
                            "id": "#main/AnnotateVDJResults/evalueVgene"
                        },
                        {
                            "source": "#main/VDJ_GatherIGCalls/gatheredCalls",
                            "id": "#main/AnnotateVDJResults/igCalls"
                        },
                        {
                            "source": "#main/AnnotateReads/metadata",
                            "id": "#main/AnnotateVDJResults/metadata"
                        },
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/AnnotateVDJResults/putativeCells"
                        },
                        {
                            "source": "#main/VDJ_GatherTCRCalls/gatheredCalls",
                            "id": "#main/AnnotateVDJResults/tcrCalls"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/AnnotateVDJResults/vdjVersion"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/BundleLogs/logs_dir"
                    ],
                    "run": "#BundleLogs.cwl",
                    "id": "#main/BundleLogs",
                    "in": [
                        {
                            "source": [
                                "#main/AnnotateReads/output",
                                "#main/AnnotateR1/output",
                                "#main/AnnotateR2/output",
                                "#main/CheckReference/output",
                                "#main/GetDataTable/output",
                                "#main/Metrics/output",
                                "#main/AddtoBam/output",
                                "#main/AnnotateMolecules/output",
                                "#main/QualityFilter/output",
                                "#main/CheckFastqs/log",
                                "#main/SplitAndSubsample/log",
                                "#main/MergeBAM/log",
                                "#main/Dense_to_Sparse_Datatable/output",
                                "#main/Dense_to_Sparse_Datatable_Unfiltered/output",
                                "#main/IndexBAM/log",
                                "#main/CellClassifier/log",
                                "#main/Bundle_Filter_Metrics/output",
                                "#main/Bundle_R2_Quality_Metrics/output"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/BundleLogs/log_files"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Bundle_Filter_Metrics/Bundle_Metrics",
                        "#main/Bundle_Filter_Metrics/output"
                    ],
                    "run": "#BundleMetrics.cwl",
                    "id": "#main/Bundle_Filter_Metrics",
                    "in": [
                        {
                            "source": "#main/QualityFilter/Filter_Metrics",
                            "id": "#main/Bundle_Filter_Metrics/Metrics"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Bundle_R2_Quality_Metrics/Bundle_Metrics",
                        "#main/Bundle_R2_Quality_Metrics/output"
                    ],
                    "run": "#BundleMetrics.cwl",
                    "id": "#main/Bundle_R2_Quality_Metrics",
                    "in": [
                        {
                            "source": "#main/AnnotateR2/R2_Quality_Metrics",
                            "id": "#main/Bundle_R2_Quality_Metrics/Metrics"
                        }
                    ]
                },
                {
                    "run": "#Cell_Classifier.cwl",
                    "out": [
                        "#main/CellClassifier/cellTypePredictions",
                        "#main/CellClassifier/log"
                    ],
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/CellClassifier",
                    "in": [
                        {
                            "source": "#main/FindDataTableForCellClassifier/molsPerCellMatrixForCellClassifier",
                            "id": "#main/CellClassifier/molsPerCellMatrix"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/CheckFastqs/SubsampleSeed",
                        "#main/CheckFastqs/SubsamplingRatio",
                        "#main/CheckFastqs/FilesToSkipSplitAndSubsample",
                        "#main/CheckFastqs/FastqReadPairs",
                        "#main/CheckFastqs/log"
                    ],
                    "run": "#CheckFastqs.cwl",
                    "id": "#main/CheckFastqs",
                    "in": [
                        {
                            "source": "#main/Internal_Settings/MinChunkSize",
                            "id": "#main/CheckFastqs/MinChunkSize"
                        },
                        {
                            "source": "#main/Reads",
                            "id": "#main/CheckFastqs/Reads"
                        },
                        {
                            "source": "#main/Subsample_Settings/Subsample_Reads",
                            "id": "#main/CheckFastqs/Subsample"
                        },
                        {
                            "source": "#main/Subsample_Settings/Subsample_Seed",
                            "id": "#main/CheckFastqs/UserInputSubsampleSeed"
                        }
                    ]
                },
                {
                    "run": "#CheckReference.cwl",
                    "out": [
                        "#main/CheckReference/Index",
                        "#main/CheckReference/Extra_Seqs",
                        "#main/CheckReference/Reference_Panel_Names",
                        "#main/CheckReference/Full_Genes",
                        "#main/CheckReference/output",
                        "#main/CheckReference/Transcript_Length",
                        "#main/CheckReference/GTF",
                        "#main/CheckReference/Target_Gene_Mapping"
                    ],
                    "requirements": [
                        {
                            "ramMin": 10000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/CheckReference",
                    "in": [
                        {
                            "source": "#main/AbSeq_Reference",
                            "id": "#main/CheckReference/AbSeq_Reference"
                        },
                        {
                            "source": "#main/Internal_Settings/Label_Version",
                            "id": "#main/CheckReference/Label_Version"
                        },
                        {
                            "source": "#main/Internal_Settings/Putative_Cell_Call",
                            "id": "#main/CheckReference/Putative_Cell_Call"
                        },
                        {
                            "source": [
                                "#main/Transcriptome_Annotation",
                                "#main/Reference_Genome"
                            ],
                            "id": "#main/CheckReference/Reference"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Sample_Tags_Version",
                            "id": "#main/CheckReference/Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Supplemental_Reference",
                            "id": "#main/CheckReference/Supplemental_Reference"
                        },
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/CheckReference/VDJ_Version"
                        }
                    ]
                },
                {
                    "run": "#DensetoSparse.cwl",
                    "scatter": [
                        "#main/Dense_to_Sparse_Datatable/Dense_Data_Table"
                    ],
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_File/Cell_Order",
                            "id": "#main/Dense_to_Sparse_Datatable/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Dense_Data_Tables",
                            "id": "#main/Dense_to_Sparse_Datatable/Dense_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Gene_List",
                            "id": "#main/Dense_to_Sparse_Datatable/Gene_List"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/Dense_to_Sparse_Datatable",
                    "out": [
                        "#main/Dense_to_Sparse_Datatable/Data_Tables",
                        "#main/Dense_to_Sparse_Datatable/output"
                    ]
                },
                {
                    "run": "#DensetoSparse.cwl",
                    "scatter": [
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/Dense_Data_Table"
                    ],
                    "in": [
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Cell_Order"
                        },
                        {
                            "source": "#main/GetDataTable/Dense_Data_Tables_Unfiltered",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Dense_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Gene_List",
                            "id": "#main/Dense_to_Sparse_Datatable_Unfiltered/Gene_List"
                        }
                    ],
                    "requirements": [
                        {
                            "ramMin": 4000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/Dense_to_Sparse_Datatable_Unfiltered",
                    "out": [
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/Data_Tables",
                        "#main/Dense_to_Sparse_Datatable_Unfiltered/output"
                    ]
                },
                {
                    "out": [
                        "#main/Dense_to_Sparse_File/Cell_Order"
                    ],
                    "run": "#DensetoSparseFile.cwl",
                    "id": "#main/Dense_to_Sparse_File",
                    "in": [
                        {
                            "source": "#main/GetDataTable/Cell_Order",
                            "id": "#main/Dense_to_Sparse_File/GDT_cell_order"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/FindDataTableForCellClassifier/molsPerCellMatrixForCellClassifier"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": {
                                    "items": "File",
                                    "type": "array"
                                },
                                "id": "#main/FindDataTableForCellClassifier/6a1da587-7c26-4715-a5ab-d88bc9457216/dataTables"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "type": "File",
                                "id": "#main/FindDataTableForCellClassifier/6a1da587-7c26-4715-a5ab-d88bc9457216/molsPerCellMatrixForCellClassifier"
                            }
                        ],
                        "id": "#main/FindDataTableForCellClassifier/6a1da587-7c26-4715-a5ab-d88bc9457216",
                        "expression": "${\n  for (var i = 0; i < inputs.dataTables.length; i++) {\n    var dataTable = inputs.dataTables[i];\n    if (dataTable.basename.indexOf(\"_RSEC_MolsPerCell.csv\") >= 0) {\n      return({molsPerCellMatrixForCellClassifier: dataTable});\n    }\n  }\n  return({molsPerCellMatrixForCellClassifier: null});\n}",
                        "class": "ExpressionTool"
                    },
                    "id": "#main/FindDataTableForCellClassifier",
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_Datatable/Data_Tables",
                            "id": "#main/FindDataTableForCellClassifier/dataTables"
                        }
                    ]
                },
                {
                    "run": "#GetDataTable.cwl",
                    "out": [
                        "#main/GetDataTable/Tag_Calls",
                        "#main/GetDataTable/Molecular_Annotation",
                        "#main/GetDataTable/Corrected_Molecular_Annotation",
                        "#main/GetDataTable/Tag_Annotation",
                        "#main/GetDataTable/Annot_Files",
                        "#main/GetDataTable/Cell_Label_Filter",
                        "#main/GetDataTable/Dense_Data_Tables",
                        "#main/GetDataTable/Dense_Data_Tables_Unfiltered",
                        "#main/GetDataTable/Expression_Data",
                        "#main/GetDataTable/Expression_Data_Unfiltered",
                        "#main/GetDataTable/UMI_Adjusted_Stats",
                        "#main/GetDataTable/UMI_Adjusted_CellLabel_Stats",
                        "#main/GetDataTable/Putative_Cells_Origin",
                        "#main/GetDataTable/Trueno_out",
                        "#main/GetDataTable/output",
                        "#main/GetDataTable/Cell_Order",
                        "#main/GetDataTable/Gene_List"
                    ],
                    "requirements": [
                        {
                            "ramMin": 64000,
                            "class": "ResourceRequirement"
                        }
                    ],
                    "id": "#main/GetDataTable",
                    "in": [
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Basic_Algo_Only",
                            "id": "#main/GetDataTable/Basic_Algo_Only"
                        },
                        {
                            "source": "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                            "id": "#main/GetDataTable/Exact_Cell_Count"
                        },
                        {
                            "source": "#main/CheckReference/Full_Genes",
                            "id": "#main/GetDataTable/Full_Genes"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Gene_Status_List",
                            "id": "#main/GetDataTable/Gene_Status_List"
                        },
                        {
                            "source": "#main/AnnotateMolecules/Mol_Annot_List",
                            "id": "#main/GetDataTable/Molecule_Annotation_List"
                        },
                        {
                            "source": "#main/Internal_Settings/Putative_Cell_Call",
                            "id": "#main/GetDataTable/Putative_Cell_Call"
                        },
                        {
                            "source": "#main/AnnotateReads/Seq_Metrics",
                            "id": "#main/GetDataTable/Seq_Metrics"
                        },
                        {
                            "source": "#main/Multiplexing_Settings/Tag_Sample_Names",
                            "id": "#main/GetDataTable/Tag_Names"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/IndexBAM/Index",
                        "#main/IndexBAM/log"
                    ],
                    "run": "#IndexBAM.cwl",
                    "id": "#main/IndexBAM",
                    "in": [
                        {
                            "source": "#main/MergeBAM/Final_Bam",
                            "id": "#main/IndexBAM/BamFile"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Internal_Settings/Label_Version",
                        "#main/Internal_Settings/Read_Filter_Off",
                        "#main/Internal_Settings/Barcode_Num",
                        "#main/Internal_Settings/Seq_Run",
                        "#main/Internal_Settings/AbSeq_UMI",
                        "#main/Internal_Settings/Putative_Cell_Call",
                        "#main/Internal_Settings/Use_DBEC",
                        "#main/Internal_Settings/Extra_Seqs",
                        "#main/Internal_Settings/MinChunkSize",
                        "#main/Internal_Settings/NumRecordsPerSplit",
                        "#main/Internal_Settings/Target_analysis",
                        "#main/Internal_Settings/Subsample_Tags",
                        "#main/Internal_Settings/VDJ_VGene_Evalue",
                        "#main/Internal_Settings/VDJ_JGene_Evalue"
                    ],
                    "in": [],
                    "run": "#InternalSettings.cwl",
                    "id": "#main/Internal_Settings",
                    "label": "Internal Settings"
                },
                {
                    "out": [
                        "#main/MergeBAM/Final_Bam",
                        "#main/MergeBAM/log"
                    ],
                    "run": "#MergeBAM.cwl",
                    "id": "#main/MergeBAM",
                    "in": [
                        {
                            "source": "#main/AddtoBam/Annotated_Bam",
                            "id": "#main/MergeBAM/BamFiles"
                        },
                        {
                            "source": "#main/AnnotateReads/Is_Trueno",
                            "id": "#main/MergeBAM/Is_Trueno"
                        },
                        {
                            "source": "#main/AnnotateReads/Sample_Name",
                            "id": "#main/MergeBAM/Sample_Name"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Metrics/Metrics_Summary",
                        "#main/Metrics/Metrics_Archive",
                        "#main/Metrics/output"
                    ],
                    "run": "#Metrics.cwl",
                    "id": "#main/Metrics",
                    "in": [
                        {
                            "source": "#main/GetDataTable/Annot_Files",
                            "id": "#main/Metrics/Annot_Files"
                        },
                        {
                            "source": "#main/AnnotateReads/Seq_Metrics",
                            "id": "#main/Metrics/Seq_Metrics"
                        },
                        {
                            "source": "#main/Internal_Settings/Seq_Run",
                            "id": "#main/Metrics/Seq_Run"
                        },
                        {
                            "source": "#main/GetDataTable/UMI_Adjusted_CellLabel_Stats",
                            "id": "#main/Metrics/UMI_Adjusted_Stats"
                        },
                        {
                            "source": "#main/AnnotateVDJResults/vdjMetricsJson",
                            "id": "#main/Metrics/vdjMetricsJson"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Multiplexing_Settings/Tag_Sample_Names",
                        "#main/Multiplexing_Settings/Sample_Tags_Version"
                    ],
                    "in": [
                        {
                            "source": "#main/Sample_Tags_Version",
                            "id": "#main/Multiplexing_Settings/_Sample_Tags_Version"
                        },
                        {
                            "source": "#main/Tag_Names",
                            "id": "#main/Multiplexing_Settings/_Tag_Sample_Names"
                        }
                    ],
                    "run": "#MultiplexingSettings.cwl",
                    "id": "#main/Multiplexing_Settings",
                    "label": "Multiplexing Settings"
                },
                {
                    "out": [
                        "#main/PairReadFiles/ReadPairs"
                    ],
                    "run": "#PairReadFiles.cwl",
                    "id": "#main/PairReadFiles",
                    "in": [
                        {
                            "source": "#main/CheckFastqs/FastqReadPairs",
                            "id": "#main/PairReadFiles/FastqReadPairs"
                        },
                        {
                            "source": "#main/SplitAndSubsample/SplitAndSubsampledFastqs",
                            "id": "#main/PairReadFiles/Reads"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Putative_Cell_Calling_Settings/Exact_Cell_Count",
                        "#main/Putative_Cell_Calling_Settings/Basic_Algo_Only"
                    ],
                    "in": [
                        {
                            "source": "#main/Basic_Algo_Only",
                            "id": "#main/Putative_Cell_Calling_Settings/_Basic_Algo_Only"
                        },
                        {
                            "source": "#main/Exact_Cell_Count",
                            "id": "#main/Putative_Cell_Calling_Settings/_Exact_Cell_Count"
                        }
                    ],
                    "run": "#PutativeCellSettings.cwl",
                    "id": "#main/Putative_Cell_Calling_Settings",
                    "label": "Putative Cell Calling Settings"
                },
                {
                    "run": "#QualityFilter.cwl",
                    "scatter": [
                        "#main/QualityFilter/Split_Read_Pairs"
                    ],
                    "in": [
                        {
                            "source": "#main/Internal_Settings/Label_Version",
                            "id": "#main/QualityFilter/Label_Version"
                        },
                        {
                            "source": "#main/Internal_Settings/Read_Filter_Off",
                            "id": "#main/QualityFilter/Read_Filter_Off"
                        },
                        {
                            "source": "#main/PairReadFiles/ReadPairs",
                            "id": "#main/QualityFilter/Split_Read_Pairs"
                        }
                    ],
                    "scatterMethod": "dotproduct",
                    "id": "#main/QualityFilter",
                    "out": [
                        "#main/QualityFilter/Filter_Metrics",
                        "#main/QualityFilter/R1",
                        "#main/QualityFilter/R2",
                        "#main/QualityFilter/output"
                    ]
                },
                {
                    "out": [
                        "#main/SplitAndSubsample/SplitAndSubsampledFastqs",
                        "#main/SplitAndSubsample/log"
                    ],
                    "run": "#SplitAndSubsample.cwl",
                    "id": "#main/SplitAndSubsample",
                    "in": [
                        {
                            "source": "#main/Reads",
                            "id": "#main/SplitAndSubsample/Fastqs"
                        },
                        {
                            "source": "#main/CheckFastqs/FilesToSkipSplitAndSubsample",
                            "id": "#main/SplitAndSubsample/FilesToSkipSplitAndSubsample"
                        },
                        {
                            "source": "#main/Internal_Settings/NumRecordsPerSplit",
                            "id": "#main/SplitAndSubsample/NumRecordsPerSplit"
                        },
                        {
                            "source": "#main/CheckFastqs/SubsamplingRatio",
                            "id": "#main/SplitAndSubsample/SubsampleRatio"
                        },
                        {
                            "source": "#main/CheckFastqs/SubsampleSeed",
                            "id": "#main/SplitAndSubsample/SubsampleSeed"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/Subsample_Settings/Subsample_Reads",
                        "#main/Subsample_Settings/Subsample_Seed"
                    ],
                    "in": [
                        {
                            "source": "#main/Subsample",
                            "id": "#main/Subsample_Settings/_Subsample_Reads"
                        },
                        {
                            "source": "#main/Subsample_seed",
                            "id": "#main/Subsample_Settings/_Subsample_Seed"
                        }
                    ],
                    "run": "#SubsampleSettings.cwl",
                    "id": "#main/Subsample_Settings",
                    "label": "Subsample Settings"
                },
                {
                    "out": [
                        "#main/Uncompress_Datatables/Uncompressed_Data_Tables",
                        "#main/Uncompress_Datatables/Uncompressed_Expression_Matrix"
                    ],
                    "run": "#UncompressDatatables.cwl",
                    "id": "#main/Uncompress_Datatables",
                    "in": [
                        {
                            "source": "#main/Dense_to_Sparse_Datatable/Data_Tables",
                            "id": "#main/Uncompress_Datatables/Compressed_Data_Table"
                        },
                        {
                            "source": "#main/GetDataTable/Expression_Data",
                            "id": "#main/Uncompress_Datatables/Compressed_Expression_Matrix"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_GatherIGCalls/gatheredCalls"
                    ],
                    "run": "#VDJ_GatherCalls.cwl",
                    "id": "#main/VDJ_GatherIGCalls",
                    "in": [
                        {
                            "source": "#main/VDJ_ig/igCalls",
                            "id": "#main/VDJ_GatherIGCalls/theCalls"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_GatherTCRCalls/gatheredCalls"
                    ],
                    "run": "#VDJ_GatherCalls.cwl",
                    "id": "#main/VDJ_GatherTCRCalls",
                    "in": [
                        {
                            "source": "#main/VDJ_tcr/tcrCalls",
                            "id": "#main/VDJ_GatherTCRCalls/theCalls"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_Settings/VDJ_Version"
                    ],
                    "in": [],
                    "run": "#VDJ_Settings.cwl",
                    "id": "#main/VDJ_Settings",
                    "label": "VDJ Settings"
                },
                {
                    "out": [
                        "#main/VDJ_SplitValidReadsIg/SplitFastaList",
                        "#main/VDJ_SplitValidReadsIg/numFiles",
                        "#main/VDJ_SplitValidReadsIg/log"
                    ],
                    "run": "#VDJ_SplitValidReads.cwl",
                    "id": "#main/VDJ_SplitValidReadsIg",
                    "in": [
                        {
                            "source": "#main/AnnotateReads/validIgReads",
                            "id": "#main/VDJ_SplitValidReadsIg/validReads"
                        }
                    ]
                },
                {
                    "out": [
                        "#main/VDJ_SplitValidReadsTcr/SplitFastaList",
                        "#main/VDJ_SplitValidReadsTcr/numFiles",
                        "#main/VDJ_SplitValidReadsTcr/log"
                    ],
                    "run": "#VDJ_SplitValidReads.cwl",
                    "id": "#main/VDJ_SplitValidReadsTcr",
                    "in": [
                        {
                            "source": "#main/AnnotateReads/validTcrReads",
                            "id": "#main/VDJ_SplitValidReadsTcr/validReads"
                        }
                    ]
                },
                {
                    "id": "#main/VDJ_ig",
                    "out": [
                        "#main/VDJ_ig/igCalls"
                    ],
                    "run": "#VDJ_ig.cwl",
                    "scatter": [
                        "#main/VDJ_ig/validReadsIg"
                    ],
                    "in": [
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_ig/VDJ_Version"
                        },
                        {
                            "source": "#main/VDJ_SplitValidReadsIg/numFiles",
                            "id": "#main/VDJ_ig/numFiles"
                        },
                        {
                            "source": "#main/VDJ_SplitValidReadsIg/SplitFastaList",
                            "id": "#main/VDJ_ig/validReadsIg"
                        }
                    ]
                },
                {
                    "id": "#main/VDJ_tcr",
                    "out": [
                        "#main/VDJ_tcr/tcrCalls"
                    ],
                    "run": "#VDJ_tcr.cwl",
                    "scatter": [
                        "#main/VDJ_tcr/validReadsTcr"
                    ],
                    "in": [
                        {
                            "source": "#main/VDJ_Settings/VDJ_Version",
                            "id": "#main/VDJ_tcr/VDJ_Version"
                        },
                        {
                            "source": "#main/VDJ_SplitValidReadsTcr/numFiles",
                            "id": "#main/VDJ_tcr/numFiles"
                        },
                        {
                            "source": "#main/VDJ_SplitValidReadsTcr/SplitFastaList",
                            "id": "#main/VDJ_tcr/validReadsTcr"
                        }
                    ]
                }
            ],
            "outputs": [
                {
                    "outputSource": "#main/GetDataTable/Cell_Label_Filter",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Cell_Label_Filter",
                    "label": "Cell Label Filter"
                },
                {
                    "outputSource": "#main/Uncompress_Datatables/Uncompressed_Data_Tables",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Data_Tables",
                    "label": "Data Tables"
                },
                {
                    "outputSource": "#main/Dense_to_Sparse_Datatable_Unfiltered/Data_Tables",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Data_Tables_Unfiltered",
                    "label": "Unfiltered Data Tables"
                },
                {
                    "outputSource": "#main/Uncompress_Datatables/Uncompressed_Expression_Matrix",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Expression_Data",
                    "label": "Expression Matrix"
                },
                {
                    "outputSource": "#main/GetDataTable/Expression_Data_Unfiltered",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Expression_Data_Unfiltered",
                    "label": "Unfiltered Expression Matrix"
                },
                {
                    "outputSource": "#main/MergeBAM/Final_Bam",
                    "type": "File",
                    "id": "#main/Final_Bam",
                    "label": "Final BAM File"
                },
                {
                    "outputSource": "#main/IndexBAM/Index",
                    "type": "File",
                    "id": "#main/Final_Bam_Index",
                    "label": "Bam Index"
                },
                {
                    "outputSource": "#main/CellClassifier/cellTypePredictions",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/ImmuneCellClassification(Experimental)"
                },
                {
                    "outputSource": "#main/BundleLogs/logs_dir",
                    "type": "Directory",
                    "id": "#main/Logs",
                    "label": "Pipeline Logs"
                },
                {
                    "outputSource": "#main/Metrics/Metrics_Summary",
                    "type": "File",
                    "id": "#main/Metrics_Summary",
                    "label": "Metrics Summary"
                },
                {
                    "outputSource": "#main/GetDataTable/Trueno_out",
                    "type": [
                        "null",
                        {
                            "items": "File",
                            "type": "array"
                        }
                    ],
                    "id": "#main/Multiplex"
                },
                {
                    "outputSource": "#main/GetDataTable/Putative_Cells_Origin",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/Putative_Cells_Origin",
                    "label": "Putative Cells Origin"
                },
                {
                    "outputSource": "#main/GetDataTable/UMI_Adjusted_Stats",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#main/UMI_Adjusted_Stats",
                    "label": "UMI Adjusted Statistics"
                }
            ],
            "id": "#main",
            "class": "Workflow"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#MergeBAM.cwl/BamFiles"
                },
                {
                    "type": "boolean",
                    "id": "#MergeBAM.cwl/Is_Trueno"
                },
                {
                    "type": "string",
                    "id": "#MergeBAM.cwl/Sample_Name"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "samtools_merge.log",
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_final.BAM"
                    },
                    "type": "File",
                    "id": "#MergeBAM.cwl/Final_Bam"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#MergeBAM.cwl/log"
                }
            ],
            "baseCommand": [
                "samtools",
                "merge"
            ],
            "id": "#MergeBAM.cwl",
            "arguments": [
                {
                    "prefix": "-@",
                    "valueFrom": "$(runtime.cores)"
                },
                {
                    "position": 0,
                    "valueFrom": "${\n    if (inputs.Is_Trueno) {\n        return \"Combined_\" + inputs.Sample_Name + \"_final.BAM\"\n    } else {\n        return inputs.Sample_Name + \"_final.BAM\"\n    }\n}"
                }
            ],
            "class": "CommandLineTool",
            "hints": [
                {
                    "coresMin": 4,
                    "class": "ResourceRequirement"
                }
            ]
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--annot-files"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Annot_Files"
                },
                {
                    "inputBinding": {
                        "prefix": "--seq-stats"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Seq_Metrics"
                },
                {
                    "inputBinding": {
                        "prefix": "--seq-run"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#Metrics.cwl/Seq_Run"
                },
                {
                    "inputBinding": {
                        "prefix": "--umi-adjusted-stats"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Metrics.cwl/UMI_Adjusted_Stats"
                },
                {
                    "inputBinding": {
                        "prefix": "--vdj-metrics-fp"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#Metrics.cwl/vdjMetricsJson"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "internal-metrics-archive.tar.gz"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Metrics_Archive"
                },
                {
                    "outputBinding": {
                        "glob": "*_Metrics_Summary.csv"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/Metrics_Summary"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#Metrics.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_metrics.py"
            ],
            "class": "CommandLineTool",
            "id": "#Metrics.cwl"
        },
        {
            "inputs": [
                {
                    "default": "Targeted",
                    "type": "string",
                    "id": "#MultiplexingSettings.cwl/Assay"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#MultiplexingSettings.cwl/_Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/_Tag_Sample_Names"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#MultiplexingSettings.cwl/Sample_Tags_Version"
                },
                {
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#MultiplexingSettings.cwl/Tag_Sample_Names"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  var enumifiedSampleTagsVersion = null;\n  if (inputs._Sample_Tags_Version) {\n  var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();\n  if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')\n  {\n    enumifiedSampleTagsVersion = 'hs';\n  }\n  else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')\n  {\n    enumifiedSampleTagsVersion = 'mm';\n  }\n  else if (_Sample_Tags_Version === 'no multiplexing')\n  {\n    enumifiedSampleTagsVersion = null;\n  }\n  else\n  {\n    throw new Error(\"Cannot parse Sample Tag Version: \" + inputs._Sample_Tags_Version);\n  }\n  }\n  return ({\n  Tag_Sample_Names: inputs._Tag_Sample_Names,\n  Sample_Tags_Version: enumifiedSampleTagsVersion\n  });\n}",
            "id": "#MultiplexingSettings.cwl"
        },
        {
            "inputs": [
                {
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/filename"
                                },
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/readFlag"
                                },
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/readPairId"
                                },
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/FastqReadPairs/library"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#PairReadFiles.cwl/FastqReadPairs"
                },
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#PairReadFiles.cwl/Reads"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "PairReadFiles takes an array of split files and pairs them, such that an R1 file is transferred to the QualityFilter with its corresponding R2 file.\nThe original FASTQ files are paired in CheckFastqs and then split and sub-sampled in SplitAndSubsample. The pairing information is taken from CheckFastqs.\n",
            "id": "#PairReadFiles.cwl",
            "outputs": [
                {
                    "type": {
                        "items": {
                            "fields": [
                                {
                                    "type": "File",
                                    "name": "#PairReadFiles.cwl/ReadPairs/R1"
                                },
                                {
                                    "type": "File",
                                    "name": "#PairReadFiles.cwl/ReadPairs/R2"
                                },
                                {
                                    "type": "int",
                                    "name": "#PairReadFiles.cwl/ReadPairs/readPairId"
                                },
                                {
                                    "type": "string",
                                    "name": "#PairReadFiles.cwl/ReadPairs/library"
                                }
                            ],
                            "type": "record"
                        },
                        "type": "array"
                    },
                    "id": "#PairReadFiles.cwl/ReadPairs"
                }
            ],
            "expression": "${\n  // use the CheckFastqs read pairing information to create a dictionary\n  // using the original fastq file name without the extension as the key\n  var fastqReadPairs = {}\n  for (var i = 0; i < inputs.FastqReadPairs.length; i++) {\n    var fileDict = inputs.FastqReadPairs[i];\n    var filename = fileDict[\"filename\"];\n\n    if (!fastqReadPairs[filename]) {\n      fastqReadPairs[filename] = {\n        readPairId: null,\n        readFlag: null,\n        library: null,\n      };\n    }\n\n    fastqReadPairs[filename].readPairId = fileDict[\"readPairId\"]\n    fastqReadPairs[filename].readFlag = fileDict[\"readFlag\"]\n    fastqReadPairs[filename].library = fileDict[\"library\"]\n  }\n\n  // now loop through the input read files which could\n  // be the original fastq files if no sub-sampling has\n  // been done, or the sub-sampled fastq files\n  var readPairs = {}\n  for (var i = 0; i < inputs.Reads.length; i++) {\n\n    // Get the fastq file\n    var f = inputs.Reads[i];\n\n    // Split on the dash to get the name of the original file\n    // and the chunk id (if it exists)\n    // We would like to ignore the case of the .fastq.gz or .fq.gz\n    // at the end of the file. JS allows one to ignore the case of\n    // an entire RegEx by adding an 'i' to the end of the pattern.\n    // Unfortunately, JS does not allow one to ignore the case for\n    // a specific RegEx group. This is one way to get around that:\n    var groups = f.basename.match(/^(.*?)(-[0-9]*)?(\\.([fF][aA][sS][tT][qQ]|[fF][qQ])\\.[gG][zZ])$/);\n\n    // If the RegEx fails, create an error\n    if (groups === undefined || groups === null) {\n      throw new Error(\"The RegEx for the fastq file name '\" + f.basename + \"' is failing.\");\n    }\n\n    // Get the base name, chunk id, and file extension\n    // The base name without the chunk id and file\n    // extension is the key from CheckFastqs\n    // The chunk id is used later to create a new unique\n    // read pair id for all fastq files (sub-sampled or not)\n    var basename = groups[1];\n    var orgChunkId = groups[2];\n    // if there is no chunk id, use an arbitrary number\n    var chunkId = 9999;\n    if (orgChunkId) {\n      // slice off the '-' and cast to an integer\n      chunkId = parseInt(orgChunkId.slice(1));\n    }\n    // double check that we have a chunk id\n    if (chunkId === undefined || chunkId === null) {\n      throw new Error(\"The fastq file sub-sampling id could not be determined!\");\n    }\n    var fileExt = groups[3];\n\n    // The basename without the chunk id and file extension\n    // should match the original file name from CheckFastqs\n    // The original file name from CheckFastqs is the key for\n    // the dictionary containing the original unique pair id\n    var filename = basename;\n    var fileDict = fastqReadPairs[filename];\n\n    // If the fileDict for this filename is not found, then try to use\n    // the original filename without the file extension as the key\n    if (fileDict === undefined || fileDict === null) {\n      // If the original filename ends in (-[0-9]*)\n      // and no sub-sampling occurs, then try to use the\n      // original filename without the extension as the key\n      var groups = f.basename.match(/^(.*?)(\\.([fF][aA][sS][tT][qQ]|[fF][qQ])\\.[gG][zZ])$/);\n\n      // If the RegEx fails, create an error\n      if (groups === undefined || groups === null) {\n        throw new Error(\"The RegEx for the fastq file name '\" + f.basename + \"' is failing.\");\n      }\n\n      // Get the base name and file extension\n      // The base name without the file extension\n      // is the key from CheckFastqs\n      var basename = groups[1];\n      var fileExt = groups[2];\n\n      var fileDict = fastqReadPairs[basename];\n\n      // If the fileDict for this filename is still not found,\n      // then the filenames are in an unexpected format and\n      // the RegEx above needs to be modified to create\n      // filenames formatted in the same way as CheckFastqs\n      // Create an error\n      if (fileDict === undefined || fileDict === null) {\n        throw new Error(\"Cannot find the fastq read pair information for '\" + filename + \"'.\");\n      }\n    }\n\n    // Get the pairing information from CheckFastqs\n    var readPairId = fileDict[\"readPairId\"];\n    var library = fileDict[\"library\"];\n    var flag = fileDict[\"readFlag\"];\n\n    // Add the chunkId to create a new unique read pair id\n    // for each file (sub-sampled or not)\n    var chunkReadPairId = readPairId + \"_\" + chunkId;\n\n    // Create a dictionary for each pair of files\n    if (!readPairs[chunkReadPairId]) {\n      readPairs[chunkReadPairId] = {\n        R1: null,\n        R2: null,\n        library: library,\n        readPairId: null,\n      };\n    }\n    // add in the R1 and R2 files, depending on the flag\n    if (flag === \"R1\") {\n      readPairs[chunkReadPairId].R1 = f\n    } else if (flag === \"R2\") {\n      readPairs[chunkReadPairId].R2 = f\n    }\n  }\n  // we are not interested in the read pair ids in readPairs\n  // flatten into an array of objects\n  var readPairsList = [];\n  var i = 1;\n  for (var key in readPairs) {\n    if (readPairs.hasOwnProperty(key)) {\n      var readPair = readPairs[key];\n      readPair.readPairId = i;\n      readPairsList.push(readPair);\n      i++;\n    }\n  }\n  // pass this array to the record array named \"ReadPairs\" on the CWL layer\n  return {ReadPairs: readPairsList}\n}",
            "class": "ExpressionTool"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Basic_Algo_Only"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/_Exact_Cell_Count"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#PutativeCellSettings.cwl/Basic_Algo_Only"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#PutativeCellSettings.cwl/Exact_Cell_Count"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  if (inputs._Exact_Cell_Count) {\n    if (inputs._Exact_Cell_Count < 1) {\n      throw(\"Illogical value for exact cell count: \" + inputs._Exact_Cell_Count);\n    }\n  }\n  return ({\n    Exact_Cell_Count: inputs._Exact_Cell_Count,\n    Basic_Algo_Only: inputs._Basic_Algo_Only,\n  });\n}",
            "id": "#PutativeCellSettings.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "prefix": "--label-version"
                    },
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#QualityFilter.cwl/Label_Version"
                },
                {
                    "inputBinding": {
                        "prefix": "--read-filter-off"
                    },
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#QualityFilter.cwl/Read_Filter_Off"
                },
                {
                    "type": {
                        "fields": [
                            {
                                "inputBinding": {
                                    "prefix": "--r1"
                                },
                                "type": "File",
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/R1"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--r2"
                                },
                                "type": "File",
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/R2"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--read-pair-id"
                                },
                                "type": "int",
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/readPairId"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--library"
                                },
                                "type": "string",
                                "name": "#QualityFilter.cwl/Split_Read_Pairs/library"
                            }
                        ],
                        "type": "record"
                    },
                    "id": "#QualityFilter.cwl/Split_Read_Pairs"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*read_quality.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#QualityFilter.cwl/Filter_Metrics"
                },
                {
                    "outputBinding": {
                        "glob": "*_R1_.fastq.gz"
                    },
                    "type": "File",
                    "id": "#QualityFilter.cwl/R1"
                },
                {
                    "outputBinding": {
                        "glob": "*_R2_.fastq.gz"
                    },
                    "type": "File",
                    "id": "#QualityFilter.cwl/R2"
                },
                {
                    "outputBinding": {
                        "glob": "*.log"
                    },
                    "type": "File",
                    "id": "#QualityFilter.cwl/output"
                }
            ],
            "baseCommand": [
                "mist_quality_filter.py"
            ],
            "class": "CommandLineTool",
            "id": "#QualityFilter.cwl"
        },
        {
            "inputs": [
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#SplitAndSubsample.cwl/Fastqs"
                },
                {
                    "type": {
                        "items": "string",
                        "type": "array"
                    },
                    "id": "#SplitAndSubsample.cwl/FilesToSkipSplitAndSubsample"
                },
                {
                    "type": [
                        "null",
                        "long"
                    ],
                    "id": "#SplitAndSubsample.cwl/NumRecordsPerSplit"
                },
                {
                    "type": "float",
                    "id": "#SplitAndSubsample.cwl/SubsampleRatio"
                },
                {
                    "type": "int",
                    "id": "#SplitAndSubsample.cwl/SubsampleSeed"
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "SplitAndSubsample splits, subsamples and formats read files to be deposited in QualityFilter.\n",
            "id": "#SplitAndSubsample.cwl",
            "steps": [
                {
                    "doc": "After scattering \"SplitAndSubsample\" on a File array, the output of each node is also an array. Thus, we are left with a nestled list. This JS expression flattens this list to deal with the split reads in PairReadFiles.cwl",
                    "out": [
                        "#SplitAndSubsample.cwl/FlattenOutput/SplitFastqList"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": {
                                    "items": {
                                        "items": "File",
                                        "type": "array"
                                    },
                                    "type": "array"
                                },
                                "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output/nestledSplitFastqList"
                            }
                        ],
                        "outputs": [
                            {
                                "type": {
                                    "items": "File",
                                    "type": "array"
                                },
                                "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output/SplitFastqList"
                            }
                        ],
                        "class": "ExpressionTool",
                        "expression": "${\n  return {SplitFastqList: [].concat.apply([], inputs.nestledSplitFastqList)}\n}\n",
                        "id": "#SplitAndSubsample.cwl/FlattenOutput/flatten_output"
                    },
                    "id": "#SplitAndSubsample.cwl/FlattenOutput",
                    "in": [
                        {
                            "source": "#SplitAndSubsample.cwl/SplitAndSubsample/SplitAndSubsampledFastqs",
                            "id": "#SplitAndSubsample.cwl/FlattenOutput/nestledSplitFastqList"
                        }
                    ]
                },
                {
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "inputBinding": {
                                    "prefix": "--fastq-file-path"
                                },
                                "type": "File",
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/Fastq"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--files-to-skip-split-and-subsample",
                                    "itemSeparator": ","
                                },
                                "type": {
                                    "items": "string",
                                    "type": "array"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/FilesToSkipSplitAndSubsample"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--num-records"
                                },
                                "type": [
                                    "null",
                                    "long"
                                ],
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/NumRecordsPerSplit"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--subsample-ratio"
                                },
                                "type": "float",
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SubsampleRatio"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--subsample-seed"
                                },
                                "type": "int",
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SubsampleSeed"
                            }
                        ],
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.9.1",
                                "class": "DockerRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "*.fastq.gz",
                                    "outputEval": "${ if (self.length === 0) { return [inputs.Fastq]; } else { return self; } }"
                                },
                                "type": {
                                    "items": "File",
                                    "type": "array"
                                },
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/SplitAndSubsampledFastqs"
                            },
                            {
                                "outputBinding": {
                                    "glob": "*.log"
                                },
                                "type": "File",
                                "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq/log"
                            }
                        ],
                        "baseCommand": [
                            "mist_split_fastq.py"
                        ],
                        "class": "CommandLineTool",
                        "id": "#SplitAndSubsample.cwl/SplitAndSubsample/split_fastq"
                    },
                    "doc": "Allocate one docker/python process per file to do the actual file splitting.",
                    "scatter": [
                        "#SplitAndSubsample.cwl/SplitAndSubsample/Fastq"
                    ],
                    "in": [
                        {
                            "source": "#SplitAndSubsample.cwl/Fastqs",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/Fastq"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/FilesToSkipSplitAndSubsample",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/FilesToSkipSplitAndSubsample"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/NumRecordsPerSplit",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/NumRecordsPerSplit"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/SubsampleRatio",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/SubsampleRatio"
                        },
                        {
                            "source": "#SplitAndSubsample.cwl/SubsampleSeed",
                            "id": "#SplitAndSubsample.cwl/SplitAndSubsample/SubsampleSeed"
                        }
                    ],
                    "id": "#SplitAndSubsample.cwl/SplitAndSubsample",
                    "out": [
                        "#SplitAndSubsample.cwl/SplitAndSubsample/SplitAndSubsampledFastqs",
                        "#SplitAndSubsample.cwl/SplitAndSubsample/log"
                    ]
                }
            ],
            "outputs": [
                {
                    "outputSource": "#SplitAndSubsample.cwl/FlattenOutput/SplitFastqList",
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#SplitAndSubsample.cwl/SplitAndSubsampledFastqs"
                },
                {
                    "outputSource": "#SplitAndSubsample.cwl/SplitAndSubsample/log",
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#SplitAndSubsample.cwl/log"
                }
            ],
            "class": "Workflow"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#SubsampleSettings.cwl/_Subsample_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#SubsampleSettings.cwl/_Subsample_Seed"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#SubsampleSettings.cwl/Subsample_Reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#SubsampleSettings.cwl/Subsample_Seed"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  var subsamplingOutputs = {\n    Subsample_Reads: inputs._Subsample_Reads,\n    Subsample_Seed: inputs._Subsample_Seed\n  }\n  return subsamplingOutputs;\n}",
            "id": "#SubsampleSettings.cwl"
        },
        {
            "inputs": [
                {
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#UncompressDatatables.cwl/Compressed_Data_Table"
                },
                {
                    "type": "File",
                    "id": "#UncompressDatatables.cwl/Compressed_Expression_Matrix"
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ],
            "outputs": [
                {
                    "outputSource": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompressed_File",
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#UncompressDatatables.cwl/Uncompressed_Data_Tables"
                },
                {
                    "outputSource": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompressed_File",
                    "type": "File",
                    "id": "#UncompressDatatables.cwl/Uncompressed_Expression_Matrix"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "id": "#UncompressDatatables.cwl/Uncompress_Datatable",
                    "out": [
                        "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompressed_File"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "inputBinding": {
                                    "position": 1
                                },
                                "type": "File",
                                "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner/Compressed_File"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "stdout": "$(inputs.Compressed_File.nameroot)",
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "$(inputs.Compressed_File.nameroot)"
                                },
                                "type": "File",
                                "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner/Uncompressed_File"
                            }
                        ],
                        "baseCommand": [
                            "gunzip"
                        ],
                        "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Uncompress_Datatable_Inner",
                        "arguments": [
                            {
                                "position": 0,
                                "valueFrom": "-c"
                            }
                        ],
                        "class": "CommandLineTool",
                        "hints": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.9.1",
                                "class": "DockerRequirement"
                            }
                        ]
                    },
                    "scatter": [
                        "#UncompressDatatables.cwl/Uncompress_Datatable/Compressed_File"
                    ],
                    "in": [
                        {
                            "source": "#UncompressDatatables.cwl/Compressed_Data_Table",
                            "id": "#UncompressDatatables.cwl/Uncompress_Datatable/Compressed_File"
                        }
                    ]
                },
                {
                    "out": [
                        "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompressed_File"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "inputBinding": {
                                    "position": 1
                                },
                                "type": "File",
                                "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner/Compressed_File"
                            }
                        ],
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "stdout": "$(inputs.Compressed_File.nameroot)",
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "$(inputs.Compressed_File.nameroot)"
                                },
                                "type": "File",
                                "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner/Uncompressed_File"
                            }
                        ],
                        "baseCommand": [
                            "gunzip"
                        ],
                        "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Uncompress_Expression_Matrix_Inner",
                        "arguments": [
                            {
                                "position": 0,
                                "valueFrom": "-c"
                            }
                        ],
                        "class": "CommandLineTool",
                        "hints": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.9.1",
                                "class": "DockerRequirement"
                            }
                        ]
                    },
                    "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix",
                    "in": [
                        {
                            "source": "#UncompressDatatables.cwl/Compressed_Expression_Matrix",
                            "id": "#UncompressDatatables.cwl/Uncompress_Expression_Matrix/Compressed_File"
                        }
                    ]
                }
            ],
            "id": "#UncompressDatatables.cwl"
        },
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--sample-name"
                    },
                    "type": "string",
                    "id": "#VDJ_Annotate_Molecules.cwl/Sample_Name"
                },
                {
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--cell-type-mapping-fp"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/cellTypeMapping"
                },
                {
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--ignore",
                        "itemSeparator": ","
                    },
                    "type": [
                        "null",
                        {
                            "items": "string",
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/chainsToIgnore"
                },
                {
                    "inputBinding": {
                        "position": 8,
                        "prefix": "--e-value-for-j"
                    },
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/evalueJgene"
                },
                {
                    "inputBinding": {
                        "position": 7,
                        "prefix": "--e-value-for-v"
                    },
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/evalueVgene"
                },
                {
                    "inputBinding": {
                        "position": 5
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/igCalls"
                },
                {
                    "inputBinding": {
                        "position": 9,
                        "prefix": "--metadata-fp"
                    },
                    "type": "File",
                    "id": "#VDJ_Annotate_Molecules.cwl/metadata"
                },
                {
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--putative-cells-json-fp"
                    },
                    "type": "File",
                    "id": "#VDJ_Annotate_Molecules.cwl/putativeCells"
                },
                {
                    "inputBinding": {
                        "position": 6
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/tcrCalls"
                },
                {
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--vdj-version"
                    },
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjVersion"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "bdgenomics/rhapsody:1.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "*_VDJ_perCellChain_unfiltered.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjCellChainDatatableUnfiltered"
                },
                {
                    "doc": "VDJ data per cell, with distribution based error correction",
                    "outputBinding": {
                        "glob": "*_VDJ_perCell.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjCellsDatatable"
                },
                {
                    "doc": "VDJ data per cell, cell type error correction",
                    "outputBinding": {
                        "glob": "*_VDJ_perCell_cellType_corrected.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjCellsDatatableCellCorrected"
                },
                {
                    "doc": "VDJ data per cell, DBEC and cell type error correction",
                    "outputBinding": {
                        "glob": "*_VDJ_perCell_DBEC_cellType_corrected.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjCellsDatatableDBECCellCorrected"
                },
                {
                    "doc": "VDJ data per cell, including non-putative cells, no error correction applied",
                    "outputBinding": {
                        "glob": "*_VDJ_perCell_unfiltered.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjCellsDatatableUnfiltered"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_readsInvalid.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjInvalidReadsDatatable"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_metrics.csv"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjMetricsCsv"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_metrics.json"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjMetricsJson"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_molecules_per_cell_and_chain_summary_boxplot.png"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjReadsAndMoleculesPerCellFigure"
                },
                {
                    "outputBinding": {
                        "glob": "*_DBEC_cutoff.png"
                    },
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjReadsPerCellByChainTypeFigure"
                },
                {
                    "outputBinding": {
                        "glob": "*_VDJ_readsValid.csv.gz"
                    },
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Molecules.cwl/vdjValidReadsDatatable"
                }
            ],
            "baseCommand": [
                "mist_annotate_molecules_vdj.py"
            ],
            "id": "#VDJ_Annotate_Molecules.cwl",
            "class": "CommandLineTool",
            "hints": [
                {
                    "ramMin": 64000,
                    "class": "ResourceRequirement"
                }
            ]
        },
        {
            "inputs": [
                {
                    "doc": ".fasta.gz",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Reads.cwl/Cdr3QueryFasta"
                },
                {
                    "type": "string",
                    "id": "#VDJ_Annotate_Reads.cwl/vdjType"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Annotate_Reads.cwl/vdjVersion"
                }
            ],
            "outputs": [
                {
                    "outputSource": "#VDJ_Annotate_Reads.cwl/PrunePyIR/PrunedPyIROutput",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_Annotate_Reads.cwl/Cdr3Call"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "out": [
                        "#VDJ_Annotate_Reads.cwl/CallCdr3/Cdr3Call"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_Annotate_Reads.cwl/CallCdr3/CallCdr3Inner/Cdr3QueryFasta"
                            },
                            {
                                "type": "string",
                                "id": "#VDJ_Annotate_Reads.cwl/CallCdr3/CallCdr3Inner/vdjType"
                            },
                            {
                                "type": [
                                    "null",
                                    "string"
                                ],
                                "id": "#VDJ_Annotate_Reads.cwl/CallCdr3/CallCdr3Inner/vdjVersion"
                            }
                        ],
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.9.1",
                                "class": "DockerRequirement"
                            },
                            {
                                "class": "InlineJavascriptRequirement"
                            },
                            {
                                "class": "ShellCommandRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "*.json.gz",
                                    "outputEval": "${\n  if (inputs.vdjVersion && inputs.Cdr3QueryFasta && self.size == 0) {\n    throw(\"No outputs from PyIR detected!\");\n  } else {\n    return(self);\n  }\n}"
                                },
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_Annotate_Reads.cwl/CallCdr3/CallCdr3Inner/Cdr3Call"
                            }
                        ],
                        "baseCommand": [
                            "mist_pyirWrapper.py"
                        ],
                        "class": "CommandLineTool",
                        "arguments": [
                            {
                                "prefix": "-r",
                                "valueFrom": "$(inputs.vdjType)"
                            },
                            {
                                "prefix": "--strand",
                                "valueFrom": "plus"
                            },
                            {
                                "prefix": "--database",
                                "valueFrom": "/mist/pyir_data"
                            },
                            {
                                "prefix": "-f",
                                "valueFrom": "json"
                            },
                            {
                                "prefix": "-m",
                                "valueFrom": "1"
                            },
                            {
                                "prefix": "-s",
                                "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"human\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"mouse\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                            },
                            {
                                "prefix": "-o",
                                "valueFrom": "${\n  if(inputs.Cdr3QueryFasta){\n    return(inputs.Cdr3QueryFasta.nameroot.split(\".\")[0]);\n  } else {\n    return(\"NA\");\n  }\n}"
                            },
                            {
                                "prefix": "-winput",
                                "valueFrom": "${\n  if (!inputs.vdjType) {\n    return(\"~/deliberatelyNotAQueryFastaToInduceFailure.fasta\");\n  } else {\n    return(inputs.Cdr3QueryFasta);\n  }\n}"
                            },
                            {
                                "shellQuote": false,
                                "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                            }
                        ],
                        "id": "#VDJ_Annotate_Reads.cwl/CallCdr3/CallCdr3Inner"
                    },
                    "id": "#VDJ_Annotate_Reads.cwl/CallCdr3",
                    "in": [
                        {
                            "source": "#VDJ_Annotate_Reads.cwl/CallConstantRegion/ConstantRegionCall",
                            "id": "#VDJ_Annotate_Reads.cwl/CallCdr3/Cdr3QueryFasta"
                        },
                        {
                            "source": "#VDJ_Annotate_Reads.cwl/vdjType",
                            "id": "#VDJ_Annotate_Reads.cwl/CallCdr3/vdjType"
                        },
                        {
                            "source": "#VDJ_Annotate_Reads.cwl/vdjVersion",
                            "id": "#VDJ_Annotate_Reads.cwl/CallCdr3/vdjVersion"
                        }
                    ]
                },
                {
                    "out": [
                        "#VDJ_Annotate_Reads.cwl/CallConstantRegion/ConstantRegionCall"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_Annotate_Reads.cwl/CallConstantRegion/CallConstantRegionInner/Cdr3QueryFasta"
                            },
                            {
                                "type": [
                                    "null",
                                    "string"
                                ],
                                "id": "#VDJ_Annotate_Reads.cwl/CallConstantRegion/CallConstantRegionInner/vdjVersion"
                            }
                        ],
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.9.1",
                                "class": "DockerRequirement"
                            },
                            {
                                "class": "InlineJavascriptRequirement"
                            },
                            {
                                "class": "ShellCommandRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "*_constant_region_called.fasta.gz",
                                    "outputEval": "${\n  if (!inputs.vdjVersion) {\n    return(null);\n  } else {\n    return(self);\n  }\n}"
                                },
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_Annotate_Reads.cwl/CallConstantRegion/CallConstantRegionInner/ConstantRegionCall"
                            }
                        ],
                        "baseCommand": [
                            "bowtie2"
                        ],
                        "class": "CommandLineTool",
                        "arguments": [
                            "--quiet",
                            "--no-head",
                            "--local",
                            "-p",
                            "1",
                            "-L",
                            "10",
                            "-N",
                            "1",
                            "--ma",
                            "4",
                            {
                                "prefix": "-f",
                                "valueFrom": "$(inputs.Cdr3QueryFasta)"
                            },
                            {
                                "prefix": "-x",
                                "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/humanVDJCidx\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/mouseVDJCidx\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                            },
                            {
                                "shellQuote": false,
                                "valueFrom": "|"
                            },
                            "awk",
                            {
                                "shellQuote": false,
                                "valueFrom": "${\n  if(inputs.Cdr3QueryFasta) {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10 |  \\\" gzip >> \" + inputs.Cdr3QueryFasta.nameroot.split(\".\")[0] + \"_constant_region_called.fasta.gz\\\"}\\'\");\n  } else {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10}\\'\");\n  }\n}"
                            },
                            {
                                "shellQuote": false,
                                "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                            }
                        ],
                        "id": "#VDJ_Annotate_Reads.cwl/CallConstantRegion/CallConstantRegionInner"
                    },
                    "id": "#VDJ_Annotate_Reads.cwl/CallConstantRegion",
                    "in": [
                        {
                            "source": "#VDJ_Annotate_Reads.cwl/Cdr3QueryFasta",
                            "id": "#VDJ_Annotate_Reads.cwl/CallConstantRegion/Cdr3QueryFasta"
                        },
                        {
                            "source": "#VDJ_Annotate_Reads.cwl/vdjVersion",
                            "id": "#VDJ_Annotate_Reads.cwl/CallConstantRegion/vdjVersion"
                        }
                    ]
                },
                {
                    "out": [
                        "#VDJ_Annotate_Reads.cwl/PrunePyIR/PrunedPyIROutput"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "inputBinding": {
                                    "position": 0
                                },
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_Annotate_Reads.cwl/PrunePyIR/PrunePyIRInner/PyIROutput"
                            }
                        ],
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.9.1",
                                "class": "DockerRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "*_pruned.csv.gz"
                                },
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_Annotate_Reads.cwl/PrunePyIR/PrunePyIRInner/PrunedPyIROutput"
                            }
                        ],
                        "baseCommand": [
                            "mist_prune_pyir.py"
                        ],
                        "class": "CommandLineTool",
                        "id": "#VDJ_Annotate_Reads.cwl/PrunePyIR/PrunePyIRInner"
                    },
                    "id": "#VDJ_Annotate_Reads.cwl/PrunePyIR",
                    "in": [
                        {
                            "source": "#VDJ_Annotate_Reads.cwl/CallCdr3/Cdr3Call",
                            "id": "#VDJ_Annotate_Reads.cwl/PrunePyIR/PyIROutput"
                        }
                    ]
                }
            ],
            "id": "#VDJ_Annotate_Reads.cwl",
            "hints": [
                {
                    "ramMax": 2000,
                    "class": "ResourceRequirement",
                    "coresMax": 1
                }
            ]
        },
        {
            "inputs": [
                {
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_GatherCalls.cwl/theCalls"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "VDJ_GatherCalls collect the outputs from the multi-processed VDJ step into one file.\n",
            "id": "#VDJ_GatherCalls.cwl",
            "steps": [
                {
                    "out": [
                        "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gatheredCalls"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "type": [
                                    {
                                        "items": [
                                            "null",
                                            "File"
                                        ],
                                        "type": "array"
                                    }
                                ],
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR/theCalls"
                            }
                        ],
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.9.1",
                                "class": "DockerRequirement"
                            },
                            {
                                "class": "InlineJavascriptRequirement"
                            },
                            {
                                "class": "ShellCommandRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "*_constant_region_called_pruned.csv.gz",
                                    "outputEval": "${\n  if (self.size == 0) {\n    throw(\"No outputs from PyIR detected in VDJ_GatherCalls!\");\n  } else {\n    return(self);\n  }\n}"
                                },
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR/gatheredCalls"
                            }
                        ],
                        "class": "CommandLineTool",
                        "arguments": [
                            {
                                "shellQuote": false,
                                "valueFrom": "${\n  if (!inputs.theCalls[0] ) {\n    return (\"echo \\\"No outputs from PyIR detected in VDJ_GatherCalls\\\"\")\n  }\n  var inputFiles = \"\"\n  if (!inputs.theCalls[0].path.split(\"_PrunePyIR\")[1]){\n    inputFiles = \"zcat\"\n    for (var i = 0; i < inputs.theCalls.length; i++) {\n      inputFiles += \" \" + inputs.theCalls[i].path\n    }\n    inputFiles += \" | \"\n  } else {\n    inputFiles = \"zcat \" + inputs.theCalls[0].path.split(\"VDJ\")[0] + \"*\" + inputs.theCalls[0].path.split(\"_PrunePyIR\")[1].split(\"_Number_\")[0] + \"_Number_*.csv.gz | \"\n  }\n  var outputFileName = \"\\\"gzip > \" + inputs.theCalls[0].nameroot.split(\"_Number_\")[0] + \"_constant_region_called_pruned.csv.gz\" + \"\\\"\"\n  var awkCommand =  \"awk \\'NR==1{F=$1;print | \" + outputFileName + \" } $1!=F { print | \" + outputFileName + \" }\\' \"\n  var outputCommand = inputFiles + awkCommand\n  return (outputCommand)\n}"
                            }
                        ],
                        "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gather_PyIR"
                    },
                    "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls",
                    "in": [
                        {
                            "source": "#VDJ_GatherCalls.cwl/theCalls",
                            "id": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/theCalls"
                        }
                    ]
                }
            ],
            "outputs": [
                {
                    "outputSource": "#VDJ_GatherCalls.cwl/VDJ_GatherCalls/gatheredCalls",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_GatherCalls.cwl/gatheredCalls"
                }
            ],
            "class": "Workflow"
        },
        {
            "inputs": [],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_JGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_VGene_Evalue"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_Settings.cwl/VDJ_Version"
                }
            ],
            "class": "ExpressionTool",
            "expression": "${\n  var vdjVersion = null;\n  if (!inputs._VDJ_Version) {\n    vdjVersion = null;}\n  else {\n    var _VDJ_Version = inputs._VDJ_Version.toLowerCase();\n    if (_VDJ_Version === \"human\" || _VDJ_Version === \"hs\" || _VDJ_Version === \"human vdj - bcr and tcr\") {\n      vdjVersion = \"human\";\n    } else if (_VDJ_Version === \"humanbcr\" || _VDJ_Version === \"human vdj - bcr only\") {\n      vdjVersion = \"humanBCR\";\n    } else if (_VDJ_Version === \"humantcr\" || _VDJ_Version === \"human vdj - tcr only\") {\n      vdjVersion = \"humanTCR\";\n    } else if (_VDJ_Version === \"mouse\" || _VDJ_Version === \"mm\" || _VDJ_Version === \"mouse vdj - bcr and tcr\") {\n      vdjVersion = \"mouse\";\n    } else if (_VDJ_Version === \"mousebcr\" || _VDJ_Version === \"mouse vdj - bcr only\") {\n      vdjVersion = \"mouseBCR\";\n    } else if (_VDJ_Version === \"humantcr\" || _VDJ_Version === \"mouse vdj - tcr only\") {\n      vdjVersion = \"mouseTCR\";\n    } else {\n      vdjVersion = inputs._VDJ_Version;\n    }\n  }\n\n  return ({\n  VDJ_Version: vdjVersion,\n  })\n}",
            "id": "#VDJ_Settings.cwl"
        },
        {
            "inputs": [
                {
                    "default": 36,
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#VDJ_SplitValidReads.cwl/num_fasta"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_SplitValidReads.cwl/validReads"
                }
            ],
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "doc": "VDJ_SplitValidReads splits fasta files to be multi-processed in the VDJ step.\n",
            "id": "#VDJ_SplitValidReads.cwl",
            "steps": [
                {
                    "out": [
                        "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/fastaList",
                        "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/numFiles",
                        "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/log"
                    ],
                    "run": {
                        "cwlVersion": "v1.0",
                        "inputs": [
                            {
                                "inputBinding": {
                                    "prefix": "--num-fasta"
                                },
                                "type": [
                                    "null",
                                    "int"
                                ],
                                "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/split_fasta/num_fasta"
                            },
                            {
                                "inputBinding": {
                                    "prefix": "--fasta-file-path"
                                },
                                "type": [
                                    "null",
                                    "File"
                                ],
                                "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/split_fasta/validReads"
                            }
                        ],
                        "requirements": [
                            {
                                "dockerPull": "bdgenomics/rhapsody:1.9.1",
                                "class": "DockerRequirement"
                            }
                        ],
                        "outputs": [
                            {
                                "outputBinding": {
                                    "glob": "*_split.fasta.gz",
                                    "outputEval": "${ if (self.length === 0) { return [inputs.validReads]; } else { return self; } }"
                                },
                                "type": [
                                    {
                                        "items": [
                                            "null",
                                            "File"
                                        ],
                                        "type": "array"
                                    }
                                ],
                                "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/split_fasta/fastaList"
                            },
                            {
                                "outputBinding": {
                                    "glob": "*.log"
                                },
                                "type": {
                                    "items": "File",
                                    "type": "array"
                                },
                                "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/split_fasta/log"
                            },
                            {
                                "outputBinding": {
                                    "glob": "*_split.fasta.gz",
                                    "outputEval": "${ return(parseInt(self.length)); }"
                                },
                                "type": "int",
                                "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/split_fasta/numFiles"
                            }
                        ],
                        "baseCommand": [
                            "mist_split_fasta.py"
                        ],
                        "class": "CommandLineTool",
                        "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/split_fasta"
                    },
                    "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads",
                    "in": [
                        {
                            "source": "#VDJ_SplitValidReads.cwl/num_fasta",
                            "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/num_fasta"
                        },
                        {
                            "source": "#VDJ_SplitValidReads.cwl/validReads",
                            "id": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/validReads"
                        }
                    ]
                }
            ],
            "outputs": [
                {
                    "outputSource": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/fastaList",
                    "type": [
                        {
                            "items": [
                                "null",
                                "File"
                            ],
                            "type": "array"
                        }
                    ],
                    "id": "#VDJ_SplitValidReads.cwl/SplitFastaList"
                },
                {
                    "outputSource": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/log",
                    "type": {
                        "items": "File",
                        "type": "array"
                    },
                    "id": "#VDJ_SplitValidReads.cwl/log"
                },
                {
                    "outputSource": "#VDJ_SplitValidReads.cwl/VDJ_SplitValidReads/numFiles",
                    "type": "int",
                    "id": "#VDJ_SplitValidReads.cwl/numFiles"
                }
            ],
            "class": "Workflow"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_ig.cwl/VDJ_Version"
                },
                {
                    "type": "int",
                    "id": "#VDJ_ig.cwl/numFiles"
                },
                {
                    "doc": ".fasta.gz",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_ig.cwl/validReadsIg"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "outputs": [
                {
                    "outputSource": "#VDJ_ig.cwl/CallCdr3ForIgs/Cdr3Call",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_ig.cwl/igCalls"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "out": [
                        "#VDJ_ig.cwl/CallCdr3ForIgs/Cdr3Call"
                    ],
                    "in": [
                        {
                            "source": "#VDJ_ig.cwl/validReadsIg",
                            "id": "#VDJ_ig.cwl/CallCdr3ForIgs/Cdr3QueryFasta"
                        },
                        {
                            "valueFrom": "Ig",
                            "id": "#VDJ_ig.cwl/CallCdr3ForIgs/vdjType"
                        },
                        {
                            "source": "#VDJ_ig.cwl/VDJ_Version",
                            "id": "#VDJ_ig.cwl/CallCdr3ForIgs/vdjVersion"
                        }
                    ],
                    "run": "#VDJ_Annotate_Reads.cwl",
                    "id": "#VDJ_ig.cwl/CallCdr3ForIgs",
                    "hints": [
                        {
                            "coresMin": "$(inputs.numFiles)",
                            "class": "ResourceRequirement"
                        }
                    ]
                }
            ],
            "id": "#VDJ_ig.cwl"
        },
        {
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#VDJ_tcr.cwl/VDJ_Version"
                },
                {
                    "type": "int",
                    "id": "#VDJ_tcr.cwl/numFiles"
                },
                {
                    "doc": ".fasta.gz",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_tcr.cwl/validReadsTcr"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "outputs": [
                {
                    "outputSource": "#VDJ_tcr.cwl/CallCdr3ForTcrs/Cdr3Call",
                    "type": [
                        "null",
                        "File"
                    ],
                    "id": "#VDJ_tcr.cwl/tcrCalls"
                }
            ],
            "class": "Workflow",
            "steps": [
                {
                    "out": [
                        "#VDJ_tcr.cwl/CallCdr3ForTcrs/Cdr3Call"
                    ],
                    "in": [
                        {
                            "source": "#VDJ_tcr.cwl/validReadsTcr",
                            "id": "#VDJ_tcr.cwl/CallCdr3ForTcrs/Cdr3QueryFasta"
                        },
                        {
                            "valueFrom": "TCR",
                            "id": "#VDJ_tcr.cwl/CallCdr3ForTcrs/vdjType"
                        },
                        {
                            "source": "#VDJ_tcr.cwl/VDJ_Version",
                            "id": "#VDJ_tcr.cwl/CallCdr3ForTcrs/vdjVersion"
                        }
                    ],
                    "run": "#VDJ_Annotate_Reads.cwl",
                    "id": "#VDJ_tcr.cwl/CallCdr3ForTcrs",
                    "hints": [
                        {
                            "coresMin": "$(inputs.numFiles)",
                            "class": "ResourceRequirement"
                        }
                    ]
                }
            ],
            "id": "#VDJ_tcr.cwl"
        }
    ],
    "$namespaces": {
        "sbg": "https://sevenbridges.com#",
        "arv": "http://arvados.org/cwl#"
    }
}