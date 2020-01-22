cwlVersion: v1.0
class: CommandLineTool
label: ccube
baseCommand: ["Rscript", "/home/pipeline/run_analysis_ccube.R"]
requirements:
  - class: DockerRequirement
    dockerPull: smcheteval/ccube:0.1

inputs:
  input_vcf:
    type: File
    inputBinding:
      position: 1

  battenberg_file:
    type: File
    inputBinding:
      position: 2

outputs:
  cellularity:
    type: File
    outputBinding:
      glob: 1A.txt

  population:
    type: File
    outputBinding:
      glob: 1B.txt

  proportion:
    type: File
    outputBinding:
      glob: 1C.txt

  cluster_assignment:
    type: File
    outputBinding:
      glob: 2A.txt

  cocluster_assignment:
    type: File
    outputBinding:
      glob: 2B.txt
