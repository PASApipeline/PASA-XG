version 1.0

# PASA-XG in one task: the pipeline's own checkpoints make a rerun of the same task
# resume, and every deliverable is collected in outputs/.
workflow pasa_xg {
  input {
    File genome
    Array[File] reads = []
    File? bam
    File? bam_index
    String input_type
    File? annotation
    Int cpu = 16
    Int memory_gb = 64
    Int disk_gb = 200
    String extra_args = ""
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/pasa-xg/pasa-xg:latest"
  }
  call run_pasa_xg {
    input: genome = genome, reads = reads, bam = bam, bam_index = bam_index,
           input_type = input_type, annotation = annotation, cpu = cpu,
           memory_gb = memory_gb, disk_gb = disk_gb, extra_args = extra_args, docker = docker
  }
  output {
    Array[File] outputs = run_pasa_xg.outputs
  }
}

task run_pasa_xg {
  input {
    File genome
    Array[File] reads
    File? bam
    File? bam_index
    String input_type
    File? annotation
    Int cpu
    Int memory_gb
    Int disk_gb
    String extra_args
    String docker
  }
  command <<<
    set -euo pipefail
    src="~{if defined(bam) then "--bam " + bam else "--reads " + sep(" ", reads)}"
    pasa-xg run --genome ~{genome} $src --input_type ~{input_type} \
      ~{if defined(annotation) then "--annotation " + annotation else ""} \
      --cpu ~{cpu} --out run ~{extra_args}
  >>>
  output {
    Array[File] outputs = glob("run/outputs/*")
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} SSD"
  }
}
