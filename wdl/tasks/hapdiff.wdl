version 1.0

task hapdiff_t {
  input {
    File ctgsPat
    File ctgsMat
    File reference
    File? vntrAnnotations
    Int minSvSize = 25
    Int threads = 32
    Int memSizeGb = 128
    Int diskSizeGb = 256
  }

  String trfString = if defined(vntrAnnotations) then "--tandem-repeats ~{vntrAnnotations}" else ""
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    hapdiff.py --reference ~{reference} ~{trfString} --pat ~{ctgsPat} --mat ~{ctgsMat} --out-dir hapdiff -t ~{threads} --sv-size ~{minSvSize} 2>&1 | tee hapdiff.log
  >>>

  output {
    File hapdiffUnphasedVcf = "hapdiff/hapdiff_unphased.vcf.gz"
    File hapdiffPhasedVcf = "hapdiff/hapdiff_phased.vcf.gz"
    File hapdiffLog = "hapdiff.log"
  }

  runtime {
    preemptible: 1
    docker: "mkolmogo/hapdiff:0.7"
    cpu: threads
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}
