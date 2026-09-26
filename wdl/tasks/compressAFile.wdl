version 1.0

workflow compressionWDL
{
	input {
		File inputFile
		String sample_name
	}

	call compressAFile {
		input: 
                inputFile = inputFile,
                sample_name = sample_name
	}

	output {
		File outCompressedFile = compressAFile.outputFile

	}
}

task compressAFile {


    input {
        File inputFile
        String sample_name
        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = round(size(inputFile, "GB")) + 40
        Int preemptible = 1
    }

    String input_basename = basename(inputFile)

    command <<<

        set -eux -o pipefail

        bgzip -@ ~{threads} ~{inputFile} > ~{sample_name}.~{input_basename}.gz

  >>>

  output {
    File outputFile = "~{sample_name}.~{input_basename}.gz"

  }

  runtime {
    docker: "meredith705/modkit:latest"
    cpu: threads
	memory: memSizeGB + " GB"
	disks: "local-disk " + diskSizeGB + " SSD"
    preemptible: preemptible
  }
}


