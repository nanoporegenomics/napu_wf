version 1.0

workflow compressionWDL
{
	input {
		File inputFile
	}

	call compressAFile {
		input: 
                inputFile = inputFile,
	}

	output {
		File outCompressedFile = compressAFile.outputFile

	}
}


task compressAFile {


    input {
        File inputFile
        Int memSizeGB = 40
        Int threads = 12
        Int diskSizeGB = round(size(inputFile, "GB")) + 40
        Int preemptible = 1
    }


    command <<<

        set -eux -o pipefail

        bgzip -@ ~{threads} ~{inputFile}

  >>>

  output {
    File outputFile = "~{inputFile}.gz"

  }

  runtime {
    docker: "meredith705/modkit:latest"
    cpu: threads
	memory: memSizeGB + " GB"
	disks: "local-disk " + diskSizeGB + " SSD"
    preemptible: preemptible
  }
}


