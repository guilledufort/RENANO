# RENANO FASTQ
## A reference-based compressor for nanopore FASTQ files
#### Documentation: https://www.biorxiv.org/content/10.1101/2021.03.26.437155v1
## Description
RENANO is a FASTQ lossless reference-based compression algorithm especially designed for nanopore sequencing FASTQ files. 

## Install from source code

### Download repository
```bash
git clone https://github.com/guilledufort/RENANO.git
```

### Requirements
0. g++ ( >= 4.8.1)
1. OpenMP library

### Install 

The following instructions will create the renano executable in the directory *renano*.
To compile renano you need to have the g++ compiler and the OpenMP library for multithreading. 

On Linux (Ubuntu or CentOS) g++ usually comes installed by default, but if not run the following:
```bash
sudo apt update
sudo apt-get install g++
```

On macOS, install GCC compiler since Clang has issues with OpenMP library:
- Install HomeBrew (https://brew.sh/)
- Install GCC (this step will be faster if Xcode command line tools are already installed using ```xcode-select --install```):
```bash
brew update
brew install gcc@9
```

The g++ installer also installs the OpenMP library, so no further steps are needed.
To check if the g++ compiler is properly installed in your system run:

On Linux
```bash
g++ --version
```
On MacOS:
```bash
g++-9 --version
```
The output should be the description of the installed software.

To compile renano run:
```bash
cd RENANO/renano
make
```

# USAGE
Run the renano executable ```/PATH/TO/renano``` with the options below:
```console 
COMPRESSION:
> Without reference:
	renano [options] [input_file [output_file]]

> With reference:
	renano [options] -r [ref_file [paf_file]] [input_file [output_file]]

> With reference and making decompression independent of the reference:
	renano [options] -s [ref_file [paf_file]] [input_file [output_file]]

COMPRESSION OPTIONS:
	-c             To use MAX COMPRESION MODE. Default is FAST MODE.

	-k <length>    Base call sequence context length. Default is 7 (max 13).

	-l <lenght>    Length of the DNA neighborhood sequence. Default is 6.

	-t <num>       Maximum number of threads allowed to use by the compressor. Default is 8.

DECOMPRESSION:
> Without reference:
	renano -d [options] foo.renano foo.fastq

> With reference:
	renano -d [options] -r [ref_file] foo.renano foo.fastq

DECOMPRESSION OPTIONS:
	-t <num>       Maximum number of threads allowed to use by the decompressor. Default is 8.

```

## Credits
The methods used for encoding the reads identifiers, and to model frequency counters,
are the ones proposed by James Bonefield in FQZComp, with some modifications.
The range coder is derived from Eugene Shelwien.
The kseq library used to parse FASTA files is authored by Heng Li.
