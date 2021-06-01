# RENANO FASTQ
## A reference-based compressor for nanopore FASTQ files
#### Documentation (preprint): https://www.biorxiv.org/content/10.1101/2021.03.26.437155v1
## Description
RENANO is a reference-based lossless FASTQ data compressor, specifically tailored to compress FASTQ files generated with nanopore sequencing technologies.
RENANO improves on its state of the art predecessor [ENANO](https://github.com/guilledufort/EnanoFASTQ/blob/master/README.md), by providing a more efficient base call sequence compression component.
Two compression algorithms are introduced, corresponding to the following scenarios: (1) a reference genome is available without cost to both the compressor and the decompressor;  (2) the reference genome is available only on the compressor side, and a compacted version of the reference is included in the compressed file.

## Install with Conda
To install directly from source, follow the instructions in the next section.

RENANO is available on conda via the bioconda channel. See [this](https://bioconda.github.io/user/install.html) page for installation instructions for conda. Once conda is installed, do the following to install renano.
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install renano
```
Note that if renano is installed this way, it should be invoked with the command `renano` rather than `./renano`. The bioconda [help page](https://bioconda.github.io/user/install.html) shows the commands if you wish to install renano in an environment.

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

## Datasets information

To test our compressor we ran experiments on the following datasets. The full information of the datasets is in our publication.

| Dataset | Num. of files | size (GB) | Description | Link |
|------|------|------|------|------|
*hss* | 1 | 268 | Human GM12878 Utah/Ceph cell line | https://github.com/nanopore-wgs-consortium/NA12878 |
*bra\** | 18 | 46 | Doubled haploid canola (Brassica napus L.) | https://www.nature.com/articles/s41598-019-45131-0#data-availability |
*sor\** | 4 | 133 | Sorghum bicolor Tx430 | https://www.nature.com/articles/s41467-018-07271-1#data-availability |
*fly\** | 1 | 17 | Drosophila ananassae | https://www.g3journal.org/content/8/10/3131#sec-1 |
*yst\** | 5 | 6 | Saccharomyces cerevisiae S288C | https://academic.oup.com/gigascience/article/6/2/giw018/2865217 |
*mic\** | 1 | 12 | Microbial community (metagenomic) | https://www.nature.com/articles/s41598-020-61989-x |

\*Datasets that require the SRA toolkit to be downloaded. 

### Downloading the datasets and the reference genomes

To download a dataset you have to run the *download_script.sh* of the specific dataset.
For example, to download *sor* run:
```bash
cd RENANO
dataset/sor/download_script.sh
```

To download the reference genome of a dataset (except for dataset *mic*) you have to run the *download_gene.sh* script of the specific dataset.
For example, to download *sor* reference genome run:
```bash
cd RENANO
dataset/sor/download_gene.sh
```

The scripts use the command *wget* to perform the download. 
To install *wget* on macOS run:
 ```bash
brew install wget
```
To install *wget* on Ubuntu or CentOS run:
 ```bash
sudo apt-get install wget
```

Some datasets require the SRA toolkit (2.9.6-1 release) to be downloaded. To install the SRA toolkit you can follow the instructions here https://ncbi.github.io/sra-tools/install_config.html, and place the toolkit's root-folder under the RENANO directory, or you can run one of the scripts we provide. There is a different script for each OS, so you have to choose the one corresponding to your OS.
For example, to install the SRA toolkit on macOS you can run:
 ```bash
cd RENANO
./install_SRA_mac.sh
```

In the case of the metagenomic dataset *mic* we do not have a reference genome sequence in advance. In this sense, in the next section we propose pipeline of operations to consruct a reference genome sequence that represents the most prevalent organisms in a dataset, which we utilized with dataset *mic*.

### Constructing a reference sequence for a metagenomic (or contaminated) dataset

To construct a referernce genome sequence for a metagenomic (or contaminated) dataset, such as *mic*, we propose a pipeline of operations and provide a series of scripts that facilitate its excecution. Let *[DS]* be the name of the dataset for which we want to construct the reference, a dataset folder with the name *[DS]*, containing the desired set FASTQ files, must be previously created in directory *RENANO/datasets/[DS]*.

Note: Steps 1 to 3 are optional. If Kraken2 is already installed, the user can create its own Kraken2 report, and then use it as input in Step 4.

- Step 1: Install Kraken2 into RENANO's root folder
```bash
cd RENANO
./install_kraken2.sh
```

- Step 2: Download and install Kraken2 pre-compiled mini database (8GB):
```bash
cd RENANO
./install_kraken2_stdDB.sh
```

- Step 3: Run Kraken2 on de desired dataset. The script receives two arguments: the dataset name *[DS]*, and the number of threads to be used. The script creates a Kraken2 report in folder *RENANO/datasets/[DS]* with name *[DS].report*. The report contains a statistical report of the species that were detected by Kraken2 in the dataset.
```bash
cd RENANO
./run_kraken2.sh -d [DS] -t [NUM_THREADS] 
```

- Step 4: Construct a representative reference sequence for dataset *[DS]* composed of the most prevalent species genomes reported in the Kraken2 report. The script receives one obligatory argument, and three optional. The obligatory argument is the dataset name *[DS]*. The script also receives the optional parameter *-r* with the name of a Kraken2 report file of the *[DS]* dataset; the default value is *RENANO/datasets/[DS].report*. The optional argument *-d* is a float number that represents the minimum percentage a species has to cover in the dataset, in terms of reads, to be included in the representative reference sequence; the default value is 0.3%. Lastly, argument *-m* is a positive integer that represents the maximum number of species that we want add to the representative reference sequence; the default value is 20. The script outputs a MULTI-FASTA file to *RENANO/datasets/[DS].fna*, which is constructed as the concatenation of the most prevalent species detected in the Kraken2 report.

```bash
cd RENANO
python create_ref.py -d [DS] -r [report_file_name] -c [MIN_COVERAGE] -m [MAX_SPECIES] 
```

## Alignment information

To obtain alignment information in [PAF format](https://lh3.github.io/minimap2/minimap2.html) for each FASTQ file we recommend using the tool [Minimap2](https://github.com/lh3/minimap2).

To install Minimap2 in the ROOT folder run:
 ```bash
cd RENANO
./install_minimap2.sh
```

To align a specific FASTQ file against a reference genome using Minimap2 run:
 ```bash
cd RENANO
minimap2/minimap2 -x map-ont --secondary=no --cs [ref_file] [fastq_file] > [paf_file]
```

To align all the files of a dataset against the corresponding reference genome use *run_minimap.sh* script.
For this script to work, both the dataset and the corresponding reference genome have to be previously downloaded, following the instructions in the previous section.
For example, to align the files in *sor* run:
 ```bash
cd RENANO
./run_minimap.sh sor
```

## Examples
We add an *example* folder with test files to run simple use examples the tool.
If installed using conda, use the command `renano` instead of `renano/renano`.
### Compress using RENANO with reference
To run the compressor with 8 threads on the example file:
```bash
cd RENANO
renano/renano -t 8 -r example/yst_genome.fna example/SAMPLE.paf example/SAMPLE.fastq example/SAMPLE.renano
```
### Decompress using RENANO with reference
To decompress with 8 threads the example compressed file:
```bash
cd RENANO
renano/renano -t 8 -d -r example/yst_genome.fna example/SAMPLE.renano example/SAMPLE_dec.fastq
```

### Check if decoding is successful
The output has to be empty.
```bash
cmp example/SAMPLE.fastq example/SAMPLE_dec.fastq
```


## Credits
The methods used for encoding the reads identifiers, and to model frequency counters,
are the ones proposed by James Bonfield in FQZComp, with some modifications.
The range coder is derived from Eugene Shelwien.
The kseq library used to parse FASTA files is authored by Heng Li.
