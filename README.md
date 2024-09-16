```
______  _    _                ______
| ___ \| |  | |               | ___ \
| |_/ /| |  | | __ ___   _____| |_/ /
| ___ \| |/\| |/ _` \ \ / / _ \    /
| |_/ /\  /\  / (_| |\ V /  __/ |\ \
\____/  \/  \/ \__,_| \_/ \___\_| \_|

```
## Description
This repo contains the BWaveR software implementation.

## Dependencies

+ zlib

It is suggested to create a conda environment to easily handle the dependencies.

## Building

To build the BWT and the MAP executables, type the following command in the terminal:
```
make BWT MAP
```


## Testing

### 1. BWT and Suffix Array computation
```
./BWT <ref> <kmer> <f_bwt>
```

<b>Inputs</b>:

+ *ref*  =  path to the reference FASTA file (uncompressed or gzipped);

+ *kmer*  =  length of kmers employed for the BWT computation;

+ *f_bwt*  =  name of the output file containing the BWT and the Suffix Array of the reference.


<b>Examples</b>:

Compute BWT and Suffix Array of E.Coli reference, using kmers of length 51, and store them in a file named *ecoli.bwt*:
```
./BWT ./test_data/EColi.fas.gz 51 ./test_data/ecoli.bwt
```



<b>Output</b>:

+ *f_bwt* = intermediate text file containing the BWT and the Suffix Array of the reference



### 2. Data Encoding and Read Mapping 

```
./MAP <reads> <f_bwt> <f_output> <block_dim> <superblock>
```


<b>Inputs</b>:

+ *reads*  =  path to the reads FASTQ file (uncompressed or gzipped);

+ *f_bwt*  =  name of the file containing the BWT and the Suffix Array of the reference;

+ *f_output*  =  name of the output file containing the mapping results;

+ *block_dim*  =  dimension of blocks employed for building the RRR sequences (max = 15, bigger is unuseful for real genomes);

+ *superblock*  =  grouping factor for blocks in the RRR sequences (currently, max = 100 in MAP_FAST);


<b>Examples</b>:

Map the read library *ecoli1000* to the E.Coli reference, using blocks of 8 bits and a grouping factor of 50. Store the mapping results in a file named *ecoli.res*:
```
./MAP ./test_data/ecoli1000.fq.gz ./test_data/ecoli.bwt ./test_data/ecoli.res 8 50
```


<b>Output</b>:

+ *f_output*  =  human-readable file containing the mapping results, using the following syntax:

```
> <read_name> | [<read_length> bp] | <strand> <N_occurrences>:
	<space-separated list of N positions in the reference>
```
The strand is represented through a <b>+</b> or a <b>-</b> symbol:

<b>+</b> means the read has occurrence(s) on the original strand represented in the reference sequence;

<b>-</b> means the read has occurrence(s) on the complementary strand.

*Example:*

```
> r373/1 | [35 bp] | + 1:
	39353

> r376/1 | [35 bp] | + 2:
	2727106	3424705
> r376/1 | [35 bp] | - 5:
	697882	472945	431543	604067	4413850

> r383/1 | [35 bp] | - 1:
	1491916
```
