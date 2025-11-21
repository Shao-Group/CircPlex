# Introduction

CircPlex is a tool for reconstrcuting circular RNA sequence from Rolling Circle Amplification (RCA) long reads. The data and scripts for reproducing the results of the paper are available [here](https://github.com/Shao-Group/CircPlex-test).

# Dependency

CircPlex requires C++11, KMC (version 3.2.4), EquiRep (version 1.0.0) and Minimap2 (version 2.24) or later to compile and run successfully. Make sure your compiler supports C++11 by checking its version.

# Installation

Clone the git repository of CircPlex using the command:

```
git clone https://github.com/Shao-Group/CircPlex.git
```

<!-- Or download the source code of latest EquiRep from [here](https://github.com/Shao-Group/EquiRep/releases/download/v1.0.0/CircPlex-1.0.0.tar.gz). -->

<!-- Use the following commands to build EquiRep:
```
cd XXX/src
./configure
make
```
The executable file `XXX` will appear at `src/XXX`. -->

CircPLex uses KMC, EquiRep and Minimap2 in its pipeline.

# Install KMC

If KMC has not been downloaded/installed, install KMC
from (https://github.com/refresh-bio/KMC).

# Install EquiRep

If EquiRep has not been downloaded/installed, install EquiRep
[(license)](https://github.com/Shao-Group/EquiRep/blob/master/LICENSE) from (https://github.com/Shao-Group/EquiRep).

# Install Minimap2

If Minimap2 has not been downloaded/installed, install Minimap2
[(license)](https://github.com/lh3/minimap2/blob/master/LICENSE.txt) from (https://github.com/lh3/minimap2).

# Usage

CircPlex processes an input FASTA/FASTQ file and generates an output FASTA file with the circular RNA seqeunces.

The usage of CircPlex is:
```
./CircPlex.sh <input_fasta_file> <output_file_prefix> <kmc-path> <equirep-path> <minimap2-path>
```
Arguments:

`<input_file>` - Path to the input FASTA/FASTQ file.

`<output_file_prefix>` - Prefix for the output FASTA file.

`<kmc-path>` - Path to KMC excutables directory.

`<equirep-path>` - Path to EquiRep executable directory.

`<minimap2-path>` - Path to Minimap2 executable directory.

# Running CircPlex on a small example

A small example of input data `input.fasta` is available in the `example` directory.

Commands to enter `example` directory and run CircPlex using `input.fasta` as input:
```
cd CircPlex
./CircPlex.sh input.fasta ./example/output <kmc-path> <equirep-path> <minimap2-path>
```

An output file named `output_circRNA_seqs.fasta` will appear in the `example` directory.

