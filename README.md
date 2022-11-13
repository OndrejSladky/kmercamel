## Algorithms

| Algorithm                     | Implemented | Basic optimizations | Support for complements |
|-------------------------------|-------------|---------------------|-------------------------|
| GREEDY with AC                | YES         | NO                  | NO                      |
| GREEDY with hashing           | YES         | YES                 | YES                     |
| Pseudosimplitigs with AC      | YES         | NO                  | NO                      |
| Pseudosimplitigs with hashing | YES         | YES                 | YES                     |


## How to install

First clone the repo:

```
git clone https://github.com/GordonHoklinder/k-mers
```

Compile the program by running `make`.


## How to run

The program has the following arguments:

- `-p path_to_fasta` - the path to fasta file. This is a required argument.
- `-k value_of_k` - the size of one k-mer. This is a required argument.
- `-a algorithm` - the algortihm which should be run. Either `greedy` or `greedyAC` for global GREEDY, `pseudosimplitigs` or `pseudosimplitigsAC` for greedily computing pseudosimplitigs.
The versions with AC use Aho-Corasick automaton. Default `greedy`.
- `-d value_of_d` - d_max used in pseudosimplitigs. Default 5.
- `-s` - if provided do not print the resulting superstring, but rather stats about it. If not print the an output in fasta file format with one record.
The name includes the statistics about the superstring and the sequence is the superstring - capital letters indicate that at given position, a k-mer starts.
- `-c` - treat k-mer and its reverse complement as equal.
- `-h` - print help.

For example:

```
./kmers ./spneumoniae.fa -a pseudosimplitigsAC -k 12 -d 7
```

runs the pseudosimplitigs on the streptococcus fasta file with `k=12` and `d=7`.

## Converting k-mer set superstring representation to traditional one

Run `./convert_superstring`. This runs a Python script which inputs the superstring masked representaion and outputs the SPSS representation.

## How to verify

First, install [jellyfish](https://github.com/gmarcais/Jellyfish).

You can verify all the algortihms for `4 < k < 32` on a given fasta file by running:

```
./verify path_to_fasta
```
