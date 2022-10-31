## Algorithms

| Algorithm                     | Implemented | Basic optimizations | Support for complements |
|-------------------------------|-------------|---------------------|-------------------------|
| GREEDY with AC                | YES         | NO                  | NO                      |
| GREEDY with hashing           | NO          | NO                  | NO                      |
| Pseudosimplitigs with AC      | YES         | NO                  | NO                      |
| Pseudosimplitigs with hashing | YES         | NO                  | NO                      |


## How to install

First clone the repo:

```
git clone https://github.com/GordonHoklinder/k-mers
```

Compile the program by running `make`.


## How to run

The program has the following arguments:

- `-p path_to_fasta` - the path to fasta file. This is a required argument.
- `-k value_of_k` - the size of one k-mer. Default 13.
- `-a algorithm` - the algortihm which should be run. Either `greedy` for global GREEDY, `pseudosimplitigs` or `pseudosimplitigsAC` for greedily computing pseudosimplitigs.
	- The latter uses Aho-Corasick automaton. Default `greedy`.
- `-d value_of_d` - d_max used in pseudosimplitigs. Default 5.
- `-s` - if provided do not print the resulting superstring, but rather stats about it. If not print the superstring - capital letters indicate that at given position, a k-mer starts.
- `-h` - print help.

For example:

```
./kmers ./spneumoniae.fa -a pseudosimplitigsAC -k 12 -d 7
```

runs the pseudosimplitigs on the streptococcus fasta file with `k=12` and `d=7`.

