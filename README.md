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

As the first argument pass the path to the fasta file.

Additionally you can pass other arguments:

- `-k` -- the size of one k-mer. Default 13.
- `-a` -- the algortihm which should be run. Either `greedy` for global GREEDY, `pseudosimplitigs` or `pseudosimplitigsAC` for greedily computing pseudosimplitigs.
- The latter uses Aho-Corasick automaton. Default `greedy`.
- `-d` -- d_max used in pseudosimplitigs. Default 5.

For example:

```
./kmers ./spneumoniae.fa -a pseudosimplitigsAC -k 12 -d 7
```

runs on the pseudosimplitigs on the streptococcus fasta file with `k=12` and `d=7`.

