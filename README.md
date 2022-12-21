# KmerCamel 🐫
KmerCamel provides implementation of three algorithms for efficiently representing set of k-mers as a masked superstring.

- Pseudosimplitigs
- Global GREEDY algorithm
- Streaming algorithm

The first two come in two different implementations:
- Encoding the k-mers as integers and using fast prefix/suffix equality checks.
- Using the Aho-Corasick automaton.

The versions with k-mer encoding are usually faster, but the versions with Aho-Corasick automaton
support arbitrarily large k-mers.

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
- `-a algorithm` - the algortihm which should be run. Either `streaming` for Streaming algorithm, `greedy` or `greedyAC` for global GREEDY, `pseudosimplitigs` or `pseudosimplitigsAC` for greedily computing pseudosimplitigs.
The versions with AC use Aho-Corasick automaton. Default `greedy`.
- `-d value_of_d` - d_max used in pseudosimplitigs. Default 5. Increasing `d` beyond `k` has no effect.
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

## How to test


For integration tests you'll have to install [jellyfish (v2)](https://github.com/gmarcais/Jellyfish)
and add it to PATH.

You can verify all the algortihms for `1 < k < 32` on a given fasta file by running:

```
make verify
```

For unittest, install googletest module by running:

```
git submodule init
git submodule update
```

You can then run the cpp unittest by `make cpptest`.

Similarly testing the convert script can be done via `make converttest`.

To run all the test, simply run `make test`.



