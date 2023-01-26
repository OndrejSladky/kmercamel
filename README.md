# KmerCamel 🐫
KmerCamel provides implementation of algorithms for efficiently representing set of k-mers as a masked superstring.
Note that this is an experimental program, focused mainly on barterial genomes and pangenomes.

The implemented algorithms are following:
- Local greedy algorithm
	
Constructs the superstring by locally finding and appending an unused k-mer with the largest overlap.

- Global GREEDY algorithm

Constructs the superstring by running the GREEDY algorithm to find the approximated shortest common superstring.

They come in two different implementations (their results may differ due to the differences in used data structures):
- Encoding the k-mers as integers and using fast prefix/suffix equality checks.
- Using the Aho-Corasick automaton.

In general, the hashing algorithms are faster than the 
Note also that all the algorithms have high memory requirements at this point.
This is especially true of the implementations with the Aho-Corasick automaton which are only experimental at this stage.

The following table summarizes the current state of development of these algorithms:

| Algorithm                     | Support for complements | Optimizations       | 
|-------------------------------|-------------------------|---------------------|
| Global GREEDY with AC         | YES                     | None yet            |
| Global GREEDY with hashing    | YES                     | Basic               |
| Local greedy with AC          | YES                     | None yet            |
| Local greedy with hashing     | YES                     | None yet            |

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
- `-a algorithm` - the algortihm which should be run. Either `greedy` or `greedyAC` for global GREEDY, `pseudosimplitigs` or `pseudosimplitigsAC` for local greedy.
The versions with AC use Aho-Corasick automaton. Default `greedy`.
- `-d value_of_d` - d_max used in local greedy. Default 5. Increasing `d` beyond `k` has no effect.
The name includes the statistics about the superstring and the sequence is the superstring - capital letters indicate that at given position, a k-mer starts.
- `-c` - treat k-mer and its reverse complement as equal.
- `-h` - print help.

For example:

```
./kmers ./spneumoniae.fa -a pseudosimplitigsAC -k 12 -d 7
```

runs the local greedy on the streptococcus fasta file with `k=12` and `d=7`.

## Converting k-mer set superstring representation to the traditional one

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

## Contact

You can contact the developer at [ondra.sladky@gmail.com](mailto:ondra.sladky@gmail.com).

