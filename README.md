# KmerCamel 🐫
KmerCamel provides implementations of algorithms for efficiently representing a set of k-mers as a masked superstring.
Note that this is an experimental program, focused mainly on bacterial genomes and pangenomes.

The implemented algorithms are the following:
- Local Greedy algorithm	 - constructs the superstring by locally finding and appending an unused k-mer with the largest overlap.

- Global Greedy algorithm - constructs the superstring by running the GREEDY algorithm to find the approximated shortest common superstring.

They come in two different implementations (their results may differ due to the differences in used data structures):
- Encoding the k-mers as integers and using fast prefix/suffix equality checks.
- Using the Aho-Corasick automaton.

In general, the hashing algorithms are faster than those using the Aho-Corasick automaton.
Note also that all the algorithms have high memory requirements at this point.
This is especially true of the implementations with the Aho-Corasick automaton which are only experimental at this stage.

The following table summarizes the current state of development of these algorithms:

| Algorithm                     | Support for complements | Optimizations       | 
|-------------------------------|-------------------------|---------------------|
| Global Greedy with AC         | YES                     | None yet            |
| Global Greedy with hashing    | YES                     | Basic               |
| Local Greedy with AC          | YES                     | None yet            |
| Local Greedy with hashing     | YES                     | None yet            |

## How to install

First clone the repo and its dependency:

```
git clone --recursive https://github.com/OndrejSladky/kmercamel
```

Compile the program by running `make`.


## How to run

The program has the following arguments:

- `-p path_to_fasta` - the path to fasta file. This is a required argument.
- `-k value_of_k` - the size of one k-mer. This is a required argument.
- `-a algorithm` - the algorithm which should be run. Either `global` or `globalAC` for Global Greedy, `local` or `localAC` for Local Greedy.
The versions with AC use Aho-Corasick automaton. Default `global`.
- `-d value_of_d` - d_max used in Local Greedy. Default 5. Increasing `d` beyond `k` has no effect.
- `-c` - treat k-mer and its reverse complement as equal.
- `-h` - print help.


The output contains the resulting superstring - capital letters indicate that at given position, a k-mer starts.

For example:

```
./kmercamel -p ./spneumoniae.fa -a localAC -k 12 -d 7
```

runs the Local Greedy on the streptococcus fasta file with `k=12` and `d=7`.

Alternatively, if your operating system supports it, you can run `./🐫` instead of `./kmercamel`.

## Converting k-mer set superstring representation to the traditional one

Run `./convert_superstring.py < input.fa`. This runs a Python script which inputs the superstring masked representation and outputs the SPSS representation.

## How to test


For integration tests you'll have to install [jellyfish (v2)](https://github.com/gmarcais/Jellyfish)
and add it to PATH.

You can verify all the algorithms for `1 < k < 32` on a given fasta file by running:

```
make verify
```

You can run the c++ unittests by `make cpptest`.

Similarly testing the convert script can be done via `make converttest`.

To run all the test, simply run `make test`.

## Contact

You can contact the developer at [sladky@iuuk.mff.cuni.cz](mailto:sladky@iuuk.mff.cuni.cz).

