# KmerCamel🐫
KmerCamel🐫 provides implementations of algorithms for efficiently representing a set of k-mers as a [masked superstring](https://doi.org/10.1101/2023.02.01.526717), based on the following paper:

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Masked superstrings as a unified framework for textual *k*-mer set representations. *bioRxiv* 2023.02.01.526717, 2023.
[https://doi.org/10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

See [supplementary materials](https://github.com/karel-brinda/masked-superstrings-supplement) of the aforementioned paper for experimental results with KmerCamel🐫.

The computation of masked superstring using KmerCamel🐫 is done in two steps -
first a superstring is computed and then its mask can be optimized.

The implemented superstring algorithms are the following:
- Local Greedy algorithm	 - constructs the superstring by locally finding and appending an unused k-mer with the largest overlap.
- Global Greedy algorithm - constructs the superstring by merging two k-mers with the largest overlap

They come in two different implementations (their results may differ due to the differences in used data structures):
- Encoding the k-mers as integers and using fast prefix/suffix equality checks.
- Using the Aho-Corasick automaton.

Note that at this point only the implementations with hash table are optimized and that the Aho-Corasick automaton
based versions of the algorithms are only experimental.

The hashing based implementations of KmerCamel🐫 support $k$-mer with $k$ at most 63,

All algorithms can be used to either work in the unidirectional model or in the bidirectional model
(i.e. treat $k$-mer and its reverse complement as the same; in this case either of them appears in the result).

The implemented mask optimization algorithms are the following:
- Minimize the number of 1s in the mask.
- Maximize the number of 1s in the mask.
- Minimize the number of runs of 1s in the mask.

## How to install

First clone the repo and its dependency:

```
git clone --recursive https://github.com/OndrejSladky/kmercamel
```

Compile the program by running `make`.


### Dependencies

The program requires the GLPK library. This is usually available by default. If not it can be installed via:

```
sudo apt-get install libglpk-dev
```

on Ubuntu or

```
brew install glpk
```

on macOS.

## How to run

The program has the following arguments:

- `-p path_to_fasta` - the path to fasta file (can be `gzip`ed). This is a required argument.
- `-k value_of_k` - the size of one k-mer (up to 63). This is a required argument.
- `-a algorithm` - the algorithm which should be run. Either `global` or `globalAC` for Global Greedy, `local` or `localAC` for Local Greedy.
The versions with AC use Aho-Corasick automaton. Default `global`.
- `-o output_path` - the path to output file. If not specified, output is printed to stdout.
- `-d value_of_d` - d_max used in Local Greedy. Default 5. Increasing `d` beyond `k` has no effect.
- `-c` - treat k-mer and its reverse complement as equal.
- `-m` - turn off memory optimizations for `global`.
- `-h` - print help.
- `-h` - print version.


The output contains the resulting superstring - capital letters indicate that at given position, a k-mer starts.

For example:

```
./kmercamel -p ./spneumoniae.fa -a local -k 31 -d 5 -c
```

runs the Local Greedy in the bidirectional model on the streptococcus fasta file with `k=31` and `d=5`.

Alternatively, if your operating system supports it, you can run `./🐫` instead of `./kmercamel`.

### Mask optimization

For mask optimization, run the subcommand `optimize` with the following arguments:

- `p path_to_fasta` - the path to fasta file (can be `gzip`ed). This is a required argument.
- `k k_value` - the size of one k-mer. This is a required argument.
- `a algorithm` - the algorithm for mask optimization. Either `ones` for maximizing the number of 1s, `runs` for minimizing the number of runs of 1s, `runsapprox` for approximately minimizing the number of runs of 1s, or `zeros` for maximizing the number of 0s. Default `ones`.
- `o output_path` - the path to output file. If not specified, output is printed to stdout.
- `c` - treat k-mer and its reverse complement as equal.
- `h` - print help.
- `v` - print version.

For example:

```
./kmercamel optimize -p ./global-spneumoniae.fa -k 31 -a runs -c
```

minimizes the number of runs of 1s in the mask of the superstring computed by Global Greedy in the bidirectional model on the streptococcus fasta file with `k=12`.

### Turn off memory optimizations for Global

In order to reduce the memory footprint of hash-table based Global Greedy,
it uses several optimizations that reduce memory but increase the running time.
If memory is not the bottleneck, turning off the memory optimizations might be desirable.
This in practice about halves the running time.


## Converting k-mer set superstring representation to the (r)SPSS representations

Run `./convert_superstring.py < input.fa`. This runs a Python script which inputs the superstring masked representation and outputs the (r)SPSS representation.

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

