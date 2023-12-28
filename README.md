# KmerCamel🐫
KmerCamel🐫 provides implementations of algorithms for efficiently representing a set of k-mers as a [masked superstring](https://doi.org/10.1101/2023.02.01.526717), based on the following paper:

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Masked superstrings as a unified framework for textual *k*-mer set representations. *bioRxiv* 2023.02.01.526717, 2023.
[https://doi.org/10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

See [supplementary materials](https://github.com/karel-brinda/masked-superstrings-supplement) of the aforementioned paper for experimental results with KmerCamel🐫.

The implemented algorithms are the following:
- Local Greedy algorithm	 - constructs the superstring by locally finding and appending an unused k-mer with the largest overlap.
- Global Greedy algorithm - constructs the superstring by merging two k-mers with the largest overlap

They come in two different implementations (their results may differ due to the differences in used data structures):
- Encoding the k-mers as integers and using fast prefix/suffix equality checks.
- Using the Aho-Corasick automaton.

Note that at this point only the implementations with hash table are optimized and that the Aho-Corasick automaton
based versions of the algorithms are only experimental.

The hashing based implementations of the default KmerCamel🐫 (`./kmercamel`) support $k$-mer with $k$ at most 31,
whereas the larger KmerCamel🐫 (`./kmercamel-large`) supports $k$-mers with $k$ at most 63 (at the cost of slight slowdown).

All algorithms can be used to either work in the unidirectional model or in the bidirectional model
(i.e. treat $k$-mer and its reverse complement as the same; in this case either of them appears in the result).

### Mask optimization

Both the local and global greedy algorithms compute a superstring of the k-mers together with a default mask that minimizes the number of 1s.
As noted in the paper, there are possibly many masks for a superstring, and various objectives for mask optimization can be considered,
depending on the application.
We refer to the [supplementary repository of the paper](https://github.com/karel-brinda/masked-superstrings-supplement/tree/main/experiments/08_optimize_masks)
for Python scripts that take a masked superstring as input and output a masked superstring with mask optimized to the number of 1s (min or max)
or to the number of runs of 1s (min).

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
- `-o output_path` - the path to output file. If not specified, output is printed to stdout.
- `-d value_of_d` - d_max used in Local Greedy. Default 5. Increasing `d` beyond `k` has no effect.
- `-c` - treat k-mer and its reverse complement as equal.
- `-m` - turn off memory optimizations for `global`.
- `-h` - print help.
- `-h` - print version.


The output contains the resulting superstring - capital letters indicate that at given position, a k-mer starts.

For example:

```
./kmercamel -p ./spneumoniae.fa -a local -k 12 -d 7 -c
```

runs the Local Greedy in the bidirectional model on the streptococcus fasta file with `k=12` and `d=7`.

Alternatively, if your operating system supports it, you can run `./🐫` instead of `./kmercamel`.

Currently, KmerCamel🐫 does not support gziped files as an input.
A possible workaround is to use `gzcat` and process substitution.

```
./kmercamel -k 13 -p <(gzcat fasta_file.fa.gz)
```

Note: on some systems you might need to use the name `zcat` instead of `gzcat`.

### Large $k$-mers

The default version of KmerCamel🐫 does not support $k > 31$. For those values, use the large KmerCamel🐫,
which supports $k < 64$.

For example:

```
./kmercamel-large -p ./spneumoniae.fa -a global -k 63 -c
```

Note: for smaller $k$ it is recommended to use default KmerCamel🐫 as it is faster.

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

