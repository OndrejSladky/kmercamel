# KmerCamel🐫
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/kmercamel.svg?style=flag&label=BioConda%20Install)](https://anaconda.org/bioconda/kmercamel)
[![KmerCamel test](https://github.com/OndrejSladky/kmercamel/actions/workflows/ci.yml/badge.svg)](https://github.com/OndrejSladky/kmercamel/actions/)

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Prerequisites](#prerequisites)
* [Getting started](#getting-started)
  * [Installation](#installation)
  * [Compression for k-mer set storage](#compression-for-k-mer-set-storage)
  * [k-mer set indexing](#k-mer-set-indexing)
  * [Time and memory](#time-and-memory)
* [Detailed instructions](#detailed-instructions)
* [How it works](#how-it-works)
* [How to test](#how-to-test)
* [Issues](#issues)
* [Changelog](#changelog)
* [Licence](#licence)
* [Contact](#contact)

<!-- vim-markdown-toc -->

## Introduction

KmerCamel🐫 is a tool for efficiently representing a set of k-mers by a [masked superstring](https://doi.org/10.1101/2023.02.01.526717).

It is based on the following paper:

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Masked superstrings as a unified framework for textual *k*-mer set representations. *bioRxiv* 2023.02.01.526717, 2023.
[https://doi.org/10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

See [supplementary materials](https://github.com/karel-brinda/masked-superstrings-supplement) of the aforementioned paper for experimental results with KmerCamel🐫.

The computation of masked superstring using KmerCamel🐫 is done in two steps -
first a superstring is computed with its default mask and then its mask can be optimized.

KmerCamel🐫 implements in the `compute` subcommand the BIGREEDY algorithm for masked superstring computation that operates in two different regimes:
- Either, it reads an input FASTA file (optionally `gzip`ed), retrieves the associated $k$-mers, and computes masked superstrings by first internally computing simplitigs and then merging these on the largest overlap.
- Or, providing the `-S` param, it inputs a FASTA file with repetition-free SPSS, such as unitigs or simplitigs (optionally `gzip`ed), with each unitig/simplitig on a separate line. The algorithm then skips the initial compute-heavy part of computing simplitigs and proceeds by merging the provided unitigs/simplitigs.

KmerCamel🐫 supports all $k < 128$ and in the first regime, it enables $k$-mer filtering by minimum threshold frequency (`-z` param).
In both regimes, KmerCamel🐫 outputs a fasta file with a single record - a masked-cased superstring, which is in the nucleotide alphabet with case of the letters determining the mask symbols. The default masks are min-ones (minimize the number of ones) and max-one masks can be computed with the `-M` param, or the `maskopt` subcommand.



## Prerequisites

* GCC
* Zlib
* GLPK (can be installed via `apt-get install libglpk-dev` on Ubuntu or `brew install glpk` on macOS)

## Getting started

### Installation

Download and compile KmerCamel🐫 by running the following commands:

```
git clone --recursive https://github.com/OndrejSladky/kmercamel
cd kmercamel && make
```

Alternatively, you can install KmerCamel from [bioconda](https://bioconda.github.io/):
```
   conda install bioconda::kmercamel
```

### Compression for k-mer set storage

```
kmercamel compute -k 31 -o ms.msfa yourfile.fa         # Compute MS with the default mask
kmercamel ms2mssep -m mask.m -s superstring.s ms.msfa  # Extract superstring and mask
bzip2 --best mask.m
xz -T1 -9 superstring.s
```

For a super-efficient compression of the superstring (often <2 bits / bp), you use some of the specialized tools based on statistical compression such as [GeCo3](https://github.com/cobilab/geco3) or [Jarvis3](https://github.com/cobilab/jarvis3).

If the masked superstrings are to be computed from simplitigs/unitigs, change the first line to `kmercamel compute -k 31 -o ms.msfa -S simplitigs.fa`.

### k-mer set indexing

Example with [FMSI](https://github.com/OndrejSladky/fmsi/activity?ref=main):
```
kmercamel compute -k 31 -o /dev/null -M mas-opt.msfa yourfile.fa   # Compute MS and the max-one mask
fmsi index -p ms-opt.msfa                                          # Create a k-mer index
```

### Time and memory

To compute masked superstrings from scratch usually takes about 0.2-1.0s / 1M k-mers, but depends on the exact dataset; in particular about 30min to compute masked superstrings for the human genome. The memory consumption on human genome is about 36 GB.
Unless the provided unitigs/simplitigs are nearly isolated $k$-mers, it is faster to compute masked superstrings from unitigs/simplitigs. Computation from simplitigs is faster proportionally to the cummulative length of the representations.

## Detailed instructions

Examples of computing masked superstrings (`compute` subcommand):
```
kmercamel compute -k 31 -o ms.msfa yourfile[.fa|.fa.gz]              # From a (gziped) fasta file, use "-" for stdin
kmercamel compute -k 31 -o ms.msfa -S simplitigs.fa                  # Faster computation from simplitigs/unitigs.
kmercamel compute -k 31 -o ms.msfa -z 2 yourfile.fa                  # Represent only k-mers appearing at least z=2 times
kmercamel compute -k 31 -o ms.msfa -u yourfile.fa                    # Treat k-mer and its reverse complement as distinct
kmercamel compute -k 31 -o ms.msfa -M ms-max-one.msfa yourfile.fa    # Also store MS with maximum ones
```
If the input file are simplitigs (or eulertigs), the execution can be significantly speeded up by adding the `-S` flag.
However, note that if `-S` is used with matchtigs (SPSS with repetitions), it may result it unnecessarily long outputs. The output will still be correct, but the default masks are not guaranteed to be min-one.

Examples of optimizing masks:
```
kmercamel maskopt -t max-one -o ms-opt.msfa -k 31 ms.msfa    # Maximize the number of 1s in the mask
kmercamel maskopt -t min-one -o ms-opt.msfa -k 31 ms.msfa    # Minimize the number of 1s in the mask
```

Format conversions:
```
kmercamel mssep2ms -m dataset.m -s dataset.s -o dataset.msfa  # M and S -> mask-cased MS in msfa
kmercamel ms2mssep -m dataset.m -s dataset.s dataset.msfa     # Mask-cased MS -> M and S
kmercamel spss2ms -k 31 -o dataset.msfa dataset.rspss         # rSPSS/general fasta to its corresponding MS
kmercamel ms2spss -k 31 -o dataset.rspss dataset.msfa         # Splitting MS in msfa into rSPSS in fa
```

Compute lower bound on the minimum possible superstring length of a k-mer set:
```
./kmercamel lowerbound -k 31 yourfile.fa        # Print the lower bound.
./kmercamel lowerbound -k 31 -S simplitigs.fa   # Computation of lower bound faster directly from maximal simplitigs.
./kmercamel lowerbound -k 31 -z 2 yourfile.fa   # Filter k-mer with fewer occurrences than 2
```

To view all options for a particular subcommand, run `kmercamel <subcommand> -h`.


## Advanced and experimental usage

KmerCamel🐫 supports two other algorithms for MS computation, _local greedy_ and _streaming_.

```
kmercamel compute -k 31 -o ms.msfa -a [streaming|local-greedy] yourfile.fa  # Use a different algorithm instead of BIGREEDY (`greedy`)
```

Additionally, KmerCamel🐫 experimentally implements BIGREEDY and local greedy algorithms in their Aho-Corasick automaton versions (`greedy-ac`, `local-greedy-ac`).
Note that they are much slower than the original versions, but they can handle arbitrarily large $k$s.


KmerCamel🐫 also supports the option to minimize the number of runs of ones in the mask.

```
kmercamel maskopt -t min-run -o ms-opt.msfa -k 31 ms.msfa         # Optimal potentially slow algorithm for minimizing the number of runs of consecutive 1s in the mask
kmercamel maskopt -t approx-min-run -o ms-opt.msfa -k 31 ms.msfa  # Inexact fast heuristic for minimizing the number of runs of consecutive 1s in the mask
```

## How it works

For details about the algorithms and their implementation, see the [Code README](./src/README.md).

## How to test

To ensure correctness of the results, KmerCamel🐫 has two levels of tests - unit tests and file-specific integration tests.

For integration tests  install [jellyfish (v2)](https://github.com/gmarcais/Jellyfish)
and add it to PATH.

You can verify all the algorithms for `1 < k < 128` on a *S. pneumoniae* by running `make verify`.
To run it on another dataset, see the [verification script](./verify.py).

You can run the C++ unittests by `make cpptest`.

To run all the test, simply run `make test`.

## Issues

Please use [Github issues](https://github.com/OndrejSladky/kmercamel/issues).


## Changelog

See [Releases](https://github.com/OndrejSladky/kmercamel/releases).


## Licence

[MIT](https://github.com/OndrejSladky/kmercamel/blob/master/LICENSE.txt)


## Contact

[Ondrej Sladky](https://iuuk.mff.cuni.cz/~sladky/) \<ondra.sladky@gmail.com\>

