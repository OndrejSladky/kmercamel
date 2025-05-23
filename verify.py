#!/usr/bin/env python3
import subprocess
import sys
import os
import argparse


def verify_instance(fasta_path: str, k: int, algorithm: str, complements: bool, masked_superstring: str, assume_simplitigs: bool) -> bool:
    """
    Check if running superstring algorithm on given fasta file produces the same set of k-mers as the original one.
    """
    with open("./bin/kmercamel.fa", "w") as k_mers:
        args = ["./kmercamel"] + (["maskopt"] if masked_superstring else (["compute"] + ([] if algorithm != "global" else ["-M", "./bin/kmercamel-maxone.fa"]))) +["-k", f"{k}", "-t" if masked_superstring else "-a", algorithm] + ([] if complements else ["-u"]) + (["-S"] if assume_simplitigs else []) + [(masked_superstring if masked_superstring != "" else fasta_path)]
        subprocess.run(args, stdout=k_mers, stderr=subprocess.DEVNULL)
    
    for s in ["", "-maxone"]:
        with open(f"./bin/converted{s}.fa", "w") as converted:
            subprocess.run(["./kmercamel", "ms2spss", "-k", f"{k}", f"./bin/kmercamel{s}.fa"], stdout=converted, stderr=subprocess.DEVNULL)

    # in result; in result with maxone, in original sequence; in merged file; in merge file for maxone
    stats = [{}, {}, {}, {}, {}]
    runs = [
        (0, "./bin/converted.fa", "converted", complements),
        (2, fasta_path, "original", complements),
    ]
    if algorithm == "global":
        runs.append((1, "./bin/converted-maxone.fa", "converted-maxone", complements))
    for i, path, result, pass_complements in runs:
        args = ["jellyfish", "count", "-m", f"{k}", "-s", "100M", "-o", f"./bin/{result}.jf", path]
        if pass_complements:
            args.insert(2, "-C")
        subprocess.run(args)
        with open(f"./bin/{result}_stats.txt", "w") as f:
            subprocess.run(["jellyfish", "stats", f"./bin/{result}.jf"], stdout=f)
        with open(f"./bin/{result}_stats.txt", "r") as f:
            for _ in range(4):
                key, value = f.readline().split()
                stats[i][key] = value
    # Count k-mers on merged file.
    res = True
    for offset, s, text in [(0, "", "default"), (1, "-maxone", "maxone")][:(2 if algorithm == "global" else 1)]:
        subprocess.run(["jellyfish", "merge", "-o", f"./bin/merged{s}.jf", f"./bin/converted{s}.jf", "./bin/original.jf"])
        with open(f"./bin/merged_stats{s}.txt", "w") as f:
            subprocess.run(["jellyfish", "stats", f"./bin/merged{s}.jf"], stdout=f)
        with open(f"./bin/merged_stats{s}.txt", "r") as f:
            for _ in range(4):
                key, value = f.readline().split()
                stats[3 + offset][key] = value
        distinct_key = "Distinct:"
        total_key = "Total:"
        if stats[offset][distinct_key] != stats[2][distinct_key] or stats[offset][distinct_key] != stats[3 + offset][distinct_key]:
            print("F")
            print(f"Failed: k={k} M={text}: expected orginal_distinct_count={stats[2][distinct_key]}, result_distinct_count={stats[offset][distinct_key]} and merged_distinct_count={stats[3 + offset][distinct_key]} to be equal.")
            res = False
        elif offset == 0 and complements and stats[0][distinct_key] != stats[0][total_key] and masked_superstring == "":
            print("W")
            print(f"Warning: k={k} M={text}: number of masked k-mers={stats[0][total_key]} is not minimal possible (minimum is {stats[0][distinct_key]}).")
        else:
            print(".", end="")
            sys.stdout.flush()
    return res


def main():
    # Initialize.
    if not os.path.exists("bin"):
        os.makedirs("bin")

    parser = argparse.ArgumentParser("check whether ./kmercamel (optimize) outputs a superstring which contains the same set of k-mers"
                                     "as the original sequence")
    parser.add_argument("path", help="path to the fasta file on which ./kmers is verified")
    parser.add_argument("--quick", help="if set do not check for full range of k", action="store_true")
    parser.add_argument("--simplitigs", help="if set assume inputs are simplitigs and verify global with this assumtion", action="store_true")
    parser.add_argument("--superstring_path", help="the path to the masked superstring to test masks")
    parser.add_argument("--k", help="the value of k for mask verification")
    args = parser.parse_args()

    success = True
    if args.simplitigs:
        print(f"Testing global from simplitigs:")
        k = int(args.k)
        success &= verify_instance(args.path, k, "global", True, "", True)
        print("")
    elif args.superstring_path is None:
        # Do the tests on superstring algoritms.
        for a in ["global", "local", "globalAC", "localAC", "streaming"]:
            print(f"Testing {a}:")
            for complements in [True, False]:
                for k in ( ([5, 8, 12] if a not in ["local", "global", "streaming"] else [5, 8, 12, 17, 31, 32, 51, 63, 127]) if args.quick else range(2, 128)):
                    success &= verify_instance(args.path, k, a, complements, "", False)
                print("")
    else:
        k = int(args.k)
        # Do the tests on mask optimization.
        for a in ["runs", "ones", "zeros", "runsapprox"]:
            print(f"Testing {a}:")
            success &= verify_instance(args.path, k, a, True, args.superstring_path, False)
            print("")

    # Print status.
    if not success:
        print("Tests failed")
        exit(1)
    print("OK")

main()
