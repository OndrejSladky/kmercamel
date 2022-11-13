#!/usr/bin/env python3
import subprocess
import sys
import os


def print_help():
    print("verify accepts a single argument - a path to the fasta file.")
    print("verify then checks whether ./kmers outputs a superstring which contains the same set of k-mers as"
          "the original sequence on this fasta file")
    print("example usage: ./verify.py ./spneumoniae.fa")


def verify_instance(fasta_path: str, k: int, algorithm: str, complements: bool) -> bool:
    """
    Check if running superstring algorithm on given fasta file produces the same set of k-mers as the original one.
    """
    with open("./bin/kmers.txt", "w") as k_mers:
        args = ["./kmers", "-p", fasta_path, "-k", f"{k}", "-a", algorithm]
        if complements:
            args.append("-c")
        subprocess.run(args, stdout=k_mers)
    with open("./bin/kmers.txt", "r") as k_mers:
        with open("./bin/converted.fa", "w") as converted:
            subprocess.run(["./convert_superstring.py"], stdin=k_mers, stdout=converted)
    # in result; in original sequence; in result without complements; in original without complements; in merged file
    stats = [{}, {}, {}, {}]
    runs = [
        (0, "./bin/converted.fa", "converted", complements),
        (1, fasta_path, "original", complements),
        (2, "./bin/converted.fa", "non-canonical", False),
    ]
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
    subprocess.run(["jellyfish", "merge", "-o", f"./bin/merged.jf", "./bin/converted.jf", "./bin/original.jf"])
    with open(f"./bin/merged_stats.txt", "w") as f:
        subprocess.run(["jellyfish", "stats", f"./bin/merged.jf"], stdout=f)
    with open(f"./bin/merged_stats.txt", "r") as f:
        for _ in range(4):
            key, value = f.readline().split()
            stats[3][key] = value
    distinct_key = "Distinct:"
    if stats[0][distinct_key] != stats[1][distinct_key] or stats[0][distinct_key] != stats[3][distinct_key]:
        print("F")
        print(f"Failed: k={k}: expected orginal_distinct_count={stats[1][distinct_key]}, result_distinct_count={stats[0][distinct_key]} and merged_distinct_count={stats[3][distinct_key]} to be equal.")
        return False
    elif complements and stats[0][distinct_key] != stats[2][distinct_key]:
        print("W")
        print(f"Warning: k={k}: number of k-mers={stats[0][distinct_key]} is not equal to number of canonical k-mers={stats[2][distinct_key]} in the superstring.")
    else:
        print(".", end="")
        sys.stdout.flush()
    return True


# Initialize.
if not os.path.exists("bin"):
    os.makedirs("bin")
if len(sys.argv) != 2:
    print_help()
    exit(1)
path = sys.argv[1]

# Do the tests.
success = True
print(f"Testing with support for reverse complements:")
for a in ["greedy", "pseudosimplitigs"]:
    print(f"Testing {a}:")
    for k in range(5, 32):
        success &= verify_instance(path, k, a, True)
    print("")
print(f"Testing without support for reverse complements:")
for a in ["greedy", "greedyAC", "pseudosimplitigs", "pseudosimplitigsAC"]:
    print(f"Testing {a}:")
    for k in range(5, 32):
        success &= verify_instance(path, k, a, False)
    print("")

# Print status.
if not success:
    print("Tests failed")
    exit(0)
print("OK")
