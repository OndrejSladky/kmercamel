#!/usr/bin/env python3
import subprocess
import sys
import os


def verify_instance(fasta_path: str, k: int, algorithm: str):
    with open("./bin/kmers.txt", "w") as k_mers:
        subprocess.run(["./kmers", "-p", fasta_path, "-k", f"{k}", "-a", algorithm], stdout=k_mers)
    with open("./bin/kmers.txt", "r") as k_mers:
        with open("./bin/converted.fa", "w") as converted:
            subprocess.run(["./convert_superstring.py"], stdin=k_mers, stdout=converted)
    stats = [{}, {}]
    for i, path, result in [(0,"./bin/converted.fa", "converted"), (1, fasta_path, "original")]:
        subprocess.run(["jellyfish", "count", "-m", f"{k}", "-s", "100M", "-o", f"./bin/{result}.jf", path ])
        with open(f"./bin/{result}_stats.txt", "w") as f:
            subprocess.run(["jellyfish", "stats", f"./bin/{result}.jf"], stdout=f)
        with open(f"./bin/{result}_stats.txt", "r") as f:
            for _ in range(4):
                key, value = f.readline().split()
                stats[i][key] = value
    if (stats[0]["Distinct:"] != stats[1]["Distinct:"]):
        print("F")
        print(f"k={k}: expected {stats[1]['Distinct:']}, got {stats[0]['Distinct:']} distinct k-mers")
    else:
        print(".",end="")
        sys.stdout.flush()
    print(f"k={k}, a={algorithm}: expected {stats[1]['Distinct:']}, got {stats[0]['Distinct:']} distinct k-mers", file=sys.stderr)


if not os.path.exists("bin"):
    os.makedirs("bin")
path = sys.argv[1]
for a in ["greedy", "greedyAC", "pseudosimplitigs", "pseudosimplitigsAC"]:
    print(f"Testing {a}:")
    for k in range(5, 32):
        verify_instance(path, k, a)


