#!/usr/bin/env python3
import subprocess
import sys
import os


def verify_instance(fasta_path: str, k: int, algorithm: str, complements: bool):
    with open("./bin/kmers.txt", "w") as k_mers:
        args = ["./kmers", "-p", fasta_path, "-k", f"{k}", "-a", algorithm]
        if complements:
            args.append("-c")
        subprocess.run(args, stdout=k_mers)
    with open("./bin/kmers.txt", "r") as k_mers:
        with open("./bin/converted.fa", "w") as converted:
            subprocess.run(["./convert_superstring.py"], stdin=k_mers, stdout=converted)
    stats = [{}, {}, {}]
    for i, path, result in [(0,"./bin/converted.fa", "converted"), (1, fasta_path, "original")]:
        args = ["jellyfish", "count", "-m", f"{k}", "-s", "100M", "-o", f"./bin/{result}.jf", path]
        if complements:
            args.insert(2, "-C")
        subprocess.run(args)
        with open(f"./bin/{result}_stats.txt", "w") as f:
            subprocess.run(["jellyfish", "stats", f"./bin/{result}.jf"], stdout=f)
        with open(f"./bin/{result}_stats.txt", "r") as f:
            for _ in range(4):
                key, value = f.readline().split()
                stats[i][key] = value
    subprocess.run(["jellyfish", "merge", "-o", f"./bin/merged.jf", "./bin/converted.jf", "./bin/original.jf"])
    with open(f"./bin/merged_stats.txt", "w") as f:
        subprocess.run(["jellyfish", "stats", f"./bin/merged.jf"], stdout=f)
    with open(f"./bin/merged_stats.txt", "r") as f:
        for _ in range(4):
            key, value = f.readline().split()
            stats[2][key] = value

    if (stats[0]["Distinct:"] != stats[1]["Distinct:"] or stats[0]["Distinct:"] != stats[2]["Distinct:"]):
        print("F")
        print(f"k={k}: expected orginal_distinct_count={stats[1]['Distinct:']}, result_distinct_count={stats[0]['Distinct:']} and merged_distinct_count={stats[2]['Distinct:']} to be equal.")
    else:
        print(".",end="")
        sys.stdout.flush()
    print(f"k={k}, a={algorithm}: orig={stats[1]['Distinct:']} result={stats[0]['Distinct:']} merged={stats[2]['Distinct:']}", file=sys.stderr)


if not os.path.exists("bin"):
    os.makedirs("bin")
path = sys.argv[1]
print(f"Testing with support for reverse complements:")
for a in ["pseudosimplitigs"]:
    print(f"Testing {a}:")
    for k in range(5, 32):
        verify_instance(path, k, a, True)
print(f"Testing without support for reverse complements:")
for a in ["greedy", "greedyAC", "pseudosimplitigs", "pseudosimplitigsAC"]:
    print(f"Testing {a}:")
    for k in range(5, 32):
        verify_instance(path, k, a, False)


