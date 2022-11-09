#!/usr/bin/env python3

from typing import List


"""
Split the given superstring into string where at every position (except for the last k-1) begins a k-mer.
"""
def split_superstring(superstring: str, k: int) -> List[str]:
    ret: List[str] = []
    next_start = 0
    was_present = False
    for i, c in enumerate(superstring):
        if c.islower():
            if was_present:
                ret.append(superstring[next_start : i + (k - 1)].upper())
            next_start = i + 1
            was_present = False
        else:
            was_present = True
    return ret



"""
Determine the k given a masked superstring.
"""
def get_k(superstring: str) -> int:
    ret = 1
    while ret <= len(superstring):
        # If the k-th position from the end is the first with capital letter, it means there starts the last k-mer of size k.
        if superstring[-ret].isupper():
            return ret
        ret += 1
    raise ValueError("Given superstring contains no valid k-mer.")



"""
Read the superstring from standard input and print the SPSS representation.
"""
def main():
    # Ignore fasta record header.
    input()
    superstring = input().strip()
    k = get_k(superstring)
    split = split_superstring(superstring, k)
    for i, s in enumerate(split):
        print(f">{i+1}")
        print(s)


if __name__ == "__main__":
    main()

