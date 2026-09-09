#!/usr/bin/env python3
"""CLI regression tests for compute -S -A; run after make kmercamel."""
import pathlib
import subprocess
import tempfile
import unittest

BINARY = pathlib.Path(__file__).resolve().parents[1] / 'kmercamel'


def canonical(sequence):
    return min(sequence, sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1])


class AlignmentTest(unittest.TestCase):
    def test_reconstruction(self):
        cases = [
            (3, ['TAC', 'ACG', 'GGC']),  # Full and shorter overlaps.
            (3, ['TACG']),              # No joins.
            (3, ['AAA', 'CCC']),        # Zero overlap.
            (3, ['TAC', 'ACG', 'CGG']), # Consecutive full overlaps.
            (4, ['ACAA', 'ATTT', 'AACA']),
            (1, ['A', 'C', 'G', 'T']),
            (31, ['A' * 31, 'A' * 30 + 'C']),
            (63, ['A' * 63, 'A' * 62 + 'C']),
            (127, ['A' * 127, 'A' * 126 + 'C']),
        ]
        for k, sequences in cases:
            for complements in (False, True):
                with self.subTest(k=k, sequences=sequences, complements=complements), tempfile.TemporaryDirectory() as directory:
                    directory = pathlib.Path(directory)
                    fasta = directory / 'input.fa'
                    alignment = directory / 'alignment.txt'
                    fasta.write_text(''.join(f'>{i}\n{s}\n' for i, s in enumerate(sequences)))
                    args = [str(BINARY), 'compute', '-k', str(k), '-S', '-A', str(alignment)]
                    if not complements:
                        args.append('-u')
                    result = subprocess.run(args + [str(fasta)], capture_output=True, text=True, check=True)
                    ms = ''.join(result.stdout.splitlines()[1:])
                    mask = alignment.read_text().strip()
                    self.assertEqual(len(ms), len(mask))
                    self.assertLessEqual(set(mask), {'0', '1'})
                    recovered = []
                    start = None
                    for i in range(len(ms) + 1):
                        active = i < len(ms) and ms[i].isupper()
                        if start is not None and (not active or mask[i] == '1'):
                            recovered.append(ms[start:i + k - 1].upper())
                            start = None
                        if active and start is None:
                            start = i
                    normalize = canonical if complements else lambda s: s
                    self.assertCountEqual(map(normalize, sequences), map(normalize, recovered))
                    # Alignment output must not change the computed superstring.
                    plain = subprocess.run([a for a in args if a not in ('-A', str(alignment))] + [str(fasta)], capture_output=True, text=True, check=True)
                    self.assertEqual(result.stdout, plain.stdout)

    def test_requires_simplitigs(self):
        with tempfile.TemporaryDirectory() as directory:
            alignment = pathlib.Path(directory) / 'alignment.txt'
            result = subprocess.run([str(BINARY), 'compute', '-k', '3', '-A', str(alignment), 'unused.fa'], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("requires '-S'", result.stderr)
            self.assertFalse(alignment.exists())

    def test_unwritable_alignment(self):
        result = subprocess.run([str(BINARY), 'compute', '-k', '3', '-S', '-A', '/dev/null/alignment.txt', 'unused.fa'], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Cannot open alignment output', result.stderr)


if __name__ == '__main__':
    unittest.main()
