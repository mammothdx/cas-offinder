#!/usr/bin/env python3

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path


IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

RC = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass
class ExpectedRow:
    line_no: int
    row_id: str
    bulge_type: str
    direction: str
    bulge_size: int
    index: int
    seq_rna: str
    seq_dna: str


def read_fasta(path: Path) -> str:
    seq_parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    return "".join(seq_parts)


def reverse_complement(seq: str) -> str:
    return seq.translate(RC)[::-1]


def load_case_inputs(case_dir: Path) -> tuple[Path, str]:
    input_path = case_dir / "input.txt"
    with input_path.open() as handle:
        genome_ref = handle.readline().strip()
        top_row = handle.readline().split()[0].strip().upper()
    genome_path = Path(genome_ref)
    if not genome_path.is_absolute():
        candidates = [
            (case_dir / genome_path).resolve(),
            (case_dir.parent.parent / genome_path).resolve(),
        ]
        for candidate in candidates:
            if candidate.exists():
                genome_path = candidate
                break
        else:
            genome_path = candidates[0]
    return genome_path, top_row


def parse_expected(case_dir: Path) -> list[ExpectedRow]:
    rows = []
    expected_path = case_dir / "expected.txt"
    with expected_path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 7:
                raise ValueError(f"{expected_path}:{line_no}: expected 7 columns, got {len(parts)}")
            rows.append(
                ExpectedRow(
                    line_no=line_no,
                    row_id=parts[0],
                    bulge_type=parts[1],
                    direction=parts[2],
                    bulge_size=int(parts[3]),
                    index=int(parts[4]),
                    seq_rna=parts[5],
                    seq_dna=parts[6],
                )
            )
    return rows


def iupac_contains(pattern: str, base: str) -> bool:
    return base.upper() in IUPAC[pattern.upper()]


def validate_bulge(row: ExpectedRow) -> list[str]:
    errors: list[str] = []
    dna_gaps = row.seq_dna.count("-")
    rna_gaps = row.seq_rna.count("-")
    if row.bulge_type == "X":
        if row.bulge_size != 0:
            errors.append("BulgeType X must have BulgeSize 0")
        if dna_gaps != 0 or rna_gaps != 0:
            errors.append("BulgeType X must not contain '-' in SeqDNA or SeqRNA")
    elif row.bulge_type == "DNA":
        if rna_gaps != row.bulge_size:
            errors.append(f"DNA bulge must have {row.bulge_size} '-' in SeqRNA, found {rna_gaps}")
        if dna_gaps != 0:
            errors.append("DNA bulge must have 0 '-' in SeqDNA")
    elif row.bulge_type == "RNA":
        if dna_gaps != row.bulge_size:
            errors.append(f"RNA bulge must have {row.bulge_size} '-' in SeqDNA, found {dna_gaps}")
        if rna_gaps != 0:
            errors.append("RNA bulge must have 0 '-' in SeqRNA")
    else:
        errors.append(f"Unknown BulgeType {row.bulge_type}")
    return errors


def validate_extracted_sequence(case_dir: Path, genome: str, row: ExpectedRow) -> list[str]:
    errors: list[str] = []
    dna_len = len(row.seq_dna.replace("-", ""))
    extracted = genome[row.index : row.index + dna_len]
    if len(extracted) != dna_len:
        return [f"Genome slice out of bounds at index {row.index} for length {dna_len}"]
    if row.direction == "-":
        extracted = reverse_complement(extracted)

    ext_idx = 0
    for col_idx, dna_char in enumerate(row.seq_dna):
        if dna_char == "-":
            continue
        if ext_idx >= len(extracted):
            errors.append(f"SeqDNA consumes beyond extracted sequence at column {col_idx}")
            break
        if dna_char.upper() != extracted[ext_idx]:
            errors.append(
                "SeqDNA mismatch at column "
                f"{col_idx}: expected genome base {extracted[ext_idx]}, found {dna_char}"
            )
        ext_idx += 1
    if ext_idx != len(extracted):
        errors.append(f"SeqDNA consumed {ext_idx} extracted bases, expected {len(extracted)}")
    return errors


def validate_alignment(case_dir: Path, top_row: str, row: ExpectedRow) -> list[str]:
    errors: list[str] = []
    if len(row.seq_rna) != len(row.seq_dna):
        return [f"SeqRNA length {len(row.seq_rna)} != SeqDNA length {len(row.seq_dna)}"]

    top_idx = 0
    for col_idx, (rna_char, dna_char) in enumerate(zip(row.seq_rna, row.seq_dna)):
        pattern_char = None
        if rna_char != "-":
            if top_idx >= len(top_row):
                errors.append(
                    f"Input top row exhausted at column {col_idx} while SeqRNA still has base {rna_char}"
                )
                break
            pattern_char = top_row[top_idx]
            top_idx += 1

        if rna_char != "-" and dna_char != "-":
            dna_base = dna_char.upper()
            if not iupac_contains(pattern_char, dna_base):
                errors.append(
                    f"Input top row mismatch at column {col_idx}: {dna_base} not allowed by {pattern_char}"
                )

            if dna_char.isupper():
                if not iupac_contains(rna_char, dna_base):
                    errors.append(
                        f"Uppercase match mismatch at column {col_idx}: {dna_base} not allowed by {rna_char}"
                    )
            else:
                if iupac_contains(rna_char, dna_base):
                    errors.append(
                        f"Lowercase mismatch at column {col_idx}: {dna_base} still allowed by {rna_char}"
                    )

    if top_idx != len(top_row):
        errors.append(f"Input top row consumed {top_idx} bases, expected {len(top_row)}")
    return errors


def validate_case(case_dir: Path) -> list[str]:
    genome_path, top_row = load_case_inputs(case_dir)
    genome = read_fasta(genome_path)
    errors: list[str] = []
    for row in parse_expected(case_dir):
        row_errors = []
        row_errors.extend(validate_extracted_sequence(case_dir, genome, row))
        row_errors.extend(validate_alignment(case_dir, top_row, row))
        row_errors.extend(validate_bulge(row))
        for error in row_errors:
            errors.append(
                f"{case_dir.name}/expected.txt:{row.line_no} {row.row_id} {row.bulge_type} "
                f"{row.direction} {row.bulge_size} {row.index}: {error}"
            )
    return errors


def main() -> int:
    root = Path(__file__).resolve().parent
    case_dirs = sorted(
        path for path in root.iterdir() if path.is_dir() and (path / "input.txt").exists() and (path / "expected.txt").exists()
    )
    all_errors: list[str] = []
    for case_dir in case_dirs:
        all_errors.extend(validate_case(case_dir))

    if all_errors:
        for error in all_errors:
            print(error)
        print(f"\nValidation failed with {len(all_errors)} error(s) across {len(case_dirs)} fixture case(s).")
        return 1

    print(f"Validated {len(case_dirs)} fixture case(s) with no errors.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
