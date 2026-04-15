#!/usr/bin/env python3
"""
Brute-force cas-offinder validator.

Independently computes every valid off-target match for a fixture case,
then verifies the tool's actual output matches exactly (ignoring row order).

Usage:
    python3 bruteforce_check.py <case_dir> <actual_output_file>

Exit 0 if actual == brute-force (set equality), exit 1 otherwise.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# IUPAC helpers
# ---------------------------------------------------------------------------

IUPAC: dict[str, set[str]] = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

# Reverse-complement: only handle ACGT (genome bases are always ACGT).
# Guide/pattern chars use IUPAC but are never RC'd in the output logic;
# only the genome (DNA) sequence is reverse-complemented.
RC = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(RC)[::-1]


def iupac_match(pattern_char: str, base: str) -> bool:
    return base.upper() in IUPAC.get(pattern_char.upper(), set())


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> str:
    parts: list[str] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts.append(line.upper())
    return "".join(parts)


@dataclass
class CaseInput:
    genome_path: Path
    pam_pattern: str
    dna_bulge: int
    rna_bulge: int
    queries: list[tuple[str, int, str]]  # (sequence, threshold, id)


def parse_input(case_dir: Path) -> CaseInput:
    input_path = case_dir / "input.txt"
    with input_path.open() as f:
        genome_ref = f.readline().strip()
        pam_parts = f.readline().strip().split()
        pam_pattern = pam_parts[0].upper()
        dna_bulge = int(pam_parts[1]) if len(pam_parts) > 1 else 0
        rna_bulge = int(pam_parts[2]) if len(pam_parts) > 2 else 0

        queries: list[tuple[str, int, str]] = []
        for line_no, line in enumerate(f):
            line = line.strip()
            if not line:
                break
            parts = line.split()
            seq = parts[0].upper()
            threshold = int(parts[1])
            qid = parts[2] if len(parts) > 2 else str(line_no)
            queries.append((seq, threshold, qid))

    genome_path = Path(genome_ref)
    if not genome_path.is_absolute():
        for candidate in [
            (case_dir / genome_path).resolve(),
            (case_dir.parent.parent / genome_path).resolve(),
        ]:
            if candidate.exists():
                genome_path = candidate
                break
        else:
            genome_path = (case_dir / genome_ref).resolve()

    return CaseInput(genome_path, pam_pattern, dna_bulge, rna_bulge, queries)


# ---------------------------------------------------------------------------
# Compare-pattern generation  (mirrors C++ parseInput)
# ---------------------------------------------------------------------------

@dataclass
class CompareEntry:
    compare: str       # padded compare string, len == patternlen
    bulge_label: str   # "0" for X, "1"/"2" for DNA, or deleted-bases string for RNA
    bulge_index: int   # position in the query where bulge occurs


def generate_compares(
    query: str,
    dna_bulge: int,
    rna_bulge: int,
    is_reversed_pam: bool,
) -> list[CompareEntry]:
    entries: list[CompareEntry] = []

    first_non_n = next((i for i, c in enumerate(query) if c != "N"), len(query))
    last_non_n = len(query) - 1 - next(
        (i for i, c in enumerate(reversed(query)) if c != "N"), len(query)
    )

    if dna_bulge == 0:
        entries.append(CompareEntry(query, "0", 0))
    else:
        if is_reversed_pam:
            tmp = query + "N" * dna_bulge
        else:
            tmp = "N" * dna_bulge + query

        for bulge_sz in range(1, dna_bulge + 1):
            pre_n = dna_bulge - bulge_sz
            for j in range(first_non_n, last_non_n + 2):
                if is_reversed_pam:
                    comp = (
                        query[:j]
                        + "N" * bulge_sz
                        + query[j:]
                        + "N" * pre_n
                    )
                else:
                    comp = (
                        "N" * pre_n
                        + query[:j]
                        + "N" * bulge_sz
                        + query[j:]
                    )
                bl = "0" if comp == tmp else str(bulge_sz)
                entries.append(CompareEntry(comp, bl, j))

    for bulge_sz in range(1, rna_bulge + 1):
        pre_n = dna_bulge + bulge_sz
        for j in range(first_non_n, last_non_n + 1 - bulge_sz):
            deleted = query[j : j + bulge_sz]
            if is_reversed_pam:
                comp = query[:j] + query[j + bulge_sz :] + "N" * pre_n
            else:
                comp = "N" * pre_n + query[:j] + query[j + bulge_sz :]
            entries.append(CompareEntry(comp, deleted, j))

    return entries


# ---------------------------------------------------------------------------
# Output-row generation  (mirrors C++ compareAll)
# ---------------------------------------------------------------------------

def trim_search_padding(seq: str, trim: int, trim_right: bool) -> str:
    if trim >= len(seq):
        return ""
    return seq[: len(seq) - trim] if trim_right else seq[trim:]


def output_alignment_start_offset(
    trim: int, trim_right: bool, direction: str
) -> int:
    if direction == "-":
        return trim if trim_right else 0
    return 0 if trim_right else trim


def indicate_mismatches(seq: str, compare: str) -> str:
    result = list(seq)
    for k, (s, c) in enumerate(zip(seq, compare)):
        if c != "N" and not iupac_match(c, s):
            result[k] = s.lower()
    return "".join(result)


# A result row: the 7-tuple that check-fixtures.sh extracts and compares.
ResultRow = tuple[str, str, str, int, int, str, str]
# (id, bulge_type, direction, bulge_size, location, seq_rna, seq_dna)


def make_row(
    genome: str,
    loci: int,
    direction: str,
    mm_count: int,
    compare: str,
    entry: CompareEntry,
    is_reversed_pam: bool,
    dna_bulge: int,
    patternlen: int,
    qid: str,
) -> ResultRow | None:
    """Build an output row exactly as C++ compareAll does, then validate it."""

    strbuf = genome[loci : loci + patternlen]
    if direction == "-":
        strbuf = reverse_complement(strbuf)
    strbuf = indicate_mismatches(strbuf, compare)

    trim_right = is_reversed_pam

    try:
        bulge_size = int(entry.bulge_label)
        is_numeric = True
    except ValueError:
        is_numeric = False
        bulge_size = len(entry.bulge_label)

    bi = entry.bulge_index

    if is_numeric:
        trim = dna_bulge - bulge_size
        offset = output_alignment_start_offset(trim, trim_right, direction)
        sr = trim_search_padding(compare, trim, is_reversed_pam)
        sr = sr[:bi] + "-" * bulge_size + sr[bi + bulge_size :]
        sd = trim_search_padding(strbuf, trim, trim_right)
        bt = "X" if bulge_size == 0 else "DNA"
    else:
        trim = dna_bulge + bulge_size
        offset = output_alignment_start_offset(trim, trim_right, direction)
        sr = trim_search_padding(compare, trim, is_reversed_pam)
        sr = sr[:bi] + entry.bulge_label + sr[bi:]
        sd = trim_search_padding(strbuf, trim, trim_right)
        sd = sd[:bi] + "-" * bulge_size + sd[bi:]
        bt = "RNA"

    location = loci + offset

    # --- self-checks (catches bugs in this script itself) ---

    if len(sr) != len(sd):
        return None  # alignment length mismatch

    # Lowercase count must equal mismatch count
    lowercase_count = sum(1 for c in sd if c.islower())
    if lowercase_count != mm_count:
        return None

    return (qid, bt, direction, bulge_size, location, sr, sd)


# ---------------------------------------------------------------------------
# PAM validation at the alignment level
# ---------------------------------------------------------------------------

def validate_pam(pam_pattern: str, seq_rna: str, seq_dna: str) -> bool:
    """
    Walk through SeqRNA / SeqDNA, consuming one PAM-pattern position per
    non-gap RNA character.  At every position where the PAM pattern is
    non-N and both RNA and DNA are non-gap, verify the DNA base is
    allowed by the PAM.
    """
    top_idx = 0
    for rna_ch, dna_ch in zip(seq_rna, seq_dna):
        if rna_ch == "-":
            continue  # gap in RNA  (DNA bulge) -- no PAM position consumed
        if top_idx >= len(pam_pattern):
            return False
        pat_ch = pam_pattern[top_idx]
        top_idx += 1
        if dna_ch == "-":
            continue  # gap in DNA  (RNA bulge) -- PAM position consumed but no DNA base
        dna_base = dna_ch.upper()
        if pat_ch != "N" and not iupac_match(pat_ch, dna_base):
            return False
    return top_idx == len(pam_pattern)


# ---------------------------------------------------------------------------
# Extracted-sequence validation
# ---------------------------------------------------------------------------

def validate_extracted(genome: str, location: int, direction: str,
                       seq_dna: str) -> bool:
    """Verify SeqDNA bases match the genome at the reported location."""
    dna_bases = seq_dna.replace("-", "")
    dna_len = len(dna_bases)
    if location < 0 or location + dna_len > len(genome):
        return False
    extracted = genome[location : location + dna_len]
    if direction == "-":
        extracted = reverse_complement(extracted)
    for ext_ch, dna_ch in zip(extracted, dna_bases):
        if dna_ch.upper() != ext_ch:
            return False
    return True


# ---------------------------------------------------------------------------
# Brute-force search
# ---------------------------------------------------------------------------

def count_mismatches(seq: str, compare: str) -> int:
    mm = 0
    for s, c in zip(seq, compare):
        if c != "N" and not iupac_match(c, s):
            mm += 1
    return mm


def brute_force(
    genome: str,
    pam_pattern: str,
    dna_bulge: int,
    rna_bulge: int,
    queries: list[tuple[str, int, str]],
) -> set[ResultRow]:
    is_reversed_pam = pam_pattern[0] != "N"
    patternlen = len(pam_pattern) + dna_bulge
    results: set[ResultRow] = set()

    for query, threshold, qid in queries:
        entries = generate_compares(query, dna_bulge, rna_bulge, is_reversed_pam)

        for entry in entries:
            assert len(entry.compare) == patternlen, (
                f"compare len {len(entry.compare)} != patternlen {patternlen}"
            )
            comp_rc = reverse_complement(entry.compare)

            for i in range(len(genome) - patternlen + 1):
                subseq = genome[i : i + patternlen]

                for direction, comp in [("+", entry.compare), ("-", comp_rc)]:
                    mm = count_mismatches(subseq, comp)
                    if mm > threshold:
                        continue

                    row = make_row(
                        genome, i, direction, mm, entry.compare,
                        entry, is_reversed_pam, dna_bulge, patternlen, qid,
                    )
                    if row is None:
                        continue

                    _, bt, d, bs, loc, sr, sd = row

                    # PAM validation at the alignment level
                    if not validate_pam(pam_pattern, sr, sd):
                        continue

                    # Genome extraction validation
                    if not validate_extracted(genome, loc, d, sd):
                        continue

                    results.add(row)

    return results


# ---------------------------------------------------------------------------
# Parse actual tool output
# ---------------------------------------------------------------------------

def parse_actual(output_file: Path) -> set[ResultRow]:
    results: set[ResultRow] = set()
    with output_file.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            qid = parts[0]
            bt = parts[1]
            sr = parts[2]
            sd = parts[3]
            loc = int(parts[5])
            direction = parts[6]
            mm = int(parts[7])
            bs = int(parts[8])

            # Validate: lowercase count == reported mismatches
            lowercase_count = sum(1 for c in sd if c.islower())
            if lowercase_count != mm:
                print(
                    f"WARNING: row {qid} {bt} {direction} {bs} {loc}: "
                    f"lowercase count {lowercase_count} != reported mm {mm}",
                    file=sys.stderr,
                )

            results.add((qid, bt, direction, bs, loc, sr, sd))
    return results


# ---------------------------------------------------------------------------
# Validate actual output rows
# ---------------------------------------------------------------------------

def validate_actual_rows(
    output_file: Path,
    genome: str,
    pam_pattern: str,
    case_name: str,
) -> list[str]:
    """Run structural validations on every actual output row."""
    errors: list[str] = []
    with output_file.open() as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            qid, bt, sr, sd = parts[0], parts[1], parts[2], parts[3]
            loc, direction = int(parts[5]), parts[6]
            mm, bs = int(parts[7]), int(parts[8])
            prefix = f"{case_name} line {line_no} ({qid} {bt} {direction} bs={bs} loc={loc})"

            # 1) SeqRNA and SeqDNA must be same length
            if len(sr) != len(sd):
                errors.append(f"{prefix}: SeqRNA len {len(sr)} != SeqDNA len {len(sd)}")
                continue

            # 2) Bulge structure
            dna_gaps = sd.count("-")
            rna_gaps = sr.count("-")
            if bt == "X":
                if bs != 0:
                    errors.append(f"{prefix}: X must have bulge_size 0")
                if dna_gaps or rna_gaps:
                    errors.append(f"{prefix}: X must have no gaps")
            elif bt == "DNA":
                if rna_gaps != bs:
                    errors.append(f"{prefix}: DNA bulge needs {bs} gaps in RNA, found {rna_gaps}")
                if dna_gaps:
                    errors.append(f"{prefix}: DNA bulge must have 0 gaps in DNA")
            elif bt == "RNA":
                if dna_gaps != bs:
                    errors.append(f"{prefix}: RNA bulge needs {bs} gaps in DNA, found {dna_gaps}")
                if rna_gaps:
                    errors.append(f"{prefix}: RNA bulge must have 0 gaps in RNA")

            # 3) Lowercase count == mismatch count
            lowercase_count = sum(1 for c in sd if c.islower())
            if lowercase_count != mm:
                errors.append(
                    f"{prefix}: lowercase bases {lowercase_count} != mismatches {mm}"
                )

            # 4) Extracted genome sequence matches SeqDNA
            if not validate_extracted(genome, loc, direction, sd):
                errors.append(f"{prefix}: SeqDNA does not match genome at location")

            # 5) PAM validation at alignment level
            if not validate_pam(pam_pattern, sr, sd):
                errors.append(f"{prefix}: PAM mismatch in alignment")

            # 6) Alignment IUPAC check (each non-gap position)
            top_idx = 0
            for col, (rch, dch) in enumerate(zip(sr, sd)):
                if rch == "-":
                    continue
                if top_idx >= len(pam_pattern):
                    errors.append(f"{prefix}: col {col} ran past PAM pattern")
                    break
                top_idx += 1
                if dch == "-":
                    continue
                db = dch.upper()
                if dch.isupper() and not iupac_match(rch, db):
                    errors.append(
                        f"{prefix}: col {col} uppercase {dch} not matched by RNA {rch}"
                    )
                if dch.islower() and iupac_match(rch, db):
                    errors.append(
                        f"{prefix}: col {col} lowercase {dch} still matched by RNA {rch}"
                    )

    return errors


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <case_dir> <actual_output_file>", file=sys.stderr)
        return 1

    case_dir = Path(sys.argv[1])
    actual_file = Path(sys.argv[2])

    ci = parse_input(case_dir)
    genome = read_fasta(ci.genome_path)

    # --- Phase 1: validate every actual row structurally ---
    row_errors = validate_actual_rows(
        actual_file, genome, ci.pam_pattern, case_dir.name
    )
    if row_errors:
        for e in row_errors:
            print(e, file=sys.stderr)
        print(
            f"\n{case_dir.name}: {len(row_errors)} structural error(s) in actual output",
            file=sys.stderr,
        )
        return 1

    # --- Phase 2: brute-force compute and set-compare ---
    bf_results = brute_force(
        genome, ci.pam_pattern, ci.dna_bulge, ci.rna_bulge, ci.queries
    )
    actual_results = parse_actual(actual_file)

    missing = bf_results - actual_results
    extra = actual_results - bf_results

    ok = True
    if missing:
        ok = False
        for row in sorted(missing):
            print(
                f"MISSING from actual: {' '.join(str(x) for x in row)}",
                file=sys.stderr,
            )
    if extra:
        ok = False
        for row in sorted(extra):
            print(
                f"EXTRA in actual: {' '.join(str(x) for x in row)}",
                file=sys.stderr,
            )

    if not ok:
        print(
            f"\n{case_dir.name}: FAIL "
            f"({len(missing)} missing, {len(extra)} extra, "
            f"{len(bf_results)} expected, {len(actual_results)} actual)",
            file=sys.stderr,
        )
        return 1

    print(
        f"{case_dir.name}: OK "
        f"({len(actual_results)} actual rows match "
        f"{len(bf_results)} brute-force results)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
