# Fixture Suite

This folder contains synthetic Cas-OFFinder fixtures that exercise:

- 3' PAM only
- 5' PAM only
- PAM on both sides
- N-only padding on both sides
- forward and reverse strand hits
- DNA bulges up to size 2
- RNA bulges up to size 2
- mixed DNA/RNA bulge cases in the same run

## Layout

Each orientation fixture contains:

- `genome.fa`: a synthetic single-chromosome FASTA
- `input.txt`: a Cas-OFFinder input file that points directly to the fixture's `genome.fa` file path
- `expected.txt`: required `(Id, BulgeType, Direction, BulgeSize, Index, SeqRNA, SeqDNA)` rows

Only the first sequence row in `input.txt` carries the PAM. The guide rows below it mask PAM positions with `N`, which matches how we want the fixtures to exercise the search logic rather than baking PAM bases into every guide.

Each fixture `genome.fa` must contain exactly one FASTA target. The validator concatenates all sequence lines within that single target and does not treat multiple FASTA records as separate chromosomes.

## Running

Build `cas-offinder`, then run:

```bash
./fixtures/check-fixtures.sh ./cas-offinder
```

To force the OpenCL backend selector, pass `C` for CPU or `G` for GPU:

```bash
./fixtures/check-fixtures.sh ./cas-offinder C
./fixtures/check-fixtures.sh ./cas-offinder G
```

The checker asserts that the expected rows are present with the exact emitted `Index`, `SeqRNA`, and `SeqDNA` alignments. It does not fail if the current build emits additional rows, which is useful while v3 output semantics are still being stabilized.

## GitHub Actions GPU lane

The GPU matrix lane in `.github/workflows/fixtures.yml` expects the repository or organization
variable `GH_GPU_RUNNER_LABEL` to contain the GitHub larger-runner label.

Use the runner label copied from the repository's **Actions > Runners** page.
