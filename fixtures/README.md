# Fixture Suite

This folder contains synthetic Cas-OFFinder fixtures that exercise:

- 3' PAM only
- 5' PAM only
- PAM on both sides
- forward and reverse strand hits
- DNA bulges up to size 2
- RNA bulges up to size 2
- mixed DNA/RNA bulge cases in the same run

## Layout

Each orientation fixture contains:

- `genome.fa`: a synthetic single-chromosome FASTA
- `input.txt`: a Cas-OFFinder input file that points at the fixture directory
- `expected.txt`: required `(Id, BulgeType, Direction, BulgeSize)` tuples

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

The checker only asserts that the expected rows are present. It does not fail if the current build emits additional rows, which is useful while v3 output semantics are still being stabilized.

## GitHub Actions GPU lane

The GPU matrix lane in `.github/workflows/fixtures.yml` expects the repository or organization
variable `GH_GPU_RUNNER_LABEL` to contain the GitHub larger-runner label.

Use the runner label copied from the repository's **Actions > Runners** page.
