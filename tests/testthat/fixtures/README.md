# Test Fixtures

These files are deliberately tiny offline fixtures for package tests.

- `sequencing_summary_dorado.tsv`: Dorado/MinKNOW-like sequencing summary rows.
- `sequencing_summary_dorado_reduced.tsv`: reduced Dorado `summary` rows
  without `passes_filtering`.
- `sequencing_summary_minknow_pod5.tsv`: current MinKNOW-style POD5-era
  sequencing summary rows with additional online-run columns.
- `barcoding_summary.tsv`: matching barcode assignments for the summary rows.
- `reads.fastq`: three nanopore-like reads with Phred qualities.
- `reads.fasta`: the same three read sequences in FASTA format.
- `pod5_manifest.tsv`: POD5 metadata manifest rows for tests that need raw-data
  provenance without committing binary POD5 files.

Do not add downloaded POD5 files or large derived artifacts here. Real ONT
POD5 examples must remain opt-in and cached outside the repository.
