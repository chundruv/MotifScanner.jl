# MotifScanner

PWM scanning of DNA sequences and of variant REF/ALT alleles (VCF), with per-motif calibrated p-values.

## Scoring VCF variants

```bash
julia --project -t 8 scripts/score_vcf_motifs.jl \
    --vcf dnvs.vcf.gz --fasta GRCh38_no_alt.fa \
    --motifs H14CORE_transfac_format.txt.gz \
    --out dnvs_h14.tsv
```

`--motifs` takes a TRANSFAC file from JASPAR (`JASPAR2026_CORE_non-redundant_pfms_transfac.txt.gz`) or HOCOMOCO v14 (`H14CORE_transfac_format.txt.gz`); plain or gzipped.

| Option | Default | Meaning |
|---|---|---|
| `--taxa` | `vertebrates` | JASPAR `tax_group` filter (comma-separated; `all` disables). Ignored for HOCOMOCO, which is human-only |
| `--tfs` | none | file with one TF gene symbol per line (e.g. TFs expressed in fetal neural cell types); dimers are kept if any partner is listed |
| `--pthresh` | `1e-4` | per-placement, single-strand null p for a site; sets `Call` and `Total*Motifs` |
| `--maxpval` | = `--pthresh` | write a row only if min(RefPval, AltPval) ≤ this. Use `1` to write every motif × variant |
| `--mindlog10p` | `0` | minimum \|Log10P_RefAlt\| |
| `--minprmax`, `--mindeltapr` | `-Inf`, `0` | legacy fraction-of-max filters (not recommended) |
| `--flank` | `40` | bp of reference on each side of the allele |
| `--pseudocount` | `0.01` | added to frequencies before log-odds (unchanged from earlier versions) |
| `--resolution` | `0.01` | bits per integer score step for the null distribution |

### Outputs

- `<out>`: one row per variant × motif passing the filters. Legacy columns are unchanged: `RefMaxScore`, `AltMaxScore`, `RefPrMax`, `AltPrMax`, `LR_RefAlt`, `PR_RefAlt`, … Calibrated columns:
  - `RefPval`, `AltPval`: null p of the best placement on each allele (zero-order human background, single strand, per placement).
  - `RefProbNorm`, `AltProbNorm`: 1 − p, i.e. the CDF value, as in motifDiff "probNorm".
  - `ProbNorm_AltRef`: alt − ref. Negative means loss (motifDiff sign).
  - `Log10P_RefAlt`: log10 AltPval − log10 RefPval. Positive means loss.
  - `Call`: `loss` (ref site, no alt site), `gain`, `retained` (site on both alleles) or `none`, at `--pthresh`.
  - `ThresholdScore`, `PThresh`, `MotifLen`, `RefLen`, `AltLen`, `RefNValid`, `AltNValid`.
- `<out>.motifs.tsv`: motif metadata and calibration: threshold score, `ThresholdPrMax` (its fraction of the max score) and `BestKmerP`. Motifs with `BestKmerP > pthresh` (very short motifs) can never be called at that threshold.
- `<out>.variants.tsv`: every parsed variant, with `RefMatch`, `RowsWritten`, `Calls_loss`, `Calls_gain`. Use it as the denominator: a variant with no rows was scored and has no site.

### Notes

- **Calibration.** `PrMax` (score / max score) is not comparable across motifs. The genome hit rate at PrMax ≥ 0.85 spans more than 1,000-fold across JASPAR vertebrate motifs, and 28% (JASPAR) / 47% (HOCOMOCO) of motifs never reach it in 400 kb of GRCh38. Each motif therefore gets its exact null score distribution by DP over discretised log-odds (`calibrate`). On 400 kb of GRCh38, the median per-placement hit rate at p ≤ 1e-4 is 0.95–1.5 × 10⁻⁴ in every width bin up to 20 bp. Motifs over 20 bp run ~2–3× higher because genomic sequence is not i.i.d.
- **REF check.** Variants whose REF does not match the FASTA are flagged (`RefMatch = false`), counted in the log, and not scored.
- **N bases.** Placements containing a non-ACGT base are excluded; previously such a base scored 0.
- **Indels.** The longer allele is scanned over more placements, so its maximum is biased upwards: deletions lean towards `loss`, insertions towards `gain`. Compare loss/gain only within variant class and indel length (e.g. case vs control within strata); never pool with SNVs.
- `scripts/score_vcf_jaspar.jl` is the previous positional-argument scorer. It applies no taxon filter and uses the uncalibrated PrMax filters.

## Library use

```julia
using MotifScanner
meta, motifs = loadtransfac("H14CORE_transfac_format.txt.gz"; taxa=["vertebrates"], tfs=nothing)
motifs = calibrate(motifs; pthresh=1e-4)
seqs = loadrefseqs(loadvcf("dnvs.vcf"), "GRCh38.fa")
scores = motifscanall(seqs, motifs; pthresh=1e-4, maxpval=1e-4)
```

Tests: `julia --project test/runtests.jl`.
