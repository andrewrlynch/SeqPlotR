# SeqPlotR vignette example data from ROSE simulated scenarios

Plain tab-separated files (read with `utils::read.delim`, then build `GRanges`). Chromosome names are UCSC style ("chr2", not "2"). Positions are 1-based, closed intervals. No VCF parser or ROSE dependency is needed to use them. Regenerate with `export_simdata.R` (needs ROSE loaded via devtools::load_all).

## Provenance
Source: /Users/andrew/code/ROSE/inst/extdata/<scenario>/ (PURPLE-style simulated output). Parsed with `ROSE::read_purple(dir, sample)`, which reads `<sc>.purple.cnv.somatic.tsv`, `.purple.purity.tsv`, `.purple.sv.vcf` and `.purple.somatic.vcf.gz`. `baf` in segments.tsv comes from the raw cnv table; `truth_*` SNV columns and events.tsv come from `<sc>.truth.snvs.tsv` and `<sc>.truth.events.tsv` (simulation ground truth, not observable in real data).

| dir | ROSE scenario |
|-----|---------------|
| bfb | sim_bfb_amplicon |
| wgd | sim_wgd_losses |
| ct | sim_chromothripsis |
| ct_bfb | sim_ct_bfb |

## Files (per directory)
- segments.tsv: chrom, start, end, cn (total CN), cn_major, cn_minor, baf (PURPLE BAF = major/cn), support_start, support_end (PURPLE segment boundary support: TELOMERE, BND, ...).
- svs.tsv: one row per junction. chrom1, pos1, strand1, chrom2, pos2, strand2 (NA for single breakends; none occur in these four), sv_type, jcn (PURPLE junction copy number), alt_reads (summed over both breakends), id. Mates are ordered by position (pos1 <= pos2). strand is ROSE's breakend orientation and follows the seq_recon convention: DUP is -/+, DEL is +/-, INV is +/+ or -/-. sv_type is ROSE's classification (INV, DUP, DEL; BND/INS/SGL would appear if present). The sv_type values are not recomputed.
- snvs.tsv: chrom, pos, ref, alt, vaf (alt_reads/depth), alt_reads, depth, sbs6 (C>A, C>G, C>T, T>A, T>C, T>G, pyrimidine-collapsed by reverse complement), clonal (ROSE flag), truth_mult (true copies carrying the SNV), truth_tau (true time, 0-1), id. ROSE's PURPLE input provides no observed multiplicity, so truth_mult is the simulation truth.
- purity.tsv: purity, ploidy (average CN), wgd (0/1 from PURPLE field).
- events.tsv (small, 2-5 rows): id, type (SV, GAIN, LOSS, WGD), footprint (bfb, chromothripsis, gain, loss, wgd), time (relative time 0-1, early to late), chrom, start, end (affected region, NA for WGD), haplotype (A/B), cluster (SV cluster id, NA if none).

## Profile
Sizes: bfb 124K, ct 128K, ct_bfb 96K, wgd 96K (total of 5 files each, all under 200 KB). SNV tables are whole.

### bfb (sim_bfb_amplicon)
- chr2 only, 1-242,193,529. purity 0.8, ploidy 2.16, no WGD.
- 12 segments, CN 1-9. Staircase amplicon at chr2:179,998,001-200,000,000: CN 2 -> 4 -> 6 -> 7 -> 8 -> 9 (peak at 185,998,001-191,998,000) -> 8 -> 7 -> 5 -> 3 -> 2, then CN 1 (cn_minor 0, LOH) from 200,000,001 to the end. Segments outside are CN 2.
- 5 SVs, all INV, all within the amplicon (179,998,001-200,000,000, a 20 Mb span), jcn 1-2: fold-back inversions (-/- at the left, +/+ at the right half).
- 1,671 SNVs, VAF 0-0.88 (median 0.38), truth_mult 1 for 1,625; the amplified region carries mult 2-8 SNVs. Events: 5 BFB cycles, times 0.3-0.7.
- Windows: amplicon `GRanges("chr2", IRanges(1.75e8, 2.05e8))`; whole chromosome for context.

### wgd (sim_wgd_losses)
- chr1 (248,956,422), chr2 (242,193,529), chr3 (198,295,559). purity 0.7, ploidy 3.42, wgd = 1.
- 5 segments, CN 2-6: chr1 CN 3 (1-121.7 Mb, 2+1) then 4 (2+2); chr2 CN 2 (2+0, whole chromosome LOH after WGD); chr3 CN 4 (1-94 Mb) then 6 (94 Mb-end, 4+2).
- 0 SVs (svs.tsv is header only; do not use for link examples).
- 1,236 SNVs, VAF 0-0.91 (median 0.25), truth_mult 1 (891), 2 (331), 4 (14): the multiplicity-2 and 4 classes mark pre-WGD SNVs. chr1 486, chr2 256, chr3 494.
- Events: loss chr2 (t 0.2), gain chr3 94 Mb-end (0.35), WGD (0.5), loss chr1 1-121.7 Mb (0.7).
- Windows: all three chromosomes (multi-window), or `GRanges("chr3", IRanges(1, 198295559))`.

### ct (sim_chromothripsis)
- chr3 only, 1-198,295,559. purity 0.7, ploidy 2.82, no WGD (background CN 3 after a gain).
- 15 segments, CN 1-3. Oscillating CN 3/1 between chr3:100,000,001 and 190,000,000 (CN 1 pieces at 100.0-104.2 Mb, 118.2-122.4, 134.0-142.8, 166.0-166.5 Mb).
- 10 SVs: 4 INV, 3 DUP, 3 DEL, all intra-chr3, 100,000,000-190,000,001, all jcn 2.
- 1,759 SNVs, VAF 0-0.76 (median 0.30), truth_mult 1 (1,234) or 2 (525; pre-gain). Events: chromothripsis (t 0.5), gain (0.75).
- Window: `GRanges("chr3", IRanges(9.5e7, 1.95e8))`.

### ct_bfb (sim_ct_bfb)
- chr4 only, 1-190,214,555. purity 0.8, ploidy 1.63, no WGD.
- 22 segments, CN 1-6. Oscillating CN 1/2 across chr4:60-122 Mb (with LOH, cn_minor 0, in the CN 1 pieces), plus two small BFB staircases at 86.64-89.10 Mb (CN 3-4-5-6) and 98.07-101.61 Mb (CN 4-5-6); CN 1 from 122 Mb to the end.
- 11 SVs: 8 INV, 2 DEL, 1 DUP, all chr4, 60,000,000-122,000,000, jcn 1-5 (one fold-back INV with jcn 5 at 89,098,897 -> 101,612,874).
- 1,265 SNVs, VAF 0-0.88 (median 0.43), truth_mult 1 (1,243), 2-5 (22). Events: chromothripsis (t 0.3) then four BFB cycles (0.45-0.8).
- Windows: `GRanges("chr4", IRanges(5.8e7, 1.24e8))` for the whole rearranged region; zoom `GRanges("chr4", IRanges(8.6e7, 1.02e8))` for the two BFB staircases.

## Notes for authors
- SNV VAF is raw alt/depth (depth around 60-90), so noisy; in amplified or WGD segments VAF clusters near mult/cn x purity-adjusted levels.
- jcn can be used as a link height or linewidth. Coordinates are simulated, round numbers near the amplicon edges.
- Build GRanges with `GRanges(df$chrom, IRanges(df$start, df$end), cn = df$cn, ...)`. For BEDPE SVs pass the data.frame directly to seq_arch/seq_recon/seq_string with map(x0 = pos1, x1 = pos2, chrom0 = chrom1, chrom1 = chrom2, strand0 = strand1, strand1 = strand2).
