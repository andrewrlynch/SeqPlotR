# Links: Arches, SV Reconstruction, and Cross-Track Connectors

This vignette walks through the **link** elements — drawables that
connect two genomic loci. The arch family stays within one track:

- `seq_arc` — Bezier arch between two loci, no stems.
- `seq_arch` — same, with vertical stems and partial-window stubs.
- `seq_recon` — `seq_arch` plus automatic SV classification by strand
  pair and chromosome.

The cross-track family connects two *different* tracks:

- `seq_string` — smooth Bezier curve whose shape (C vs. S) is inferred
  from the strand pair.
- `seq_synteny` — filled trapezoid between homologous / syntenic blocks
  in two tracks.
- `seq_zoom` — four-corner polygon projecting a region from one track
  (e.g. an overview) onto another (e.g. a detail view).

``` r

library(SeqPlotR)
#> 
#> Attaching package: 'SeqPlotR'
#> The following object is masked from 'package:base':
#> 
#>     %||%
library(GenomicRanges)
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo

win <- GRanges("chr1", IRanges(1e6, 3e6))
```

Most examples use simulated structural variants, copy-number segments
and SNVs shipped with the package. A small helper reads one table of a
scenario, and the junctions of two scenarios are loaded once here: a
simulated chromothripsis on chr3 (`ct`) and a simulated chromothripsis
followed by breakage-fusion-bridge cycles on chr4 (`ct_bfb`).

``` r

sim <- function(sc, f) {
  utils::read.delim(system.file("extdata", "sim", sc, paste0(f, ".tsv"),
                                package = "SeqPlotR"))
}

ct_sv  <- sim("ct", "svs")
bfb_sv <- sim("ct_bfb", "svs")
ct_win  <- GRanges("chr3", IRanges(9.5e7, 1.95e8))
bfb_win <- GRanges("chr4", IRanges(5.8e7, 1.24e8))
```

### Anatomy of a `seq_link`

Every link has two anchors. SeqPlotR encodes them in a single `data`
argument with paired column names — there is no `data2`/`mapping2`.
`data` may be:

1.  A `data.frame` (BEDPE-like). All anchor fields come from named
    columns.
2.  A `GRanges` whose `mcols` carry the anchor-1 fields. The GRanges
    contributes anchor 0 via its `seqnames` / `start` / `strand`; anchor
    1 is read from `mcols` columns named in
    [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md).
3.  A plain `GRanges` (single locus). Anchor 0 is `start`, anchor 1 is
    `end` of the same range — used for compact within-locus arcs.

The [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md)
vocabulary is **always explicit**:

| field | meaning | required? |
|----|----|----|
| `x0`, `x1` | genomic position of each anchor | always |
| `chrom0`, `chrom1` | chromosome / seqname of each anchor | for data.frame; auto for GRanges |
| `strand0`, `strand1` | strand of each anchor | for `seq_recon`; optional otherwise |
| `y0`, `y1` | data-scale baselines | optional (default `0`) |
| `height` | arch peak height in data-scale units | optional |

GRanges accessors `seqnames` and `strand` are injected as
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) specials
alongside `start` / `end` / `width` / `mid`, so a typical
GRanges-with-mcols call can use bare names:
`map(chrom0 = seqnames, strand0 = strand, ...)`.

#### Where links live

- **Inside a `seq_track`** (added via `%+%`): both `t0` and `t1` are
  locked to that track. Use this for within-track arcs that summarise
  intra-region relationships.
- **At the `seq_plot` level** (added via `%+%`): with both `t0` and `t1`
  unset, the link attaches to the most recently added track (just like
  any other element). With `t0` / `t1` set, the link is **deferred** —
  stored on the plot and drawn last, after every track has been laid
  out.

### `seq_arc` — within-track arch, no stems

`seq_arc` draws a single Bezier arch between two loci. The simplest form
encodes the loci as `start` / `end` of one GRanges row. Here each
junction of a simulated chromothripsis-plus-BFB event on chr4 becomes
one range, and its junction copy number sets the arch height.

``` r

arc_gr <- GRanges(bfb_sv$chrom1, IRanges(bfb_sv$pos1, bfb_sv$pos2),
                  jcn = bfb_sv$jcn)

seq_plot() %|%
  seq_track(track_id = "A",
            data = arc_gr,
            mapping = map(x = start, y = jcn),
            windows = bfb_win) %+%
  seq_arc(map(x0 = start, x1 = end, height = jcn),
          aesthetics = aes(color = "#4385BE", linewidth = 1.2)) -> p
p$plot()
```

![](links_files/figure-html/arc-simple-1.png)

`map(height = jcn)` makes each arch’s peak as tall as the junction copy
number; the peak heights are interpreted in the track’s `yscale`.

#### Curve and orientation

`aes(curve = ...)` controls how high the arch bulges:

- `"length"` (default) — bulge scales with the genomic span of the link.
- `"equal"` — fixed bulge of 0.2 data-units regardless of span.
- a numeric value — used directly as the bulge offset.

`aes(orientation = "+")` arches up (default), `"-"` arches down.

``` r

seq_plot() %|%
  seq_track(track_id = "A",
            data = arc_gr,
            mapping = map(x = start, y = jcn),
            windows = bfb_win) %+%
  seq_arc(map(x0 = start, x1 = end, height = jcn),
          aesthetics = aes(color = "#AF3029",
                           linewidth = 1.5,
                           curve = "equal")) -> p
p$plot()
```

![](links_files/figure-html/arc-curve-1.png)

### `seq_arch` — arch with stems and stubs

`seq_arch` adds vertical stems from each anchor’s baseline (`y0` / `y1`)
up to the arch endpoints. When one anchor falls outside every visible
window, `seq_arch` instead draws a half-height **stub** at the visible
end with a small angled hook and a label naming the partner chromosome.

#### A BEDPE-style data.frame

The junction table of a simulated chromothripsis with BFB cycles on chr4
is already BEDPE-like. Mapping the junction copy number to the arch
height makes the high-copy fold-back inversion stand out. Heights are in
track-relative units, so jcn is rescaled to 0.2-1 first. The 2 kb
fold-back inversions are too narrow to show an arch and appear as single
stems.

``` r

bfb_sv$h <- 0.2 + 0.8 * bfb_sv$jcn / max(bfb_sv$jcn)

seq_plot() %|%
  seq_track(track_id = "A", windows = bfb_win) %+%
  seq_arch(data = bfb_sv,
           mapping = map(x0 = pos1, x1 = pos2,
                         chrom0 = chrom1, chrom1 = chrom2,
                         height = h),
           aesthetics = aes(color = "#205EA6", linewidth = 1.2)) -> p
p$plot()
```

![](links_files/figure-html/arch-bedpe-1.png)

#### Stubs for half-visible links

When an anchor lies outside the visible window, `seq_arch` draws a stub
labelled with the partner chromosome. Toggle with
`aes(plotStubs = TRUE/FALSE)`; control hook geometry with
`aes(stubAngle = 45, stubLength = 0.02)`.

``` r

# Two in-window pairs and three half-visible (anchor 1 off the right edge,
# or onto chr2).
stub_df <- data.frame(
  chr1 = "chr1",
  start1 = c(1.15e6, 1.40e6, 1.85e6, 2.20e6, 2.60e6),
  chr2 = c("chr1", "chr1", "chr1", "chr2", "chr1"),
  start2 = c(1.55e6, 2.05e6, 4.50e6, 5.00e5, 4.20e6),
  height = c(0.7, 0.5, 0.8, 0.6, 0.9),
  stringsAsFactors = FALSE
)

seq_plot() %|%
  seq_track(track_id = "A", windows = win) %+%
  seq_arch(data = stub_df,
           mapping = map(x0 = start1, x1 = start2,
                         chrom0 = chr1, chrom1 = chr2,
                         height = height),
           aesthetics = aes(color   = "#66800B",
                            linewidth = 1.1,
                            stubLength = 0.04)) -> p
p$plot()
```

![](links_files/figure-html/arch-stubs-1.png)

### `seq_recon` — SV reconstruction by strand pair

`seq_recon` extends `seq_arch` and colours each arch by the strand pair
on its two anchors:

| strand pair        | class                  | tier | default colour          |
|--------------------|------------------------|------|-------------------------|
| `+/+`              | head-to-head inversion | low  | `flexoki_palette(9)[3]` |
| `-/-`              | tail-to-tail inversion | low  | `flexoki_palette(9)[4]` |
| `-/+`              | tandem duplication     | mid  | `flexoki_palette(9)[1]` |
| `+/-`              | deletion               | mid  | `flexoki_palette(9)[2]` |
| different `chrom`s | translocation          | high | `flexoki_palette(9)[9]` |

Both `strand0` and `strand1` must be supplied via
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) — `seq_recon`
errors at `prep()` time if either is missing.

#### Junctions of a chromothripsis

A simulated chromothripsis on chr3 has four inversions, three
duplications and three deletions. The simulated strands already follow
the `seq_recon` convention, so the table maps directly:

``` r

ct_sv[, c("chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2", "sv_type")]
#>    chrom1      pos1 strand1 chrom2      pos2 strand2 sv_type
#> 1    chr3 100000000       +   chr3 176732097       +     INV
#> 2    chr3 104223287       -   chr3 134031414       +     DUP
#> 3    chr3 118243812       +   chr3 156228134       +     INV
#> 4    chr3 122383459       -   chr3 132882005       -     INV
#> 5    chr3 132882004       +   chr3 156228135       -     DEL
#> 6    chr3 142835471       -   chr3 190000000       +     DUP
#> 7    chr3 165998161       +   chr3 166468693       -     DEL
#> 8    chr3 169361302       +   chr3 176732098       -     DEL
#> 9    chr3 169361303       -   chr3 187244210       +     DUP
#> 10   chr3 187244211       -   chr3 190000001       -     INV
```

Plot with the canonical recon mapping:

``` r

recon_map <- map(x0 = pos1, x1 = pos2,
                 chrom0 = chrom1, chrom1 = chrom2,
                 strand0 = strand1, strand1 = strand2)

seq_plot() %|%
  seq_track(track_id = "SV", windows = ct_win) %+%
  seq_recon(data = ct_sv, mapping = recon_map) -> p
p$plot()
```

![](links_files/figure-html/recon-plot-1.png)

The horizontal guide lines mark the inversion and Dup/Del tiers, with
HH/TT and DEL/DUP labels on the left margin. This scenario has no
translocations, so the TRA tier stays empty and is not drawn.

#### Adding a translocation

To show the translocation tier, one hand-made row joins a chr3
breakpoint to chr5. It is not part of the simulated data.

``` r

tra_row <- data.frame(chrom1 = "chr3", pos1 = 1.45e8, strand1 = "+",
                      chrom2 = "chr5", pos2 = 3.0e7, strand2 = "-")
sv_tra <- rbind(ct_sv[, names(tra_row)], tra_row)

seq_plot() %|%
  seq_track(track_id = "SV", windows = ct_win) %+%
  seq_recon(data = sv_tra, mapping = recon_map) -> p
p$plot()
```

![](links_files/figure-html/recon-tra-1.png)

#### Recoloring classes

Override any default colour via
`aes(h2hColor = ..., t2tColor = ..., dupColor = ..., delColor = ..., transColor = ...)`.
A high-contrast palette on the chromothripsis junctions plus the
hand-made translocation:

``` r

seq_plot() %|%
  seq_track(track_id = "SV", windows = ct_win) %+%
  seq_recon(data = sv_tra, mapping = recon_map,
            aesthetics = aes(h2hColor   = "#AF3029",
                             t2tColor   = "#205EA6",
                             dupColor   = "#66800B",
                             delColor   = "#BC5215",
                             transColor = "#5E409D")) -> p
p$plot()
```

![](links_files/figure-html/recon-recolor-1.png)

#### Subsetting the rendered tiers

`drawClasses` controls which tiers are drawn (and which guides / labels
appear). Drop the translocation tier to focus on the intra-chromosome
calls:

``` r

seq_plot() %|%
  seq_track(track_id = "SV", windows = ct_win) %+%
  seq_recon(data = sv_tra, mapping = recon_map,
            drawClasses = c("Inversion", "Dup/Del")) -> p
p$plot()
```

![](links_files/figure-html/recon-subset-1.png)

## Cross-track connectors

Cross-track links are added at the **plot level** — after both
referenced tracks have been defined — with `t0` / `t1` naming the two
tracks explicitly. Anchors are still encoded in a single `data` object
via the BEDPE-style
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) vocabulary.

A recurring subtlety: when a cross-track link carries its own `data`,
keep track-level `mapping` minimal. Anything in a track’s mapping is
merged into the link’s mapping before evaluation, so expressions like
`map(x = start, y = logR)` on the track will try to find `logR` in the
link’s BEDPE table and error out. The idiomatic fix is to attach the
per-track `data` and `mapping` to the elements that use them, not to the
track itself.

### `seq_string` — smooth curves between two tracks

`seq_string` draws a cubic Bezier curve between anchors in `t0` and
`t1`. The `strand0` / `strand1` fields drive automatic curve-shape
inference:

- matching strands (`+/+`, `-/-`) → **C** curve
- opposing strands (`+/-`, `-/+`) → **S** curve

Force a specific shape with `aes(type = "c")` or `aes(type = "s")`.

#### SV breakpoints across copy-number and SNV tracks

A simulated chromothripsis with BFB cycles on chr4, shown as a
copy-number segment track above an SNV allele-fraction track. Each SV
joins its first breakpoint on the copy-number step to its second
breakpoint in the SNV track. Per-link `y0` / `y1` values anchor each
string to a data-scale y in its track: `y0` sits on the bottom of the
copy-number axis and `y1` on the top of the VAF axis, so the curves run
through the gap between the tracks.

``` r

bfb_seg <- sim("ct_bfb", "segments")
bfb_seg <- bfb_seg[bfb_seg$end > 5.8e7 & bfb_seg$start < 1.24e8, ]
seg_gr  <- GRanges(bfb_seg$chrom, IRanges(bfb_seg$start, bfb_seg$end),
                   cn = bfb_seg$cn)

snv <- sim("ct_bfb", "snvs")
snv <- snv[snv$pos > 5.8e7 & snv$pos < 1.24e8, ]
snv_gr <- GRanges(snv$chrom, IRanges(snv$pos, width = 1), vaf = snv$vaf)

bfb_sv$cn0  <- 1
bfb_sv$vaf1 <- 0.8
```

``` r

seq_plot() %|%
  seq_track(track_id = "CN", windows = bfb_win, track_height = 1.2) %+%
    seq_segment(data = seg_gr,
                mapping = map(x = start, x_end = end, y = cn, y_end = cn),
                aesthetics = aes(color = "#205EA6", linewidth = 2)) %__%
  seq_track(track_id = "SNV", windows = bfb_win, track_height = 1.2) %+%
    seq_point(data = snv_gr, mapping = map(x = start, y = vaf),
              aesthetics = aes(color = "#66800B", size = 0.25,
                               alpha = 0.5)) %+%
  seq_string(data = bfb_sv,
             map(x0 = pos1, x1 = pos2, chrom0 = chrom1, chrom1 = chrom2,
                 strand0 = strand1, strand1 = strand2,
                 y0 = cn0, y1 = vaf1),
             t0 = "CN", t1 = "SNV",
             aesthetics = aes(color = "#AF3029", linewidth = 1.1,
                              alpha = 0.75)) -> p
p$plot()
```

![](links_files/figure-html/string-cross-track-1.png)

Each string lands at `(x0, y0)` in the CN track and `(x1, y1)` in the
SNV track. Strands drive the curve shape: deletions and duplications
(opposing strands) bend as S-curves, inversions (matching strands) bend
as C-curves.

#### Forcing a shape

Drop the strand fields (or set `aes(type = "c")`) to render every string
as a C-curve regardless of strand:

``` r

seq_plot() %|%
  seq_track(track_id = "CN", windows = bfb_win, track_height = 1.2) %+%
    seq_segment(data = seg_gr,
                mapping = map(x = start, x_end = end, y = cn, y_end = cn),
                aesthetics = aes(color = "#205EA6", linewidth = 2)) %__%
  seq_track(track_id = "SNV", windows = bfb_win, track_height = 1.2) %+%
    seq_point(data = snv_gr, mapping = map(x = start, y = vaf),
              aesthetics = aes(color = "#66800B", size = 0.25,
                               alpha = 0.5)) %+%
  seq_string(data = bfb_sv,
             map(x0 = pos1, x1 = pos2, chrom0 = chrom1, chrom1 = chrom2,
                 y0 = cn0, y1 = vaf1),
             t0 = "CN", t1 = "SNV",
             aesthetics = aes(color = "#BC5215", linewidth = 1.1,
                              alpha = 0.75, type = "c")) -> p
p$plot()
```

![](links_files/figure-html/string-force-c-1.png)

### `seq_synteny` — filled trapezoids between tracks

`seq_synteny` draws a filled quadrilateral connecting a **region** in
`t0` to a homologous region in `t1`. Both edges of each region come from
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md):

- `x0`, `x0_end` — left/right genomic edge in track `t0`
- `x1`, `x1_end` — left/right genomic edge in track `t1`

When `y0` / `y1` are absent from
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md), the polygon
attaches to the bottom edge of `t0`’s inner panel and the top edge of
`t1`’s inner panel — so it fills the gap between the two tracks exactly.

#### Syntenic blocks between two assemblies

Two gene-model tracks representing homologous segments in two
assemblies, with a trapezoid per block coloured by the mapped `fill`
column:

``` r

genes_top <- GRanges("chr1",
  IRanges(start = c(1.20e6, 1.70e6, 2.30e6),
          width = c(2.0e5,  1.5e5,  2.5e5)),
  gene_id    = c("G1", "G2", "G3"),
  feature    = "exon",
  strand_col = c("+", "-", "+"),
  color      = c("#205EA6", "#AF3029", "#66800B")
)
genes_bot <- GRanges("chr1",
  IRanges(start = c(1.25e6, 1.60e6, 2.40e6),
          width = c(1.8e5,  2.0e5,  2.2e5)),
  gene_id    = c("G1p", "G2p", "G3p"),
  feature    = "exon",
  strand_col = c("+", "-", "+"),
  color      = c("#205EA6", "#AF3029", "#66800B")
)

syn_df <- data.frame(
  c0 = "chr1",
  p0 = c(1.20e6, 1.70e6, 2.30e6),
  p0e = c(1.40e6, 1.85e6, 2.55e6),
  c1 = "chr1",
  p1 = c(1.25e6, 1.60e6, 2.40e6),
  p1e = c(1.43e6, 1.80e6, 2.62e6),
  col = c("#205EA6", "#AF3029", "#66800B"),
  stringsAsFactors = FALSE
)

seq_plot() %|%
  seq_track(track_id = "AsmA", windows = win, track_height = 0.9) %+%
    seq_gene(data = genes_top,
             mapping = map(group = gene_id, type = feature,
                           strand = strand_col, color = color)) %__%
  seq_track(track_id = "AsmB", windows = win, track_height = 0.9) %+%
    seq_gene(data = genes_bot,
             mapping = map(group = gene_id, type = feature,
                           strand = strand_col, color = color)) %+%
  seq_synteny(data = syn_df,
              map(x0 = p0, x0_end = p0e, x1 = p1, x1_end = p1e,
                  chrom0 = c0, chrom1 = c1, fill = col),
              t0 = "AsmA", t1 = "AsmB",
              aesthetics = aes(alpha = 0.35, linewidth = 0.3)) -> p
p$plot()
```

![](links_files/figure-html/synteny-1.png)

The trapezoids start at `AsmA`’s bottom inner edge and end at `AsmB`’s
top inner edge, with left / right corners driven by the `x0` / `x0_end`
and `x1` / `x1_end` fields. The per-link `fill` column is picked up from
the mapping; the `aes(alpha = ...)` applies globally to the fill.

### `seq_zoom` — overview-to-detail projection

`seq_zoom` connects the same genomic region as it appears in two tracks
that span different window ranges — a chromosome-scale overview stacked
above a zoomed detail panel, for example. Only `x0` / `x0_end` are
required; `x1` / `x1_end` default to the same region (the most common
case: “show *this* region in both tracks”).

The polygon auto-detects which track is above and attaches to its bottom
edge; the lower track gets the top edge. A small stem offset
(`aes(stemOffset = 0.01)`) nudges the polygon clear of the track
borders.

#### Chromosome-scale overview into a zoomed detail

``` r

# Overview spans 5 Mb; detail zooms into 1.6 – 2.4 Mb.
overview_win <- GRanges("chr1", IRanges(1,      5e6))
detail_win   <- GRanges("chr1", IRanges(1.6e6,  2.4e6))

xs_o  <- seq(1,      5e6,   length.out = 160)
xs_d  <- seq(1.6e6,  2.4e6, length.out = 120)
sig_o <- GRanges("chr1", IRanges(xs_o, width = 1),
                 score = 0.5 + 0.3 * sin((xs_o - 1) / 6e5))
sig_d <- GRanges("chr1", IRanges(xs_d, width = 1),
                 score = 0.5 + 0.3 * sin((xs_d - 1) / 6e5) +
                         rnorm(length(xs_d), 0, 0.03))

zoom_region <- GRanges("chr1", IRanges(1.6e6, 2.4e6))

seq_plot() %|%
  seq_track(track_id = "Overview", windows = overview_win,
            track_height = 0.5) %+%
    seq_area(data = sig_o, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#4385BE", alpha = 0.55)) %__%
  seq_track(track_id = "Detail", windows = detail_win,
            track_height = 1.3) %+%
    seq_area(data = sig_d, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#205EA6", alpha = 0.8)) %+%
  seq_zoom(data = zoom_region,
           map(x0 = start, x0_end = end),
           t0 = "Overview", t1 = "Detail",
           aesthetics = aes(fill = "#878580", alpha = 0.2,
                            color = "#878580", linewidth = 0.4)) -> p
p$plot()
```

![](links_files/figure-html/zoom-1.png)

The grey quadrilateral fans from the narrow slice of the overview track
out to the full width of the detail panel. Multiple zoom polygons are
fine — drop more rows into `zoom_region` to connect more than one region
at once.

#### Projecting different regions onto the two tracks

Supply separate `x1` / `x1_end` fields when the region in `t1` isn’t
identical to the one in `t0` — useful when connecting paralogous or
translocated segments:

``` r

zoom_df <- data.frame(
  c0   = "chr1", p0  = 1.60e6, p0e  = 2.00e6,
  c1   = "chr1", p1  = 1.80e6, p1e  = 2.30e6,
  stringsAsFactors = FALSE
)

seq_plot() %|%
  seq_track(track_id = "Overview", windows = overview_win,
            track_height = 0.5) %+%
    seq_area(data = sig_o, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#4385BE", alpha = 0.55)) %__%
  seq_track(track_id = "Detail", windows = detail_win,
            track_height = 1.3) %+%
    seq_area(data = sig_d, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#205EA6", alpha = 0.8)) %+%
  seq_zoom(data = zoom_df,
           map(x0 = p0, x0_end = p0e, x1 = p1, x1_end = p1e,
               chrom0 = c0, chrom1 = c1),
           t0 = "Overview", t1 = "Detail",
           aesthetics = aes(fill = "#AF3029", alpha = 0.18,
                            color = "#AF3029", linewidth = 0.4)) -> p
p$plot()
```

![](links_files/figure-html/zoom-distinct-regions-1.png)

#### Fanning out to two detail panels with a patchwork layout

Zoom polygons compose naturally with `seq_plot(layout = "...")`. Drop a
wide overview into the top row and pair it with two independent detail
panels below, then emit one `seq_zoom` per detail region. Because
plot-level links validate that `t0` and `t1` refer to tracks already
added to the plot, the zooms go at the end of the chain.

``` r

fan_layout <- "
AAAA
BBCC
"
seq_preview_layout(layout = fan_layout)
```

![](links_files/figure-html/zoom-fanout-preview-1.png)

``` r

# Overview spans chr1:1-5 Mb; B zooms into 1.2-1.8 Mb, C into 3.2-4.0 Mb.
overview_win <- GRanges("chr1", IRanges(1,      5e6))
detail_b_win <- GRanges("chr1", IRanges(1.2e6,  1.8e6))
detail_c_win <- GRanges("chr1", IRanges(3.2e6,  4.0e6))

xs_a  <- seq(1,      5e6,   length.out = 200)
xs_b  <- seq(1.2e6,  1.8e6, length.out =  80)
xs_c  <- seq(3.2e6,  4.0e6, length.out =  80)
mk_sig <- function(xs, jitter = 0) {
  GRanges("chr1", IRanges(xs, width = 1),
          score = 0.5 + 0.3 * sin((xs - 1) / 5e5) +
                  rnorm(length(xs), 0, jitter))
}
sig_a <- mk_sig(xs_a, jitter = 0.00)
sig_b <- mk_sig(xs_b, jitter = 0.03)
sig_c <- mk_sig(xs_c, jitter = 0.03)

zoom_b <- GRanges("chr1", IRanges(1.2e6, 1.8e6))
zoom_c <- GRanges("chr1", IRanges(3.2e6, 4.0e6))

seq_plot(layout = fan_layout) %+%
  seq_track(track_id = "A", windows = overview_win) %+%
    seq_area(data = sig_a, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#4385BE", alpha = 0.55)) %+%
  seq_track(track_id = "B", windows = detail_b_win) %+%
    seq_area(data = sig_b, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#66800B", alpha = 0.8)) %+%
  seq_track(track_id = "C", windows = detail_c_win) %+%
    seq_area(data = sig_c, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#AF3029", alpha = 0.8)) %+%
  seq_zoom(data = zoom_b, map(x0 = start, x0_end = end),
           t0 = "A", t1 = "B",
           aesthetics = aes(fill = "#66800B", alpha = 0.22,
                            color = "#66800B", linewidth = 0.5)) %+%
  seq_zoom(data = zoom_c, map(x0 = start, x0_end = end),
           t0 = "A", t1 = "C",
           aesthetics = aes(fill = "#AF3029", alpha = 0.22,
                            color = "#AF3029", linewidth = 0.5)) -> p
p$plot()
```

![](links_files/figure-html/zoom-fanout-1.png)

Each zoom polygon inherits the color of its destination track, so the
left and right fan-outs are visually distinguishable at a glance. The
same idiom scales to more panels — add more letters to the layout
string’s second row and a matching `seq_zoom` per region.

### `seq_highlight` — multi-track region annotation

`seq_highlight` draws a continuous filled band that passes through every
track from `t0` down to `t1` (inclusive) at the same genomic
coordinates; widths in NPC follow each track’s own scale, so the band
naturally fans or compresses across tracks with different windows.
Unlike `seq_zoom` (which is a polygon between a pair of tracks),
`seq_highlight` is a single shape that crosses many panels, intended for
highlighting one or more loci of interest across a stack — the typical
ChIP-seq / ATAC-seq use case.

#### Stacked tracks with shared windows

``` r

win <- GRanges("chr1", IRanges(1, 1e6))
xs  <- seq(1, 1e6, length.out = 400)
mk_sig <- function(seed, freq, jitter = 0.03) {
  set.seed(seed)
  GRanges("chr1", IRanges(xs, width = 1),
          score = 0.5 + 0.3 * sin((xs - 1) / freq) +
                  rnorm(length(xs), 0, jitter))
}
sig_a <- mk_sig(1, 1.5e5)
sig_b <- mk_sig(2, 2.0e5)
sig_c <- mk_sig(3, 2.5e5)

# Two genomic regions to flag across all three tracks.
hl <- GRanges("chr1", IRanges(c(2.0e5, 6.5e5),
                              c(2.6e5, 7.4e5)))

seq_plot() %|%
  seq_track(track_id = "A", windows = win) %+%
    seq_area(data = sig_a, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#4385BE", alpha = 0.7)) %__%
  seq_track(track_id = "B", windows = win) %+%
    seq_area(data = sig_b, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#205EA6", alpha = 0.7)) %__%
  seq_track(track_id = "C", windows = win) %+%
    seq_area(data = sig_c, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#24837B", alpha = 0.7)) %+%
  seq_highlight(data = hl, map(x0 = start, x0_end = end),
                t0 = "A", t1 = "C",
                aesthetics = aes(fill = "#AF3029", alpha = 0.18)) -> p
p$plot()
```

![](links_files/figure-html/highlight-stack-1.png)

Each highlight region produces a single band that spans all three
tracks. Per-track rectangles share their inner panel y-range, and
trapezoids in the gaps between tracks bridge each pair.

#### Different scales: bands fan and compress

When the tracks span different genomic ranges, the highlight follows
each track’s own scale, producing a fanning shape across the scale
change. Pair an overview with a zoomed detail, drop a single
`seq_highlight`, and the band traces the zoom relationship visually:

``` r

overview_win <- GRanges("chr1", IRanges(1,      5e6))
detail_win   <- GRanges("chr1", IRanges(1.6e6,  2.4e6))

xs_o  <- seq(1,      5e6,   length.out = 200)
xs_d  <- seq(1.6e6,  2.4e6, length.out = 120)
sig_o <- GRanges("chr1", IRanges(xs_o, width = 1),
                 score = 0.5 + 0.3 * sin((xs_o - 1) / 6e5))
sig_d <- GRanges("chr1", IRanges(xs_d, width = 1),
                 score = 0.5 + 0.3 * sin((xs_d - 1) / 6e5))

# A single 1.8 - 2.2 Mb region of interest.
hl <- GRanges("chr1", IRanges(1.8e6, 2.2e6))

seq_plot() %|%
  seq_track(track_id = "Overview", windows = overview_win,
            track_height = 0.5) %+%
    seq_area(data = sig_o, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#4385BE", alpha = 0.55)) %__%
  seq_track(track_id = "Detail", windows = detail_win,
            track_height = 1.3) %+%
    seq_area(data = sig_d, mapping = map(x = start, y = score),
             aesthetics = aes(fill = "#205EA6", alpha = 0.8)) %+%
  seq_highlight(data = hl, map(x0 = start, x0_end = end),
                t0 = "Overview", t1 = "Detail",
                aesthetics = aes(fill = "#878580", alpha = 0.22,
                                 color = "#878580", linewidth = 0.4)) -> p
p$plot()
```

![](links_files/figure-html/highlight-fan-1.png)

The band is narrow inside the overview (a 5 Mb window) and wide inside
the zoomed detail (an 800 kb window covering the same region), with a
trapezoid bridging the two.
