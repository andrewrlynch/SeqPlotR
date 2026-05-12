# Add a track, element, or plot-level feature to a SeqPlot or SeqTrack

The single composition operator for SeqPlotR. Dispatches on the class of
the right-hand side:

## Usage

``` r
e1 %+% e2
```

## Arguments

- e1:

  A `SeqPlot` or `SeqTrack` object.

- e2:

  A `SeqTrack`, `SeqElement`, `SeqLink`, or `SeqAnnotation`.

## Value

`e1`, invisibly modified in place (R6 reference semantics).

## Details

- `SeqTrack` -\> added to the plot layout

- `SeqElement` -\> added to the current (last) track

- `SeqLink` -\> stored in `plot$plot_links` (plot-level, deferred)

- `SeqAnnotation` -\> stored in `plot$plot_annotations` (plot-level,
  deferred)

When the LHS is a `SeqTrack`, only `SeqElement` and `SeqLink` are
accepted. Links added via `%+%` on a `SeqTrack` have their `t0` and `t1`
locked to the parent track's `track_id` — no override is permitted.
