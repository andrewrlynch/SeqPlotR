# Auto-print a SeqPlot

Renders the plot to the current graphics device. Matches ggplot2's
convention that a bare `seq_plot() %+% ...` expression in the console or
a knitr chunk draws itself.

## Usage

``` r
# S3 method for class 'SeqPlot'
print(x, ...)
```

## Arguments

- x:

  A `SeqPlot` (i.e. a `SeqPlotR6` instance).

- ...:

  Ignored.

## Value

The plot, invisibly.
