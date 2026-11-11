# SeqAnnotation R6 base class

Internal R6 generator for plot-level annotations added to a `SeqPlot`
with `%+%`. Subclasses implement `prep(layout)` and `draw()`;
`allowance` is the bottom space (in lines) reserved for them by
`layoutGrid()`.

## Public fields

- `allowance`:

  Bottom space reserved, in lines.

- `coordCanvas`:

  Canvas positions set by `prep()`.

## Methods

### Public methods

- [`SeqAnnotation$prep()`](#method-SeqAnnotation-prep)

- [`SeqAnnotation$draw()`](#method-SeqAnnotation-draw)

- [`SeqAnnotation$clone()`](#method-SeqAnnotation-clone)

------------------------------------------------------------------------

### `SeqAnnotation$prep()`

Override in subclasses. Default implementation errors.

#### Usage

    SeqAnnotation$prep(layout)

#### Arguments

- `layout`:

  The plot's `layout` list.

------------------------------------------------------------------------

### `SeqAnnotation$draw()`

Override in subclasses. Default implementation errors.

#### Usage

    SeqAnnotation$draw()

------------------------------------------------------------------------

### `SeqAnnotation$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqAnnotation$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
