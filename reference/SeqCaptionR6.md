# SeqCaption R6 class

Internal R6 generator backing
[`seq_caption()`](http://andrewlynch.io/SeqPlotR/reference/seq_caption.md).

## Super class

`SeqAnnotation` -\> `SeqCaption`

## Public fields

- `text`:

  Caption text.

- `banner`:

  Optional banner text.

- `aesthetics`:

  `SeqAes` with caption and banner styling.

## Methods

### Public methods

- [`SeqCaption$new()`](#method-SeqCaption-initialize)

- [`SeqCaption$prep()`](#method-SeqCaption-prep)

- [`SeqCaption$draw()`](#method-SeqCaption-draw)

- [`SeqCaption$clone()`](#method-SeqCaption-clone)

------------------------------------------------------------------------

### `SeqCaption$new()`

Construct a SeqCaptionR6.

#### Usage

    SeqCaption$new(text, banner = NULL, allowance = 2.8, aesthetics = aes())

#### Arguments

- `text`:

  Caption text.

- `banner`:

  Optional banner text.

- `allowance`:

  Bottom space reserved, in lines.

- `aesthetics`:

  Optional `SeqAes`.

------------------------------------------------------------------------

### `SeqCaption$prep()`

Place the caption at the centre top of the bottom canvas margin.

#### Usage

    SeqCaption$prep(layout)

#### Arguments

- `layout`:

  The plot's `layout` list (uses `canvasMarginBounds`).

------------------------------------------------------------------------

### `SeqCaption$draw()`

Draw the caption line and the banner under it.

#### Usage

    SeqCaption$draw()

------------------------------------------------------------------------

### `SeqCaption$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SeqCaption$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
