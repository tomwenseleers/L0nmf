# Learn and reuse a common channel aggregation

Cluster outcome channels using size-adjusted Poisson distances and Ward
linkage. The resulting map can be reused across samples or learned from
reference spectra supplied as rows with `normalize=TRUE`.

## Usage

``` r
channel_groups(
  y,
  k = 32L,
  channels = colnames(y),
  normalize = FALSE,
  ncores = NULL
)

aggregate_channels(y, groups, channels = colnames(y))
```

## Arguments

- y:

  Numeric matrix with channels in columns.

- k:

  Maximum number of channel groups.

- channels:

  Unique channel identifiers, normally m/z values.

- normalize:

  Give every input row unit total, useful for spectral libraries.

- ncores:

  Native thread count, defaulting to physical cores.

- groups:

  A map returned by `channel_groups()`.

## Value

`channel_groups()` returns a compact reusable map accepted as
`precluster` by
[`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
and
[`coupledNMF()`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md).
`aggregate_channels()` returns a matrix whose columns sum the mapped
channels. Unseen channels are collected in the map's reserved group.
Save maps with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

## Examples

``` r
y <- matrix(rpois(100*12, 3), 100, 12)
groups <- channel_groups(y, k=4L, ncores=1L)
z <- aggregate_channels(y, groups)
stopifnot(sum(y) == sum(z))
```
