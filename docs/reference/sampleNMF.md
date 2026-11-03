# Extract one sample's factorization from a coupled fit

Reconstruct a standard \[splogconNMF()\] object on a sample's original
grid. Its profiles include alignment and exposure, and have unit column
sums; the matching spectra carry their original-grid amplitudes.
\`compound_id\` maps its local components to the coupled consensus
compounds.

## Usage

``` r
sampleNMF(object, sample = 1L)
```

## Arguments

- object:

  A \[coupledNMF()\] fit.

- sample:

  One sample index or name.

## Value

A \`splogconNMF\` object without an attached observed or fitted matrix.
