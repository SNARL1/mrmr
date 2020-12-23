
# mrmr: mark recapture miscellany in R

![example branch
parameter](https://github.com/SNARL1/mrmr/workflows/R-CMD-check/badge.svg)

This package automates common data processing steps for mark recapture
data, with an emphasis on data collected by folks at the Sierra Nevada
Aquatic Research Lab (SNARL). Many of the assumptions made in the data
processing and modeling are pinned to the type of capture-recapture data
collected by the SNARL team, and as a result, this package should not be
thought of as a general purpose mark-recapture modeling toolbox (for
that, see the RMark, unmarked, multimark or Rcapture packages).

## Installation

``` r
remotes::install_github("SNARL1/mrmr")
```

## Usage

This package consists of three core functions to use consecutively:

1.  `clean_data()` ingests capture-recapture, introduction, and survey
    data.
2.  `fit_model()` fits a Jolly-Seber mark recapture model to quickly
    estimate demographic parameters and population abundance.
3.  `plot_model()` generates plots of abundance and recruitment over
    time, or of the survival of introduced cohorts.

For an in-depth look at how to use these functions, see the vignette [An
Introduction to the mrmr
package](https://snarl1.github.io/mrmr/articles/intro-to-mrmr.html).

## Background

Customized mark recapture modeling can be time and resource intensive,
especially when there are unique features of a study system that prevent
or complicate the immediate use of a general purpose package. For
example, the model described in Joseph and Knapp 2018 **Ecosphere**
(<https://doi.org/10.1002/ecs2.2499>) - which accounts for introductions
into a population, and continuous, time-varying, incompletely observed
individual-level covariates on detection and survival - takes
considerable effort to develop and deploy. There is a need in this
particular research group for models that can be used in near-real time
as data are available to inform conservation decisions. This package
fills this need, with the primary goal of estimating abundance and
recruitment through time in montane amphibian populations.

## License

This package is licensed under the MIT license.

Copyright 2019 Maxwell B. Joseph

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
“Software”), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
