TnSeq
-----

# Overview

An implementation of the algorithm described in
"A Zero-Inflated Poisson Model for Insertion Tolerance Analysis of Genes Based on Tn-Seq Data"
Liu F., Wang C., Wu Z., Zhang Q., and Liu P.

# Installation

There are two likely methods of installing this package:

* If one has Hadley's devtools installed, then run the following from R:

> install_github("abelew/TnSeq")

Conversely one may download the archive file, unzip it and run
make and make install

## Usage

TLDR: Look in vignettes/

Take a look at the .Rmd vignettes for the current state of the run-ability of this code.

The files 'case_study_cov.Rmd' and 'case_study_nocov.Rmd' are the originals.

em_nocovariance.Rmd is the result of a rewrite of the code and includes a couple
simplistic plots of the result.

em_oo.Rmd is using a second rewrite of the code to use S3 objects in order to make
packaging the data required for the EM algorithm a bit more flexible.  For the moment
it likely does not work, but once it does it should make it much easier to add in
the covariance model as well as some other old code from the 'essentials' package
which I could never get to work.
