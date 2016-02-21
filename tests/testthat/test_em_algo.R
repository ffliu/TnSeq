library(testthat)
library(TnSeq)
context("Test usability of TnSeq")

datafile <- system.file("mtb_example.tsv", package="TnSeq")
tnseq_df <- read.table(datafile, header=TRUE)
previous_result_file <- system.file("testrun_nocov.rda", package="TnSeq")
load(previous_result_file)

tnseq_result <- perform_tnseq_em(tnseq_df)
previous_result <- as.data.frame(testrun_nocov)

tnseq_ess <- as.numeric(tnseq_result$ess)
previous_ess <- as.numeric(previous_result$ess)
test_that("Is the known 'ess' column the same as the new one?", {
    expect_equal(previous_ess, tnseq_ess)
})

tnseq_sd <- as.numeric(tnseq_result$sd)
previous_sd <- as.numeric(previous_result$sd)
test_that("Is the known 'sd' column the same as the new one?", {
    expect_equal(previous_sd, tnseq_sd)
})

tnseq_p <- tnseq_result$pvalue
previous_p <- previous_result$pvalue
test_that("Is the known 'pvalue' column the same as the new one?", {
    expect_equal(previous_p, tnseq_p)
})

tnseq_q <- tnseq_result$qvalue
previous_q <- previous_result$qvalue
test_that("Is the known 'qvalue' column the same as the new one?", {
    expect_equal(previous_q, tnseq_q)
})


