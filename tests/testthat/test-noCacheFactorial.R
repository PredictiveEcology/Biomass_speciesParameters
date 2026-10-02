## Reading the factorial .rds files in Init is not Cache()d: the cache copy would be the same .rds
## (measured 2026-10-02 on cohortDataFactorial_medium.rds, 6 GB in memory: readRDS of the input 34.3 s,
## of the cache copy 32.9 s), so a hit saves about a second while a miss, on every new cachePath,
## spends ~6 min writing another 1.6 GB.
test_that("Init does not Cache() the factorial table reads", {
  f <- testthat::test_path("..", "..", "Biomass_speciesParameters.R")
  skip_if_not(file.exists(f), "module source not available")
  env <- new.env()
  exprs <- parse(f, keep.source = FALSE)
  for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("Init")))
    eval(e, env)
  initBody <- paste(deparse(body(env$Init)), collapse = "\n")
  expect_false(grepl("prepInputs_cohortDataFactorial", initBody, fixed = TRUE))
  expect_false(grepl("prepInputs_speciesTableFactorial", initBody, fixed = TRUE))
})
