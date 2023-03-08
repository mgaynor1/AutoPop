test_that("gametes samples is correct", {
 expect_warning(cytotype_repro_mate(cgen =c(100, 100, 100),
                      b = 0.02,
                      cc = 0.052,
                      gnum.vec=c(10, 10, 10),
                      s = 0.1,
                      mc = 1)
  , "NaNs produced")

})

