test_that("pipeline runs end-to-end", {
  # tiny toy df
  tmpdir <- tempdir()
  df <- data.frame(
    gene  = c("G1","G1","G2"),
    Start = c(100,150,200),
    End   = c(200,180,300),
    index = c("chr1_100_200_G1/ENSG1",
              "chr1_150_180_G1/ENSG1",
              "chr1_200_300_G2/ENSG2"),
    Cell1_IRratio = c(0.3, 0.2, 0.5),
    Cell2_IRratio = c(0.4, 0.3, 0.6)
  )

  # run nested_intron
  df2 <- nested_intron(df)
  expect_true(all(c("Nested","Intron_cluster") %in% names(df2)))

  # run irtranscript (short nperms for speed)
  res <- irtranscript(df2, score = c("Cell1_IRratio","Cell2_IRratio"),
                      save_dir= tmpdir,nperms = 100)
  # 1. File should exist
  outfile <- file.path(tmpdir, "Cell1_IRratio_irTranscript.txt")
  expect_true(file.exists(outfile))

  # 2. File should have the expected columns
  df <- read.delim(outfile)
  expect_true(all(c("gene","obs.IRratio","z_score") %in% names(df)))


})
