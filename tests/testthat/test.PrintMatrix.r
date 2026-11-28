test_that("PrintMatrix writes a matrix and a list of matrices",
{
  m <- matrix(c(1, 2.5, 3, 4), 2, dimnames = list(c("a", "b"), c("u", "v")))
  tmp <- withr::local_tempfile(fileext = ".csv")
  PrintMatrix(m, output.file = tmp)
  expect_equal(readLines(tmp), c('"","u","v"', '"a",1,3', '"b",2.5,4'))

  m.list <- list(first = m, second = m * 2)
  expect_null(PrintMatrix(m.list, output.file = tmp))
  expect_equal(readLines(tmp), c("first", "1,3", "2.5,4", "", "second", "2,6", "5,8", ""))
  # output is not left diverted, and nothing is dumped in the working directory
  expect_equal(sink.number(), 0)
  expect_false(file.exists("dumpdata.R"))

  PrintMatrix(unname(m.list), output.file = tmp)
  expect_equal(readLines(tmp)[c(1, 5)], c("1", "2"))
})

test_that("PrintMatrix closes the sink on error",
{
  tmp <- withr::local_tempfile()
  expect_error(PrintMatrix(list(a = 1, b = new.env()), output.file = tmp))
  expect_equal(sink.number(), 0)
})
