context("test-confusion.R")

test_that("false_negative", {
  expect_false(false_negative(predicted = NULL, actual = NULL))
  expect_false(false_negative(predicted = c(1L), actual = NULL))
  expect_true(false_negative(predicted = NULL, actual = c(1L)))
  expect_false(false_negative(predicted = c(1L), actual = c(1L)))
  expect_false(false_negative(predicted = c(1L, 2L), actual = c(1L)))
  expect_true(false_negative(predicted = c(1L), actual = c(1L, 2L)))
})

test_that("false_positive", {
  expect_false(false_positive(predicted = NULL, actual = NULL))
  expect_true(false_positive(predicted = c(1L), actual = NULL))
  expect_false(false_positive(predicted = NULL, actual = c(1L)))
  expect_false(false_positive(predicted = c(1L), actual = c(1L)))
  expect_true(false_positive(predicted = c(1L, 2L), actual = c(1L)))
  expect_false(false_positive(predicted = c(1L), actual = c(1L, 2L)))
})
