test_that("Equal class weight", {
  groups=list(as.factor(c(1,0)), as.factor(c(0,1)))
  expect_equal(importance.idx(groups), 0.5)
})

test_that("Best-Case class weight", {
  groups=list(as.factor(c(1,0)), as.factor(c(1,1)))
  expect_equal(importance.idx(groups), 0.25)

  groups=list(as.factor(c(0,0)), as.factor(c(1,1)))
  expect_equal(importance.idx(groups), 0)
})

test_that("Unequal class weight", {
  groups=list(as.factor(c(0,0, 1, 1)), as.factor(c(0, 1, 1, 1)))
  expect_equal(importance.idx(groups), 0.4375)
})

test_that("String Names", {
  groups=list(as.factor(c("red", "red", "blue", "blue")), as.factor(c("red", "blue", "blue", "blue")))
  expect_equal(importance.idx(groups), 0.4375)
})
