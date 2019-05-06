test_that("Equal class weight", {
  groups=list(left=list(Y=as.factor(c(1,0))), right=list(Y=as.factor(c(0,1))))
  expect_equal(impurity.idx(groups), 0.5)
})

test_that("Best-Case class weight", {
  groups=list(left=list(Y=as.factor(c(1,0))), right=list(Y=as.factor(c(1,1))))
  expect_equal(impurity.idx(groups), 0.25)

  groups=list(left=list(Y=as.factor(c(0,0))), right=list(Y=as.factor(c(1,1))))
  expect_equal(impurity.idx(groups), 0)
})

test_that("Unequal class weight", {
  groups=list(left=list(Y=as.factor(c(0,0,1,1))), right=list(Y=as.factor(c(0,1,1,1))))
  expect_equal(impurity.idx(groups), 0.4375)
})

test_that("String Names", {
  groups=list(left=list(Y=as.factor(c("red", "blue"))),
              right=list(Y=as.factor(c("blue", "blue"))))
  expect_equal(impurity.idx(groups), 0.25)
})
