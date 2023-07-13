library(indicspecies)

data(wetland)
groups <- c(rep(1, 17), rep(2, 14), rep(3,10))

wetkm <- kmeans(wetland, centers=3)
groupskm <- wetkm$cluster

test_that("Can run strassoc",{
  for(func in c("r", "r.g", "IndVal", "IndVal.g", "A", "A.g", "B")) expect_type(strassoc(wetland, cluster=groups, func=func), "double")   
})
test_that("Can run strassoc with bootstrapping",{
  expect_type(strassoc(wetland, cluster=groups, func="A.g", nboot.ci = 99), "list")   
})