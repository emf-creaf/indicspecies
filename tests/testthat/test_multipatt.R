library(indicspecies)

data(wetland)
groups <- c(rep(1, 17), rep(2, 14), rep(3,10))

wetkm <- kmeans(wetland, centers=3)
groupskm <- wetkm$cluster

test_that("Can run multipatt with group combinations",{
  expect_s3_class(multipatt(wetland, groups, control = how(nperm=999), func = "IndVal.g"), "multipatt")   
  expect_s3_class(multipatt(wetland, groups, control = how(nperm=999), func = "IndVal"), "multipatt")   
  expect_s3_class(multipatt(wetland, groups, control = how(nperm=999), func = "r"), "multipatt")   
  expect_s3_class(multipatt(wetland, groups, control = how(nperm=999), func = "r.g"), "multipatt")   
})

test_that("Can run multipatt without group combinations",{
  expect_s3_class(multipatt(wetland, groups, control = how(nperm=999), func = "IndVal.g", duleg = TRUE), "multipatt")   
  expect_s3_class(multipatt(wetland, groups, control = how(nperm=999), func = "IndVal", duleg = TRUE), "multipatt")   
  expect_s3_class(multipatt(wetland, groups, control = how(nperm=999), func = "r", duleg = TRUE), "multipatt")   
  expect_s3_class(multipatt(wetland, groups, control = how(nperm=999), func = "r.g", duleg = TRUE), "multipatt")   
})

test_that("Coverage can be calculated", {
  indval <- multipatt(wetland, groups, control = how(nperm=999))
  expect_type(coverage(wetland, indval), "double")
})