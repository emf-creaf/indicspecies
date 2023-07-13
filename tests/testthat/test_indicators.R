library(indicspecies)

data(wetland)
groups <- c(rep(1, 17), rep(2, 14), rep(3,10))

wetkm <- kmeans(wetland, centers=3)
groupskm <- wetkm$cluster

test_that("Can run indicators",{
  expect_s3_class(ind1 <- indicators(X=wetland, cluster=groups, group=2, 
                             max.order = 3, verbose=FALSE, 
                             At=0.5, Bt=0.2), "indicators")   
  expect_s3_class(ind2 <- indicators(X=wetland, cluster=groups, group=2, 
                             max.order = 3, verbose=FALSE, 
                             At=0.5, Bt=0.2, nboot.ci = 100), "indicators") 
  expect_no_error(print(ind1))
  expect_no_error(print(ind2))
})
