library(indicspecies)

data(wetland)
groups <- c(rep(1, 17), rep(2, 14), rep(3,10))

wetkm <- kmeans(wetland, centers=3)
groupskm <- wetkm$cluster

test_that("Can run signassoc",{
  expect_type(signassoc(wetland, cluster=groups,  alternative = "two.sided", 
                        control = how(nperm=199)), "double")   
  expect_type(signassoc(wetland, cluster=groups,  alternative = "greater", 
                        control = how(nperm=199)), "double")  
  expect_type(signassoc(wetland, cluster=groups,  alternative = "less", 
                        control = how(nperm=199)), "double")  
})
