test_that("bulk native projection preserves storage, inputs and thread parity", {
  H <- Matrix::sparseMatrix(i=c(2:6,3:8),j=c(rep(1,5),rep(2,6)),
    x=c(1,5,1,4,1,2,1,3,8,1,2),dims=c(9L,2L))
  values <- H@x
  project <- getFromNamespace("slcn_project_design_cpp", "logconNMF")
  local <- getFromNamespace("spolo_lc", "logconNMF")
  a <- project(H,1L);b <- project(H,2L)
  expect_identical(a,b);expect_identical(H@x,values)
  expect_identical(a@i,H@i);expect_identical(a@p,H@p)
  for(j in 1:2) {
    at <- seq.int(H@p[j]+1L,H@p[j+1L])
    expect_identical(a@x[at],local(H@x[at])$h)
  }
})
