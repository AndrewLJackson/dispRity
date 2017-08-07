context("model.test")

test_that("select.model.list [internal]", {
    data(disparity)
    full_list <- select.model.list(disparity, observed = FALSE)
    expect_is(full_list, "list")
    expect_equal(names(full_list), c("central_tendency", "variance", "sample_size", "subsamples"))
    expect_equal(as.vector(unlist(lapply(full_list, class))), c(rep("numeric", 2), "integer", "numeric"))
    expect_equal(as.vector(unlist(lapply(full_list, length))),rep(7, 4))

    ## Seed
    # set.seed(123) ; data(BeckLee_mat99) ; data(BeckLee_ages) ; data(BeckLee_tree)

    # ## Continuous time slices
    # continous_data <- time.subsamples(BeckLee_mat99, BeckLee_tree, method = "continuous", time = rev(seq(from = 0, to = 120, by = 5)), model = "gradual")

    # ## Bootstrapping the data
    # data_bootstrapped <- boot.matrix(continous_data)

    # ## Getting the true median
    # disp1_true_median <- dispRity(continous_data, metric = median)
    # ## Measuring disparity as the median of the variance in each dimension
    # disp2_median_var <- dispRity(data_bootstrapped, metric = c(variances, median))

    # disp1 <- select.model.list(disp1_true_median, observed = TRUE)
    # disp2 <- select.model.list(disp2_median_var, observed = FALSE)
})

test_that("pooled.variance [internal]", {
    data(disparity)
    full_list <- select.model.list(disparity, observed = FALSE)

    test1 <- pooled.variance(full_list, rescale.variance = TRUE)
    expect_is(test1, "list")
    expect_equal(names(test1), c("central_tendency", "variance", "sample_size", "subsamples"))
    test2 <- pooled.variance(full_list, rescale.variance = FALSE)
    expect_is(test2, "numeric")
    expect_equal(length(test2), 1)
})

test_that("BM.parameters [internal]", {
    dummy_list <- list("central_tendency" = rep(1, 5),
                       "variances" = rep(1, 5), 
                       "sample_size" = rep(10, 5),
                       "subsamples" = seq(1:5))
    test <- BM.parameters(dummy_list)

    expect_is(test, "array")
    expect_equal(test[1], c("anc_state" = 1))
    expect_equal(test[2], c("sigma_squared" = -0.2))
})