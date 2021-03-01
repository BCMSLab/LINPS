# runs before build database
context("Testing the output in the results/ directory.")

test_that("Correct directory structure exist.", {
    
    expect_true(dir.exists('../results/'))
})
