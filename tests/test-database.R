# runs after build database
context("Test the database file")

# setwd('tests/')

test_that("Database file exists.", {
    expect_true(file.exists('../results/LINPS.sqlite'))
})

test_that("Can connect to the database.", {
    con <- DBI::dbConnect(RSQLite::SQLite(),
                          dbname = "../results/LINPS.sqlite")
    
    expect_true(DBI::dbIsValid(con))
    
    DBI::dbDisconnect(con)
})

test_that("Contain tables in the correct format.", {
    con <- DBI::dbConnect(RSQLite::SQLite(),
                          dbname = "../results/LINPS.sqlite")
    expect_equal(DBI::dbListTables(con),
                 c('bif', 'diff_expr', 'models', 'networks', 'nodes', 'perturbations'))

    DBI::dbDisconnect(con)
})
