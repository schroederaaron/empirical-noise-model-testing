# Loaded by testthat before the test files. Run from the directory holding
# external/docker_r_libs:  Rscript -e 'GOF_ROOT <- "Distribution_GOF"; testthat::test_dir(file.path(GOF_ROOT, "tests/testthat"))'
if (!exists("GOF_ROOT")) GOF_ROOT <- normalizePath(Sys.getenv("GOF_ROOT", file.path("..", "..")))
if (!exists("gof_all")) source(file.path(GOF_ROOT, "R", "setup.R"))
