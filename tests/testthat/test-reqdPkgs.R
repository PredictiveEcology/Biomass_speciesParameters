## With `options(spades.reqdPkgsAttach = FALSE)` (SpaDES.core >= 3.2.1.9030) a module sees only the
## packages in its own `reqdPkgs`: nothing another module lists is attached for it. A function the module
## calls without `pkg::` then has to come from one of its own reqdPkgs (or base R, or SpaDES.core), and a
## `pkg::` call needs `pkg` in reqdPkgs so that it is installed. `dplyr::summarise(N = n())` broke that way:
## dplyr was not listed, so `n()` was found only while another module had attached dplyr.

moduleFiles <- function() {
  root <- file.path(modulePath, moduleName)
  c(file.path(root, paste0(moduleName, ".R")),
    list.files(file.path(root, "R"), pattern = "\\.[Rr]$", full.names = TRUE))
}

reqdPkgNames <- function() {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  unique(Require::extractPkgName(unlist(md$reqdPkgs)))
}

## base R's packages are always available; a module need not list them
basePkgs <- setdiff(rownames(utils::installed.packages(priority = "base")), "tcltk")

test_that("every pkg:: in the module names a package in reqdPkgs", {
  pd <- do.call(rbind, lapply(moduleFiles(), function(f) utils::getParseData(parse(f, keep.source = TRUE))))
  used <- unique(pd$text[pd$token == "SYMBOL_PACKAGE"])
  missing <- setdiff(used, c(reqdPkgNames(), basePkgs))
  expect_identical(missing, character(0), info = paste(missing, collapse = ", "))
})

test_that("every function the module calls without pkg:: comes from its own reqdPkgs, base R or the module", {
  skip_if_not_installed("codetools")
  ## only the function definitions: the module file also holds defineModule(sim, ...), which needs a sim
  env <- new.env()
  for (f in moduleFiles()) for (e in as.list(parse(f, keep.source = FALSE)))
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function")))
      assign(as.character(e[[2]]), eval(e[[3]], env), envir = env)
  funs <- Filter(is.function, mget(ls(env, all.names = TRUE), envir = env))
  pkgs <- reqdPkgNames()
  pkgs <- pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  ## attaching a package also attaches its Depends; mirror SpaDES.core's imports
  deps <- unlist(lapply(pkgs, function(p) {
    d <- utils::packageDescription(p)$Depends
    if (is.null(d)) character() else trimws(sub("\\(.*", "", strsplit(d, ",")[[1]]))
  }))
  pkgs <- setdiff(unique(c(pkgs, deps, "SpaDES.core")), "R")
  pkgs <- pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  known <- unique(c(ls(env, all.names = TRUE),
                    unlist(lapply(pkgs, getNamespaceExports)),
                    unlist(lapply(basePkgs, function(p) tryCatch(getNamespaceExports(p), error = function(e) NULL))),
                    ls(baseenv(), all.names = TRUE)))
  called <- unique(unlist(lapply(funs, function(f) codetools::findGlobals(f, merge = FALSE)$functions)))
  ## data.table reads `.()` (= list()) and `J()` itself inside DT[...]; they are not exported functions
  dtSpecials <- c(".", "J")
  unresolved <- sort(setdiff(called, c(known, dtSpecials)))
  expect_identical(unresolved, character(0), info = paste(unresolved, collapse = ", "))
})
