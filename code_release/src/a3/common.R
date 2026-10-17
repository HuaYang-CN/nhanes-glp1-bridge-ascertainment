S <- file.path(OUTPUT_DIR,"a3")
P <- PROJECT_ROOT; B1 <- file.path(OUTPUT_DIR,"b1")
for (pkg in c("survey","readr","haven")) if (!requireNamespace(pkg,quietly=TRUE)) stop("Required existing package unavailable: ",pkg)
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=FALSE,digits=17)
out <- function(x, name, folder="diagnostics") {
 p <- file.path(S,folder,name)
 if(file.exists(p)) stop("Refuse output overwrite: ",p)
 write_csv(x,p,na="NA"); invisible(x)
}
dat <- readRDS(file.path(B1,"diagnostics/rebuilt/bridge_analysis_ready.rds"))
primary <- dat$base_eligible %in% TRUE & !is.na(dat$primary_observable)
parentmask <- is.finite(dat$WTPH2YR) & dat$WTPH2YR>0 & is.finite(dat$SDMVSTRA) & is.finite(dat$SDMVPSU)
PRIMARY_N <- sum(primary)
dat$.primary <- primary
make_parent <- function() svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTPH2YR,nest=TRUE,data=dat[parentmask,])
checks <- list()
check <- function(id, value, expected=TRUE, tolerance=0, reason="") {
 ok <- if(is.numeric(value)&&is.numeric(expected)) is.finite(value)&&abs(value-expected)<=tolerance else identical(value,expected)
 checks[[length(checks)+1]] <<- data.frame(check_id=id,actual=as.character(value),expected=as.character(expected),tolerance=tolerance,status=if(ok) "PASS" else "FAIL",reason=reason)
 invisible(ok)
}
assert_checks <- function() if(any(vapply(checks,function(x)x$status=="FAIL",logical(1)))) stop("HOLD_DATA_OR_METHODS: validation failed")
record_environment <- function() {
 writeLines(c(capture.output(sessionInfo()),paste(".libPaths:",paste(.libPaths(),collapse=" | ")),paste("lonely:",getOption("survey.lonely.psu")),paste("adjust domain lonely:",getOption("survey.adjust.domain.lonely"))),file.path(S,"logs/sessionInfo.txt"))
 out(data.frame(package=c("survey","readr","haven"),version=vapply(c("survey","readr","haven"),function(p)as.character(packageVersion(p)),character(1))),"package_versions.csv","provenance")
}
