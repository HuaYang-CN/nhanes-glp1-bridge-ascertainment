# Portable process-local runtime. Scientific definitions live in the modules.
cfg <- read.csv(CONFIG_FILE, stringsAsFactors=FALSE, check.names=FALSE, fileEncoding="UTF-8")
val <- function(k,optional=FALSE) {
 x<-cfg$value[cfg$key==k]
 if(length(x)!=1L || (!optional && !nzchar(x))) stop("Missing/empty configuration: ",k)
 x
}
PROJECT_ROOT<-normalizePath(val("PROJECT_ROOT"),winslash="/",mustWork=TRUE)
RAW_DIR<-normalizePath(val("RAW_DATA_DIR"),winslash="/",mustWork=TRUE)
OUTPUT_DIR<-normalizePath(val("OUTPUT_DIR"),winslash="/",mustWork=TRUE)
lib<-val("R_LIBRARY_DIR",TRUE)
if(nzchar(lib)) .libPaths(c(normalizePath(lib,winslash="/",mustWork=TRUE),.libPaths()))
if(as.character(getRversion())!="4.5.1") stop("This candidate is validated with R 4.5.1 only")
expected_versions<-c(survey="4.4-2",haven="2.5.5",dplyr="1.1.4",readr="2.1.5",tibble="3.3.0")
for(pkg in names(expected_versions)) {
 if(!requireNamespace(pkg,quietly=TRUE)) stop("Install outside this pipeline is required for missing dependency: ",pkg)
 actual<-packageDescription(pkg)$Version
 if(actual!=expected_versions[[pkg]]) stop("Unsupported package version ",pkg,": ",actual,"; expected ",expected_versions[[pkg]])
}
suppressPackageStartupMessages({library(haven);library(dplyr);library(readr);library(survey);library(tibble)})
options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=FALSE,digits=17)
locale<-Sys.setlocale("LC_CTYPE","English_United States.utf8")
if(!identical(locale,"English_United States.utf8")) stop("Verified Windows UTF-8 locale unavailable")
inside<-function(p,root) identical(tolower(p),tolower(root)) || startsWith(tolower(p),paste0(tolower(root),"/"))
safe_output<-function(p) {
 if(!is.character(p)||length(p)!=1L||!nzchar(p)) stop("Non-path writer destination")
 if(file.exists(p)) stop("Refuse existing output: ",p)
 ancestor<-dirname(p)
 while(!dir.exists(ancestor)) ancestor<-dirname(ancestor)
 resolved<-normalizePath(ancestor,winslash="/",mustWork=TRUE)
 if(!inside(resolved,OUTPUT_DIR)||inside(resolved,RAW_DIR)||inside(resolved,RELEASE_ROOT)) stop("Writer escaped OUTPUT_DIR: ",p)
 if(any(strsplit(gsub("\\\\","/",p),"/",fixed=TRUE)[[1]] %in% "..")) stop("Parent traversal in writer destination")
 dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
 p
}
write_csv<-function(x,file,...) readr::write_csv(x,safe_output(file),...)
saveRDS<-function(object,file,...) base::saveRDS(object,safe_output(file),...)
writeLines<-function(text,con,...) base::writeLines(text,safe_output(con),...)
assert_stage_path<-safe_output
STAGE_ROOT<-file.path(OUTPUT_DIR,"b1")
CODE_DIR<-file.path(RELEASE_ROOT,"src/b1")
DERIVED_DIR<-file.path(STAGE_ROOT,"diagnostics/rebuilt")
RESULT_ROOT<-file.path(STAGE_ROOT,"results/trackA/rebuilt")
OUTPUT_ROOT<-RESULT_ROOT
FIG_DIR<-file.path(RESULT_ROOT,"figures")
LOG_DIR<-file.path(OUTPUT_DIR,"logs")
INPUT_SHA256<-val("INPUT_SHA256",TRUE)
for(d in c(DERIVED_DIR,RESULT_ROOT,LOG_DIR,file.path(OUTPUT_DIR,"b2"),file.path(OUTPUT_DIR,"b2/results"),file.path(OUTPUT_DIR,"b2/diagnostics"),file.path(OUTPUT_DIR,"a3/diagnostics"),file.path(OUTPUT_DIR,"a3/logs"),file.path(OUTPUT_DIR,"a3/provenance"))) dir.create(d,recursive=TRUE,showWarnings=FALSE)
print_environment<-function(){print(sessionInfo());print(.libPaths());print(options()[c("survey.lonely.psu","survey.adjust.domain.lonely")]);print(expected_versions)}
