source(file.path(RELEASE_ROOT,"src/a3/common.R"),encoding="UTF-8")
record_environment()
fields <- c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH3","BMXBMI","LBXGH","LBXSCR","egfr","MCQ160E","MCQ160F","MCQ160B","DIQ010","DIQ160","KIQ022","WTPH2YR","SDMVSTRA","SDMVPSU","base_eligible","primary_observable","dx_mi","dx_stroke","dx_chf")
check("required_fields",all(fields %in% names(dat)));assert_checks()
out(data.frame(variable=names(dat),class=vapply(dat,function(x)paste(class(x),collapse="/"),character(1)),typeof=vapply(dat,typeof,character(1)),label=vapply(dat,function(x){z<-attr(x,"label");if(is.null(z)) "" else z},character(1))),"input_schema.csv")
check("unique_primary_IDs",length(unique(dat$SEQN[primary])),PRIMARY_N)
check("no_duplicate_SEQN",!anyDuplicated(dat$SEQN));check("no_missing_SEQN",!anyNA(dat$SEQN))
check("primary_binary_partition",sum(dat$primary_observable[primary] %in% c(0,1)),PRIMARY_N)
b1 <- readr::read_csv(file.path(B1,"diagnostics/analytic_population_membership.csv"),show_col_types=FALSE)
b1 <- b1[b1$lane=="rebuilt",]; ids1<-b1$SEQN[b1$primary]
b2 <- readr::read_csv(file.path(OUTPUT_DIR,"b2/diagnostics/population_membership_and_reasons.csv"),show_col_types=FALSE)
ids2<-b2$SEQN[b2$primary]
check("exact_membership_B1",identical(sort(dat$SEQN[primary]),sort(ids1)))
check("exact_membership_B2",identical(sort(dat$SEQN[primary]),sort(ids2)))
ix<-match(dat$SEQN[primary],b2$SEQN)
check("clinical_classification_B2",all(dat$primary_observable[primary]==b2$primary_status[ix]))
check("primary_valid_positive_weights_design",all(parentmask[primary]))
par<-make_parent();dom<-subset(par,.primary)
check("age_public_bounds",all(dat$RIDAGEYR[primary]>=20 & dat$RIDAGEYR[primary]<=80))
check("DIQ010_allowed_original_cohort",all(dat$DIQ010[primary] %in% c(2,3)))
for (item in list(c("MCQ160E","dx_mi"),c("MCQ160F","dx_stroke"),c("MCQ160B","dx_chf"))) {
 raw<-dat[[item[1]]];x<-rep(NA_integer_,length(raw));x[raw %in% 1]<-1L;x[raw %in% 2]<-0L
 y<-dat[[item[2]]];check(paste0("archived_yes_no_",item[1]),identical(is.na(x),is.na(y))&&all(x[!is.na(x)]==y[!is.na(y)]))
}
out(do.call(rbind,checks),"preflight_QA.csv");assert_checks()
cohort<-dat[primary,]
out(cohort[c("SEQN","primary_observable","WTPH2YR","SDMVSTRA","SDMVPSU")],"verified_primary_membership.csv")
psu<-unique(par$variables[c("SDMVSTRA","SDMVPSU")]);psu$primary_n<-vapply(seq_len(nrow(psu)),function(i)sum(cohort$SDMVSTRA==psu$SDMVSTRA[i]&cohort$SDMVPSU==psu$SDMVPSU[i]),integer(1))
out(psu,"parent_PSU_occupancy.csv")
nums<-c("RIDAGEYR","BMXBMI","LBXGH","LBXSCR","egfr")
out(do.call(rbind,lapply(nums,function(v){x<-cohort[[v]];data.frame(variable=v,observed_n=sum(!is.na(x)),missing_n=sum(is.na(x)),min=min(x,na.rm=TRUE),max=max(x,na.rm=TRUE),nonfinite_n=sum(!is.na(x)&!is.finite(x)))})),"raw_ranges_and_availability.csv")
out(do.call(rbind,lapply(c("RIAGENDR","RIDRETH3","MCQ160E","MCQ160F","MCQ160B","DIQ010","DIQ160","KIQ022"),function(v){t<-table(cohort[[v]],useNA="always");data.frame(variable=v,code=names(t),n=as.integer(t))})),"raw_coding_counts.csv")
cat("PREFLIGHT PASS; observed primary n=",sum(primary),"; parent n=",sum(parentmask),". No weighted characteristics calculated.\n")
