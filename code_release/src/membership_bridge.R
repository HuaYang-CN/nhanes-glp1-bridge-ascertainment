# Independently reuse the verified baseline assignment expressions, not old RDS/ledgers.
extract_assignments<-function(script,allowed,envir) {
 expr<-parse(file.path(RELEASE_ROOT,"src/b1",script),encoding="UTF-8");selected<-character()
 for(e in expr)if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]])) {
  key<-as.character(e[[2]]);if(key %in% allowed){eval(e,envir);selected<-c(selected,key)}
 }
 stopifnot(setequal(selected,allowed))
}
dat<-readRDS(file.path(DERIVED_DIR,"bridge_analysis_ready.rds"));e<-new.env(parent=globalenv());e$dat<-dat
extract_assignments("02_reproduce_primary_results_v4.R",c("make_design","primary_d","expanded_d","a1c_exclusion_d","primary_design","pathways","id_d","id_design","subgroup_d","subgroup_design","subgroup_levels","subgroup_var"),e)
ef<-new.env(parent=globalenv());ef$dat<-dat
extract_assignments("03_sensitivity_fasting_glucose.R",c("bmi_plus_binary","or_status","make_design","dat2","fast_d","des","des_aug_pos"),ef)
comparison<-e$id_d;positive<-comparison[comparison$primary_observable==1,]
membership<-data.frame(lane="rebuilt",SEQN=dat$SEQN,primary=dat$SEQN %in% e$primary_d$SEQN,expanded=dat$SEQN %in% e$expanded_d$SEQN,a1c_exclusion=dat$SEQN %in% e$a1c_exclusion_d$SEQN,comparison=dat$SEQN %in% comparison$SEQN,measurement_positive=dat$SEQN %in% positive$SEQN,fasting=dat$SEQN %in% ef$fast_d$SEQN)
stopifnot(nrow(positive)==sum(comparison$measurement_only==1)+sum(comparison$both_identified==1))
neither<-comparison[comparison$primary_observable==0&comparison$diagnosis_proxy==0,]
stopifnot(nrow(comparison)==sum(comparison$measurement_only==1)+sum(comparison$diagnosis_only==1)+sum(comparison$both_identified==1)+nrow(neither))
counts<-lapply(names(e$pathways),function(v){d<-e$primary_d;x<-d[[v]];data.frame(lane="rebuilt",domain=paste0("Pathway: ",v),outcome_nonmissing_n=sum(!is.na(x)),outcome_missing_n=sum(is.na(x)),positive_n=sum(x==1,na.rm=TRUE))})
write_csv(membership,file.path(STAGE_ROOT,"diagnostics/analytic_population_membership.csv"))
write_csv(bind_rows(counts),file.path(STAGE_ROOT,"diagnostics/denominators_and_numerators.csv"))
cat("Fresh independent baseline membership/pathway bridge PASS\n")
