source(file.path(RELEASE_ROOT,"src/b2/common.R"),encoding="UTF-8")
a_pred<-ifelse(is.na(dat$LBXGH),NA_integer_,as.integer(dat$LBXGH>=5.7&dat$LBXGH<6.5));a_high<-ifelse(is.na(dat$LBXGH),NA_integer_,as.integer(dat$LBXGH>=6.5))
f_pred<-ifelse(is.na(dat$LBXGLU),NA_integer_,as.integer(dat$LBXGLU>=100&dat$LBXGLU<=125));f_high<-ifelse(is.na(dat$LBXGLU),NA_integer_,as.integer(dat$LBXGLU>=126));f_cont<-ifelse(is.na(dat$LBXGLU),NA_integer_,as.integer(dat$LBXGLU>=100&dat$LBXGLU<126))
checks<-tibble(flag=c("a1c_prediabetes","a1c_diabetes_range","fpg_prediabetes","fpg_diabetes_range","no_reported_diabetes"),matches=c(equal_missing_values(a_pred,dat$a1c_prediabetes),equal_missing_values(a_high,dat$a1c_diabetes_range),equal_missing_values(f_pred,dat$fpg_prediabetes),equal_missing_values(f_high,dat$fpg_diabetes_range),identical(as.integer(dat$DIQ010 %in% c(2,3)),dat$no_reported_diabetes)))
save_csv(checks,"glycemic_flag_verification.csv","diagnostics");stopifnot(all(checks$matches))
raw<-read_xpt(file.path(RAW_DIR,"GLU_L.xpt"));stopifnot(!anyNA(raw$SEQN),!anyDuplicated(raw$SEQN))
precision<-list()
for(domain in list(list("raw_GLU_all",rep(TRUE,nrow(raw))),list("raw_GLU_positive_weight",mask_true(raw$WTSAF2YR>0)),list("raw_GLU_fast_candidate",raw$SEQN %in% dat$SEQN[fast_candidate]),list("raw_GLU_complete_fasting",raw$SEQN %in% dat$SEQN[fast]))) {
 x<-raw$LBXGLU[domain[[2]]];precision[[domain[[1]]]]<-tibble(domain=domain[[1]],n=length(x),nonmissing_n=sum(!is.na(x)),noninteger_values_n=sum(!is.na(x)&abs(x-round(x))>1e-10),fractional_125_to_126_n=sum(!is.na(x)&x>125&x<126),rule_reclassification_n=sum(!is.na(x)&x>125&x<126))
}
a_cat<-case_when(is.na(dat$LBXGH)~"Missing",dat$LBXGH<5.7~"Below5.7",dat$LBXGH<6.5~"Prediabetes range",TRUE~"Diabetes range")
f_cat<-case_when(is.na(dat$LBXGLU)~"Missing",dat$LBXGLU<100~"Below100",dat$LBXGLU<=125~"Archived prediabetes range",dat$LBXGLU<126~"Fractional125to126",TRUE~"Diabetes range")
cross<-list()
for(domain in list(list("base_phlebotomy",base),list("primary",primary),list("fasting_complete",fast))) {x<-tibble(a1c_category=a_cat[domain[[2]]],fpg_category=f_cat[domain[[2]]],reported_diabetes=dat$DIQ010[domain[[2]]]) |> count(a1c_category,fpg_category,reported_diabetes,name="n");x$domain<-domain[[1]];cross[[domain[[1]]]]<-x}
G1<-primary&(is.na(dat$LBXGH)|dat$LBXGH<6.5);G3obs<-primary&!is.na(dat$LBXGH);G1obs<-G3obs&dat$LBXGH<6.5;G2<-fast&dat$LBXGH<6.5&dat$LBXGLU<126
out<-list()
add_ph<-function(mask,scenario,scope,description) {
 row<-estimate_domain(parents$phleb,mask,dat$primary_observable,"primary_observable","WTPH2YR");row$scenario<-scenario;row$scope<-scope;row$rule<-description;row$baseline_n<-sum(primary);row$excluded_from_baseline_n<-sum(primary&!mask);row$retained_missing_HbA1c_n<-sum(mask&is.na(dat$LBXGH));row$outcome_reclassified_n<-0L;out[[length(out)+1]] <<- row
}
add_ph(primary,"G0","scenario_specific","Self-reported diabetes exclusion only")
add_ph(G1,"G1","scenario_specific","Additionally exclude HbA1c>=6.5, retain missing HbA1c")
add_ph(G3obs,"G3_observed_a1c","observed_laboratory_population","Observed HbA1c, original outcome and self-report exclusion")
add_ph(G1obs,"G1_observed_a1c","observed_laboratory_population","Observed HbA1c<6.5")
add_ph(G1,"G0","G0_G1_common_population","G0 outcome on G0/G1 intersection")
add_ph(G1,"G1","G0_G1_common_population","G1 outcome on G0/G1 intersection")
add_ph(G1obs,"G3_observed_a1c","observed_common_population","Observed original outcome on low-HbA1c intersection")
add_ph(G1obs,"G1_observed_a1c","observed_common_population","Same original outcome, observed HbA1c<6.5")
ff_cont<-make_fasting_flags(dat,f_cont);fast_cont<-fast_candidate&!is.na(dat$LBXGH)&!is.na(dat$LBXGLU)&!is.na(ff_cont$a1c)&!is.na(ff_cont$augmented);fast_rule_common<-fast&fast_cont
for(item in list(list("G0_fasting",fast,"scenario_specific"),list("G2",G2,"scenario_specific"),list("G3_observed_both_labs",fast,"observed_laboratory_population"),list("G0_fasting",G2,"G0_G2_common_population"),list("G2",G2,"G0_G2_common_population"),list("FPG_continuous_rule",fast_cont,"continuous_threshold_sensitivity"))) {
 scenario<-item[[1]];mask<-item[[2]];scope<-item[[3]]
 values<-if(scenario=="FPG_continuous_rule") list(HbA1c_only=ff_cont$a1c,HbA1c_or_FPG=ff_cont$augmented,gained_by_FPG=as.integer(ff_cont$a1c==0&ff_cont$augmented==1)) else list(HbA1c_only=dat$.fast_a1c,HbA1c_or_FPG=dat$.fast_aug,gained_by_FPG=dat$.fast_gain)
 for(v in names(values)) {row<-estimate_domain(parents$fasting,mask,values[[v]],v,"WTSAF2YR");row$scenario<-scenario;row$scope<-scope;row$baseline_n<-sum(fast);row$excluded_from_baseline_n<-sum(fast&!mask);row$outcome_reclassified_n<-if(scenario=="FPG_continuous_rule") sum(fast_rule_common&values[[v]]!=switch(v,HbA1c_only=dat$.fast_a1c,HbA1c_or_FPG=dat$.fast_aug,gained_by_FPG=dat$.fast_gain)) else 0L;row$rule<-if(scenario=="G2") "HbA1c<6.5 AND FPG<126, both observed, original outcomes" else if(scenario=="FPG_continuous_rule") "100<=FPG<126 instead of<=125" else "Original archived glycemic rules";out[[length(out)+1]]<-row}
}
transition<-tibble(SEQN=dat$SEQN,primary=primary,fasting=fast,G1=G1,G2=G2,HbA1c_high=a_high,FPG_high=f_high,HbA1c_missing=is.na(dat$LBXGH),FPG_missing=is.na(dat$LBXGLU),FPG=dat$LBXGLU,FPG_fractional_gap=!is.na(dat$LBXGLU)&dat$LBXGLU>125&dat$LBXGLU<126,original_fpg_prediabetes=f_pred,continuous_fpg_prediabetes=f_cont,original_augmented=dat$.fast_aug,continuous_augmented=ff_cont$augmented)
# Historical reference comparisons are isolated in tests, after calculations.
save_csv(bind_rows(out),"glycemic_sensitivity.csv")
save_csv(bind_rows(precision),"FPG_value_precision_and_threshold_counts.csv","diagnostics")
save_csv(bind_rows(cross),"laboratory_range_crossclassification.csv","diagnostics")
save_csv(transition |> filter(primary|fast|FPG_fractional_gap),"glycemic_record_transitions.csv","diagnostics")
save_csv(raw |> filter(!is.na(LBXGLU),LBXGLU>125,LBXGLU<126) |> select(SEQN,LBXGLU,WTSAF2YR),"raw_fractional_FPG_gap_records.csv","diagnostics")
print(bind_rows(precision));print(bind_rows(out) |> filter(scope=="scenario_specific") |> select(scenario,estimand,n,estimate,lower,upper,excluded_from_baseline_n))
