source(file.path(RELEASE_ROOT,"src/b2/common.R"),encoding="UTF-8")
spec<-list(
 list("Primary observable criteria","primary_observable",primary,"WTPH2YR","primary/Figure1A_primary_and_sensitivity.csv","Primary observable criteria",primary),
 list("Expanded clinical proxies","expanded_proxy",expanded,"WTPH2YR","primary/Figure1A_primary_and_sensitivity.csv","Expanded clinical proxies",expanded),
 list("Exclude HbA1c ≥ 6.5%","primary_observable",a1c_exclusion,"WTPH2YR","primary/Figure1A_primary_and_sensitivity.csv","Exclude HbA1c ≥ 6.5%",a1c_exclusion))
for(item in list(c("Measurement-based criteria","primary_observable"),c("Self-reported diagnosis proxy","diagnosis_proxy"),c("Identified by both","both_identified"),c("Measurement only","measurement_only"),c("Diagnosis proxy only","diagnosis_only"))) spec[[length(spec)+1]]<-list(item[1],item[2],comparison,"WTPH2YR","primary/Figure1C_identification_strategy.csv",item[1],comparison)
spec[[length(spec)+1]]<-list("Conditional measurement-only","measurement_only",cond,"WTPH2YR","ascertainment_gap/KEY_measurement_only_among_measurement_positive.csv","Measurement-only among measurement-positive adults",comparison)
for(item in list(c("HbA1c-only fasting",".fast_a1c","HbA1c-only prediabetes definition"),c("HbA1c-or-FPG fasting",".fast_aug","HbA1c or fasting-glucose prediabetes definition"))) spec[[length(spec)+1]]<-list(item[1],item[2],fast,"WTSAF2YR","sensitivity_fpg/fasting_glucose_sensitivity_main.csv",item[3],fast)
spec[[length(spec)+1]]<-list("Conditional fasting gain",".fast_gain",fast&dat$.fast_aug %in% 1,"WTSAF2YR","sensitivity_fpg/fasting_glucose_conditional_gain.csv","Added by fasting glucose among HbA1c/FPG-positive adults",fast)
for(v in names(paths)) spec[[length(spec)+1]]<-list(unname(paths[v]),v,primary,"WTPH2YR","primary/Figure1B_criterion_pathways.csv",unname(paths[v]),primary)
output<-list();occupancy<-list();lonely<-list();identities<-list()
for(s in spec) {
 label<-s[[1]];var<-s[[2]];mask<-s[[3]];w<-s[[4]];parent<-if(w=="WTPH2YR") parents$phleb else parents$fasting
 source<-ref_table(s[[5]]) |> filter(characteristic==s[[6]]);stopifnot(nrow(source)==1)
 filtered<-make_design(dat[mask_true(s[[7]]),],w);filtered$variables$.outcome<-dat[[var]][filtered$variables$.rowid]
 filtered<-filtered[mask_true(mask[filtered$variables$.rowid]),]
 m0<-svymean(~.outcome,filtered,na.rm=TRUE);p0<-as.numeric(coef(m0)[1])*100;se0<-as.numeric(SE(m0)[1])*100
 b<-estimate_domain(parent,mask,dat[[var]],label,w);r<-ratio_estimate(parent,mask,dat[[var]],label,w)
 a<-b;a$implementation<-"A0_filtered_TrackA";a$estimate<-source$estimate;a$lower<-source$lower;a$upper<-source$upper;a$se<-se0;a$ci_method<-if("ci_method" %in% names(source)) source$ci_method else "logit";a$ci_df<-degf(filtered);a$df_policy<-"archived_filtered_design";a$source_file<-file.path(B1,"results/trackA/rebuilt",s[[5]])
 b$implementation<-"B0_full_parent_domain";b$source_file<-a$source_file;r$implementation<-"B0_parent_df_ratio";r$source_file<-a$source_file
 for(key in c("estimate","se","lower","upper","ci_df")) {b[[paste0(key,"_difference_vs_A0")]]<-b[[key]]-a[[key]];r[[paste0(key,"_difference_vs_A0")]]<-r[[key]]-a[[key]];a[[paste0(key,"_difference_vs_A0")]]<-0}
 b$se_relative_change_percent<-100*(b$se-a$se)/a$se;r$se_relative_change_percent<-100*(r$se-a$se)/a$se;a$se_relative_change_percent<-0
 output[[label]]<-bind_rows(a,b,r)
 known<-mask_true(mask)&!is.na(dat[[var]]);cov<-coverage(parent,known);occ<-cov$occupancy;occ$estimand<-label;occ$weight_variable<-w;occupancy[[label]]<-occ
 identities[[label]]<-tibble(estimand=label,filtered_mean_minus_TrackA_pp=p0-source$estimate,B0_minus_A0_pp=b$estimate-a$estimate,ratio_minus_domain_point_pp=r$estimate-b$estimate,ratio_minus_domain_SE_pp=r$se-b$se,point_agreement=abs(r$estimate-b$estimate)<1e-8,SE_agreement=abs(r$se-b$se)<1e-8)
 for(policy in c("adjust_false","adjust_true","fail_false")) {
  options(survey.lonely.psu=if(policy=="fail_false") "fail" else "adjust",survey.adjust.domain.lonely=policy=="adjust_true")
  row<-estimate_domain(parent,mask,dat[[var]],label,w);row$policy<-policy;lonely[[length(lonely)+1]]<-row
 }
 options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=FALSE)
}
save_csv(bind_rows(output),"design_comparison.csv")
save_csv(bind_rows(occupancy),"design_domain_psu_occupancy.csv","diagnostics")
save_csv(bind_rows(lonely),"lonely_policy_diagnostics.csv","diagnostics")
save_csv(bind_rows(identities),"design_point_and_SE_identity_checks.csv","diagnostics")
print(bind_rows(output) |> filter(implementation=="B0_full_parent_domain") |> select(estimand,n,estimate,se,lower,upper,estimate_difference_vs_A0,se_difference_vs_A0,lower_difference_vs_A0,upper_difference_vs_A0))
