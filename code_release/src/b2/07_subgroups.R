source(file.path(RELEASE_ROOT,"src/b2/common.R"),encoding="UTF-8")
groups<-list(Age=case_when(dat$RIDAGEYR>=20&dat$RIDAGEYR<=64~"20–64 y",dat$RIDAGEYR>=65&dat$RIDAGEYR<=74~"65–74 y",dat$RIDAGEYR>=75~"≥75 y",TRUE~NA_character_),Sex=dat$sex,`Race and ethnicity`=dat$race_ethnicity,Income=dat$income_group,Coverage=dat$coverage_group)
refs<-ref_table("primary/Figure1D_subgroups.csv");primary_design<-make_design(dat[primary,],"WTPH2YR")
rows<-list();occupancy<-list();policies<-list()
for(i in seq_len(nrow(refs))) {
 ref<-refs[i,];mask<-primary&groups[[ref$section]] %in% ref$characteristic
 domain<-estimate_domain(parents$phleb,mask,dat$primary_observable,ref$characteristic,"WTPH2YR")
 ratio<-ratio_estimate(parents$phleb,mask,dat$primary_observable,ref$characteristic,"WTPH2YR")
 idx<-primary_design$variables$.rowid;design0<-primary_design;design0$variables$.num<-ifelse(mask[idx],dat$primary_observable[idx],0);design0$variables$.den<-as.numeric(mask[idx]);z0<-svyratio(~.num,~.den,design0);se0<-100*as.numeric(SE(z0)[1])
 a<-ratio;a$estimate<-ref$estimate;a$lower<-ref$lower;a$upper<-ref$upper;a$se<-se0;a$ci_method<-ref$ci_method;a$implementation<-"A0_archived_primary_ratio_parentdf15";a$ci_df<-degf(primary_design)
 ratio$implementation<-"B0_full_parent_ratio_parentdf";domain$implementation<-"B0_full_parent_domain_subgroupdf"
 for(kind in c("A0","ratio","domain")) {
  row<-switch(kind,A0=a,ratio=ratio,domain=domain);row$section<-ref$section;row$characteristic<-ref$characteristic
  row$estimate_difference_vs_A0<-row$estimate-ref$estimate;row$se_difference_vs_A0<-row$se-se0;row$lower_difference_vs_A0<-row$lower-ref$lower;row$upper_difference_vs_A0<-row$upper-ref$upper;row$CI_width_pp<-row$upper-row$lower;row$RSE_percent<-100*row$se/abs(row$estimate)
  flags<-c(if(row$n<30) "n_below30",if(row$domain_df<8) "represented_df_below8",if(is.finite(row$kish_n)&&row$kish_n<30) "Kish_n_below30",if(is.finite(row$RSE_percent)&&row$RSE_percent>30) "RSE_above30percent",if(is.finite(row$CI_width_pp)&&row$CI_width_pp>20) "CI_width_above20pp",if(row$single_psu_domain_strata>0) "domain_single_PSU_strata",if(row$status!="ESTIMATED") "CI_not_estimable")
  row$reliability_flags<-paste(flags,collapse=";");row$flag_basis<-"Prespecified verification transparency screens, not official NCHS suppression rules";rows[[length(rows)+1]]<-row
 }
 cov<-coverage(parents$phleb,mask);occ<-cov$occupancy;occ$section<-ref$section;occ$characteristic<-ref$characteristic;occupancy[[i]]<-occ
 for(policy in c("adjust_false","adjust_true","fail_false")) {
  options(survey.lonely.psu=if(policy=="fail_false") "fail" else "adjust",survey.adjust.domain.lonely=policy=="adjust_true")
  row<-estimate_domain(parents$phleb,mask,dat$primary_observable,ref$characteristic,"WTPH2YR");row$policy<-policy;row$section<-ref$section;policies[[length(policies)+1]]<-row
 }
 options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=FALSE)
}
# Static API diagnosis only. No repaired subgroup hypothesis test is evaluated.
formula<-~factor(primary_observable)+age_group
rows_expr<-formula[[2]][[2]];selected_names<-as.character(rows_expr)
source<-deparse(getFromNamespace("svychisq.survey.design","survey"))
stopifnot(any(grepl("as.character(rows)",source,fixed=TRUE)))
failure<-tibble(archived_formula=paste(deparse(formula),collapse=" "),selected_names=paste(selected_names,collapse=";"),missing_selected_names=paste(setdiff(selected_names,names(dat)),collapse=";"),cause="survey4.4-2 converts formula LHS term call to column names; factor(primary_observable) becomes factor and primary_observable, so factor is absent",evidence="Installed svychisq.survey.design source rowvar <- unique(design$variables[, as.character(rows)]); Stage02B1 errors preserved",diagnostic_action="Corrected test not executed; original P values remain NA; original manuscript had no prespecified subgroup hypothesis testing")
save_csv(bind_rows(rows),"subgroup_reliability.csv")
save_csv(bind_rows(occupancy),"subgroup_parent_PSU_occupancy.csv","diagnostics")
save_csv(bind_rows(policies),"subgroup_lonely_policy_diagnostics.csv","diagnostics")
save_csv(failure,"subgroup_formula_failure_diagnosis.csv","diagnostics")
print(bind_rows(rows) |> filter(implementation=="B0_full_parent_domain_subgroupdf") |> select(section,characteristic,n,domain_df,estimate,lower,upper,lower_difference_vs_A0,upper_difference_vs_A0,reliability_flags));print(failure)
