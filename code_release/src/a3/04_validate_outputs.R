source(file.path(RELEASE_ROOT,"src/a3/common.R"),encoding="UTF-8")
t1<-readr::read_csv(file.path(S,"diagnostics/statistics_full_precision.csv"),show_col_types=FALSE,na="NA")
comp<-readr::read_csv(file.path(S,"diagnostics/completeness_full_precision.csv"),show_col_types=FALSE,na="NA")
par<-make_parent();pv<-par$variables;occ<-unique(pv[c("SDMVSTRA","SDMVPSU")]);weights_all<-pv$WTPH2YR
independent_values <- function(r) {
 raw<-pv[[r$source_field]];cat<-r$category
 if(r$statistic=="weighted_mean") return(as.numeric(raw))
 if(r$variable=="age_category") return(if(cat=="20-64") as.numeric(raw>=20&raw<65) else if(cat=="65-74") as.numeric(raw>=65&raw<75) else as.numeric(raw>=75))
 if(r$variable=="bmi_category") return(if(cat=="<27")as.numeric(raw<27) else if(cat=="27 to <30")as.numeric(raw>=27&raw<30) else if(cat=="30 to <35")as.numeric(raw>=30&raw<35) else as.numeric(raw>=35))
 if(r$variable=="hba1c_category")return(if(cat=="<5.7")as.numeric(raw<5.7) else if(cat=="5.7 to <6.5")as.numeric(raw>=5.7&raw<6.5) else as.numeric(raw>=6.5))
 if(r$variable=="egfr_low")return(as.numeric(raw<60))
 if(r$variable %in% c("MCQ160E","MCQ160F","MCQ160B"))return(ifelse(raw %in% c(1,2),as.numeric(raw==1),NA_real_))
 as.numeric(raw==as.numeric(r$category_code))
}
manual<-list()
for(i in seq_len(nrow(t1))) {
 r<-t1[i,];x<-independent_values(r);valid<-pv$.primary & !is.na(x);n<-sum(valid)
 check(paste0(r$result_id,"_observed"),n,as.integer(r$observed_n))
 check(paste0(r$result_id,"_missing"),PRIMARY_N-n,as.integer(r$missing_n))
 check(paste0(r$result_id,"_finite"),all(is.finite(x[valid])))
 if(r$statistic=="weighted_percentage")check(paste0(r$result_id,"_category_n"),sum(x[valid]==1),as.integer(r$category_n))
 w<-weights_all;den<-sum(w[valid]);mean_manual<-sum(w[valid]*x[valid])/den
 residual<-rep(0,nrow(pv));residual[valid]<-w[valid]*(x[valid]-mean_manual)/den
 psu_sum<-vapply(seq_len(nrow(occ)),function(j)sum(residual[pv$SDMVSTRA==occ$SDMVSTRA[j]&pv$SDMVPSU==occ$SDMVPSU[j]]),numeric(1))
 variance<-sum(vapply(split(psu_sum,occ$SDMVSTRA),function(z)length(z)/(length(z)-1)*sum((z-mean(z))^2),numeric(1)))
 se_manual<-sqrt(variance)
 ratio<-par;ratio$variables$.num<-ifelse(valid,x,0);ratio$variables$.den<-as.numeric(valid)
 z<-svyratio(~.num,~.den,ratio);mult<-if(r$statistic=="weighted_percentage")100 else 1
 check(paste0(r$result_id,"_manual_mean"),mult*mean_manual,r$estimate,1e-8)
 check(paste0(r$result_id,"_manual_Taylor_SE"),mult*se_manual,r$se,1e-8)
 check(paste0(r$result_id,"_parent_ratio_mean"),mult*as.numeric(coef(z)),r$estimate,1e-8)
 check(paste0(r$result_id,"_parent_ratio_SE"),mult*as.numeric(SE(z)),r$se,1e-8)
 check(paste0(r$result_id,"_weighted_denominator"),den,r$weighted_denominator,1e-6)
 cv<-pv[valid,];df<-nrow(unique(cv[c("SDMVSTRA","SDMVPSU")]))-length(unique(cv$SDMVSTRA))
 check(paste0(r$result_id,"_represented_df"),as.numeric(df),as.numeric(r$df))
 check(paste0(r$result_id,"_nonnegative_SE"),r$se>=0)
 if(r$status=="ESTIMATED") {
   ci<-if(r$statistic=="weighted_mean")mean_manual+c(-1,1)*qt(.975,df)*se_manual else plogis(qlogis(mean_manual)+c(-1,1)*qt(.975,df)*se_manual/(mean_manual*(1-mean_manual)))
   check(paste0(r$result_id,"_CI_lower"),mult*ci[1],r$lower_ci,1e-6)
   check(paste0(r$result_id,"_CI_upper"),mult*ci[2],r$upper_ci,1e-6)
   check(paste0(r$result_id,"_CI_order"),r$lower_ci<=r$estimate && r$estimate<=r$upper_ci)
   if(r$statistic=="weighted_percentage")check(paste0(r$result_id,"_proportion_bounds"),r$lower_ci>=0&&r$upper_ci<=100&&r$estimate>=0&&r$estimate<=100)
 } else check(paste0(r$result_id,"_nonestimable_reason"),nzchar(r$not_estimable_reason)&&is.na(r$lower_ci)&&is.na(r$upper_ci))
 manual[[i]]<-data.frame(result_id=r$result_id,manual_mean=mult*mean_manual,manual_Taylor_SE=mult*se_manual,parent_ratio_mean=mult*as.numeric(coef(z)),parent_ratio_SE=mult*as.numeric(SE(z)),represented_df=df)
}
for(v in c("age_category","RIAGENDR","RIDRETH3","bmi_category","hba1c_category")) {
 r<-t1[t1$variable==v,]
 check(paste0("partition_n_",v),sum(r$category_n),r$observed_n[1])
 check(paste0("partition_weighted_",v),sum(r$estimate),100,1e-8)
}
cohort<-dat[primary,]
for(i in seq_len(nrow(comp))) {
 r<-comp[i,];x<-cohort[[r$variable]]
 questionnaire<-r$variable %in% c("MCQ160B","MCQ160E","MCQ160F","DIQ010","DIQ160","KIQ022")
 observed<-if(questionnaire)sum(x %in% (if(r$variable=="DIQ010")c(1,2,3) else c(1,2)))else sum(!is.na(x))
 check(paste0("S1_",r$variable,"_observed"),observed,as.integer(r$observed_n))
 check(paste0("S1_",r$variable,"_missing"),PRIMARY_N-observed,as.integer(r$missing_n))
 check(paste0("S1_",r$variable,"_percent"),100*(PRIMARY_N-observed)/PRIMARY_N,r$unweighted_missing_percent,1e-12)
 check(paste0("S1_",r$variable,"_NA"),sum(is.na(x)),as.integer(r$raw_NA_n))
 if(questionnaire) {
   for(pair in list(c("1","yes_n"),c("2","no_n"),c("7","refused_n"),c("9","dont_know_n")))check(paste0("S1_",r$variable,"_",pair[2]),sum(x %in% as.numeric(pair[1])),as.integer(r[[pair[2]]]))
 }
}
check("eGFR_source_creatinine_missing_mask",identical(is.na(cohort$egfr),is.na(cohort$LBXSCR)))
check("no_high_lab_primary_exclusion",sum(cohort$LBXGH>=6.5,na.rm=TRUE)>0)
for(v in c("RIDAGEYR","BMXBMI","LBXGH","LBXSCR","egfr")) {
 x<-cohort[[v]];range<-switch(v,RIDAGEYR=c(20,80),BMXBMI=c(0,100),LBXGH=c(0,30),LBXSCR=c(0,30),egfr=c(0,250))
 check(paste0("plausible_range_",v),all(x[!is.na(x)]>=range[1]&x[!is.na(x)]<=range[2]),reason="Broad validation only; no outlier removed")
}
lock<-readr::read_csv(file.path(OUTPUT_DIR,"b2/results/denominator_audit.csv"),show_col_types=FALSE)
lock<-lock[lock$record_type %in% "criterion_route",]
lock$result_id<-sprintf("C02-DENOM-%03d",seq_len(nrow(lock))+17L)
route_rows<-list()
for(item in list(c("crit_bmi35","BMXBMI","C02-DENOM-018"),c("crit_egfr","egfr","C02-DENOM-019"),c("crit_a1c_prediabetes","LBXGH","C02-DENOM-020"),c("crit_mi","MCQ160E","C02-DENOM-021"),c("crit_stroke","MCQ160F","C02-DENOM-022"),c("crit_chf_proxy","MCQ160B","C02-DENOM-023"))) {
 z<-lock[lock$result_id==item[3],];check(paste0("lock_registry_identity_",item[3]),nrow(z),1L)
 actual<-sum(is.na(cohort[[item[1]]]));check(paste0("locked_route_missing_",item[1]),as.numeric(actual),as.numeric(z$missing_n))
 raw<-comp[comp$variable==item[2],]
 route_rows[[item[1]]]<-data.frame(route=item[1],source_raw_variable=item[2],raw_missing_n=raw$missing_n,route_missing_n=actual,locked_route_missing_n=z$missing_n,locked_result_id=item[3],explanation="Route depends on BMI threshold and three-state logic; raw measurement/history NA or unknown is not a route negative. Below BMI threshold the route can be known negative despite unavailable raw component.")
}
bmi35<-t1[t1$variable=="bmi_category"&t1$category=="35 or higher",];z<-lock[lock$result_id=="C02-DENOM-018",]
check("BMI35_matches_locked_route_count",as.numeric(bmi35$category_n),as.numeric(z$positive_n))
check("BMI35_matches_locked_weighted_percent",bmi35$estimate,as.numeric(z$estimate),1e-8)
check("Mexican_demographic_valid_domain_df",as.numeric(t1$df[t1$variable=="RIDRETH3"&t1$category_code=="1"]),15)
out(do.call(rbind,manual),"independent_estimator_checks.csv")
out(do.call(rbind,route_rows),"raw_vs_route_missingness.csv")
out(do.call(rbind,checks),"independent_numerical_QA.csv");assert_checks()
cat("Independent numerical checks:",length(checks),"PASS; manual PSU Taylor variance and full-parent ratios verified.\n")
