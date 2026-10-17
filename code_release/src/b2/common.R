S<-file.path(OUTPUT_DIR,"b2")
B1<-file.path(OUTPUT_DIR,"b1");P<-PROJECT_ROOT
suppressPackageStartupMessages({library(dplyr);library(readr);library(survey);library(tibble);library(haven)})
if(!exists("output_attempt",inherits=TRUE)) output_attempt<-""
assert_output<-function(p) {
 ancestor<-dirname(p)
 while(!dir.exists(ancestor)) ancestor<-dirname(ancestor)
 resolved<-normalizePath(ancestor,winslash="/",mustWork=TRUE)
 root<-normalizePath(S,winslash="/",mustWork=TRUE)
 stopifnot(identical(tolower(resolved),tolower(root))||startsWith(tolower(resolved),paste0(tolower(root),"/")))
 dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
 if(file.exists(p)) stop("Existing output cannot be overwritten: ",p)
 p
}
save_csv<-function(x,name,folder="results") {
 p<-assert_output(file.path(S,folder,output_attempt,name));write_csv(x,p);cat("OUTPUT:",p,"\n");invisible(x)
}
save_text<-function(x,name,folder="diagnostics") {
 p<-assert_output(file.path(S,folder,output_attempt,name));writeLines(x,p,useBytes=TRUE);cat("OUTPUT:",p,"\n")
}
tri_or<-function(...) {
 m<-cbind(...);out<-rep(NA_integer_,nrow(m));positive<-rowSums(m==1,na.rm=TRUE)>0
 out[positive]<-1L;known_negative<-rowSums(is.na(m))==0&rowSums(m==0,na.rm=TRUE)==ncol(m);out[known_negative]<-0L;out
}
bmi_binary<-function(bmi,cutoff,flag) {
 out<-rep(NA_integer_,length(bmi));out[!is.na(bmi)&bmi<cutoff]<-0L
 use<-!is.na(bmi)&bmi>=cutoff&!is.na(flag);out[use]<-as.integer(flag[use]==1);out
}
egfr_formula<-function(creatinine,age,sex) {
 k<-ifelse(sex==2,0.7,ifelse(sex==1,0.9,NA_real_));alpha<-ifelse(sex==2,-0.241,ifelse(sex==1,-0.302,NA_real_));female<-ifelse(sex==2,1.012,ifelse(sex==1,1,NA_real_));ratio<-creatinine/k
 ifelse(is.na(creatinine)|is.na(age)|is.na(k),NA_real_,142*pmin(ratio,1)^alpha*pmax(ratio,1)^(-1.200)*0.9938^age*female)
}
make_fasting_flags<-function(d,fpg= d$fpg_prediabetes) {
 combined<-ifelse(is.na(d$a1c_prediabetes)&is.na(fpg),NA_integer_,as.integer((!is.na(d$a1c_prediabetes)&d$a1c_prediabetes==1)|(!is.na(fpg)&fpg==1)))
 a_route<-bmi_binary(d$BMXBMI,27,d$a1c_prediabetes);aug_route<-bmi_binary(d$BMXBMI,27,combined)
 list(a1c=tri_or(d$crit_bmi35,d$crit_egfr,a_route,d$crit_mi,d$crit_stroke),augmented=tri_or(d$crit_bmi35,d$crit_egfr,aug_route,d$crit_mi,d$crit_stroke))
}
dat<-readRDS(file.path(B1,"diagnostics/rebuilt/bridge_analysis_ready.rds"));dat$.rowid<-seq_len(nrow(dat))
stopifnot(!anyNA(dat$SEQN),!anyDuplicated(dat$SEQN))
required<-c("WTPH2YR","WTSAF2YR","SDMVSTRA","SDMVPSU","base_eligible","primary_observable","diagnosis_proxy","expanded_proxy","LBXGH","LBXGLU","DIQ010","DIQ160","egfr")
stopifnot(all(required %in% names(dat)))
f<-make_fasting_flags(dat);dat$.fast_a1c<-f$a1c;dat$.fast_aug<-f$augmented;dat$.fast_gain<-as.integer(f$a1c==0&f$augmented==1)
mask_true<-function(x) x %in% TRUE
base<-mask_true(dat$base_eligible);primary<-base&!is.na(dat$primary_observable);expanded<-base&!is.na(dat$expanded_proxy)
comparison<-primary&!is.na(dat$diagnosis_proxy);cond<-comparison&dat$primary_observable %in% 1
a1c_exclusion<-primary&(is.na(dat$LBXGH)|dat$LBXGH<6.5)
fast_candidate<-mask_true(dat$RIDAGEYR>=20&dat$medicare==1&dat$no_reported_diabetes==1&!is.na(dat$BMXBMI)&dat$WTSAF2YR>0)
fast<-fast_candidate&!is.na(dat$LBXGH)&!is.na(dat$LBXGLU)&!is.na(dat$.fast_a1c)&!is.na(dat$.fast_aug)
make_design<-function(d,w) svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=as.formula(paste0("~",w)),nest=TRUE,data=d)
parents<-list()
for(item in list(c("phleb","WTPH2YR"),c("fasting","WTSAF2YR"))) {
 w<-dat[[item[2]]];positive<-!is.na(w)&is.finite(w)&w>0
 if(any(!is.na(w)&!is.finite(w))) stop("Nonfinite source weight")
 if(any(positive&(is.na(dat$SDMVSTRA)|is.na(dat$SDMVPSU)))) stop("Positive weight with invalid design key")
 parents[[item[1]]]<-make_design(dat[positive,],item[2])
}
cat("INPUT:",file.path(B1,"diagnostics/rebuilt/bridge_analysis_ready.rds"),"\nParent n:",nrow(parents$phleb$variables),nrow(parents$fasting$variables),"\n")
paths<-c(crit_bmi35="BMI ≥ 35",crit_egfr="BMI ≥ 30 + eGFR <60 mL/min/1.73 m²",crit_a1c_prediabetes="BMI ≥ 27 + prediabetes",crit_mi="BMI ≥ 27 + prior MI",crit_stroke="BMI ≥ 27 + prior stroke",crit_chf_proxy="BMI ≥ 30 + heart failure proxy",crit_bp_proxy="BMI ≥ 30 + uncontrolled BP proxy")
coverage<-function(parent,mask) {
 vars<-parent$variables;in_domain<-mask_true(mask[vars$.rowid]);v<-vars[in_domain,]
 allpsu<-unique(vars[c("SDMVSTRA","SDMVPSU")]);allpsu$n_domain<-vapply(seq_len(nrow(allpsu)),function(i) sum(in_domain&vars$SDMVSTRA==allpsu$SDMVSTRA[i]&vars$SDMVPSU==allpsu$SDMVPSU[i]),integer(1))
 counts<-table(v$SDMVSTRA,interaction(v$SDMVSTRA,v$SDMVPSU,drop=TRUE))
 strata<-length(unique(v$SDMVSTRA));psus<-nrow(unique(v[c("SDMVSTRA","SDMVPSU")]))
 dompsu<-allpsu[allpsu$n_domain>0,];bystrata<-table(dompsu$SDMVSTRA)
 list(n=nrow(v),strata=strata,psus=psus,represented_df=psus-strata,parent_df=degf(parent),zero_domain_psus=sum(allpsu$n_domain==0),zero_domain_strata=sum(!unique(vars$SDMVSTRA) %in% v$SDMVSTRA),single_psu_domain_strata=sum(bystrata==1),occupancy=allpsu)
}
estimate_domain<-function(parent,mask,outcome,label,weight,policy="domain_df",ci="logit") {
 mask<-mask_true(mask);known<-mask&!is.na(outcome);cov<-coverage(parent,known)
 idx<-parent$variables$.rowid;des<-parent;des$variables$.outcome<-outcome[idx];des<-des[known[idx],]
 df<-if(policy=="parent_df") degf(parent) else cov$represented_df
 row<-tibble(estimand=label,weight_variable=weight,universe_n=sum(mask),n=cov$n,missing_outcome_n=sum(mask&is.na(outcome)),positive_n=sum(known&outcome==1),contributing_strata=cov$strata,contributing_psus=cov$psus,domain_df=cov$represented_df,parent_df=cov$parent_df,ci_df=df,zero_domain_psus=cov$zero_domain_psus,zero_domain_strata=cov$zero_domain_strata,single_psu_domain_strata=cov$single_psu_domain_strata,estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_,ci_method="not_estimable",df_policy=policy,lonely_psu=getOption("survey.lonely.psu"),adjust_domain_lonely=getOption("survey.adjust.domain.lonely"),status="NOT_ESTIMABLE",error="",kish_n=NA_real_,weighted_denominator=NA_real_)
 if(cov$n==0) {row$error<-"Empty domain";return(row)}
 ans<-tryCatch({
  m<-svymean(~.outcome,des,na.rm=TRUE);est<-as.numeric(coef(m)[1]);se<-as.numeric(SE(m)[1]);row$estimate<-100*est;row$se<-100*se
  ww<-as.numeric(weights(des));row$kish_n<-sum(ww)^2/sum(ww^2);row$weighted_denominator<-sum(ww)
  if(df<=0) {row$error<-"Point/SE available, CI not estimable: nonpositive represented-domain df";return(row)}
  success<-FALSE
  if(ci=="ratio_logit") {
   if(is.finite(est)&&is.finite(se)&&est>0&&est<1) {q<-qt(.975,df);lo<-plogis(qlogis(est)-q*se/(est*(1-est)));hi<-plogis(qlogis(est)+q*se/(est*(1-est)));row$lower<-100*lo;row$upper<-100*hi;row$ci_method<-"domain_ratio_logit";success<-TRUE}
  } else {
   for(method in c("logit","beta")) {
    z<-try(svyciprop(~.outcome,des,method=method,df=df,na.rm=TRUE),silent=TRUE)
    if(!inherits(z,"try-error")) {limits<-try(confint(z),silent=TRUE);if(!inherits(limits,"try-error")&&all(is.finite(limits))) {row$lower<-100*as.numeric(limits[1]);row$upper<-100*as.numeric(limits[2]);row$ci_method<-method;success<-TRUE;break}}
   }
  }
  if(!success && is.finite(est)&&is.finite(se)) {q<-qt(.975,df);row$lower<-100*max(0,est-q*se);row$upper<-100*min(1,est+q*se);row$ci_method<-"wald_fallback";success<-TRUE}
  row$status<-if(success&&all(is.finite(c(row$estimate,row$se,row$lower,row$upper)))) "ESTIMATED" else "NOT_ESTIMABLE";row
 },error=function(e) {row$error<-conditionMessage(e);row})
 ans
}
ratio_estimate<-function(parent,mask,outcome,label,weight) {
 row<-estimate_domain(parent,mask,outcome,label,weight,policy="parent_df",ci="ratio_logit")
 idx<-parent$variables$.rowid;known<-mask_true(mask)&!is.na(outcome);des<-parent
 des$variables$.num<-ifelse(known[idx],outcome[idx],0);des$variables$.den<-as.numeric(known[idx])
 z<-tryCatch(svyratio(~.num,~.den,des),error=function(e)e)
 if(inherits(z,"error")) {row$status<-"NOT_ESTIMABLE";row$error<-conditionMessage(z);return(row)}
 est<-as.numeric(coef(z)[1]);se<-as.numeric(SE(z)[1]);row$estimate<-100*est;row$se<-100*se
 if(is.finite(est)&&is.finite(se)&&est>0&&est<1&&row$ci_df>0) {q<-qt(.975,row$ci_df);row$lower<-100*plogis(qlogis(est)-q*se/(est*(1-est)));row$upper<-100*plogis(qlogis(est)+q*se/(est*(1-est)));row$ci_method<-"full_parent_ratio_logit";row$status<-"ESTIMATED"}
 row
}
ref_table<-function(file) read_csv(file.path(B1,"results/trackA/rebuilt",file),show_col_types=FALSE)
add_cells<-function(measure,proxy) list(measurement=measure,proxy=proxy,both=as.integer(measure==1&proxy==1),measurement_only=as.integer(measure==1&proxy==0),proxy_only=as.integer(measure==0&proxy==1),neither=as.integer(measure==0&proxy==0))
scenario_estimates<-function(mask,measure,proxy,scenario,scope="scenario_specific",weight="WTPH2YR") {
 cells<-add_cells(measure,proxy);eligible<-mask_true(mask);known<-eligible&!is.na(measure)&!is.na(proxy);parent<-if(weight=="WTPH2YR") parents$phleb else parents$fasting
 out<-list();nmeas<-sum(known&measure==1);nproxy<-sum(known&proxy==1)
 for(v in names(cells)) {
  row<-estimate_domain(parent,known,cells[[v]],v,weight);row$scenario<-scenario;row$scope<-scope;row$candidate_n<-sum(eligible);row$comparison_n<-sum(known);row$measurement_positive_n<-nmeas;row$proxy_positive_n<-nproxy;out[[v]]<-row
 }
 row<-estimate_domain(parent,known&measure %in% 1,cells$measurement_only,"conditional_measurement_only",weight);row$scenario<-scenario;row$scope<-scope;row$candidate_n<-sum(eligible);row$comparison_n<-sum(known);row$measurement_positive_n<-nmeas;row$proxy_positive_n<-nproxy;out$conditional<-row
 bind_rows(out)
}
equal_missing_values<-function(a,b,tolerance=1e-10) identical(is.na(a),is.na(b))&&all(abs(a[!is.na(a)]-b[!is.na(b)])<=tolerance)
