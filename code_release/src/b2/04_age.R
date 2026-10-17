source(file.path(RELEASE_ROOT,"src/b2/common.R"),encoding="UTF-8")
top<-dat$RIDAGEYR %in% 80;assumed_age<-dat$RIDAGEYR;assumed_age[top]<-85
e0<-egfr_formula(dat$LBXSCR,dat$RIDAGEYR,dat$RIAGENDR);stopifnot(equal_missing_values(e0,dat$egfr))
e1<-egfr_formula(dat$LBXSCR,assumed_age,dat$RIAGENDR)
k0<-bmi_binary(dat$BMXBMI,30,ifelse(is.na(e0),NA_integer_,as.integer(e0<60)));k1<-bmi_binary(dat$BMXBMI,30,ifelse(is.na(e1),NA_integer_,as.integer(e1<60)))
m1<-tri_or(dat$crit_bmi35,k1,dat$crit_a1c_prediabetes,dat$crit_mi,dat$crit_stroke);m0<-dat$primary_observable
stopifnot(equal_missing_values(k0,dat$crit_egfr))
newprimary<-base&!is.na(m1);commonprimary<-primary&newprimary;newcomparison<-newprimary&!is.na(dat$diagnosis_proxy);commoncomparison<-comparison&newcomparison
trans<-tibble(SEQN=dat$SEQN,topcoded=top,base=base,original_primary_domain=primary,original_comparison_domain=comparison,eGFR_E0=e0,eGFR_E1=e1,eGFR_numeric_changed=!is.na(e0)&!is.na(e1)&abs(e1-e0)>1e-10,crossed_egfr60=!is.na(e0)&!is.na(e1)&(e0<60)!=(e1<60),route_E0=k0,route_E1=k1,route_changed=!is.na(k0)&!is.na(k1)&k0!=k1,primary_E0=m0,primary_E1=m1,composite_changed=!is.na(m0)&!is.na(m1)&m0!=m1,unknownness_changed=is.na(m0)!=is.na(m1),measurement_only_E0=as.integer(m0==1&dat$diagnosis_proxy==0),measurement_only_E1=as.integer(m1==1&dat$diagnosis_proxy==0),other_primary_route_positive=tri_or(dat$crit_bmi35,dat$crit_a1c_prediabetes,dat$crit_mi,dat$crit_stroke)==1)
scope_counts<-list()
for(domain in list(list("all_public",rep(TRUE,nrow(dat))),list("base_eligible",base),list("original_primary",primary),list("original_comparison",comparison))) {
 x<-trans[domain[[2]],];scope_counts[[domain[[1]]]]<-tibble(domain=domain[[1]],universe_n=nrow(x),topcoded_n=sum(x$topcoded),egfr_numeric_changed_n=sum(x$eGFR_numeric_changed),egfr60_crossing_n=sum(x$crossed_egfr60),route_changed_n=sum(x$route_changed),composite_changed_n=sum(x$composite_changed),unknownness_changed_n=sum(x$unknownness_changed),measurement_only_changed_n=sum(!is.na(x$measurement_only_E0)&!is.na(x$measurement_only_E1)&x$measurement_only_E0!=x$measurement_only_E1),route_change_already_positive_elsewhere_n=sum(x$route_changed&x$other_primary_route_positive,na.rm=TRUE))
}
dist<-list()
for(domain in list(list("all_public_topcoded",top),list("primary_topcoded",top&primary))) for(scenario in c("E0","E1")) {
 x<-if(scenario=="E0") e0[domain[[2]]] else e1[domain[[2]]];known<-x[!is.na(x)]
 q<-if(length(known)) quantile(known,c(0,.25,.5,.75,1)) else rep(NA_real_,5)
 dist[[length(dist)+1]]<-tibble(domain=domain[[1]],scenario=scenario,n=length(x),nonmissing_n=length(known),mean=mean(known),sd=sd(known),min=q[1],q25=q[2],median=q[3],q75=q[4],max=q[5])
}
out<-list()
for(scenario in c("E0","E1")) {
 m<-if(scenario=="E0") m0 else m1;km<-if(scenario=="E0") k0 else k1;pm<-if(scenario=="E0") primary else newprimary
 for(scope in c("scenario_specific","common_population")) {
  mask<-if(scope=="scenario_specific") pm else commonprimary
  row<-estimate_domain(parents$phleb,mask,m,"primary_observable","WTPH2YR");row$scenario<-scenario;row$scope<-scope;row$assumption<-if(scenario=="E0") "Public age80 retained" else "Illustrative age85 for topcoded group, not recovered age";out[[length(out)+1]]<-row
  row<-estimate_domain(parents$phleb,mask,km,"BMI_eGFR_route","WTPH2YR");row$scenario<-scenario;row$scope<-scope;out[[length(out)+1]]<-row
  comparisonmask<-if(scope=="common_population") commoncomparison else pm&!is.na(dat$diagnosis_proxy)
  tab<-scenario_estimates(comparisonmask,m,dat$diagnosis_proxy,scenario,scope);out[[length(out)+1]]<-tab
 }
}
save_csv(bind_rows(out),"age_eGFR_estimates.csv")
save_csv(bind_rows(scope_counts),"age_reclassification_counts.csv","diagnostics")
save_csv(bind_rows(dist),"age_eGFR_distributions.csv","diagnostics")
save_csv(trans |> filter(top|route_changed|composite_changed),"age_record_transitions.csv","diagnostics")
print(bind_rows(scope_counts));print(bind_rows(out) |> filter(scope=="scenario_specific",estimand %in% c("primary_observable","conditional_measurement_only")) |> select(scenario,estimand,n,estimate,lower,upper))
