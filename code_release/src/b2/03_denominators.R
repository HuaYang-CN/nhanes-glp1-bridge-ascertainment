source(file.path(RELEASE_ROOT,"src/b2/common.R"),encoding="UTF-8")
cohorts<-list(all_public=rep(TRUE,nrow(dat)),base_eligible=base,primary=primary,expanded=expanded,a1c_exclusion=a1c_exclusion,proxy_known=base&!is.na(dat$diagnosis_proxy),comparison=comparison,conditional=cond,fasting_candidate=fast_candidate,fasting_complete=fast)
old<-read_csv(file.path(B1,"diagnostics/analytic_population_membership.csv"),show_col_types=FALSE) |> filter(lane=="rebuilt")
idx<-match(dat$SEQN,old$SEQN);stopifnot(!anyNA(idx))
membershipcheck<-bind_rows(lapply(c("primary","expanded","a1c_exclusion","comparison","measurement_positive","fasting"),function(v) {current<-switch(v,measurement_positive=cond,fasting=fast,cohorts[[v]]);tibble(domain=v,current_n=sum(current),Stage02B1_n=sum(old[[v]]),changed_members=sum(current!=old[[v]][idx]),status=if(all(current==old[[v]][idx])) "MATCH" else "DISCREPANCY")}))
save_csv(membershipcheck,"independent_membership_verification.csv","diagnostics")
stopifnot(all(membershipcheck$changed_members==0))
flow_masks<-list(All_public=rep(TRUE,nrow(dat)),Age20=mask_true(dat$RIDAGEYR>=20),Medicare=mask_true(dat$RIDAGEYR>=20&dat$medicare==1),Positive_phlebotomy_weight=mask_true(dat$RIDAGEYR>=20&dat$medicare==1&dat$WTPH2YR>0),No_reported_diabetes=mask_true(dat$RIDAGEYR>=20&dat$medicare==1&dat$WTPH2YR>0&dat$no_reported_diabetes==1),BMI_observed=base,Primary_known=primary,Proxy_comparison_known=comparison)
flow<-list();previous<-rep(TRUE,nrow(dat))
for(name in names(flow_masks)) {now<-flow_masks[[name]];flow[[name]]<-tibble(record_type="selection_flow",domain=name,source_population_n=nrow(dat),eligible_population_n=sum(base),universe_n=sum(previous),nonmissing_n=sum(now),missing_n=sum(previous&!now),population_removed_n=sum(previous&!now),reason=paste("Sequential loss at",name));previous<-now}
out<-list();summary<-list()
cohort_outcomes<-list(base_eligible=dat$primary_observable,primary=dat$primary_observable,expanded=dat$expanded_proxy,a1c_exclusion=dat$primary_observable,proxy_known=dat$diagnosis_proxy,comparison=dat$primary_observable,conditional=dat$measurement_only,fasting_candidate=dat$.fast_a1c,fasting_complete=dat$.fast_a1c)
for(name in names(cohort_outcomes)) {
 mask<-cohorts[[name]];value<-cohort_outcomes[[name]];w<-if(startsWith(name,"fasting")) "WTSAF2YR" else "WTPH2YR";parent<-if(w=="WTPH2YR") parents$phleb else parents$fasting
 row<-estimate_domain(parent,mask,value,name,w);row$record_type<-"cohort_denominator";row$domain<-name;row$source_population_n<-nrow(dat);row$eligible_population_n<-if(w=="WTPH2YR") sum(base) else sum(fast_candidate);row$ascertainable_population_n<-sum(mask&!is.na(value));row$conditional_domain_n<-if(name=="conditional") sum(cond) else NA_integer_;row$nonmissing_n<-row$n;row$missing_n<-row$missing_outcome_n;row$missing_percent<-100*row$missing_outcome_n/row$universe_n;out[[length(out)+1]]<-row
}
for(v in names(paths)) {
 row<-estimate_domain(parents$phleb,primary,dat[[v]],v,"WTPH2YR");row$record_type<-"criterion_route";row$domain<-v;row$source_population_n<-nrow(dat);row$eligible_population_n<-sum(base);row$ascertainable_population_n<-sum(primary);row$nonmissing_n<-row$n;row$missing_n<-row$missing_outcome_n;row$missing_percent<-100*row$missing_n/row$universe_n;row$reason<-switch(v,crit_bmi35="BMI known in primary universe",crit_egfr="Unknown eGFR at BMI>=30",crit_a1c_prediabetes="Unknown HbA1c at BMI>=27",crit_mi="Unknown MI history at BMI>=27",crit_stroke="Unknown stroke history at BMI>=27",crit_chf_proxy="Unknown CHF history at BMI>=30",crit_bp_proxy="Unknown treatment or measured BP status at BMI>=30");out[[length(out)+1]]<-row
}
overlap<-list();names_cohort<-setdiff(names(cohorts),"all_public")
for(i in seq_along(names_cohort)) for(j in seq_along(names_cohort)) if(j>i) {a<-names_cohort[i];b<-names_cohort[j];ma<-cohorts[[a]];mb<-cohorts[[b]];overlap[[length(overlap)+1]]<-tibble(cohort_a=a,cohort_b=b,n_a=sum(ma),n_b=sum(mb),intersection_n=sum(ma&mb),a_only_n=sum(ma&!mb),b_only_n=sum(!ma&mb))}
changes<-tibble(SEQN=dat$SEQN,base=base,primary=primary,expanded=expanded,a1c_exclusion=a1c_exclusion,comparison=comparison,fasting_candidate=fast_candidate,fasting_complete=fast,primary_status=dat$primary_observable,expanded_status=dat$expanded_proxy,LBXGH=dat$LBXGH,LBXGLU=dat$LBXGLU,expanded_unknown=is.na(dat$expanded_proxy),HbA1c_high=dat$LBXGH>=6.5,proxy_unknown=is.na(dat$diagnosis_proxy),missing_HbA1c=is.na(dat$LBXGH),missing_FPG=is.na(dat$LBXGLU))
addedexpanded<-base&!primary&expanded;lostexpanded<-primary&!expanded;both<-primary&expanded
transitions<-tibble(analysis=c("expanded_vs_primary","HbA1c_exclusion_vs_primary","comparison_vs_primary","fasting_complete_vs_positive_weight"),source_n=c(sum(primary),sum(primary),sum(primary),sum(fast_candidate)),target_n=c(sum(expanded),sum(a1c_exclusion),sum(comparison),sum(fast)),added_n=c(sum(addedexpanded),0,0,0),excluded_n=c(sum(lostexpanded),sum(primary&!a1c_exclusion),sum(primary&!comparison),sum(fast_candidate&!fast)),same_person_common_n=c(sum(both),sum(a1c_exclusion),sum(comparison),sum(fast)),known_endpoint_reclassification_n=c(sum(both&dat$primary_observable!=dat$expanded_proxy),0,0,0),reason=c("Unknown expanded status may remove primary negatives; extra route positivity may resolve unknown primary", "Observed HbA1c>=6.5 excluded; missing HbA1c retained", "Unknown diagnosis-proxy composite removed", "Missing lab or unknown paired definition after positive fasting weight"),original_estimates_same_people=c(FALSE,FALSE,FALSE,FALSE))
common_rows<-list()
for(v in c("primary_observable","expanded_proxy")) {row<-estimate_domain(parents$phleb,both,dat[[v]],v,"WTPH2YR");row$record_type<-"common_primary_expanded";row$domain<-v;row$paired_comparison<-TRUE;common_rows[[v]]<-row}
counts <- read_csv(file.path(B1,"diagnostics/denominators_and_numerators.csv"),show_col_types=FALSE) |> filter(lane=="rebuilt",grepl("^Pathway:",domain))
routes<-bind_rows(out) |> filter(record_type=="criterion_route")
check<-tibble(route=routes$domain,current_nonmissing=routes$n,Stage02B1_nonmissing=counts$outcome_nonmissing_n[match(paste0("Pathway: ",routes$domain),counts$domain)],current_missing=routes$missing_n)
check$status<-ifelse(check$current_nonmissing==check$Stage02B1_nonmissing,"MATCH","DISCREPANCY")
save_csv(bind_rows(unname(c(flow,out,common_rows))),"denominator_audit.csv")
save_csv(bind_rows(overlap),"cohort_overlap_counts.csv","diagnostics")
save_csv(transitions,"sensitivity_population_transitions.csv","diagnostics")
save_csv(changes |> filter(base|fast_candidate),"population_membership_and_reasons.csv","diagnostics")
save_csv(check,"independent_pathway_count_verification.csv","diagnostics")
print(transitions);print(routes |> select(domain,universe_n,n,missing_n,missing_percent,estimate,lower,upper))
