source(file.path(RELEASE_ROOT,"src/b2/common.R"),encoding="UTF-8")
raw<-read_xpt(file.path(RAW_DIR,"DIQ_L.xpt"));stopifnot(!anyNA(raw$SEQN),!anyDuplicated(raw$SEQN))
diq_present<-dat$SEQN %in% raw$SEQN
code_category<-function(x) case_when(x==1~"Yes",x==2~"No",x==3~"Borderline",x==7~"Refused",x==9~"Don't know",is.na(x)~"Missing",TRUE~"Unexpected code")
diabetes<-code_category(dat$DIQ010);pred<-code_category(dat$DIQ160)
expected_skip<-(!is.na(dat$DIQ010)&dat$DIQ010==1)|(!is.na(dat$RIDAGEYR)&dat$RIDAGEYR<12)|(!is.na(dat$RIDAGEYR)&dat$RIDAGEYR>=12&dat$DIQ010 %in% 3)
routing<-case_when(!diq_present~"DIQ component record absent",!is.na(dat$DIQ160)&expected_skip~"Observed response under expected skip",dat$DIQ160==1~"Observed Yes",dat$DIQ160==2~"Observed No",dat$DIQ160==7~"Observed Refused",dat$DIQ160==9~"Observed Don't know",is.na(dat$DIQ160)&dat$DIQ010 %in% 1~"Structural skip: reported diabetes",is.na(dat$DIQ160)&dat$RIDAGEYR<12~"Structural skip: age below12",is.na(dat$DIQ160)&dat$RIDAGEYR>=12&dat$DIQ010 %in% 3~"Structural skip: Borderline",is.na(dat$DIQ160)&(is.na(dat$DIQ010)|is.na(dat$RIDAGEYR))~"Missing: unresolved routing",is.na(dat$DIQ160)~"Other missing: expected question",TRUE~"Unexpected response code")
rows<-list()
for(domain in list(list("all_public_DEMO",rep(TRUE,nrow(dat))),list("age20_Medicare",mask_true(dat$RIDAGEYR>=20&dat$medicare==1)),list("base_eligible",base),list("primary_known",primary),list("original_common",comparison),list("primary_proxy_unknown",primary&is.na(dat$diagnosis_proxy)))) {
 mask<-domain[[2]];x<-tibble(DIQ010_category=diabetes[mask],DIQ160_category=pred[mask],routing=routing[mask],proxy_unknown=is.na(dat$diagnosis_proxy[mask])) |> count(DIQ010_category,DIQ160_category,routing,proxy_unknown,name="n");x$record_type<-"routing_frequency";x$population<-domain[[1]];rows[[length(rows)+1]]<-x
}
history_D2<-dat$dx_prediabetes
assumed_borderline<-dat$RIDAGEYR>=12&dat$DIQ010 %in% 3
history_D2[mask_true(assumed_borderline)]<-1L
proxy_D2<-tri_or(dat$crit_bmi35,dat$crit_dx_kidney,bmi_binary(dat$BMXBMI,27,history_D2),dat$crit_mi,dat$crit_stroke)
scenarios<-list(D0=list(mask=primary,proxy=dat$diagnosis_proxy,assumption="Original DIQ160-only glycemic history"),D1=list(mask=primary&dat$DIQ010 %in% 2,proxy=dat$diagnosis_proxy,assumption="Exclude Borderline respondents"),D2=list(mask=primary,proxy=proxy_D2,assumption="DIQ010 Borderline is assumed possible history, raw DIQ160 unchanged"),D3=list(mask=primary,proxy=dat$diagnosis_proxy,assumption="Structural DIQ160 remains unknown, explicit ascertainability"))
known<-lapply(scenarios,function(z) z$mask&!is.na(z$proxy));all_common<-Reduce(`&`,known);d0d2_common<-known$D0&known$D2
out<-list()
for(name in names(scenarios)) {
 z<-scenarios[[name]]
 for(scope in c("scenario_specific","all_scenarios_common",if(name %in% c("D0","D2")) "D0_D2_common" else character())) {
  mask<-switch(scope,scenario_specific=z$mask,all_scenarios_common=all_common,D0_D2_common=d0d2_common)
  tab<-scenario_estimates(mask,dat$primary_observable,z$proxy,name,scope);tab$analytical_assumption<-z$assumption;tab$record_type<-"scenario_estimate";out[[length(out)+1]]<-tab
 }
}
changes<-tibble(SEQN=dat$SEQN,routing=routing,original_proxy=dat$diagnosis_proxy,D2_proxy=proxy_D2,primary_universe=primary,newly_ascertainable=primary&is.na(dat$diagnosis_proxy)&!is.na(proxy_D2),lost_ascertainability=primary&!is.na(dat$diagnosis_proxy)&is.na(proxy_D2),known_proxy_value_changed=primary&!is.na(dat$diagnosis_proxy)&!is.na(proxy_D2)&dat$diagnosis_proxy!=proxy_D2,assumed_borderline=mask_true(assumed_borderline),raw_DIQ160=dat$DIQ160)
save_csv(bind_rows(rows),"diq_routing_frequencies.csv")
save_csv(bind_rows(out),"diq_scenario_estimates.csv")
save_csv(changes |> filter(primary_universe),"diq_record_transitions.csv","diagnostics")
summary<-changes |> filter(primary_universe) |> group_by(routing) |> summarise(n=n(),proxy_unknown=sum(is.na(original_proxy)),newly_ascertainable=sum(newly_ascertainable),lost_ascertainability=sum(lost_ascertainability),known_proxy_changed=sum(known_proxy_value_changed),.groups="drop")
save_csv(summary,"diq_proxy_missingness_attribution.csv","diagnostics")
raw_value_check<-equal_missing_values(dat$DIQ160[match(raw$SEQN,dat$SEQN)],raw$DIQ160,tolerance=0)
save_csv(tibble(variable="DIQ160",raw_aligned_values_and_missingness_exact=raw_value_check,n_compared=nrow(raw),note="Subsetting a numeric R vector drops label attributes; compare values/missingness independently of display attributes"),"DIQ_raw_input_value_integrity.csv","diagnostics")
stopifnot(raw_value_check)
stopifnot(identical(scenarios$D0$proxy,scenarios$D3$proxy))
print(bind_rows(out) |> filter(estimand=="conditional_measurement_only") |> select(scenario,scope,candidate_n,comparison_n,measurement_positive_n,positive_n,estimate,lower,upper))
