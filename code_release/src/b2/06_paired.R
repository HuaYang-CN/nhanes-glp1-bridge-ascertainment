source(file.path(RELEASE_ROOT,"src/b2/common.R"),encoding="UTF-8")
out<-list();identities<-list()
for(spec in list(list("Measurement_minus_proxy",comparison,dat$primary_observable,dat$diagnosis_proxy,"WTPH2YR"),list("HbA1c_FPG_minus_HbA1c",fast,dat$.fast_aug,dat$.fast_a1c,"WTSAF2YR"))) {
 name<-spec[[1]];mask<-spec[[2]];a<-spec[[3]];b<-spec[[4]];w<-spec[[5]];parent<-if(w=="WTPH2YR") parents$phleb else parents$fasting
 known<-mask_true(mask)&!is.na(a)&!is.na(b);stopifnot(sum(known)==sum(mask))
 des<-parent;idx<-des$variables$.rowid;des$variables$.a<-a[idx];des$variables$.b<-b[idx];des$variables$.delta<-a[idx]-b[idx];des<-des[known[idx],]
 cov<-coverage(parent,known);df<-cov$represented_df
 result<-tryCatch({
  joint<-svymean(~.a+.b,des);delta<-svymean(~.delta,des);est<-as.numeric(coef(delta)[1]);se<-as.numeric(SE(delta)[1]);v<-vcov(joint);joint_se<-sqrt(v[1,1]+v[2,2]-2*v[1,2]);a_est<-as.numeric(coef(joint)[1]);b_est<-as.numeric(coef(joint)[2])
  forward<-as.integer(a==1&b==0);reverse<-as.integer(a==0&b==1);fo<-estimate_domain(parent,known,forward,"forward_discordance",w);re<-estimate_domain(parent,known,reverse,"reverse_discordance",w)
  valid<-is.finite(se)&&se>0&&df>0;ci<-if(valid) est+c(-1,1)*qt(.975,df)*se else c(NA_real_,NA_real_);p<-if(valid) 2*pt(abs(est/se),df,lower.tail=FALSE) else NA_real_
  row<-tibble(contrast=name,track="TrackB_PAIR",weight_variable=w,n=sum(known),estimate_a=100*a_est,estimate_b=100*b_est,difference_pp=100*est,se_pp=100*se,lower_pp=100*ci[1],upper_pp=100*ci[2],design_df=df,parent_df=degf(parent),contributing_strata=cov$strata,contributing_psus=cov$psus,p_value=p,p_value_method="Two-sided design-based t, paired indicator difference",ci_method="Taylor paired mean difference, t interval",inference_scope="Exploratory Stage02B2, not original prespecified testing",forward_discordance_n=sum(known&forward==1),reverse_discordance_n=sum(known&reverse==1),forward_weighted_percent=fo$estimate,reverse_weighted_percent=re$estimate,nominal_alpha=.05,status=if(valid) "ESTIMATED" else "NOT_ESTIMABLE",error=if(valid) "" else "Nonfinite/zero SE or nonpositive df",lonely_psu=getOption("survey.lonely.psu"),adjust_domain_lonely=getOption("survey.adjust.domain.lonely"),displayed_input_subtraction_pp=round(100*a_est,1)-round(100*b_est,1),rounded_full_precision_difference_pp=round(100*est,1))
  check<-tibble(contrast=name,delta_minus_marginal_difference_pp=100*(est-(a_est-b_est)),delta_minus_discordance_difference_pp=100*est-(fo$estimate-re$estimate),paired_SE_minus_joint_covariance_SE_pp=100*(se-joint_se),same_people=sum(known)==sum(mask),same_weight_and_design=TRUE,identity_pass=abs(est-(a_est-b_est))<=1e-10&&abs(100*est-(fo$estimate-re$estimate))<=1e-8&&abs(se-joint_se)<=1e-10)
  identities[[name]]<-check;row
 },error=function(e) tibble(contrast=name,track="TrackB_PAIR",weight_variable=w,n=sum(known),status="NOT_ESTIMABLE",error=conditionMessage(e)))
 out[[name]]<-result
}
tab<-bind_rows(out);tab$p_value_Holm_two_contrasts<-p.adjust(tab$p_value,method="holm");tab$nominal_inference<-ifelse(is.na(tab$p_value),"NOT_ESTIMABLE",ifelse(tab$p_value<.05,"Nominal P below0.05","Nominal P at least0.05"));tab$Holm_inference<-ifelse(is.na(tab$p_value_Holm_two_contrasts),"NOT_ESTIMABLE",ifelse(tab$p_value_Holm_two_contrasts<.05,"Holm P below0.05","Holm P at least0.05"))
save_csv(tab,"paired_contrasts.csv");save_csv(bind_rows(identities),"paired_covariance_and_identity_checks.csv","diagnostics")
print(tab |> select(contrast,n,difference_pp,se_pp,lower_pp,upper_pp,design_df,p_value,p_value_Holm_two_contrasts));print(bind_rows(identities))
