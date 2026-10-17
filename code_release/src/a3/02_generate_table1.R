source(file.path(RELEASE_ROOT,"src/a3/common.R"),encoding="UTF-8")
if(!file.exists(file.path(RELEASE_ROOT,"schemas/table1_variables.csv"))) stop("Missing package variable specification")
spec <- readr::read_csv(file.path(RELEASE_ROOT,"schemas/table1_variables.csv"),show_col_types=FALSE,na=c("NA"))
dat$age_category <- cut(dat$RIDAGEYR,c(20,65,75,Inf),right=FALSE,labels=c("20-64","65-74","75 or older"))
dat$bmi_category <- cut(dat$BMXBMI,c(-Inf,27,30,35,Inf),right=FALSE,labels=c("<27","27 to <30","30 to <35","35 or higher"))
dat$hba1c_category <- cut(dat$LBXGH,c(-Inf,5.7,6.5,Inf),right=FALSE,labels=c("<5.7","5.7 to <6.5","6.5 or higher"))
dat$egfr_low <- as.integer(dat$egfr<60)
par <- make_parent(); dom<-subset(par,.primary)
source_rows<-list();occ_rows<-list();binary_rows<-list()
for(i in which(spec$summary_method!="completeness")) {
 s<-spec[i,];v<-s$variable;raw<-dom$variables[[v]]
 if(s$summary_method=="continuous") {
   valid<-!is.na(raw);check(paste0("finite_",v),all(is.finite(raw[valid])))
   levels<-NA_character_;codes<-NA_character_;labels<-NA_character_
 } else {
   encoded<-strsplit(s$categorical_levels,"|",fixed=TRUE)[[1]]
   if(any(grepl("=",encoded,fixed=TRUE))) {
     pieces<-strsplit(encoded,"=",fixed=TRUE);codes<-vapply(pieces,`[`,character(1),1);labels<-vapply(pieces,`[`,character(1),2)
   } else {codes<-encoded;labels<-encoded}
   if(v %in% c("MCQ160B","MCQ160E","MCQ160F")) {check(paste0("source_codes_",v),all(is.na(raw)|raw %in% c(1,2,7,9)));valid<-raw %in% c(1,2);codes<-"1";labels<-"Yes"}
   else {valid<-!is.na(raw);check(paste0("source_codes_",v),all(as.character(raw[valid]) %in% codes))}
   if(v=="egfr_low") {codes<-"1";labels<-"eGFR below 60"}
 }
 d<-dom[valid,];n<-nrow(d$variables);miss<-PRIMARY_N-n
 strata<-length(unique(d$variables$SDMVSTRA));psus<-nrow(unique(d$variables[c("SDMVSTRA","SDMVPSU")]))
 df<-psus-strata
 check(paste0("observed_missing_",v),n+miss,PRIMARY_N)
 basepsu<-unique(par$variables[c("SDMVSTRA","SDMVPSU")]);basepsu$variable<-v
 basepsu$valid_n<-vapply(seq_len(nrow(basepsu)),function(j)sum(d$variables$SDMVSTRA==basepsu$SDMVSTRA[j]&d$variables$SDMVPSU==basepsu$SDMVPSU[j]),integer(1))
 occ_rows[[v]]<-basepsu
 single<-sum(table(basepsu$SDMVSTRA[basepsu$valid_n>0])==1)
 cats<-if(s$summary_method=="continuous") "" else labels
 for(k in seq_along(cats)) {
   x<-if(s$summary_method=="continuous") as.numeric(d$variables[[v]]) else as.numeric(as.character(d$variables[[v]])==codes[k])
   d$variables$.value<-x
   ans<-svymean(~.value,d);m<-as.numeric(coef(ans)[1]);se<-as.numeric(SE(ans)[1]);mult<-if(s$summary_method=="continuous") 1 else 100
   limits<-c(NA_real_,NA_real_);status<-"NOT_ESTIMABLE";reason<-"";method<-if(s$summary_method=="continuous") "Taylor mean; Student t" else "Taylor proportion; design-based logit"
   if(df<=0) reason<-"Nonpositive represented valid-observation-domain df"
   else if(s$summary_method!="continuous" && (m<=0||m>=1)) reason<-"Boundary proportion; logit CI not estimable"
   else {
     limits<-tryCatch(if(s$summary_method=="continuous") m+c(-1,1)*qt(.975,df)*se else as.numeric(confint(svyciprop(~.value,d,method="logit",df=df))),error=function(e){reason<<-conditionMessage(e);c(NA_real_,NA_real_)})
     if(all(is.finite(limits))) status<-"ESTIMATED" else if(reason=="") reason<-"Nonfinite CI"
   }
   flags<-c(if(s$summary_method!="continuous"&&sum(x==1)<30) "category_n_below30",if(s$summary_method!="continuous"&&m>0&&se/m>.3) "RSE_above30percent",if(df<8) "represented_df_below8",if(single>0) "domain_single_PSU_strata")
   source_rows[[length(source_rows)+1]]<-data.frame(result_id=sprintf("C03A-T1-%03d",length(source_rows)+1),variable=v,statistic=if(s$summary_method=="continuous") "weighted_mean" else "weighted_percentage",category=cats[k],category_code=if(s$summary_method=="continuous") "" else codes[k],source_field=s$source_field,original_label=s$original_label,definition=s$missing_rule,population="Original verified primary analytical population",universe_n=PRIMARY_N,observed_n=n,missing_n=miss,category_n=if(s$summary_method=="continuous") NA_integer_ else sum(x==1),weight="WTPH2YR",weighted_denominator=sum(weights(d)),estimate=mult*m,se=mult*se,lower_ci=mult*limits[1],upper_ci=mult*limits[2],ci_method=method,df=df,parent_df=as.numeric(degf(par)),represented_strata=strata,represented_PSUs=psus,single_PSU_domain_strata=single,unit=s$unit,display_label=s$display_label,status=status,not_estimable_reason=reason,reliability_notes=paste(flags,collapse=";"),lonely_psu="adjust",adjust_domain_lonely=FALSE,destination="tables/Stage03A_04_Table1_Preview.md",code_path="code/02_generate_table1.R",input_path=file.path(B1,"diagnostics/rebuilt/bridge_analysis_ready.rds"),input_sha256=INPUT_SHA256,documentation_url=s$documentation_url,verification_status="PENDING_INDEPENDENT_QA",analysis_timing="Newly specified Stage03A descriptive analysis")
 }
 if(s$summary_method=="yes") binary_rows[[v]]<-data.frame(variable=v,observed_n=n,yes_n=sum(x==1),no_n=sum(x==0),missing_n=miss)
 if(s$summary_method=="categorical") {
   counts<-vapply(codes,function(z)sum(as.character(raw[valid])==z),integer(1))
   check(paste0("category_partition_",v),sum(counts),n)
 }
}
out(do.call(rbind,source_rows),"statistics_full_precision.csv")
out(do.call(rbind,occ_rows),"valid_domain_PSU_occupancy.csv")
out(do.call(rbind,binary_rows),"binary_partitions.csv")
out(do.call(rbind,checks),"table1_generation_QA.csv");assert_checks()
cat("Generated",length(source_rows),"Table 1 statistical rows. NOT_ESTIMABLE:",sum(vapply(source_rows,function(x)x$status!="ESTIMATED",logical(1))),"\n")
