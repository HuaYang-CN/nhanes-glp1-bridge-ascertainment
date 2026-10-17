source(file.path(RELEASE_ROOT,"src/a3/common.R"),encoding="UTF-8")
spec<-readr::read_csv(file.path(RELEASE_ROOT,"schemas/table1_variables.csv"),show_col_types=FALSE,na="NA")
spec<-spec[spec$destination!="Table 1",];cohort<-dat[primary,];rows<-list()
for(i in seq_len(nrow(spec))) {
 s<-spec[i,];v<-s$variable;x<-cohort[[v]]
 questionnaire<-v %in% c("MCQ160E","MCQ160F","MCQ160B","DIQ010","DIQ160","KIQ022")
 valid<-if(questionnaire) x %in% (if(v=="DIQ010") c(1,2,3) else c(1,2)) else !is.na(x)
 refused<-if(questionnaire)sum(x %in% 7) else 0L;dk<-if(questionnaire)sum(x %in% 9) else 0L
 structural<-if(v=="DIQ160")sum(is.na(x)&cohort$DIQ010 %in% 3) else 0L
 lab<-if(v %in% c("LBXGH","LBXSCR","egfr")) sum(is.na(x)) else 0L
 ordinary<-sum(is.na(x))-structural-lab
 check(paste0("missing_partition_",v),sum(!valid),refused+dk+structural+lab+ordinary)
 check(paste0("completeness_balance_",v),sum(valid)+sum(!valid),PRIMARY_N)
 if(questionnaire)check(paste0("documented_codes_",v),all(is.na(x)|x %in% (if(v=="DIQ010") c(1,2,3,7,9) else c(1,2,7,9))))
 if(v=="DIQ160") {check("DIQ160_expected_skip_not_observed",all(is.na(x[cohort$DIQ010 %in% 3])))}
 parts<-c(if(structural)paste("Expected structural routing skip",structural),if(refused)paste("Refused",refused),if(dk)paste("Don't know",dk),if(lab)paste(if(v=="egfr")"Missing derived value due to unavailable creatinine" else "Missing laboratory measurement",lab),if(ordinary)paste("Item/unknown NA",ordinary))
 if(!length(parts))parts<-"No missing/invalid values"
 rows[[i]]<-data.frame(result_id=sprintf("C03A-S1-%03d",i),variable=v,source_field=s$source_field,original_label=s$original_label,display_label=s$display_label,universe_n=PRIMARY_N,observed_n=sum(valid),missing_n=sum(!valid),unweighted_missing_percent=100*sum(!valid)/PRIMARY_N,raw_NA_n=sum(is.na(x)),refused_n=refused,dont_know_n=dk,expected_structural_skip_n=structural,missing_lab_or_derived_n=lab,ordinary_item_or_unknown_NA_n=ordinary,missingness_type=paste(parts,collapse="; "),yes_n=if(questionnaire)sum(x %in% 1)else NA_integer_,no_n=if(questionnaire)sum(x %in% 2)else NA_integer_,borderline_n=if(v=="DIQ010")sum(x %in% 3)else NA_integer_,applicable_universe=s$applicable_universe,weight="WTPH2YR",percentage_method="Unweighted missing count / computed primary universe",caveat=s$interpretive_caveat,documentation_url=s$documentation_url,destination="tables/Stage03A_06_TableS1_Preview.md",code_path="code/03_generate_completeness.R",input_sha256=INPUT_SHA256,verification_status="PENDING_INDEPENDENT_QA")
}
out(do.call(rbind,rows),"completeness_full_precision.csv")
out(do.call(rbind,checks),"completeness_generation_QA.csv");assert_checks()
cat("Generated",nrow(spec),"completeness rows; questionnaire unknowns/skips remain missing.\n")
