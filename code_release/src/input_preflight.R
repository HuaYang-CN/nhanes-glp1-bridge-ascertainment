schema<-read.csv(file.path(RELEASE_ROOT,"schemas/input_variables.csv"),stringsAsFactors=FALSE,check.names=FALSE)
components<-unique(schema$component);frames<-list();checks<-list()
for(stem in components) {
 p<-file.path(RAW_DIR,paste0(stem,".xpt"))
 if(!file.exists(p))stop("Missing public NHANES component: ",stem)
 x<-haven::read_xpt(p);names(x)<-toupper(names(x));spec<-schema[schema$component==stem,]
 if(!all(spec$variable %in% names(x)))stop("Missing required variable in ",stem,": ",paste(setdiff(spec$variable,names(x)),collapse=","))
 if(anyNA(x$SEQN)||anyDuplicated(x$SEQN)||any(!is.finite(x$SEQN)|x$SEQN!=round(x$SEQN)|x$SEQN<=0))stop("Invalid/duplicate participant identifiers: ",stem)
 for(i in seq_len(nrow(spec))) {
  v<-spec$variable[i];label<-attr(x[[v]],"label");if(is.null(label))label<-""
  if(!is.numeric(x[[v]]))stop("Non-numeric required source field: ",stem,"/",v)
  if(label!=spec$label[i])stop("Source label/unit mismatch: ",stem,"/",v)
  checks[[length(checks)+1]]<-data.frame(component=stem,variable=v,label=label,n=nrow(x),missing=sum(is.na(x[[v]])),status="PASS")
 }
 for(w in intersect(c("WTPH2YR","WTSAF2YR"),names(x)))if(any(!is.na(x[[w]])&(!is.finite(x[[w]])|x[[w]]<0)))stop("Invalid weight: ",stem,"/",w)
 frames[[stem]]<-x$SEQN
}
for(stem in names(frames))if(any(!frames[[stem]] %in% frames$DEMO_L))stop("Component identifier not present in DEMO_L: ",stem)
write_csv(bind_rows(checks),file.path(OUTPUT_DIR,"preflight/input_schema_checks.csv"))
cat("Input component/ID/label/unit/type preflight PASS\n")
