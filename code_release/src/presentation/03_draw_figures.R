# Display-only graphics. Never loads participant data or executes survey estimators.
out <- normalizePath(file.path(OUTPUT_DIR,"presentation"), winslash="/", mustWork=TRUE)
library(grid)
# Library configuration is supplied by the portable process-local runtime.
stopifnot(requireNamespace("readr",quietly=TRUE))
d <- as.data.frame(readr::read_csv(file.path(out,"results/panel_plot_data.csv"),col_types=readr::cols(.default=readr::col_character()),na=character(),show_col_types=FALSE))
stopifnot(nrow(d)==35,length(unique(d$figure))==6)
stopifnot(!anyDuplicated(d$datum_id))
pal <- c(blue="#0072B2",orange="#D55E00",grey="#666666")
font <- "Times New Roman"
gptext <- function(size=8,bold=FALSE,col="black") gpar(fontfamily=font,fontsize=size,fontface=if(bold)2 else 1,col=col)
txt <- function(s,x,y,size=8,bold=FALSE,just="left",col="black") grid.text(s,x=x,y=y,just=just,gp=gptext(size,bold,col))
marker <- function(x,y,pch,col){
 gp <- gpar(col=col,fill=col,lwd=.65)
 if(pch==16)grid.circle(x,y,r=unit(.8,"mm"),gp=gp)
 else if(pch==15)grid.rect(x,y,width=unit(1.6,"mm"),height=unit(1.6,"mm"),gp=gp)
 else if(pch==17)grid.polygon(x=unit(x,"npc")+unit(c(-.9,0,.9),"mm"),y=unit(y,"npc")+unit(c(-.7,1,-.7),"mm"),gp=gp)
 else grid.polygon(x=unit(x,"npc")+unit(c(-.9,0,.9,0),"mm"),y=unit(y,"npc")+unit(c(0,.9,0,-.9),"mm"),gp=gp)
}
spec <- list(
 Figure_1=list(h=c(59,80),x=c(60,40),heads=c("Observable criteria and sensitivity definitions","Criterion-specific clinical routes")),
 Figure_2=list(h=c(66,42,42),x=c(70,50,15),heads=c("D0 identification categories","Conditional measurement-only ascertainment","Within-person prevalence difference (exploratory/post hoc)")),
 Figure_3=list(h=c(58,44),x=c(60,15),heads=c("Fasting-glucose sensitivity in the same participants","Paired incremental ascertainment (exploratory/post hoc)")),
 Supplementary_Figure_S1=list(h=172,x=NA,heads=""),
 Supplementary_Figure_S2=list(h=85,x=50,heads="Conditional measurement-only by questionnaire scenario"),
 Supplementary_Figure_S3=list(h=c(54,68),x=c(60,60),heads=c("Age-topcoding assumptions: primary criteria","Laboratory diabetes-range restriction: fasting sensitivity")))
panel <- function(fig,p,xmax,heading){
 r <- d[d$figure==fig & d$panel==p,,drop=FALSE]
 if(fig=="Supplementary_Figure_S1") return(flow(r,p))
 contrast <- all(r$unit=="percentage_points") || any(grepl("PAIR",r$result_id))
 xmin <- if(contrast) -xmax else 0
 cx0 <- .435; cx1 <- .705; nx <- function(v) cx0+(v-xmin)/(xmax-xmin)*(cx1-cx0)
 txt(p,.012,.944,12,TRUE); txt(heading,.052,.946,9,TRUE)
 if(fig=="Figure_1" & p=="B")txt("Potentially overlapping routes",.052,.883,8)
 if(fig=="Figure_2" & p=="A")txt(paste0("Common comparison population n=",unique(r$n)),.052,.872,8)
 if(fig=="Figure_2" & p=="B")txt(paste0("Measurement-positive domain; overall comparison n=",r$comparison_n),.052,.827,8)
 if(fig=="Figure_3" & p=="A")txt(paste0("Fasting population n=",unique(r$n),"; ",unique(r$weight)),.052,.858,8)
 if(fig=="Supplementary_Figure_S2")txt("Scenario populations differ; D2 assumes Borderline possible history",.052,.872,8)
 if(fig=="Supplementary_Figure_S3" & p=="A")txt("Assumed age85 is illustrative",.052,.842,8)
 if(fig=="Supplementary_Figure_S3" & p=="B")txt(paste0("G2: HbA1c <6.5% and fasting glucose <126 mg/dL; n=",unique(r$n)),.052,.858,8)
 top <- if(nrow(r)==1) .54 else if(nrow(r)==2) .62 else .745
 bot <- if(nrow(r)==1) top else .285
 ys <- seq(top,bot,length.out=nrow(r))
 axisy <- .17
 txt(if(contrast)"Difference, pp (95% CI)" else "% (95% CI)",.74,if(nrow(r)==1).72 else .827,8,TRUE)
 grid.lines(c(cx0,cx1),c(axisy,axisy),gp=gpar(lwd=.65))
 ticks <- if(contrast)c(-xmax,0,xmax) else seq(0,xmax,by=if(xmax==70)20 else 20)
 for(t in ticks){
  grid.lines(rep(nx(t),2),c(axisy,axisy-.022),gp=gpar(lwd=.65));txt(as.character(t),nx(t),axisy-.063,8,just="centre")
 }
 if(contrast)grid.lines(rep(nx(0),2),c(axisy,top+.06),gp=gpar(lty=2,col="black",lwd=.65))
 for(i in seq_len(nrow(r))){
  z <- r[i,];y <- ys[i];col <- pal[z$color]
  label <- z$label
  if(fig=="Figure_2" & p=="B")label <- "Measurement-only"
  txt(label,.052,y+.026,8)
  ntext <- paste0("n=",z$n)
  if(fig=="Figure_2" & p=="A")ntext <- paste0("Category n=",z$positive_n,"; denominator ",z$n)
  if(fig=="Figure_2" & p=="B")ntext <- paste0("Numerator ",z$positive_n,"; domain n=",z$n)
  if(fig=="Supplementary_Figure_S2")ntext <- paste0(z$positive_n,"/",z$n,"; comparison n=",z$comparison_n)
  txt(ntext,.052,y-.045,8)
  est <- as.numeric(z$estimate);lo <- as.numeric(z$lower);hi <- as.numeric(z$upper)
  stopifnot(lo>=xmin,hi<=xmax)
  grid.lines(c(nx(lo),nx(hi)),rep(y,2),gp=gpar(col=col,lwd=1.15))
  for(v in c(lo,hi))grid.lines(rep(nx(v),2),c(y-.017,y+.017),gp=gpar(col=col,lwd=.8))
  marker(nx(est),y,as.integer(z$shape),col)
  disp <- paste0(if(contrast && est>0)"+" else "",z$displayed_estimate," (",gsub(" to ","-",z$displayed_ci,fixed=TRUE),")")
  txt(disp,.74,y,8)
 }
}
flow <- function(r,p){
 txt(p,.012,.968,12,TRUE)
 counts <- r[r$plotted_field=="nonmissing_n",]
 chain <- counts[1:7,]
 yy <- seq(.92,.375,length.out=7)
 for(i in seq_len(7)){
  y <- yy[i];grid.rect(.41,y,width=.63,height=.065,gp=gpar(fill="white",col="black",lwd=.75))
  txt(chain$label[i],.41,y+.013,9,just="centre");txt(paste0("n=",chain$n[i]),.41,y-.017,9,TRUE,just="centre")
  if(i>1){
   grid.lines(c(.41,.41),c(yy[i-1]-.034,y+.034),arrow=arrow(length=unit(1.5,"mm"),type="closed"),gp=gpar(lwd=.7))
   txt(paste0("Excluded ",chain$missing_n[i]),.76,(yy[i-1]+y)/2,8)
  }
 }
 y <- .183;xs <- c(.265,.755)
 branch <- list(counts[8,],r[r$plotted_field=="n",])
 for(i in seq_len(2)){
  if(i==1)grid.lines(c(.41,.41,xs[i],xs[i]),c(.34,.285,.285,y+.052),arrow=arrow(length=unit(1.5,"mm"),type="closed"),gp=gpar(lwd=.7))
  else txt(paste0("Separate fasting domain; ",branch[[i]]$weight),xs[i],.28,8,just="centre")
  grid.rect(xs[i],y,width=.43,height=.10,gp=gpar(fill="white",col="black",lwd=.75))
  lab <- if(i==1)"D0 comparison\nstatus ascertainable" else "Fasting sensitivity\nHbA1c and glucose observed"
  txt(lab,xs[i],y+.016,9,just="centre");txt(paste0("n=",branch[[i]]$n),xs[i],y-.031,9,TRUE,just="centre")
 }
 txt(paste0("D0 unknown: ",counts$missing_n[8]),.12,.307,8,just="centre")
 txt("Analysis domains may overlap; counts are not additive",.5,.069,9,just="centre")
}
export <- function(path,width,height,draw){
 cairo_pdf(paste0(path,".pdf"),width=width/25.4,height=height/25.4,family=font,onefile=TRUE,bg="white")
 grid.newpage();draw();dev.off()
 png(paste0(path,".png"),width=width,height=height,units="mm",res=300,type="cairo",bg="white")
 grid.newpage();draw();dev.off()
 jpeg(paste0(path,".jpg"),width=width,height=height,units="mm",res=300,quality=100,type="cairo",bg="white")
 grid.newpage();draw();dev.off()
}
for(fig in names(spec)){
 s <- spec[[fig]];ps <- unique(d$panel[d$figure==fig]);h <- sum(s$h)
 folder <- if(grepl("Supplementary",fig))"supplementary" else "main"
 draw <- function(){
  ytop <- 1
  for(i in seq_along(ps)){
   hf <- s$h[i]/h;pushViewport(viewport(x=.5,y=ytop-hf/2,width=1,height=hf))
   panel(fig,ps[i],s$x[i],s$heads[i]);popViewport();ytop <- ytop-hf
  }
 }
 export(file.path(out,"figures",folder,fig),170,h,draw)
 for(i in seq_along(ps)){
  pp <- ps[i];jj <- i
  export(file.path(out,"figures",folder,"panels",paste0(fig,pp)),170,s$h[i],function()panel(fig,pp,s$x[jj],s$heads[jj]))
  # Panel PNGs retained as additional review artifacts.
 }
 cat(fig,"panels",length(ps),"height",h,"mm\n")
}
sink(file.path(out,"logs/R_sessionInfo.txt"));print(sessionInfo());print(capabilities());print(installed.packages()[,c("Package","Version")]);sink()
