# Statistical analysis: hypercapnic acidemia and 28-day mortality in septic shock.
source("R/00_setup.R")
RAW<-as.data.table(readRDS(file.path(WORK_DIR,"analytic.rds"))); A<-copy(RAW); set.seed(2026)
A[,`:=`(ev=as.integer(!is.na(death_h)&death_h<=672),tt=pmin(fu_h,672)/24,lpf=log(pmin(pmax(pf,30),700)),lac2=as.integer(!is.na(lact)&lact>2),acid=as.integer(ph<7.35))]
A[,sevz:=(sev-mean(sev,na.rm=TRUE))/sd(sev,na.rm=TRUE),by=db]; A[is.na(nvaso),nvaso:=1]
num<-c("lact","bun","creat","na","hb","wbc","plt","hr","map","sevz","lpf","hco3","alb")
for(v in num) A[,(v):=fifelse(is.na(get(v)),median(get(v),na.rm=TRUE),get(v)),by=db]
A[,lac4:=as.integer(lact>=4)]; A[,vfd:=fifelse(ev==1,0,pmax(28-pmax(fifelse(is.na(vend),0,vend),0)/24,0))]
M2<-"age+male+cad+af+ckd+dm+liver+cancer+pna+osa"
M3<-paste(M2,"+sevz+imv+niv+nvaso+lpf+lact+bun+creat+na+hb+wbc+plt+hr+map")
ff<-function(x,d) if(length(unique(d$db))>1) paste(x,"+db") else x
lg<-function(f,d,w=NULL) glm(as.formula(f),binomial,d,weights=w)
ex<-function(fit,tr){s<-summary(fit)$coef; r<-s[rownames(s)%in%tr,,drop=FALSE]; b<-r[,1]; se<-r[,2]; cbind(exp(b),exp(b-1.96*se),exp(b+1.96*se),2*pnorm(-abs(b/se)))}
exc<-function(fit,tr){s<-summary(fit)$coef; r<-s[rownames(s)%in%tr,,drop=FALSE]; b<-r[,1]; se<-if("robust se"%in%colnames(r)) r[,"robust se"] else r[,"se(coef)"]; cbind(exp(b),exp(b-1.96*se),exp(b+1.96*se),r[,"Pr(>|z|)"])}
cx<-function(f,d,w=NULL,y="Surv(tt,ev)"){f<-as.formula(paste(y,"~",f,if(length(unique(d$db))>1) "+strata(db)" else "")); coxph(f,d,weights=w,robust=!is.null(w))}
QT<-c("qlo","qchc","qhca")
# ---- early-stop gate ----
g<-ex(lg(paste("ev~q+",M3),A[copd==0]),"qhca"); cat("GATE",round(g,4),"\n"); if(g[4]>=0.05){stop("Primary association not significant in MIMIC-IV; analysis stopped.")}
an<-function(d,full){o<-list(); o$n<-nrow(d); o$t1<-d[,.(n=.N,ev=sum(ev)),keyby=q]; o$rate<-c(d28=mean(d$ev),hd=mean(d$hd),icud=mean(d$icud))
 o$desc<-d[,.(age=mean(age),male=mean(male),paco2=median(paco2),ph=mean(ph),hco3=mean(hco3),lact=median(lact),pf=median(exp(lpf)),sev=mean(sev,na.rm=T),nv2=mean(nvaso>=2),imv=mean(imv),pna=mean(pna),d28=mean(ev)),keyby=q]
 o$crude<-ex(lg("ev~q",d),QT); o$m2<-ex(lg(ff(paste("ev~q+",M2),d),d),QT); o$m3<-ex(lg(ff(paste("ev~q+",M3),d),d),QT); o$hr<-exc(cx(paste("q+",M3),d),QT)
 o$km<-1-pchisq(survdiff(Surv(tt,ev)~q,d)$chisq,3)
 h<-d[paco2>45]; h[,`:=`(dph=(7.40-ph)/0.05,hc5=(paco2-45)/5)]; o$phc<-ex(lg(ff(paste("ev~dph+hc5+",M3),h),h),c("dph","hc5"))
 hv<-d[q%in%c("chc","hca")]; hv[,x:=as.integer(q=="hca")]; o$hvc<-ex(lg(ff(paste("ev~x+",M3),hv),hv),"x")
 d[,`:=`(hc5=pmax(paco2-45,0)/5,lo5=pmax(35-paco2,0)/5)]; fi<-lg(ff(paste("ev~lo5+hc5*acid+",M3),d),d); o$int<-summary(fi)$coef["hc5:acid",4]
 o$hcA<-ex(lg(ff(paste("ev~lo5+hc5+",M3),d[acid==1]),d[acid==1]),"hc5"); o$hcN<-ex(lg(ff(paste("ev~lo5+hc5+",M3),d[acid==0]),d[acid==0]),"hc5")
 f0<-lg(ff(paste("ev~",M3),d),d); fl<-lg(ff(paste("ev~paco2+",M3),d),d); fr<-lg(ff(paste("ev~ns(paco2,4)+",M3),d),d); o$rcs<-c(anova(f0,fr,test="LRT")[2,5],anova(fl,fr,test="LRT")[2,5])
 a0<-roc(d$ev,fitted(f0),quiet=T); f1<-lg(ff(paste("ev~q+",M3),d),d); a1<-roc(d$ev,fitted(f1),quiet=T); o$inc<-c(anova(f0,f1,test="LRT")[2,5],auc(a0),auc(a1),roc.test(a0,a1)$p.value)
 o$s_hd<-ex(lg(ff(paste("hd~q+",M3),d),d),QT); o$s_icu<-ex(lg(ff(paste("icud~q+",M3),d),d),QT)
 if(all(d$db=="MIMIC-IV")) for(hz in c(90,365)){d$e2<-as.integer(!is.na(d$death_h)&d$death_h<=hz*24); d$t2<-pmin(d$fu_h,hz*24)/24; o[[paste0("s_",hz)]]<-exc(coxph(as.formula(paste("Surv(t2,e2)~q+",M3)),d),QT)}
 s<-summary(lm(as.formula(ff(paste("vfd~q+",M3),d)),d))$coef; o$s_vfd<-cbind(s[QT,1],s[QT,1]-1.96*s[QT,2],s[QT,1]+1.96*s[QT,2],s[QT,4]); o$rate["vfd"]<-median(d$vfd)
 n0<-d[imv==0]; n0[,st:=fifelse(!is.na(intub_h)&intub_h<=168,1L,fifelse(!is.na(death_h)&death_h<=168,2L,0L))]; n0[,ti:=pmin(fifelse(st==1,intub_h,fifelse(st==2,death_h,pmin(fu_h,168))),168)]
 o$rate["int7"]<-mean(n0$st==1); X<-model.matrix(as.formula(ff(paste("~q+",sub("\\+imv","",M3)),n0)),n0)[,-1]
 cr<-tryCatch(crr(n0$ti,n0$st,X),error=function(e) NULL); if(!is.null(cr)){s<-summary(cr)$coef[QT,]; o$s_int<-cbind(s[,2],exp(log(s[,2])-1.96*s[,3]),exp(log(s[,2])+1.96*s[,3]),s[,5])}
 o$s_intq<-1-pchisq(cuminc(n0$ti,n0$st,n0$q)$Tests[1,1],3)
 sv<-d[ev==0&los>0]; s<-summary(lm(as.formula(ff(paste("log(los)~q+",M3),sv)),sv))$coef; o$s_los<-cbind(exp(s[QT,1]),exp(s[QT,1]-1.96*s[QT,2]),exp(s[QT,1]+1.96*s[QT,2]),s[QT,4])
 f<-copy(full); o$em<-lapply(c(copd="copd",imv="imv",pna="pna",lac4="lac4"),function(v){fi<-lg(ff(paste0("ev~q*",v,"+",M3),f),f); p<-summary(fi)$coef[paste0("qhca:",v),4]
  r<-sapply(0:1,function(l) ex(lg(ff(paste("ev~q+",M3),f[get(v)==l]),f[get(v)==l]),"qhca")); list(r0=r[,1],r1=r[,2],ror=ex(fi,paste0("qhca:",v)),p=p,n=c(sum(f[[v]]==0),sum(f[[v]]==1)))})
 b<-d[q%in%c("n","hca")]; b[,hi:=as.integer(q=="hca")]; cv<-strsplit(gsub(" ","",ff(M3,b)),"\\+")[[1]]; fps<-as.formula(paste("hi~",paste(cv,collapse="+")))
 w1<-weightit(fps,b,method="glm",estimand="ATE"); b$ipw<-pmin(w1$weights,quantile(w1$weights,.99)); b$oww<-weightit(fps,b,method="glm",estimand="ATO")$weights
 qb<-function(wv,f="ev~hi"){b$wq<-wv; g<-suppressWarnings(glm(as.formula(f),quasibinomial,b,weights=wq)); s<-summary(g)$coef["hi",]; c(exp(s[1]),exp(s[1]-1.96*s[2]),exp(s[1]+1.96*s[2]),s[4])}
 o$hi_m3<-ex(lg(ff(paste("ev~hi+",M3),b),b),"hi"); o$iptw<-qb(b$ipw); o$ow<-qb(b$oww); o$dr<-qb(b$ipw,ff(paste("ev~hi+",M3),b)); o$iptwhr<-exc(cx("hi",b,b$ipw),"hi")
 mt<-matchit(fps,b,method="nearest",caliper=.2,std.caliper=TRUE); md<-match.data(mt); s<-summary(glm(ev~hi,binomial,md))$coef["hi",]; o$psm<-c(exp(s[1]),exp(s[1]-1.96*s[2]),exp(s[1]+1.96*s[2]),s[4],sum(md$hi))
 o$smd<-c(max(abs(summary(mt)$sum.all[-1,3])),max(abs(summary(mt)$sum.matched[-1,3])));
 wsmd<-function(w){X<-model.matrix(as.formula(paste("~",paste(cv,collapse="+"))),b)[,-1,drop=FALSE]; max(apply(X,2,function(x){m1<-weighted.mean(x[b$hi==1],w[b$hi==1]); m0<-weighted.mean(x[b$hi==0],w[b$hi==0]); s<-sqrt((var(x[b$hi==1])+var(x[b$hi==0]))/2); if(s>0) abs(m1-m0)/s else 0}))}; o$wsmd<-c(raw=wsmd(rep(1,nrow(b))),iptw=wsmd(b$ipw),ow=wsmd(b$oww)); o$ev<-evalues.OR(o$hi_m3[1],o$hi_m3[2],o$hi_m3[3],rare=mean(b$ev)<.15)[2,1]
 sen<-function(dd,f=M3,tr="qhca"){if(nrow(dd)<300||sum(dd$ev)<40||sum(dd$q=="hca")<30) return(NULL); c(ex(lg(ff(paste("ev~q+",f),dd),dd),tr),nrow(dd))}
 rq<-function(dd,cut=45,phc=7.35){dd<-copy(dd); dd[,q:=factor(fifelse(paco2<35,"lo",fifelse(paco2<=cut,"n",fifelse(ph>=phc,"chc","hca"))),levels=c("n","lo","chc","hca"))]; dd}
 S<-list(); S[["Excluding cardiac arrest"]]<-sen(d[arrest==0]); S[["No invasive ventilation before the blood gas"]]<-sen(d[imvpre==0]); S[["Excluding any chronic lung disease"]]<-sen(d[clung==0])
 S[["Excluding sleep apnea and obesity"]]<-sen(d[osa==0&obese==0]); S[["Sepsis-3 septic shock (lactate >2 mmol/L)"]]<-sen(d[lac2==1]); S[["Blood gas within 0-6 h"]]<-sen(d[habg<=6])
 S[["Excluding death within 24 h"]]<-sen(d[is.na(death_h)|death_h>24]); S[["Additional adjustment for bicarbonate"]]<-sen(d,paste(M3,"+hco3")); S[["Additional adjustment for respiratory rate and GCS"]]<-sen(copy(d)[,rr:=fifelse(is.na(rr),median(rr,na.rm=T),rr)][,gcs:=fifelse(is.na(gcs),median(gcs,na.rm=T),gcs)],paste(M3,"+rr+gcs"))
 S[["Acidemia defined as pH <7.30"]]<-sen(rq(d,45,7.30)); S[["Hypercapnia defined as PaCO2 >50 mmHg"]]<-sen(rq(d,50,7.35)); S[["Including COPD (adjusted for COPD)"]]<-sen(full,paste(M3,"+copd"))
 S[["Cox model (HR)"]]<-c(o$hr[3,],nrow(d))
 im<-copy(d); R2<-RAW[match(paste(im$db,im$id),paste(db,id))]; for(v in c("lact","alb","bun","hr","map","hco3")) set(im,j=v,value=R2[[v]]); set(im,j="lpf",value=log(pmin(pmax(R2$pf,30),700)))
 mi<-mice(im[,c("ev","q",strsplit(gsub(" ","",M3),"\\+")[[1]],"db"),with=FALSE],m=5,printFlag=FALSE,seed=1); pz<-summary(pool(with(mi,glm(as.formula(ff(paste("ev~q+",M3),d)),family=binomial)))); r<-pz[pz$term=="qhca",]
 S[["Multiple imputation (m = 5)"]]<-c(exp(r$estimate),exp(r$estimate-1.96*r$std.error),exp(r$estimate+1.96*r$std.error),r$p.value,nrow(d)); o$S<-S
 sgv<-list(`Age ≥65 years`=d$age>=65,`Male sex`=d$male==1,`Invasive ventilation`=d$imv==1,`Pneumonia`=d$pna==1,`Lactate ≥4 mmol/L`=d$lac4==1,`Two or more vasopressors`=d$nvaso>=2,`Chronic kidney disease`=d$ckd==1,`Malignancy`=d$cancer==1,`Chronic liver disease`=d$liver==1)
 o$sg<-lapply(sgv,function(z){d$z<-as.integer(z); pint<-tryCatch(summary(lg(ff(paste("ev~q*z+",M3),d),d))$coef["qhca:z",4],error=function(e) NA)
  g<-function(x) tryCatch(ex(lg(ff(paste("ev~q+",M3),x),x),"qhca"),error=function(e) matrix(NA,1,4)); list(yes=g(d[z==1]),no=g(d[z==0]),nyes=sum(d$z==1),nno=sum(d$z==0),pint=pint)})
 o}
R<-list(); R[["MIMIC-IV"]]<-an(A[copd==0],A)
saveRDS(R,file.path(OUT_DIR,"results.rds")); saveRDS(A,file.path(WORK_DIR,"analysis.rds"))
for(x in names(R)){o<-R[[x]]; cat("\n==",x,o$n,"\n"); print(round(o$m3,3)); print(round(o$hr,3)); print(round(o$phc,3)); print(round(o$hvc,3)); cat("int",o$int,"\n"); print(round(rbind(o$hcA,o$hcN),3)); print(round(o$s_vfd,3)); print(sapply(o$em,function(z) round(c(z$r0[1],z$r1[1],z$p),3))); print(round(rbind(o$iptw,o$ow,o$dr,o$psm[1:4]),3)); print(round(do.call(rbind,o$S),3))}
