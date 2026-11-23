source("R/00_setup.R")
con<-dbConnect(duckdb::duckdb(),file.path(WORK_DIR,"extract.duckdb"),read_only=TRUE)
Q<-function(s) as.data.table(dbGetQuery(con,s)); H<-function(a,b) as.numeric(difftime(a,b,units="hours")); fl<-list()
# ================= MIMIC =================
m<-Q("select * from st"); fl$m0<-nrow(m); dx<-Q("select hadm_id,icd_code from dx"); fx<-function(p) as.integer(m$hadm_id%in%dx[grepl(p,icd_code)]$hadm_id)
m[,sep:=fx("^(99591|99592|78552|A40|A41|R652)")]; m<-m[sep==1]; fl$m_sep<-nrow(m)
m[,los:=H(outtime,intime)]; m<-m[los>=24]; fl$m_los<-nrow(m)
VP<-Q("select * from inp"); VP<-merge(VP[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); VP[,`:=`(hs=H(starttime,intime),he=H(endtime,intime))]; m<-m[stay_id%in%VP[hs<=24&he>=-6]$stay_id]; fl$m_vaso<-nrow(m)
L<-Q("select stay_id,specimen_id,itemid,charttime,v,value from lab"); L<-L[stay_id%in%m$stay_id]; L<-merge(L,m[,.(stay_id,intime)],by="stay_id"); L[,h:=H(charttime,intime)]
sp<-L[itemid==52033,.(specimen_id,st=value)]; po<-L[itemid==50821,unique(specimen_id)]
g<-dcast(L[itemid%in%c(50818,50820,50821,50816,50803)&h>=-6&h<=24],stay_id+specimen_id+h~itemid,value.var="v",fun.aggregate=function(x) x[1]); setnames(g,c("50818","50820","50821","50816","50803"),c("paco2","ph","pao2","fio2l","hco3a"),skip_absent=TRUE)
g[,st:=sp$st[match(specimen_id,sp$specimen_id)]]; g<-g[(st=="ART."|(is.na(st)&specimen_id%in%po))&!is.na(paco2)&!is.na(ph)&paco2>=10&paco2<=150&ph>=6.6&ph<=7.8][order(h)][,.SD[1],by=stay_id]
m<-merge(m,g[,.(stay_id,habg=h,paco2,ph,pao2,fio2l,hco3a)],by="stay_id"); fl$m_abg<-nrow(m)
ref<-m$habg; names(ref)<-m$stay_id
gv<-function(it,lo,hi,fun="near",a=-24,z=6){x<-L[itemid%in%it&v>=lo&v<=hi]; x[,hh:=h-ref[as.character(stay_id)]]; x<-x[hh>=a&hh<=z]
 x<-if(fun=="near") x[order(abs(hh))][,.(v=v[1]),by=stay_id] else x[,.(v=if(fun=="max") max(v) else min(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(lact=gv(50813,.1,30),alb=gv(50862,1,6),bun=gv(51006,1,250),creat=gv(50912,.1,25),na=gv(50983,100,180),cl=gv(50902,60,140),hco3s=gv(50882,2,60),hb=gv(51222,3,20),wbc=gv(51301,.1,200),plt=gv(51265,1,1500),glu=gv(50931,20,1500),k=gv(50971,1.5,9),bili=gv(50885,.1,60,"max",-24,24),pltmin=gv(51265,1,1500,"min",-24,24),crmax=gv(50912,.1,25,"max",-24,24))]
m[,hco3:=fifelse(!is.na(hco3a)&hco3a>2&hco3a<60,hco3a,hco3s)]
ce<-Q("select * from ce"); ce<-ce[stay_id%in%m$stay_id]; ce<-merge(ce,m[,.(stay_id,intime)],by="stay_id"); ce[,h:=H(charttime,intime)]
cm<-function(its,lo,hi,fun=mean){x<-ce[itemid%in%its&v>=lo&v<=hi&h>=0&h<=24,.(v=fun(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(hr=cm(220045,20,250),rr=cm(c(220210,224690),3,70),map=cm(c(220052,220181,225312),20,200),mapmin=cm(c(220052,220181,225312),20,200,min),spo2=cm(220277,50,100),fio2c=cm(223835,21,100,max))]
w8<-ce[itemid%in%c(226512,224639)&v>=25&v<=300][order(h)][,.SD[1],by=stay_id]; m[,wt:=w8$v[match(stay_id,w8$stay_id)]]
gc<-ce[itemid%in%c(220739,223900,223901)&h>=0&h<=24,.(s=sum(v),n=.N),by=.(stay_id,charttime)][n==3,.(gcs=min(s)),by=stay_id]; m[,gcs:=gc$gcs[match(stay_id,gc$stay_id)]]
m[,vaso:=1L]
m[,fio2:=fifelse(!is.na(fio2l),fio2l,fifelse(!is.na(fio2c),fio2c,NA_real_))]; m[,pf:=pao2/(fio2/100)]
sc<-function(x,br) fifelse(is.na(x),0L,as.integer(findInterval(x,br)))
m[,sev:=as.numeric(fifelse(is.na(pltmin),0L,4L-sc(pltmin,c(20,50,100,150)))+sc(bili,c(1.2,2,6,12))+fifelse(vaso==1,3L,fifelse(!is.na(mapmin)&mapmin<70,1L,0L))+sc(crmax,c(1.2,2,3.5,5)))]
pv<-Q("select * from pv"); pv<-merge(pv,m[,.(stay_id,intime,habg)],by="stay_id"); pv[,`:=`(hs=H(starttime,intime)-habg,he=H(endtime,intime)-habg)]
m[,imv:=as.integer(stay_id%in%pv[itemid==225792&hs<=2&he>=-1]$stay_id)]; m[,imvpre:=as.integer(stay_id%in%pv[itemid==225792&hs<0&he>=-1]$stay_id)]
iv<-pv[itemid==225792&hs>2,.(t=min(hs)),by=stay_id]; m[,intub_h:=iv$t[match(stay_id,iv$stay_id)]]
m[,niv:=as.integer(stay_id%in%pv[itemid==225794&hs<=24]$stay_id)]
m[,t0:=pmax(habg,0)]; ve<-VP[,.(ve=max(he)),by=stay_id]; m[,vend:=ve$ve[match(stay_id,ve$stay_id)]-t0]; nv<-VP[hs<=24&he>=-6,.(nv=uniqueN(itemid)),by=stay_id]; m[,nvaso:=nv$nv[match(stay_id,nv$stay_id)]]
dth<-as.POSIXct(m$deathtime,tz="UTC"); dd<-as.POSIXct(ifelse(is.na(m$dod),NA,paste(as.character(m$dod),"23:59:00")),tz="UTC"); dth[is.na(dth)]<-dd[is.na(dth)]
m[,death_h:=pmax(H(dth,intime)-t0,0)]; m[,fu_h:=pmin(fifelse(is.na(death_h),365*24,death_h),365*24)]
m[,`:=`(hd=hospital_expire_flag,icud=as.integer(!is.na(death_h)&death_h<=los-t0+1),los=los-t0)]
m[,`:=`(copd=fx("^(49[0-2]|496|J4[0-4])"),clung=fx("^(49[0-6]|J4[0-7]|515|516|J84)"),osa=fx("^(32723|78051|78053|78057|G4733|27803|E662)"),arrest=fx("^(4275|I46)"),af=fx("^(42731|42732|I48)"),cad=fx("^(41[0-4]|I2[0-5])"),ckd=fx("^(585|N18)"),dm=fx("^(250|E1[0-4])"),liver=fx("^(571|K70|K7[2-6])"),cancer=fx("^(1[4-9][0-9]|20[0-8]|C[0-9])"),sepsis=fx("^(99591|99592|78552|A40|A41|R652)"),pna=fx("^(48[0-6]|J1[2-8])"),obese=fx("^(2780|E66)"),male=as.integer(gender=="M"),db="MIMIC-IV",id=stay_id)]
cols<-c("id","db","age","male","wt","habg","paco2","ph","pao2","fio2","pf","hco3","lact","alb","bun","creat","na","cl","hb","wbc","plt","glu","k","hr","rr","map","gcs","vaso","sev","imv","imvpre","intub_h","niv","death_h","fu_h","hd","icud","los","copd","clung","osa","arrest","af","cad","ckd","dm","liver","cancer","sepsis","pna","obese","vend","nvaso")
D<-m[,..cols]; D[,q:=factor(fifelse(paco2<35,"lo",fifelse(paco2<=45,"n",fifelse(ph>=7.35,"chc","hca"))),levels=c("n","lo","chc","hca"))]
saveRDS(D,file.path(WORK_DIR,"analytic.rds")); jsonlite::write_json(fl,file.path(OUT_DIR,"flow.json"),auto_unbox=TRUE); print(unlist(fl))
print(D[,.(n=.N,d28=mean(!is.na(death_h)&death_h<=672)),keyby=.(db,copd,q)])
print(D[,lapply(.SD,function(x) round(mean(is.na(x)),2)),by=db,.SDcols=c("pf","hco3","lact","alb","bun","hr","map","sev")])
