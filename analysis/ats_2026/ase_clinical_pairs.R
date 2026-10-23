b<-Sys.getenv('CLIF_ANALYSIS_ROOT');if(!nzchar(b))stop('Set CLIF_ANALYSIS_ROOT to the aggregate analysis workspace');out<-file.path(b,'ASE_RUSH_updated_20261022/clinical_figure');d0<-read.csv(file.path(out,'selected_monthly_counts.csv'));keys<-read.csv(file.path(out,'selected_common_pairs.csv'));groups<-c('No presumed infection','Presumed infection without ASE','ASE');margins<-list();audit<-list();effects<-list()
for(i in seq_len(nrow(keys))){
 k<-keys[i,];d<-subset(d0,organism==k$organism & antimicrobial==k$antimicrobial);site_groups<-table(d$site,factor(d$ase_group,levels=groups));complete_sites<-rownames(site_groups)[rowSums(site_groups>0)==3];d<-subset(d,site%in%complete_sites);raw<-d;d$Y<-d$n_non_susceptible;d$T<-d$n_interpretable;status<-'model_standardized';reason<-''
 repeat{n0<-nrow(d);ss<-aggregate(cbind(Y,T)~site,d,sum);ss<-subset(ss,Y>0 & T>Y)$site;d<-subset(d,site%in%ss);if(!nrow(d))break;mm<-aggregate(cbind(Y,T)~calendar_month,d,sum);mm<-subset(mm,Y>0 & T>Y)$calendar_month;d<-subset(d,calendar_month%in%mm);if(nrow(d)==n0||!nrow(d))break}
 if(nrow(d)==0 || length(unique(d$site))<5 || length(unique(d$ase_group))<3 || length(unique(d$calendar_month))<18 || sum(d$Y)<30 || sum(d$T-d$Y)<30){status<-'pooled_observed';reason<-'Insufficient within-site/month outcome variation or support for adjusted model'}
 fit<-NULL;V<-NULL;grid<-NULL
 if(status=='model_standardized'){
  d$site<-factor(d$site);d$month<-factor(d$calendar_month);d$ase_group<-factor(d$ase_group,levels=groups);S<-nlevels(d$site)
  warnings<-character();fit<-tryCatch(withCallingHandlers(glm(cbind(Y,T-Y)~site+month+ase_group,data=d,family=quasibinomial(),control=glm.control(maxit=100)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)NULL)
  if(is.null(fit)||!fit$converged||anyNA(coef(fit))||length(warnings)){status<-'pooled_observed';reason<-'Adjusted model failed or warned'}else{
   X<-model.matrix(fit);mu<-fitted(fit);B<-tryCatch(chol2inv(chol(crossprod(X,X*as.numeric(d$T*mu*(1-mu))))),error=function(e)NULL)
   if(is.null(B)){status<-'pooled_observed';reason<-'Singular adjusted information matrix'}else{
    U<-rowsum(X*as.numeric(d$Y-d$T*mu),d$site);V<-B%*%crossprod(U)%*%B*S/(S-1)*(nrow(X)-1)/(nrow(X)-ncol(X));crit<-qt(.975,S-1)
    sm<-paste(d$site,d$month,sep='|');tt<-table(sm);common<-names(tt[tt==3]);grid<-d[!duplicated(sm)&sm%in%common,]
    if(length(unique(grid$site))<5){status<-'pooled_observed';reason<-'Fewer than five sites with shared monthly support after model screening'}
   }
  }
 }
 if(status=='model_standardized'){
  weights<-1/(length(unique(grid$site))*as.numeric(table(grid$site)[as.character(grid$site)]));stopifnot(abs(sum(weights)-1)<1e-8)
  for(g in groups){nd<-grid;nd$ase_group<-factor(g,levels=groups);xp<-model.matrix(delete.response(terms(fit)),nd,contrasts.arg=fit$contrasts);prob<-plogis(as.numeric(xp%*%coef(fit)));m<-sum(weights*prob);grad<-colSums(xp*as.numeric(weights*prob*(1-prob)));se<-sqrt(as.numeric(t(grad)%*%V%*%grad))/(m*(1-m));margins[[length(margins)+1]]<-data.frame(organism=k$organism,antimicrobial=k$antimicrobial,ase_group=g,percent=100*m,lower=100*plogis(qlogis(m)-crit*se),upper=100*plogis(qlogis(m)+crit*se),estimate_type=status,sites=length(unique(grid$site)),rank=k$rank_among_comparable_pairs,n_interpretable=sum(subset(raw,ase_group==g)$n_interpretable))}
  cv<-rep(0,length(coef(fit)));names(cv)<-names(coef(fit));cv['ase_groupASE']<-1;cv['ase_groupPresumed infection without ASE']<--1;be<-sum(cv*coef(fit));se<-sqrt(as.numeric(t(cv)%*%V%*%cv));effects[[length(effects)+1]]<-data.frame(organism=k$organism,antimicrobial=k$antimicrobial,estimate=exp(be),lower=exp(be-crit*se),upper=exp(be+crit*se),p=2*pt(-abs(be/se),S-1))
 }else{
  for(g in groups){rr<-subset(raw,ase_group==g);n<-sum(rr$n_interpretable);ns<-sum(rr$n_non_susceptible);stopifnot(n>0);margins[[length(margins)+1]]<-data.frame(organism=k$organism,antimicrobial=k$antimicrobial,ase_group=g,percent=100*ns/n,lower=NA_real_,upper=NA_real_,estimate_type=status,sites=length(unique(rr$site)),rank=k$rank_among_comparable_pairs,n_interpretable=n)}
 }
 audit[[i]]<-data.frame(organism=k$organism,antimicrobial=k$antimicrobial,rank=k$rank_among_comparable_pairs,n_interpretable=k$interpretable_tests,status=status,reason=reason)
}
r<-do.call(rbind,margins);a<-do.call(rbind,audit);write.csv(r,file.path(out,'common_pair_percentages.csv'),row.names=FALSE);write.csv(a,file.path(out,'common_pair_model_audit.csv'),row.names=FALSE)
e<-do.call(rbind,effects);e$q_selected_screen<-p.adjust(e$p,'BH');write.csv(e,file.path(out,'common_pair_exploratory_contrasts.csv'),row.names=FALSE)
stopifnot(nrow(r)==3*nrow(keys),all(is.finite(r$percent)),all(r$percent>=0 & r$percent<=100));print(a,row.names=FALSE)
