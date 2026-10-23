b<-Sys.getenv('CLIF_ANALYSIS_ROOT');if(!nzchar(b))stop('Set CLIF_ANALYSIS_ROOT to the aggregate analysis workspace');src<-file.path(b,'ASE_RUSH_updated_20261022/ase');out<-src
fits<-readRDS(file.path(src,'phenotype_models.rds'));groups<-c('No presumed infection','Presumed infection without ASE','ASE');r<-list()
for(key in names(fits)[grepl('^organism_representation\\|',names(fits))]){
 fit<-fits[[key]];org<-strsplit(key,'|',fixed=TRUE)[[1]][2];d<-model.frame(fit);Y<-model.response(d)[,1];T<-rowSums(model.response(d));X<-model.matrix(fit);mu<-fitted(fit);S<-length(unique(d$site));w<-T*mu*(1-mu)
 B<-chol2inv(chol(crossprod(X,X*w)));U<-rowsum(X*as.numeric(Y-T*mu),as.character(d$site));V<-B%*%crossprod(U)%*%B*S/(S-1)*(nrow(X)-1)/(nrow(X)-ncol(X));crit<-qt(.975,S-1)
 sm<-paste(d$site,d$month,sep='|');counts<-table(sm);common<-names(counts[counts==3]);grid<-d[!duplicated(sm) & sm%in%common,];ng<-length(unique(grid$site))
 if(ng<5)next
 ww<-1/(ng*as.numeric(table(grid$site)[as.character(grid$site)]));stopifnot(abs(sum(ww)-1)<1e-8)
 for(g in groups){
  nd<-grid;nd$ase_group<-factor(g,levels=groups);xp<-model.matrix(delete.response(terms(fit)),nd,contrasts.arg=fit$contrasts);prob<-plogis(as.numeric(xp%*%coef(fit)));m<-sum(ww*prob);grad<-colSums(xp*as.numeric(ww*prob*(1-prob)));se<-sqrt(as.numeric(t(grad)%*%V%*%grad))/(m*(1-m))
  r[[length(r)+1]]<-data.frame(organism=org,ase_group=g,standardized_percent=100*m,lower=100*plogis(qlogis(m)-crit*se),upper=100*plogis(qlogis(m)+crit*se),model_sites=S,shared_support_sites=ng,shared_site_months=nrow(grid))
 }
}
r<-do.call(rbind,r);stopifnot(all(is.finite(r$standardized_percent)),all(r$lower<=r$standardized_percent),all(r$upper>=r$standardized_percent));write.csv(r,file.path(out,'standardized_organism_representation.csv'),row.names=FALSE)
# Verify that coefficient contrasts reproduce the published main comparison exactly.
old<-read.csv(file.path(src,'adjusted_phenotype_comparisons.csv'));old<-subset(old,kind=='organism_representation' & scenario=='primary' & contrast=='ASE versus presumed infection without ASE')
for(i in seq_len(nrow(old))){fit<-fits[[paste('organism_representation',old$organism[i],'NA','primary',sep='|')]];be<-coef(fit);ratio<-exp(be['ase_groupASE']-be['ase_groupPresumed infection without ASE']);stopifnot(abs(ratio-old$estimate[i])<1e-10)}
print(subset(r,organism %in% c('pseudomonas_aeruginosa','escherichia_coli','klebsiella_pneumoniae','staphylococcus_aureus')),row.names=FALSE)
