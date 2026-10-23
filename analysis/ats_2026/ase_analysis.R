suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2)})
base <- Sys.getenv('CLIF_ANALYSIS_ROOT');if(!nzchar(base))stop('Set CLIF_ANALYSIS_ROOT to the aggregate analysis workspace')
out<-file.path(base,'ASE_RUSH_updated_20261022/ase')
load_kind<-function(kind){files<-list.files(file.path(out,'inputs'),paste0('__',kind,'[.]csv$'),full.names=TRUE);bind_rows(lapply(files,function(f)mutate(read_csv(f,show_col_types=FALSE),site=sub('__.*','',basename(f)))))}
monthly<-read_csv(file.path(base,'ASE_RUSH_updated_20261022/monthly_main_analysis_input.csv'),show_col_types=FALSE)%>%filter(fluid_category=='All cultures')%>%select(site,ase_group,calendar_month,n_culture_events,n_positive_culture_events,n_icu_days)%>%mutate(calendar_month=as.Date(calendar_month))
orgs<-c('staphylococcus_aureus','escherichia_coli','klebsiella_pneumoniae','pseudomonas_aeruginosa','enterococcus_faecalis','enterococcus_faecium','stenotrophomonas_maltophilia','enterobacter_cloacae','serratia_marcescens','acinetobacter_baumannii','proteus_mirabilis','candida_albicans','candida_glabrata')
det<-load_kind('ase_monthly_organism_detection_counts')%>%filter(site!='emory',organism_category%in%orgs)%>%mutate(calendar_month=as.Date(calendar_month))%>%rename(source_days=n_icu_days)%>%inner_join(monthly,by=c('site','ase_group','calendar_month'))
stopifnot(!anyDuplicated(det[c('site','ase_group','calendar_month','organism_category')]),all(abs(det$source_days-det$n_icu_days)<1e-7),all(det$n_detection_events<=det$n_positive_culture_events))
source(file.path(out,'susceptibility_screen_functions.R'))
ast<-load_kind('ase_monthly_organism_antimicrobial_susceptibility')%>%mutate(calendar_month=as.Date(calendar_month))%>%screen_susceptibility_months()
keys<-read_csv(file.path(out,'original_77_pair_screen.csv'),show_col_types=FALSE)
ast<-semi_join(ast,keys,by=c('organism_category'='organism','antimicrobial_category'='antimicrobial'))
stopifnot(!anyDuplicated(ast[c('site','ase_group','calendar_month','organism_category','antimicrobial_category')]),all(ast$n_susceptible+ast$n_non_susceptible==ast$n_interpretable))
write_csv(ast%>%count(site,ase_group,observation_status,name='pair_months'),file.path(out,'ast_coverage_audit.csv'))
groups<-c('No presumed infection','Presumed infection without ASE','ASE');results<-list();margins<-list();fit_audit<-list();savedfits<-list()
fit_one<-function(input,y,den,proportion=FALSE,kind,org,drug=NA_character_,scenario='primary'){
 d<-input%>%filter(.data[[den]]>0,n_observed_culture_events>0)%>%mutate(Y=.data[[y]],T=.data[[den]])
 # Uninformative fixed-effect levels cannot identify a group contrast; avoid separated nuisance intercepts.
 before<-nrow(d)
 repeat{
  n0<-nrow(d)
  d<-d%>%group_by(site)%>%filter(sum(Y)>0,if(proportion)sum(T-Y)>0 else TRUE)%>%ungroup()%>%group_by(calendar_month)%>%filter(sum(Y)>0,if(proportion)sum(T-Y)>0 else TRUE)%>%ungroup()
  if(nrow(d)==n0)break
 }
 S<-n_distinct(d$site);basekey<-tibble(kind=kind,organism=org,antimicrobial=drug,scenario=scenario)
 audit<-bind_cols(basekey,tibble(input_rows=before,model_rows=nrow(d),sites=S,successes=sum(d$Y),denominator=sum(d$T),status='not_fitted'))
 if(S<5||n_distinct(d$ase_group)<3||n_distinct(d$calendar_month)<18||sum(d$Y)<30||(proportion&&sum(d$T-d$Y)<30)){
  fit_audit[[length(fit_audit)+1]]<<-mutate(audit,status='insufficient_information');return(invisible(NULL))
 }
 d<-d%>%mutate(site=factor(site),month=factor(calendar_month),ase_group=factor(ase_group,levels=groups))
 warns<-character()
 fit<-tryCatch(withCallingHandlers(glm(if(proportion)cbind(Y,T-Y)~site+month+ase_group else Y~site+month+ase_group+offset(log(T)),data=d,family=if(proportion)quasibinomial()else quasipoisson(),control=glm.control(maxit=100)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error')||!isTRUE(fit$converged)||anyNA(coef(fit))){fit_audit[[length(fit_audit)+1]]<<-mutate(audit,status='failed_or_nonconverged');return(invisible(NULL))}
 X<-model.matrix(fit);mu<-fitted(fit);w<-if(proportion)d$T*mu*(1-mu)else mu;res<-if(proportion)d$Y-d$T*mu else d$Y-mu
 B<-tryCatch(chol2inv(chol(crossprod(X,X*w))),error=function(e)NULL)
 if(is.null(B)){fit_audit[[length(fit_audit)+1]]<<-mutate(audit,status='singular_information');return(invisible(NULL))}
 U<-rowsum(X*res,as.character(d$site));V<-B%*%crossprod(U)%*%B*S/(S-1)*(nrow(X)-1)/(nrow(X)-ncol(X));crit<-qt(.975,S-1)
 pi_name<-'ase_groupPresumed infection without ASE';ase_name<-'ase_groupASE'
 for(contrast in c('ASE versus presumed infection without ASE','Presumed infection without ASE versus no presumed infection','ASE versus no presumed infection')){
  v<-setNames(rep(0,length(coef(fit))),names(coef(fit)))
  if(contrast=='ASE versus presumed infection without ASE'){v[ase_name]<-1;v[pi_name]<- -1}else if(contrast=='Presumed infection without ASE versus no presumed infection')v[pi_name]<-1 else v[ase_name]<-1
  beta<-sum(v*coef(fit));se<-sqrt(as.numeric(t(v)%*%V%*%v))
  results[[length(results)+1]]<<-bind_cols(basekey,tibble(contrast=contrast,sites=S,rows=nrow(d),estimate=exp(beta),lower=exp(beta-crit*se),upper=exp(beta+crit*se),p=2*pt(-abs(beta/se),S-1),site_list=paste(levels(d$site),collapse=';'),warnings=paste(unique(warns),collapse=';')))
 }
 if(kind=='non_susceptible_fraction'&&scenario=='primary'){
  # Standardize all three groups to the same observed site-month support; use equal site weights.
  common<-d%>%count(site,month)%>%filter(n==3)%>%select(site,month)
  grid<-semi_join(d,common,by=c('site','month'))%>%group_by(site,month)%>%slice_head(n=1)%>%ungroup()
  if(nrow(grid)&&n_distinct(grid$site)>=5){
   weights<-1/(n_distinct(grid$site)*as.numeric(table(grid$site)[as.character(grid$site)]))
   stopifnot(abs(sum(weights)-1)<1e-8)
   for(group in groups){
    nd<-grid;nd$ase_group<-factor(group,levels=groups)
    xp<-model.matrix(delete.response(terms(fit)),nd,contrasts.arg=fit$contrasts);prob<-plogis(as.numeric(xp%*%coef(fit)));m<-sum(weights*prob);gradient<-colSums(xp*as.numeric(weights*prob*(1-prob)));se_m<-sqrt(as.numeric(t(gradient)%*%V%*%gradient));se_l<-se_m/(m*(1-m))
    margins[[length(margins)+1]]<<-bind_cols(basekey,tibble(ase_group=group,standardized_fraction=m,lower=plogis(qlogis(m)-crit*se_l),upper=plogis(qlogis(m)+crit*se_l),sites=n_distinct(grid$site),site_months=nrow(grid)))
   }
  }
 }
 fit_audit[[length(fit_audit)+1]]<<-mutate(audit,status='estimated',warnings=paste(unique(warns),collapse=';'))
 savedfits[[paste(kind,org,drug,scenario,sep='|')]]<<-fit
}
for(org in orgs){
 d<-filter(det,organism_category==org)
 fit_one(d,'n_detection_events','n_icu_days',FALSE,'organism_detection_rate',org)
 fit_one(d,'n_detection_events','n_positive_culture_events',TRUE,'organism_representation',org)
 cat('Organism',org,'done\n')
}
for(i in seq_len(nrow(keys))){
 key<-keys[i,];d<-filter(ast,organism_category==key$organism,antimicrobial_category==key$antimicrobial,fraction_model_eligible)
 fit_one(d,'n_non_susceptible','n_interpretable',TRUE,'non_susceptible_fraction',key$organism,key$antimicrobial)
 fit_one(filter(d,testing_fraction>=.5),'n_non_susceptible','n_interpretable',TRUE,'non_susceptible_fraction',key$organism,key$antimicrobial,'testing_50pct')
 if(i%%10==0)cat('AST',i,'/',nrow(keys),'done\n')
}
r<-bind_rows(results)%>%group_by(kind,scenario,contrast)%>%mutate(q=p.adjust(p,'BH'))%>%ungroup()
write_csv(r,file.path(out,'adjusted_phenotype_comparisons.csv'));write_csv(bind_rows(margins),file.path(out,'standardized_non_susceptible_percentages.csv'));write_csv(bind_rows(fit_audit),file.path(out,'model_audit.csv'));saveRDS(savedfits,file.path(out,'phenotype_models.rds'))
descdet<-det%>%group_by(ase_group,organism_category)%>%summarise(detections=sum(n_detection_events),icu_days=sum(n_icu_days),positive_cultures=sum(n_positive_culture_events),.groups='drop')%>%mutate(detections_per_100_icu_days=100*detections/icu_days,percent_positive_cultures_with_organism=100*detections/positive_cultures)
write_csv(descdet,file.path(out,'organism_descriptive_by_group.csv'))
descast<-ast%>%filter(fraction_model_eligible)%>%group_by(ase_group,organism_category,antimicrobial_category)%>%summarise(non_susceptible=sum(n_non_susceptible),interpretable=sum(n_interpretable),sites=n_distinct(site),.groups='drop')%>%mutate(non_susceptible_percent=100*non_susceptible/interpretable)
write_csv(descast,file.path(out,'susceptibility_observed_by_group.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_session_info.txt'))
cat('PRIMARY ORGANISM FINDINGS\n');print(r%>%filter(scenario=='primary',kind!='non_susceptible_fraction',contrast=='ASE versus presumed infection without ASE')%>%select(kind,organism,estimate,lower,upper,q),n=30)
cat('PRIMARY AST FINDINGS\n');print(r%>%filter(scenario=='primary',kind=='non_susceptible_fraction',contrast=='ASE versus presumed infection without ASE')%>%arrange(q)%>%select(organism,antimicrobial,sites,estimate,lower,upper,q)%>%head(20),n=20)
