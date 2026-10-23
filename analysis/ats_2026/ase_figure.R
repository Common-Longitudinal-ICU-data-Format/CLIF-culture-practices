b<-Sys.getenv('CLIF_ANALYSIS_ROOT');if(!nzchar(b))stop('Set CLIF_ANALYSIS_ROOT to the aggregate analysis workspace');out<-file.path(b,'ASE_RUSH_updated_20261022/clinical_figure')
pos<-read.csv(file.path(b,'ASE_RUSH_updated_20261022/ase','standardized_organism_representation.csv'));ns<-read.csv(file.path(out,'common_pair_percentages.csv'));sel<-read.csv(file.path(out,'selected_common_pairs.csv'))
orgs<-unique(sel$organism);labels<-vapply(gsub('_',' ',orgs),function(x)paste0(toupper(substr(x,1,1)),substring(x,2)),character(1));groups<-c('No presumed infection','Presumed infection without ASE','ASE');heads<-c('Non-infected','Infected','ASE')
old_org<-c('pseudomonas_aeruginosa','escherichia_coli','klebsiella_pneumoniae','staphylococcus_aureus','candida_albicans','candida_glabrata','stenotrophomonas_maltophilia','enterococcus_faecium','enterococcus_faecalis','enterobacter_cloacae','serratia_marcescens','acinetobacter_baumannii','proteus_mirabilis');old_col<-c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999','#008B8B','#17365D','#800040','#6B8E23');cols<-setNames(old_col,old_org)[orgs]
druglabel<-function(x){x<-gsub('_','–',x);x<-paste0(toupper(substr(x,1,1)),substring(x,2));x[x=='Trimethoprim–sulfamethoxazole']<-'Trimethoprim–\nsulfamethoxazole';x[x=='Piperacillin–tazobactam']<-'Piperacillin–\ntazobactam';x}
ytop<-numeric(length(orgs));cursor<-0
for(i in seq_along(orgs)){ytop[i]<-cursor;cursor<-cursor-(2.05+1.05*sum(sel$organism==orgs[i]))}
limits<-c(cursor+.40,.62)
plot_figure<-function(){
 layout(matrix(1:4,nrow=1),widths=c(4.6,4.46,4.46,4.46));par(family='sans',fg='#26343D',col.axis='#26343D',col.lab='#26343D',cex=1,oma=rep(0,4))
 # Dedicated label column leaves all three data panels exactly the same width.
 par(mar=c(3.6,.3,3,.3));plot(NA,xlim=c(0,1),ylim=limits,axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')
 for(i in seq_along(orgs)){org<-orgs[i];y<-ytop[i];ss<-sel[sel$organism==org,];ss<-ss[order(ss$antimicrobial),]
  text(.98,y+.1,labels[i],adj=1,font=3,cex=1.38);text(.98,y-.38,'% of positive ICU culture events',adj=1,cex=.92,col='#536670');text(.98,y-.91,'Resistance',adj=1,cex=1.05,col='#536670')
  for(j in seq_len(nrow(ss)))text(.98,y-(.60+1.05*j),druglabel(ss$antimicrobial[j]),adj=1,cex=1.16)
 }
 widths<-numeric(3)
 for(g in 1:3){par(mar=c(3.6,.5,3,.7));plot(NA,xlim=c(-3,103),ylim=limits,axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i');widths[g]<-par('pin')[1];abline(v=seq(0,100,25),col='#E6EBEE');axis(1,at=seq(0,100,25),labels=paste0(seq(0,100,25),'%'),cex.axis=1.1,tck=-.012);mtext(heads[g],side=3,line=.8,font=2,cex=1.55);mtext('Percentage',side=1,line=2.3,cex=1.15)
  for(i in seq_along(orgs)){org<-orgs[i];y<-ytop[i];col<-cols[org];rect(-3,y-.62,103,y+.62,col=adjustcolor(col,alpha.f=.09),border=NA)
   draw<-function(v,yy,lo,hi,obs=FALSE,rep=FALSE){
    if(is.finite(lo)&&is.finite(hi)){stopifnot(lo<=v,v<=hi,lo>=0,hi<=100);segments(lo,yy,hi,yy,col=col,lwd=3);segments(c(lo,hi),yy-.11,c(lo,hi),yy+.11,col=col,lwd=2.5)}
    points(v,yy,pch=if(rep)23 else 21,bg=col,col='#4C5560',lwd=.9,cex=if(rep)2.4 else 2.1)
    low<-if(is.finite(lo))lo else v;high<-if(is.finite(hi))hi else v;right<-high<78
    if(!right && low<18){text(v,yy+.30,sprintf('%.1f%%%s',v,if(obs)'†'else ''),cex=1.14)}else text(if(right)high+3.1 else low-3.1,yy,sprintf('%.1f%%%s',v,if(obs)'†'else ''),adj=if(right)0 else 1,cex=1.14)
   }
   d<-subset(pos,organism==org & ase_group==groups[g]);stopifnot(nrow(d)==1);draw(d$standardized_percent,y,d$lower,d$upper,rep=TRUE)
   ss<-sel[sel$organism==org,];ss<-ss[order(ss$antimicrobial),]
   for(j in seq_len(nrow(ss))){d<-subset(ns,organism==org & antimicrobial==ss$antimicrobial[j] & ase_group==groups[g]);stopifnot(nrow(d)==1);draw(d$percent,y-(.6+1.05*j),d$lower,d$upper,d$estimate_type=='pooled_observed')}
   if(i<length(orgs))abline(h=ytop[i+1]+.72,col='#C5D0D5',lty=3)
  };box(col='#AAB9C1')
 }
 stopifnot(diff(range(widths))<1e-10)
 write.csv(data.frame(panel=heads,width_inches=widths),file.path(out,'panel_width_validation.csv'),row.names=FALSE)
}
png(file.path(out,'ASE_clinical_phenotype_panels.png'),width=18,height=19.7,units='in',res=260);plot_figure();dev.off()
jpeg(file.path(out,'ASE_clinical_phenotype_panels.jpg'),width=18,height=19.7,units='in',res=260,quality=98);plot_figure();dev.off()

cairo_pdf(file.path(out,'ASE_clinical_phenotype_panels.pdf'),width=18,height=19.7);plot_figure();dev.off()
