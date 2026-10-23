suppressPackageStartupMessages({library(readr);library(jsonlite)})
analysis_root<-Sys.getenv('CLIF_ANALYSIS_ROOT');if(!nzchar(analysis_root))stop('Set CLIF_ANALYSIS_ROOT to the aggregate analysis workspace');out<-file.path(analysis_root,'CLIF_refresh_20261022/atlas')
rd<-function(n)read_csv(file.path(out,n),show_col_types=FALSE)
org<-rd('display_organisms.csv');src<-rd('display_sources.csv');groups<-rd('display_organism_types.csv');meta<-fromJSON(file.path(out,'analysis_metadata.json'));ast<-rd('resistance_display_counts.csv');drugs<-rd('resistance_drug_tracks.csv');right_org<-rd('resistance_organisms.csv');right_groups<-rd('resistance_organism_types.csv')
source_labels<-rd('specimen_display_labels.csv');src$label<-source_labels$label[match(src$specimen,source_labels$specimen)];stopifnot(!anyNA(src$label))
ast$organism_order<-match(ast$organism,right_org$organism)-1;ast$drug_track<-match(ast$antimicrobial,drugs$antimicrobial)
stopifnot(!anyNA(ast$organism_order),!anyNA(ast$drug_track))
class_colors<-c('Beta-lactams'='#377EB8','Aminoglycosides'='#4DAF4A','Fluoroquinolones'='#E58B22','Tetracyclines'='#984EA3','Glycopeptides'='#CC6688','Folate antagonists'='#9A7438')
# Pale class backgrounds identify untested rings without altering heatmap cells.
class_tints<-setNames(vapply(class_colors,function(z){v<-col2rgb(z)/255;rgb(.14*v[1]+.86,.14*v[2]+.86,.14*v[3]+.86)},character(1)),names(class_colors))
ast_palette<-colorRampPalette(c('#FFF7BC','#FEC44F','#EC7014','#B30000','#49006A'))(101)
N<-nrow(org);gap<-2;step<-(360-gap*nrow(groups))/N;angle<-numeric(N);group_start<-group_end<-numeric(nrow(groups));cursor<-90
for(g in seq_len(nrow(groups))){ix<-which(org$group_id==groups$group_id[g]);group_start[g]<-cursor;angle[ix]<-cursor-(seq_along(ix)-.5)*step;group_end[g]<-cursor-length(ix)*step;cursor<-group_end[g]-gap}
stopifnot(abs(cursor+270)<1e-8)
build_geometry<-function(oo,gg){
 nn<-nrow(oo);sp<-(360-gap*nrow(gg))/nn;aa<-numeric(nn);starts<-ends<-numeric(nrow(gg));cur<-90
 for(j in seq_len(nrow(gg))){ii<-which(oo$group_id==gg$group_id[j]);starts[j]<-cur;aa[ii]<-cur-(seq_along(ii)-.5)*sp;ends[j]<-cur-length(ii)*sp;cur<-ends[j]-gap}
 stopifnot(abs(cur+270)<1e-8)
 list(org=oo,groups=gg,N=nn,angle=aa,group_start=starts,group_end=ends,step=sp)
}
left_geom<-build_geometry(org,groups);right_geom<-build_geometry(right_org,right_groups)

outer<-1;inner<-.52;tilewidth<-.90*step;centers<-c(0,2.8)
xy<-function(a,r,cx=0)cbind(cx+r*cos(a*pi/180),r*sin(a*pi/180))
arc<-function(a,b,r,cx,n=5)xy(seq(a,b,length.out=n),r,cx)
annulus<-function(a,b,r0,r1,col,cx=0,n=5){p<-rbind(arc(a,b,r1,cx,n),arc(b,a,r0,cx,n));polygon(p[,1],p[,2],col=col,border=NA)}
curved_label<-function(label,mid,r,available,col,cx,cex_max=2.4){
 chars<-strsplit(label,'')[[1]];cex<-cex_max;widths<-strwidth(chars,cex=cex,font=2)*1.05;span<-sum(widths)/r*180/pi
 if(span>available*.90){cex<-cex*available*.90/span;widths<-strwidth(chars,cex=cex,font=2)*1.05}
 offsets<-(cumsum(widths)-widths/2-sum(widths)/2)/r*180/pi;lower<-sin(mid*pi/180)<0;aa<-if(lower)mid+offsets else mid-offsets
 for(j in seq_along(chars)){q<-xy(aa[j],r,cx)[1,];text(q[1],q[2],chars[j],srt=aa[j]+if(lower)90 else -90,cex=cex,font=2,col=col,adj=c(.5,.5))}
}
frame<-function(cx,geom=left_geom){
 org<-geom$org;groups<-geom$groups;N<-geom$N;angle<-geom$angle;group_start<-geom$group_start;group_end<-geom$group_end
 for(k in seq_len(N)){a<-angle[k];deg<-((a+180)%%360)-180;left<-abs(deg)>90;rot<-if(left)deg+180 else deg;p<-xy(a,1.023,cx)[1,];text(p[1],p[2],org$label[k],srt=rot,adj=if(left)c(1,.5)else c(0,.5),cex=1.55,font=2,col='#18323C')}
 for(g in seq_len(nrow(groups))){
  annulus(group_end[g],group_start[g],.501,.513,groups$color[g],cx,n=max(8,groups$n_categories[g]));mid<-(group_start[g]+group_end[g])/2;available<-group_start[g]-group_end[g]
  if(groups$organism_type[g]=='Mycobacteria/AFB'){
   curved_label('Mycobacteria',mid,.479,available,groups$color[g],cx)
  }else curved_label(groups$label[g],mid,.479,available,groups$color[g],cx)
 }
}
detection_color<-function(n,source){alpha<-.30+.70*log1p(n)/log1p(meta$display_shared_maximum_cell_count);adjustcolor(source,alpha.f=alpha)}
plot_atlas<-function(mode){
 d<-rd(paste0(mode,'_display_counts.csv'));totals<-rd(paste0(mode,'_display_source_totals.csv'));whole<-meta$whole_cohort_totals[[mode]];view<-meta$display_views[[mode]]
 stopifnot(nrow(d)==view$reported_cells,nrow(org)==view$organism_categories,all(d$n_detections>0))
 draw<-function(){
  par(mar=c(0,0,0,0),family='sans',xpd=NA);plot.new();plot.window(xlim=c(-1.5,4.3),ylim=c(-1.525,1.525),asp=1,xaxs='i',yaxs='i')
  # Detection retains all five groups; resistance has independent geometry without fungi/yeast.
  dw<-(outer-inner)/nrow(src)
  for(k in seq_len(nrow(src))){r1<-outer-(k-1)*dw;annulus(0,360,r1-dw*.96,r1,'#F0F2F4',centers[1],n=500)}
  for(g in seq_len(nrow(groups)))annulus(group_end[g]-gap,group_end[g],inner,1.006,'white',centers[1])
  for(k in seq_len(nrow(d))){row<-d[k,];a<-angle[row$organism_order+1];r1<-outer-(row$track-1)*dw;annulus(a-tilewidth/2,a+tilewidth/2,r1-dw*.96,r1,detection_color(row$n_detections,src$color[row$track]),centers[1])}
  aw<-(outer-inner)/nrow(drugs)
  for(k in seq_len(nrow(drugs))){r1<-outer-(k-1)*aw;annulus(0,360,r1-aw*.96,r1,class_tints[drugs$antibiotic_class[k]],centers[2],n=500)}
  for(g in seq_len(nrow(right_geom$groups)))annulus(right_geom$group_end[g]-gap,right_geom$group_end[g],inner,1.006,'white',centers[2])
  for(k in seq_len(nrow(ast))){row<-ast[k,];a<-right_geom$angle[row$organism_order+1];r1<-outer-(row$drug_track-1)*aw;cc<-ast_palette[round(row$non_susceptible_percent)+1];annulus(a-.90*right_geom$step/2,a+.90*right_geom$step/2,r1-aw*.96,r1,cc,centers[2])}
  frame(centers[1],left_geom);frame(centers[2],right_geom)
  # Compact center keys use larger text within the reduced central opening.
  cx<-centers[1]
  text(cx-.255,.32,'Source',adj=0,cex=2.2,font=2,col='#193640');text(cx+.30,.32,'n',adj=1,cex=2.2,font=2,col='#193640')
  ys<-seq(.26,-.32,length.out=nrow(src))
  for(k in seq_len(nrow(src))){y<-ys[k];rect(cx-.33,y-.010,cx-.31,y+.010,col=src$color[k],border=NA);text(cx-.292,y,sprintf('%02d',k),adj=0,cex=1.7,col='#627782');text(cx-.235,y,src$label[k],adj=0,cex=2.15,font=2,col='#263E48');ix<-match(src$specimen[k],totals$specimen);text(cx+.30,y,format(totals$n_detections[ix],big.mark=',',trim=TRUE),adj=1,cex=2.0,col='#526C78')}
  # The detection-count bar is intentionally omitted.
  cx<-centers[2];text(cx-.26,.32,'Antibiotic',adj=0,cex=2.2,font=2,col='#193640');text(cx+.30,.32,'Tested n',adj=1,cex=2.2,font=2,col='#193640')
  # Class headers align with the outside-to-inside antibiotic ring order.
  classes<-unique(drugs$antibiotic_class);line_step<-.49/(nrow(drugs)+length(classes)-1);y<-.275
  for(cl in classes){
   rect(cx-.355,y-.006,cx-.345,y+.006,col=class_colors[cl],border=NA);text(cx-.335,y,cl,adj=0,cex=1.7,font=2,col=class_colors[cl]);y<-y-line_step
   for(k in which(drugs$antibiotic_class==cl)){text(cx-.335,y,sprintf('%02d',k),adj=0,cex=1.65,font=2,col='#49006A');text(cx-.27,y,drugs$label[k],adj=0,cex=2.05,font=2,col='#263E48');text(cx+.30,y,format(drugs$n_interpretable[k],big.mark=',',trim=TRUE),adj=1,cex=1.95,col='#526C78');y<-y-line_step}
  }
  text(cx,-.255,'Resistance (%)',cex=2.05,font=2,col='#49006A');xx<-seq(cx-.30,cx+.30,length.out=102)
  for(j in 1:101)rect(xx[j],-.32,xx[j+1],-.28,col=ast_palette[j],border=NA)
  tt<-c(0,25,50,75,100);text(cx-.30+tt/100*.60,-.35,paste0(tt,'%'),cex=1.85,col='#49006A',font=2)
  for(j in seq_along(class_tints))rect(cx-.23+(j-1)*.024,-.405,cx-.208+(j-1)*.024,-.385,col=class_tints[j],border='#AAB2BA',lwd=.5);text(cx-.065,-.395,'No tests',adj=0,cex=1.8,col='#526C78')

 }
 prefix<-paste0('organism_source_atlas_',mode,'_twin_wheels_balanced_keys')
 png(file.path(out,paste0(prefix,'.png')),width=24000,height=12621,res=252.42,type='cairo');draw();dev.off()
 cairo_pdf(file.path(out,paste0(prefix,'.pdf')),width=95.0819672131,height=50);draw();dev.off()
 svg(file.path(out,paste0(prefix,'.svg')),width=95.0819672131,height=50);draw();dev.off()
 cat(prefix,': detection organisms',N,'; resistance organisms',nrow(right_org),';',nrow(src),'source tracks;',nrow(drugs),'antibiotic tracks\n')
}
plot_atlas('all11');plot_atlas('exclude_emory')
