#CHC analysis
#SKHsu_20190412
#2019_0530_update_v2 (including expression)
#2020_0901_update_publication (chopped, including mating assays)

rm(list=ls())
library(biomaRt)
library(ggfortify)
library(agricolae)
library(limma)
library(edgeR)
library(topGO)
library(pheatmap)
library(ggplot2)
####import basic function####
cont_table=function(query,background,classifyer){
  p1=length(Reduce(intersect,list(query,background,classifyer)))
  q1=length(intersect(query,background))-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}
genomic_enrichment=function(query,background,ann_summary,windowsize,steps,chr){
  X0=0
  X=windowsize
  S=steps
  pos=c()
  counts=c()
  fold_enrichment=c()
  p.val=c()
  x=ann_summary[ann_summary[,1]%in%chr,]
  while (X0<max(x[,c(2,3)])){
    gene_within=x[x[,2]<X0+X/2&x[,2]>X0-X/2,5]
    p=c()
    fe=c()
    for (j in 1:length(query)){
      p=cbind(p,fisher.test(cont_table(query[[j]],background,gene_within),alternative = "greater")$p.value) 
      fe=cbind(fe,(length(intersect(gene_within,query[[j]]))/length(gene_within))/(length(query[[j]])/length(background)))
    }
    p.val=rbind(p.val,p)
    fold_enrichment=rbind(fold_enrichment,fe)
    colnames(p.val)=names(query)
    colnames(fold_enrichment)=names(query)
    counts=c(counts,length(x[x[,2]<X0+X/2&x[,2]>X0-X/2,5]))
    pos=c(pos,X0)
    X0=X0+S
  }
  out_p=data.frame("mid_pos"=pos,"total_gene_counts"=counts,p.val)
  out_fe=data.frame("mid_pos"=pos,"total_gene_counts"=counts,fold_enrichment)
  out=list("p"=out_p,"fe"=out_fe)
  return(out)
}
ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}
ROC_like_curve=function(a,b,bin,mycol="black",positive=T){
  lvl=seq(0,1,bin)
  if (positive==T){
    inter=c()
    for (i in lvl){
      conserved_genes=intersect(rownames(a),rownames(b))
      tmp1=a[conserved_genes,]
      tmp2=b[conserved_genes,]
      index=sum(rank(tmp1$logFC,ties.method = "random")>length(conserved_genes)*(1-i)&rank(tmp2$logFC,ties.method = "random")>length(conserved_genes)*(1-i))
      inter=c(inter,index/(length(conserved_genes)*i))
    }
    print(data.frame(r=cor(tmp1$logFC,tmp2$logFC),AUC=auc(lvl,inter,0,1,type = "spline")))
    points(lvl,inter,asp=1,type="l",lwd=3,col=mycol)
  }
  else {
    inter=c()
    for (i in lvl){
      conserved_genes=intersect(rownames(a),rownames(b))
      tmp1=a[conserved_genes,]
      tmp2=b[conserved_genes,]
      index=sum(rank(tmp1$logFC,ties.method = "random")<length(conserved_genes)*i&rank(tmp2$logFC,ties.method = "random")<length(conserved_genes)*i)
      inter=c(inter,index/(length(conserved_genes)*i))
    }
    print(data.frame(r=cor(tmp1$logFC,tmp2$logFC),AUC=auc(lvl,inter,0,1,type = "spline")))
    points(lvl,inter,asp=1,type="l",lwd=3,col=mycol)
  }
}
ROC_like_curve_GO=function(a,b,bin,mycol="black"){
  lvl=seq(0,1,bin)
  inter=c()
  for (i in lvl){
    conserved_GO=intersect(a$GO.ID,b$GO.ID)
    tmp1=a[a$GO.ID%in%conserved_GO,]
    tmp2=b[b$GO.ID%in%conserved_GO,]
    index=length(intersect(tmp1[0:(length(conserved_GO)*i),1],tmp2[0:(length(conserved_GO)*i),1]))
    inter=c(inter,index/(length(conserved_GO)*i))
  }
  print(auc(lvl,inter,0,1,type = "spline"))
  points(lvl,inter,asp=1,type="l",lwd=3,col=mycol)
}
pseudo_chr=function(x){
  for (i in 1:length(unique(x$CHR))){
    x$pseudo_CHR[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=i
  }
  return(x)
}
pseudo_pos=function(x){
  for (i in 1:length(unique(x$CHR))){
    if (i==1) x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=x$BP[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]
    else x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=x$BP[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]+max(x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i-1]])
  }
  return(x)
}
p_resolution=function(x){
  x$P=10^(-as.numeric(strsplit2(x$"-logp",split = "=")[,2]))
  return(x)
}
SNPID_gen=function(x){
  x$SNPID=paste(x$CHR,x$BP,sep = "_")
  return(x)
}
snp_identifier=function(query,snpset,threshold1,threshold2){
  subset1=snpset[snpset$CHR%in%query$Chr&snpset$CHR=="X"&snpset$P<threshold1,]
  subset2=snpset[snpset$CHR%in%query$Chr&snpset$CHR!="X"&snpset$P<threshold2,]
  subset=rbind(subset1,subset2)
  count=sum(query$start<subset$BP&query$end>subset$BP,na.rm = T)
  SNPID=as.character(na.omit(subset$SNPID[query$start<subset$BP&query$end>subset$BP]))
  
  query$SNP_counts=count
  if(length(SNPID)>0) query$SNP_ID=paste(SNPID,collapse = ",")
  else query$SNP_ID=NA
  return(query)
}
ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}
####gene expression evolution####
count_corrected=read.csv("/Volumes/Temp1/shengkai/remap_run/florida_result_fullset_v7/corrected_readcount.csv",header = T,row.names = 1)
group_gene=apply(strsplit2(colnames(count_corrected),split = "_")[,c(1,3)],1,function(x) paste(x,collapse = ""))
y_corrected=DGEList(counts=count_corrected,group = group_gene)
y_corrected=calcNormFactors(y_corrected)
LRT_res_f=read.csv("/Volumes/Temp1/shengkai/remap_run/florida_result_fullset_v7/test/LRT_res_hotevoF.csv",header = T,stringsAsFactors = F,row.names = 1)
LRT_res_m=read.csv("/Volumes/Temp1/shengkai/remap_run/florida_result_fullset_v7/test/LRT_res_hotevoM.csv",header = T,stringsAsFactors = F,row.names = 1)
LRT_res_merged=cbind(LRT_res_f,LRT_res_m)


#to do: order the genes by annotation if possible#done
GOI_CHC=c("FBgn0033246","FBgn0027571","FBgn0042627","FBgn0029975","FBgn0035471","FBgn0032394", #FAS
          "FBgn0037762","FBgn0050008","FBgn0034382","FBgn0037765", #Elogase
          "FBgn0086687","FBgn0043043","FBgn0029172", #desaturase
          "FBgn0010019","FBgn0015623","FBgn0033524","FBgn0038037","FBgn0030615", #cytochrome P450
          "FBgn0035438","FBgn0031479","FBgn0038465","FBgn0004577","FBgn0011828", #peroxidase
          "FBgn0260941","FBgn0032055","FBgn0038033","FBgn0036512","FBgn0031684","FBgn0030612",
          "FBgn0037832","FBgn0037623","FBgn0036698","FBgn0004108","FBgn0024740","FBgn0037819",
          "FBgn0031478")
GOI_CHC_name=c("ACC","FASN1","FASN2","KAR","TER","HADC1",
               "EloF","CG30008","CG18609","CG9458",
               "Desat1","Desat2","DesatF",
               "Cyp4g1","Cpr","Cyp49a1","Cyp9f2","Cyp4S3",
               "PHGPx","Prx6005","Irc","Pxd","Pxn",
               "app","Sgp","CG10097","CG16979","ND-13A","CG5599","Desi","CG9801","CG7724","Nrt","Lip2","CG14688","CG8814")

ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
GOI_OBP=intersect(read.table("/Volumes/Temp1/shengkai/CHC/OBP_FlyBase_IDs.txt",stringsAsFactors = F)[,1],rownames(LRT_res_f))
OBP_tab=ID_converter(GOI_OBP,db = ensembl,attributes = c("flybase_gene_id","external_gene_name"),filters = "flybase_gene_id")

fisher.test(cont_table(rownames(LRT_res_m)[LRT_res_f$padj<0.05&LRT_res_m$padj<0.05&LRT_res_f$logFC>0&LRT_res_m$logFC>0],
                       rownames(LRT_res_m),GOI_CHC),alternative = "greater")

png("/Volumes/Temp1/shengkai/CHC/gene_expression_CHC_GOI_both.png",height = 15,width = 15,units = "cm",res = 600,pointsize = 8)
par(mar=c(5.1,5.1,2.1,2.1))
bp=barplot(t(LRT_res_merged[GOI_CHC,c(1,6)]),beside = T,names.arg = GOI_CHC_name,xlab=expression(log[2](FC)),
           las=1,horiz = T,xlim=c(-1.2,1),font.axis=1.5,col=c("salmon","royalblue"))
points(ifelse(as.vector(t(LRT_res_merged[GOI_CHC,c(1,6)]))>0,as.vector(t(LRT_res_merged[GOI_CHC,c(1,6)]))+0.02,as.vector(t(LRT_res_merged[GOI_CHC,c(1,6)]))-0.02),
       bp,pch=ifelse(as.vector(t(LRT_res_merged[GOI_CHC,c(5,10)]))<0.05,"*",""))
dev.off()
#cor:0.50 (pearson) 0.60 (spearman)

png("/Volumes/Temp1/shengkai/CHC/FigureS1.png",height = 20,width = 20,units = "cm",res = 600,pointsize = 8)
par(mar=c(5.1,5.1,2.1,2.1),bg=NA)
bp=barplot(t(LRT_res_merged[rev(GOI_CHC),c(1,6)]),beside = T,names.arg = rev(GOI_CHC_name),xlab=expression(log[2](FC)),
           las=1,horiz = T,xlim=c(-1.2,1),font.axis=1.5,col=c("salmon","royalblue"))
points(ifelse(as.vector(t(LRT_res_merged[rev(GOI_CHC),c(1,6)]))>0,as.vector(t(LRT_res_merged[rev(GOI_CHC),c(1,6)]))+0.02,as.vector(t(LRT_res_merged[rev(GOI_CHC),c(1,6)]))-0.02),
       bp,pch=ifelse(as.vector(t(LRT_res_merged[rev(GOI_CHC),c(5,10)]))<0.05,"*",""))
dev.off()

myColor = colorRampPalette(c("blue", "white", "red"))(100)
myBreaks = c(seq(-1.2, 0, length.out=ceiling(100/2) + 1),
             seq(0.9/100, 0.9,length.out=floor(100/2)))

pheatmap(na.omit(LRT_res_merged[GOI_CHC,c(1,6)]),cellwidth = 20,cellheight = 20,
         cluster_rows = F,cluster_cols = F,labels_col = c("female","male"),
         color = myColor,breaks = myBreaks,
         filename = "/Volumes/Temp1/shengkai/CHC/Figure1c_goi.png",res=600)

bp=barplot(t(LRT_res_merged[OBP_tab[,1],c(1,6)]),beside = T,names.arg = OBP_tab[,2],xlab=expression(log[2](FC)),
           las=1,horiz = T,xlim=c(-2,2),font.axis=3,col=c("salmon","royalblue"))
points(ifelse(as.vector(t(LRT_res_merged[OBP_tab[,1],c(1,6)]))>0,as.vector(t(LRT_res_merged[OBP_tab[,1],c(1,6)]))+0.02,as.vector(t(LRT_res_merged[OBP_tab[,1],c(1,6)]))-0.02),
       bp,pch=ifelse(as.vector(t(LRT_res_merged[OBP_tab[,1],c(5,10)]))<0.05,"*",""))


#PCA
# pca_all_m=prcomp(t(log(cpm(y_corrected)[,substr(group_gene,2,2)=="m"])))
# pca_GOI_CHC_m=prcomp(t(log(cpm(y_corrected)[rownames(y_corrected)%in%GOI_CHC,substr(group_gene,2,2)=="m"])))
# plot(pca_all_m$x,col=ifelse(substr(group_gene[substr(group_gene,2,2)=="m"],1,1)=="B","forestgreen","salmon"),
#      pch=c(rep(19,5),rep(1:10,each=3)),asp=1)
# plot(pca_GOI_CHC_m$x,col=ifelse(substr(group_gene[substr(group_gene,2,2)=="m"],1,1)=="B","forestgreen","salmon"),
#      pch=c(rep(19,5),rep(1:10,each=3)),asp=1)


####gene expression divergence####
count_dat=read.csv("/Volumes/Temp1/shengkai/remap_run/monster_fullset_readcounts.csv",header = T,row.names = 1)
count_dat_use=count_dat[,strsplit2(colnames(count_dat),"_")[,3]%in%"m"&substr(colnames(count_dat),1,1)%in%c("B","H")]
count_dat_use_sort=count_dat_use[,order(colnames(count_dat_use))]
count_dat_use_sort_filtered=count_dat_use_sort[apply(cpm(count_dat_use_sort),1,function(x) !sum(x<0.1)>=1),]

repl=gl(10,3)
y_dvg=DGEList(counts = count_dat_use_sort_filtered[,-1:-5],group = repl)
y_dvg=calcNormFactors(y_dvg)
ModelDesign=model.matrix(~repl)
DGE=estimateDisp(y_dvg,design = ModelDesign,robust = T)
GLM=glmFit(DGE,design = ModelDesign)
LRT_res=glmLRT(GLM,coef = c(2:10))
dvg_res_table=LRT_res$table
dvg_res_table$padj=p.adjust(dvg_res_table$PValue,method = "BH")
write.table(dvg_res_table,"/Volumes/Temp1/shengkai/CHC/gene_expression_divergence/DE_analysis.txt",quote = F,sep = "\t")
sum(dvg_res_table$padj<0.05&apply(dvg_res_table[,1:9],1,function(x) any(abs(x)>0.5)))#3062, 1022
background=rownames(y_dvg)
sig_dvg_gene=rownames(y_dvg)[dvg_res_table$padj<0.05]
png("/Volumes/Temp1/shengkai/CHC/gene_expression_divergence/sig_dvg_gene_heatmap.png",width = 15,height = 15,units = "cm",res = 600,pointsize = 6)
pheatmap(log10(cpm(y_dvg[sig_dvg_gene,])),scale="row",show_rownames = F)
dev.off()
pheatmap(log10(cpm(count_dat_use_sort_filtered[sig_dvg_gene,])),scale="row",show_rownames = F)

pca_dvg=prcomp(t(log10(cpm(y_dvg[sig_dvg_gene,]))))
plot(pca_dvg$x,col=repl)

pheatmap(t(pca_dvg$x)[1:10,],cluster_cols = F,cluster_rows = F)

tmp=factor(as.integer(background%in%sig_dvg_gene))
names(tmp)=background
tgd=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
resTopGO.classic=runTest(tgd, algorithm = "classic", statistic = "Fisher")
resTopGO.weight01=runTest(tgd, algorithm = "weight01", statistic = "Fisher")
GO_res_dvg_gene=GenTable(tgd,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
head(GO_res_dvg_gene,20)
write.table(GO_res_dvg_gene,"/Volumes/Temp1/shengkai/CHC/gene_expression_divergence/GO_res_dvg.txt",quote = F,sep = "\t")

flyatlas=read.table("/Volumes/Temp1/shengkai/common_info/fly_atlas_enrichment.table",header = T,stringsAsFactors = F)
p.val=c()
odds=c()
exp.num=c()
tissue_enriched_ID=list()
for(i in seq(4,36,2)){
  tissue_specific_ID=flyatlas$FB[flyatlas[,i]>2]
  p=fisher.test(cont_table(sig_dvg_gene,background,tissue_specific_ID),alternative = "greater")$p.value
  o=fisher.test(cont_table(sig_dvg_gene,background,tissue_specific_ID),alternative = "greater")$estimate
  g=intersect(sig_dvg_gene,tissue_specific_ID)
  e=length(sig_dvg_gene)*sum(tissue_specific_ID%in%background)/length(background)
  p.val=rbind(p.val,p)
  odds=rbind(odds,o)
  tissue_enriched_ID[[(i-2)/2]]=g
  exp.num=rbind(exp.num,e)
}
odds=odds[,1]
exp.num=exp.num
exp.num=round(exp.num,0)
obs.num= t(sapply(tissue_enriched_ID,function(x) sapply(x,length)))
padj=p.adjust(p.val,method = "BH")
names(tissue_enriched_ID)=strsplit2(colnames(flyatlas)[seq(4,36,2)],split = "_")[,1]
names(padj)=c("Br","Hd","Cr","Mg","Hg","Tb","Ov","Ts","Ag","Tg","Cs","Sg","SmV","SmM","Fb","Ey","Hr")
names(odds)=c("Br","Hd","Cr","Mg","Hg","Tb","Ov","Ts","Ag","Tg","Cs","Sg","SmV","SmM","Fb","Ey","Hr")
names(exp.num)=c("Br","Hd","Cr","Mg","Hg","Tb","Ov","Ts","Ag","Tg","Cs","Sg","SmV","SmM","Fb","Ey","Hr")
names(obs.num)=c("Br","Hd","Cr","Mg","Hg","Tb","Ov","Ts","Ag","Tg","Cs","Sg","SmV","SmM","Fb","Ey","Hr")
barplot(-log10(padj),las=2,horiz = T)


flyatlas2=read.table("/Volumes/Temp1/shengkai/common_info/flyatlas2/version_18.05.25/flyaltlas2_log2fc.txt",header = T,stringsAsFactors = F)
tissues=read.delim("/Volumes/Temp1/shengkai/common_info/flyatlas2/version_18.05.25/Tissue.txt",header = T,stringsAsFactors = F)
p.val2=c()
odds2=c()
tissue_enriched_ID2=list()
for(i in 2:dim(flyatlas2)[2]){
  tissue_specific_ID2=flyatlas2$FBgnID[!is.na(flyatlas2[,i])&flyatlas2[,i]>1]
  p=fisher.test(cont_table(sig_dvg_gene,background,tissue_specific_ID2),alternative = "greater")$p.value
  o=fisher.test(cont_table(sig_dvg_gene,background,tissue_specific_ID2),alternative = "greater")$estimate
  g=intersect(sig_dvg_gene,tissue_specific_ID2)
  p.val2=rbind(p.val2,p)
  odds2=c(odds2,o)
  tissue_enriched_ID2[[i-1]]=g
}
padj2=p.adjust(p.val2[1:13],method = "BH")
odds2=odds2[1:13]
names(tissue_enriched_ID2)=colnames(flyatlas2)[2:14]
names(padj2)=strsplit2(colnames(flyatlas2)[2:14],split = "_")[,3]
names(odds2)=strsplit2(colnames(flyatlas2)[2:14],split = "_")[,3]
png("/Volumes/Temp1/shengkai/CHC/gene_expression_divergence/tissue_enrichment.png",height = 8.7,width = 8.7,units = "cm",res=600,pointsize = 6)
barplot(-log10(padj2),las=1,horiz = T,xlab=expression(-log[10](P[adj])))
dev.off()

fisher.test(cont_table(sig_dvg_gene,background,GOI_CHC),alternative = "greater")

GOI_ACP=read.table("/Volumes/Temp1/shengkai/CHC/ACP_flybase.txt",stringsAsFactors = F)[,1]
GOI_SFP=read.table("/Volumes/Temp1/shengkai/CHC/SFP_flybase.txt",stringsAsFactors = F)[,1]
fisher.test(cont_table(sig_dvg_gene,background,GOI_ACP),alternative = "greater")
fisher.test(cont_table(sig_dvg_gene,background,GOI_SFP),alternative = "greater")

boxplot(log10(cpm(count_dat_use_sort_filtered)[intersect(GOI_ACP,background)[8],])~c(rep(0,5),repl),
        main=intersect(GOI_ACP,background)[8],xlab="Population",ylab=expression(log[10](CPM)))
boxplot(log10(cpm(count_dat_use_sort_filtered)[intersect(GOI_SFP,background)[2],])~c(rep(0,5),repl),
        main=intersect(GOI_SFP,background)[2],xlab="Population",ylab=expression(log[10](CPM)))

pheatmap(log10(cpm(count_dat_use_sort_filtered)[intersect(GOI_SFP,sig_dvg_gene),]),scale = "row",cluster_rows = F,cluster_cols = F)

annot_col=data.frame(population=c(rep("Anc.",5),paste0("H",formatC(repl,digits = 1,flag = "0"))))
annot_colors=list(population=c(Anc.="forestgreen",H01="brown1",H02="chocolate",H03="coral",H04="orange",
                               H05="salmon",H06="orchid",H07="orangered",H08="firebrick",H09="magenta",
                               H10="darkorange"))
rownames(annot_col)=colnames(count_dat_use_sort_filtered)

repr.div.gene = intersect(sig_dvg_gene,genesInTerm(tgd,"GO:0032504")[[1]])                  
pca_goi=prcomp(t(log10(cpm(count_dat_use_sort_filtered[repr.div.gene,]))))
ve_goi=pca_goi$sdev^2/sum(pca_goi$sdev^2)
sum(ve_goi[1:9])

png("./gene_expression_divergence/reproductive_genes_PCA_ve.png",width = 8.7,height = 8.7,units = "cm",res=600,pointsize = 8)
par(mar=c(5,5,2,2))
plot(ve_goi*100,type="b",pch=19,xlab="Principal components (PCs)",ylab="Variance explained (%)")
abline(h=5,col="red",lty=2)
dev.off()

plotscore=pca_goi$x
plotscore[,2]=-pca_goi$x[,2]
plotscore[,4]=-pca_goi$x[,4]
plotscore[,6]=-pca_goi$x[,6]

pheatmap(t(plotscore)[1:9,],show_rownames = T,cluster_cols = F,cluster_rows = F,show_colnames = F,
         annotation_col = annot_col,annotation_colors = annot_colors,annotation_names_col = F)

png("/Volumes/Temp1/shengkai/CHC/figure3b_new2.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 8)
layout(matrix(3:1,3,1),heights = c(1,1,1.5))
par(mar=c(5,5,2,1))
boxplot(plotscore[,1]~annot_col$population,cex=0.75,lwd=0.75,col=NULL,
        ylab=paste("PC",1),xlab="Population",las=1,ylim=c(-1,1),border=brewer.pal(6,"Set1")[1])
for (i in 2:3){
  par(mar=c(0,5,2,1))
  boxplot(plotscore[,i]~annot_col$population,xaxt="n",cex=0.75,lwd=0.75,col=NULL,
          ylab=paste("PC",i),xlab="",las=1,ylim=c(-1,1),border=brewer.pal(6,"Set1")[i])
}
dev.off()

png("/Volumes/Temp1/shengkai/CHC/figureS4.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 8)
layout(matrix(3:1,3,1),heights = c(1,1,1.5))
par(mar=c(5,5,2,1))
boxplot(plotscore[,1]~annot_col$population,cex=0.75,lwd=0.75,col=NULL,
        ylab=paste("PC",4),xlab="Population",las=1,ylim=c(-1,1),border=brewer.pal(6,"Set1")[4])
for (i in 5:6){
  par(mar=c(0,5,2,1))
  boxplot(plotscore[,i]~annot_col$population,xaxt="n",cex=0.75,lwd=0.75,col=NULL,
          ylab=paste("PC",i),xlab="",las=1,ylim=c(-1,1),border=brewer.pal(6,"Set1")[i])
}
dev.off()


pheatmap(log10(cpm(count_dat_use_sort_filtered[repr.div.gene,])),
         scale="row",show_rownames = F,cluster_cols = F,show_colnames = F,cutree_rows = 2,
         annotation_col = annot_col,annotation_colors = annot_colors,annotation_names_col = F,
         filename = "/Volumes/Temp1/shengkai/CHC/figure3b.png",height = 8.7,width = 8.7,unit="cm",res=600,
         pointsize=8)

avg_repro_gene=apply(cpm(count_dat_use_sort_filtered[repr.div.gene,]),1,function(x) tapply(x,annot_col$population,mean))
dist_repro_genes=dist(log10(avg_repro_gene))

#1-3: 1.43
#1-6: 1.58
#3-6: 1.27

evoavg_repro_gene=apply(cpm(count_dat_use_sort_filtered[repr.div.gene,]),1,function(x) tapply(x,substr(colnames(count_dat_use_sort_filtered),1,1),mean))
logFCavg_repro_gene=log2(evoavg_repro_gene[2,]/evoavg_repro_gene[1,])
evoavg_all=apply(cpm(count_dat_use_sort_filtered[sig_dvg_gene,]),1,function(x) tapply(x,substr(colnames(count_dat_use_sort_filtered),1,1),mean))
logFCavg_all=log2(evoavg_all[2,]/evoavg_all[1,])
wilcox.test(logFCavg_repro_gene,logFCavg_all)#p<0.001
png("/Volumes/Temp1/shengkai/CHC/FigureS6.png",height = 8.7,width = 8.7,units = "cm",res=600,pointsize = 8)
par(mar=c(5,5,2,2))
plot(density(logFCavg_all),xlim=c(-1,1),ylim=c(0,3.5),main="",xlab=expression(paste("Avg. ",log[2],"FC")))
lines(density(logFCavg_repro_gene),col="red")
legend("topleft",legend = c("all diverged genes","diverged reproductive genes"),bty = "n",
       lty=1,col=c("black","red"))
dev.off()

####modularity of reproductive genes####
net_conn=read.table("/Volumes/Temp1/shengkai/common_info/flynet_supervised_0.6.txt",header = F)

set.seed(5487)
random_gene=sample(rownames(pca_goi$rotation),26)
PC1_gene=rownames(pca_goi$rotation)[abs(pca_goi$rotation[,1])>quantile(abs(pca_goi$rotation[,1]),0.9)]
PC2_gene=rownames(pca_goi$rotation)[abs(pca_goi$rotation[,2])>quantile(abs(pca_goi$rotation[,2]),0.9)]

nrow(net_conn[net_conn[,1]%in%random_gene&net_conn[,2]%in%random_gene,])
nrow(net_conn[net_conn[,1]%in%PC1_gene|net_conn[,2]%in%PC1_gene,])
nrow(net_conn[net_conn[,1]%in%PC2_gene|net_conn[,2]%in%PC2_gene,])

unique(net_conn[net_conn[,1]%in%PC1_gene,1])


####Figure1_for_V5####
#phenotypic distribution
pa=rnorm(1e5,1,1)
p1=rnorm(1e5,3,1)
p2=rnorm(1e5,2.9,1)
p3=rnorm(1e5,3.1,1)

png("/Volumes/Temp1/shengkai/CHC/V5_figure1_phenotype.png",width = 8,height = 6,units = "cm",res=600,pointsize = 8)
par(mar=c(4,4,1,2),mfrow=c(2,2))
plot(density(pa),xlab="Phenotypic value",main="",xlim=c(-3,5))
polygon(density(pa),col=alpha("forestgreen",0.2),border = "forestgreen")

plot(density(p1),xlab="Phenotypic value",main="",xlim=c(-3,5))
lines(density(pa),col="forestgreen",lty=2)
polygon(density(p1),col=alpha("salmon",0.2),border = "salmon")

plot(density(p2),xlab="Phenotypic value",main="",xlim=c(-3,5))
lines(density(pa),col="forestgreen",lty=2)
polygon(density(p2),col=alpha("violet",0.2),border = "violet")

plot(density(p3),xlab="Phenotypic value",main="",xlim=c(-3,5))
lines(density(pa),col="forestgreen",lty=2)
polygon(density(p3),col=alpha("firebrick",0.2),border = "firebrick")
dev.off()
#genetic architecture
png("/Volumes/Temp1/shengkai/CHC/V5_figure1_GA.png",width = 6,height = 6,units = "cm",res=600)
par(mar=c(1,1,1,1),mfrow=c(2,1))
plot(1,1,type="n",xlab="",ylab="",main="",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
i=1
while(i<=100){
  segments(x0 = i/100,x1 = i/100,y0 = 0,y1 = 0.2,col = rainbow(125)[i])
  i=i+1
}
text(0.5,0.9,labels = "Phenotype")
polygon(x=c(0.01,1,0.5),y=c(0.2,0.2,0.8),col=alpha("green",0.2),border = NA)
segments(x0 = c(0.01,1),x1 = c(0.5,0.5),y0 = c(0.2,0.2),y1 = c(0.8,0.8))

plot(1,1,type="n",xlab="",ylab="",main="",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
i=1
while(i<=100){
  segments(x0 = i/100,x1 = i/100,y0 = 0,y1 = 0.2,col = rainbow(125)[i])
  i=i+1
}
text(0.5,0.9,labels = "Phenotype")
polygon(x=c(0.01,0.33,0.17),y=c(0.2,0.2,0.5),col=alpha("orange",0.2),border = NA)
polygon(x=c(0.34,0.66,0.5),y=c(0.2,0.2,0.5),col=alpha("forestgreen",0.2),border = NA)
polygon(x=c(0.67,1,0.83),y=c(0.2,0.2,0.5),col=alpha("blue",0.2),border = NA)
segments(x0 = c(0.01,0.33,0.34,0.66,0.67,1,0.17,0.5,0.83),
         x1 = c(0.17,0.17,0.5,0.5,0.83,0.83,0.5,0.5,0.5),
         y0 = c(0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5),
         y1 = c(0.5,0.5,0.5,0.5,0.5,0.5,0.8,0.8,0.8))
dev.off()

#AFC
set.seed(10001)
af0=rep(0.1,100)+runif(100,-0.05,0.05)

af1_red=af0+runif(100,-0.05,0.05)
af2_red=af0+runif(100,-0.05,0.05)
af3_red=af0+runif(100,-0.05,0.05)
idx1=sample(1:100,20)
af1_red[idx1]=rep(0.8,20)+runif(20,-0.15,0.15)
idx2=sample(1:100,20)
af2_red[idx2]=rep(0.8,20)+runif(20,-0.15,0.15)
idx3=sample(1:100,20)
af3_red[idx3]=rep(0.8,20)+runif(20,-0.15,0.15)

af1_ss=af0
af2_ss=af0
af3_ss=af0
idx1_ss=sample(1:33,20)
af1_ss[idx1_ss]=rep(0.8,20)+runif(20,-0.15,0.15)
af1_ss[-idx1_ss]=rep(0.1,80)+runif(20,-0.1,0.1)
idx2_ss=sample(34:66,20)
af2_ss[idx2_ss]=rep(0.8,20)+runif(20,-0.15,0.15)
af2_ss[-idx2_ss]=rep(0.1,80)+runif(20,-0.1,0.1)
idx3_ss=sample(67:100,20)
af3_ss[idx3_ss]=rep(0.8,20)+runif(20,-0.15,0.15)
af3_ss[-idx3_ss]=rep(0.1,80)+runif(20,-0.1,0.1)

png("/Volumes/Temp1/shengkai/CHC/V5_figure1_AFC.png",width = 6,height = 9,units = "cm",pointsize = 8,res = 600)
layout(mat = matrix(1:6,3,2),widths = 1,heights = 1)
par(mar=c(4,4,1,1))
plot(rbind(af0,af1_red)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in sample(1:100,100)){
  points(rbind(af0,af1_red)[,i],lwd=1,col=alpha(rainbow(125)[i],0.6),type="b")
}
plot(rbind(af0,af2_red)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in sample(1:100,100)){
  points(rbind(af0,af2_red)[,i],lwd=1,col=alpha(rainbow(125)[i],0.6),type="b")
}
plot(rbind(af0,af3_red)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in sample(1:100,100)){
  points(rbind(af0,af3_red)[,i],lwd=1,col=alpha(rainbow(125)[i],0.6),type="b")
}
plot(rbind(af0,af1_ss)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in sample(1:100,100)){
  points(rbind(af0,af1_ss)[,i],lwd=1,col=alpha(rainbow(125)[i],0.6),type="b")
}
plot(rbind(af0,af2_ss)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in sample(1:100,100)){
  points(rbind(af0,af2_ss)[,i],lwd=1,col=alpha(rainbow(125)[i],0.6),type="b")
}
plot(rbind(af0,af3_ss)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in sample(1:100,100)){
  points(rbind(af0,af3_ss)[,i],lwd=1,col=alpha(rainbow(125)[i],0.6),type="b")
}
dev.off()

#GA+selection
png("/Volumes/Temp1/shengkai/CHC/V5_figure1_GA_S.png",width = 6,height = 9,units = "cm",pointsize = 8,res = 600)
layout(mat = matrix(1:6,3,2),widths = 1,heights = 1)
par(mar=c(0.5,0.5,0.5,0.5))
set.seed(10001)
plot(1,1,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab = "")
x.cord=rep(seq(0,1,length.out = 10),10)
y.cord=rep(seq(0,1,length.out = 10),each=10)
rf=sample(1:100,20)
points(x.cord,y.cord,col=ifelse(1:100%in%rf,"orange","grey"),pch=19)

plot(1,1,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab = "")
x.cord=rep(seq(0,1,length.out = 10),10)
y.cord=rep(seq(0,1,length.out = 10),each=10)
rf=sample(1:100,20)
points(x.cord,y.cord,col=ifelse(1:100%in%rf,"orange","grey"),pch=19)

plot(1,1,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab = "")
x.cord=rep(seq(0,1,length.out = 10),10)
y.cord=rep(seq(0,1,length.out = 10),each=10)
rf=sample(1:100,20)
points(x.cord,y.cord,col=ifelse(1:100%in%rf,"orange","grey"),pch=19)

plot(1,1,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab = "")
x.cord=c(rep(seq(0,1,length.out = 10)[1:5],10),rep(seq(0,1,length.out = 10)[6:10],10))
y.cord=c(rep(seq(0,1,length.out = 10),each=5),rep(seq(0,1,length.out = 10),each=5))
rf=sample(1:25,20)
points(x.cord,y.cord,col=ifelse(1:100%in%rf,"orange","grey"),pch=rep(c(19,17,15,8),each=25))

plot(1,1,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab = "")
x.cord=c(rep(seq(0,1,length.out = 10)[1:5],10),rep(seq(0,1,length.out = 10)[6:10],10))
y.cord=c(rep(seq(0,1,length.out = 10),each=5),rep(seq(0,1,length.out = 10),each=5))
rf=sample(26:50,20)
points(x.cord,y.cord,col=ifelse(1:100%in%rf,"orange","grey"),pch=rep(c(19,17,15,8),each=25))

plot(1,1,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab = "")
x.cord=c(rep(seq(0,1,length.out = 10)[1:5],10),rep(seq(0,1,length.out = 10)[6:10],10))
y.cord=c(rep(seq(0,1,length.out = 10),each=5),rep(seq(0,1,length.out = 10),each=5))
rf=sample(76:100,20)
points(x.cord,y.cord,col=ifelse(1:100%in%rf,"orange","grey"),pch=rep(c(19,17,15,8),each=25))
dev.off()

####Figure3a_v1####
set.seed(100)
af1=rep(0.1,8)+runif(8,-0.05,0.05)
pie_tab_af1=matrix(c(af1,1-af1),2,8,byrow = T)
colnames(pie_tab_af1)=paste("locus",1:8)
rownames(pie_tab_af1)=paste("Allele",1:2)
pie_tab_gg_af1=reshape2::melt(pie_tab_af1)

af2=af1
af2[sample(1:8,4)]=1
pie_tab_af2=matrix(c(af2,1-af2),2,8,byrow = T)
colnames(pie_tab_af2)=paste("locus",1:8)
rownames(pie_tab_af2)=paste("Allele",1:2)
pie_tab_gg_af2=reshape2::melt(pie_tab_af2)

af3=af1
af3[sample(1:8,4)]=1
pie_tab_af3=matrix(c(af3,1-af3),2,8,byrow = T)
colnames(pie_tab_af3)=paste("locus",1:8)
rownames(pie_tab_af3)=paste("Allele",1:2)
pie_tab_gg_af3=reshape2::melt(pie_tab_af3)

af4=af1
af4[1:4]=1
pie_tab_af4=matrix(c(af4,1-af4),2,8,byrow = T)
colnames(pie_tab_af4)=paste("locus",1:8)
rownames(pie_tab_af4)=paste("Allele",1:2)
pie_tab_gg_af4=reshape2::melt(pie_tab_af4)

af5=af1
af5[5:8]=1
pie_tab_af5=matrix(c(af5,1-af5),2,8,byrow = T)
colnames(pie_tab_af5)=paste("locus",1:8)
rownames(pie_tab_af5)=paste("Allele",1:2)
pie_tab_gg_af5=reshape2::melt(pie_tab_af5)


png("/Volumes/Temp1/shengkai/CHC/Figure3a_1.png",width = 6,height = 4,units = "cm",res = 600,pointsize = 6,bg="transparent")
gg=ggplot(pie_tab_gg_af1,aes(x="",y=value,fill=Var1))+
  geom_bar(width = 1,stat = "identity")+
  coord_polar("y", start=0)+
  facet_wrap(.~Var2,nrow = 2)+
  scale_fill_manual(values=c("orange","grey"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),
        panel.grid=element_blank(),axis.title = element_blank(),
        legend.title = element_blank(),strip.text.x = element_text(size = 6),
        legend.text=element_text(size=6),legend.position = "none",
        plot.background = element_rect(fill="transparent",colour=NA))
g=ggplot_gtable(ggplot_build(gg))
stript=which(grepl('strip-t', g$layout$name))
fills=brewer.pal(8,"Set2")
k=1
for (i in stript) {
  j=which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill=fills[k]
  k=k+1
}
grid.draw(g)
dev.off()

png("/Volumes/Temp1/shengkai/CHC/Figure3a_2.png",width = 6,height = 4,units = "cm",res = 600,pointsize = 6,bg="transparent")
gg=ggplot(pie_tab_gg_af2,aes(x="",y=value,fill=Var1))+
  geom_bar(width = 1,stat = "identity")+
  coord_polar("y", start=0)+
  facet_wrap(.~Var2,nrow = 2)+
  scale_fill_manual(values=c("orange","grey"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),
        panel.grid=element_blank(),axis.title = element_blank(),
        legend.title = element_blank(),strip.text.x = element_text(size = 6),
        legend.text=element_text(size=6),legend.position = "none",
        plot.background = element_rect(fill="transparent",colour=NA))
g=ggplot_gtable(ggplot_build(gg))
stript=which(grepl('strip-t', g$layout$name))
fills=brewer.pal(8,"Set2")
k=1
for (i in stript) {
  j=which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill=fills[k]
  k=k+1
}
grid.draw(g)
dev.off()

png("/Volumes/Temp1/shengkai/CHC/Figure3a_3.png",width = 6,height = 4,units = "cm",res = 600,pointsize = 6,bg="transparent")
gg=ggplot(pie_tab_gg_af3,aes(x="",y=value,fill=Var1))+
  geom_bar(width = 1,stat = "identity")+
  coord_polar("y", start=0)+
  facet_wrap(.~Var2,nrow = 2)+
  scale_fill_manual(values=c("orange","grey"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),
        panel.grid=element_blank(),axis.title = element_blank(),
        legend.title = element_blank(),strip.text.x = element_text(size = 6),
        legend.text=element_text(size=6),legend.position = "none",
        plot.background = element_rect(fill="transparent",colour=NA))
g=ggplot_gtable(ggplot_build(gg))
stript=which(grepl('strip-t', g$layout$name))
fills=brewer.pal(8,"Set2")
k=1
for (i in stript) {
  j=which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill=fills[k]
  k=k+1
}
grid.draw(g)
dev.off()

png("/Volumes/Temp1/shengkai/CHC/Figure3a_4.png",width = 6,height = 4,units = "cm",res = 600,pointsize = 6,bg="transparent")
gg=ggplot(pie_tab_gg_af4,aes(x="",y=value,fill=Var1))+
  geom_bar(width = 1,stat = "identity")+
  coord_polar("y", start=0)+
  facet_wrap(.~Var2,nrow = 2)+
  scale_fill_manual(values=c("orange","grey"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),
        panel.grid=element_blank(),axis.title = element_blank(),
        legend.title = element_blank(),strip.text.x = element_text(size = 6),
        legend.text=element_text(size=6),legend.position = "none",
        plot.background = element_rect(fill="transparent",colour=NA))
g=ggplot_gtable(ggplot_build(gg))
stript=which(grepl('strip-t', g$layout$name))
fills=brewer.pal(8,"Set2")
k=1
for (i in stript) {
  j=which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill=fills[k]
  k=k+1
}
grid.draw(g)
dev.off()

png("/Volumes/Temp1/shengkai/CHC/Figure3a_5.png",width = 6,height = 4,units = "cm",res = 600,pointsize = 6,bg="transparent")
gg=ggplot(pie_tab_gg_af5,aes(x="",y=value,fill=Var1))+
  geom_bar(width = 1,stat = "identity")+
  coord_polar("y", start=0)+
  facet_wrap(.~Var2,nrow = 2)+
  scale_fill_manual(values=c("orange","grey"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),
        panel.grid=element_blank(),axis.title = element_blank(),
        legend.title = element_blank(),strip.text.x = element_text(size = 6),
        legend.text=element_text(size=6),legend.position = "none",
        plot.background = element_rect(fill="transparent",colour=NA))
g=ggplot_gtable(ggplot_build(gg))
stript=which(grepl('strip-t', g$layout$name))
fills=brewer.pal(8,"Set2")
k=1
for (i in stript) {
  j=which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill=fills[k]
  k=k+1
}
grid.draw(g)
dev.off()

new_list=list("Sig. divergent"=sig_dvg_gene,"N.S."=setdiff(rownames(y_dvg),sig_dvg_gene))
new_list2=list("ACPs (13)"=GOI_ACP,"SFPs (21)"=GOI_SFP)
pie_tab=sapply(new_list2,function(x) sapply(new_list,function(y) sum(x%in%y)))
pie_tab_gg=reshape2::melt(t(pie_tab)/colSums(pie_tab))

png("/Volumes/Temp1/shengkai/CHC/Figure3c.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 8)
gg=ggplot(pie_tab_gg,aes(x="",y=value,fill=Var2))+
  geom_bar(width = 1,stat = "identity")+
  coord_polar("y", start=0)+
  facet_wrap(.~Var1,nrow = 1)+
  scale_fill_manual(values=c("orange","grey"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),
        panel.grid=element_blank(),axis.title = element_blank(),
        legend.title = element_blank(),strip.text.x = element_text(size = 8),
        legend.text=element_text(size=8),legend.position = "bottom")
g=ggplot_gtable(ggplot_build(gg))
stript=which(grepl('strip-t', g$layout$name))
fills=brewer.pal(12,"Paired")[c(3,5)]
k=1
for (i in stript) {
  j=which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill=fills[k]
  k=k+1
}
grid.draw(g)
dev.off()

pheatmap(t(pca_goi$x)[1:9,],
         scale="none",cluster_cols = F,cluster_rows = F,
         show_colnames = F,annotation_col = annot_col,annotation_colors = annot_colors,
         annotation_names_col = F,labels_row = paste0("PC",1:9," (",round(ve_goi[1:9],3)*100,"%)"))
         

####Figure3a_v2####
set.seed(100)
af0=rep(0.1,20)+runif(20,-0.05,0.05)

af1_red=af0+runif(20,-0.05,0.05)
af2_red=af0+runif(20,-0.05,0.05)
idx1=sample(1:20,10)
af1_red[idx1]=rep(0.9,10)+runif(10,-0.05,0.05)
idx2=sample(1:20,10)
af2_red[idx2]=rep(0.9,10)+runif(10,-0.05,0.05)

png("Figure3a_1_new.png",width = 8.7,height = 5,units = "cm",pointsize = 8,res = 600)
par(mar=c(5,4,3,1),mfrow=c(1,2))
plot(rbind(af0,af1_red)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     main="R1",xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in 1:20){
  points(rbind(af0,af1_red)[,i],lwd=2,col=alpha(ifelse(i>10,rainbow(200)[i],rainbow(200)[120+i]),0.5),type="b")
}
plot(rbind(af0,af2_red)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     main="R2",xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in 1:20){
  points(rbind(af0,af2_red)[,i],lwd=2,col=alpha(ifelse(i>10,rainbow(200)[i],rainbow(200)[120+i]),0.5),type="b")
}
dev.off()

af1_ss=af0
af2_ss=af0
af1_ss[1:10]=rep(0.9,10)+runif(10,-0.05,0.05)
af1_ss[11:20]=rep(0.1,10)+runif(10,-0.1,0.1)
af2_ss[1:10]=rep(0.1,10)+runif(10,-0.1,0.1)
af2_ss[11:20]=rep(0.9,10)+runif(10,-0.05,0.05)

png("Figure3a_2_new.png",width = 8.7,height = 5,units = "cm",pointsize = 8,res = 600)
par(mar=c(5,4,3,1),mfrow=c(1,2))
plot(rbind(af0,af1_ss)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     main="R1",xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in 1:20){
  points(rbind(af0,af1_ss)[,i],lwd=2,col=alpha(ifelse(i>10,rainbow(200)[i],rainbow(200)[120+i]),0.5),type="b")
}
plot(rbind(af0,af2_ss)[,1],xlim=c(0.5,2.5),ylim=c(0,1),xlab="Generation",ylab="Frequency",
     main="R2",xaxt="n",type="n")
axis(1,at=1:2,labels = c(0,100))
for (i in 1:20){
  points(rbind(af0,af2_ss)[,i],lwd=2,col=alpha(ifelse(i>10,rainbow(200)[i],rainbow(200)[120+i]),0.5),type="b")
}
dev.off()


####discusssion parallelism####
mean_all_rep_m=apply(cpm(y_corrected)[,substr(group_gene,2,2)=="m"],1,function(x) tapply(x,c(rep(19,5),rep(1:10,each=3)),mean))
logfc_all_rep_m=log2(t(mean_all_rep_m)/mean_all_rep_m[11,])[,1:10]
mean_GOI_CHC_rep_m=apply(cpm(y_corrected)[rownames(y_corrected)%in%GOI_CHC,substr(group_gene,2,2)=="m"],1,function(x) tapply(x,c(rep(19,5),rep(1:10,each=3)),mean))
logfc_GOI_CHC_rep_m=log2(t(mean_GOI_CHC_rep_m)/mean_GOI_CHC_rep_m[11,])[,1:10]
mean_GOI_reproduction_rep_m=apply(cpm(y_corrected)[rownames(y_corrected)%in%intersect(sig_dvg_gene,genesInTerm(tgd,"GO:0032504")[[1]]),substr(group_gene,2,2)=="m"],1,function(x) tapply(x,c(rep(19,5),rep(1:10,each=3)),mean))
logfc_GOI_reproduction_rep_m=log2(t(mean_GOI_reproduction_rep_m)/mean_GOI_reproduction_rep_m[11,])[,1:10]
mean_GOI_ATP_rep_m=apply(cpm(y_corrected)[rownames(y_corrected)%in%genesInTerm(tgd,"GO:0006099")[[1]],substr(group_gene,2,2)=="m"],1,function(x) tapply(x,c(rep(19,5),rep(1:10,each=3)),mean))
logfc_GOI_ATP_rep_m=log2(t(mean_GOI_ATP_rep_m)/mean_GOI_ATP_rep_m[11,])[,1:10]
mean_GOI_signal_rep_m=apply(cpm(y_corrected)[rownames(y_corrected)%in%genesInTerm(tgd,"GO:0007165")[[1]],substr(group_gene,2,2)=="m"],1,function(x) tapply(x,c(rep(19,5),rep(1:10,each=3)),mean))
logfc_GOI_signal_rep_m=log2(t(mean_GOI_signal_rep_m)/mean_GOI_signal_rep_m[11,])[,1:10]


png("/Volumes/Temp1/shengkai/CHC/comparison_parallelism_1.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 8)
png("/Volumes/Temp1/shengkai/CHC/comparison_parallelism_2.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 8)
plot(density(cor(logfc_all_rep_m)[lower.tri(cor(logfc_all_rep_m))]),main="", xlab="similarity",xlim=c(0.2,1))
lines(density(cor(logfc_GOI_CHC_rep_m)[lower.tri(cor(logfc_all_rep_m))]),col="red")
lines(density(cor(logfc_GOI_ATP_rep_m)[lower.tri(cor(logfc_all_rep_m))]),col="orange")
lines(density(cor(logfc_GOI_reproduction_rep_m)[lower.tri(cor(logfc_all_rep_m))]),col="blue")
lines(density(cor(logfc_GOI_signal_rep_m)[lower.tri(cor(logfc_all_rep_m))]),col="green")
dev.off()

png("/Volumes/Temp1/shengkai/CHC/Figure3c_new.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 8)
plot_dat=data.frame(h=c(1-cor(logfc_all_rep_m)[lower.tri(cor(logfc_all_rep_m))],
                        1-cor(logfc_GOI_reproduction_rep_m)[lower.tri(cor(logfc_all_rep_m))],
                        1-cor(logfc_GOI_CHC_rep_m)[lower.tri(cor(logfc_all_rep_m))]),
                    cat=as.factor(c(rep("Background",45),rep("Reproductive",45),rep("CHC-related",45))))
qplot(cat,h,data = plot_dat,geom = "violin",fill=cat,alpha=I(0.5),xlab = "",ylab="Heterogeneity across replicates")+
  theme_classic()+
  theme(axis.text.x=element_text(size=8),legend.position = "NULL")+
  geom_boxplot(width=0.1,fill="white")+
  scale_fill_manual(values=c("grey","salmon","cyan"))
dev.off()

with(plot_dat,anova(lm(h~cat)))
with(plot_dat,wilcox.test(h[cat=="Background"],h[cat=="Reproductive"]))
####quantification by openChrom####
setwd("/Volumes/Temp1/shengkai/CHC/2018_Drosophila_analyzed_v2/quantification/")
RI_calculation=function(x,alkane){
  x=as.numeric(x)
  if(x<min(alkane$RT)) return(0)
  else if(x>max(alkane$RT)) return(3001)
  else if(any(alkane$RT%in%x)) return(100*alkane$nC[alkane$RT%in%x])
  else{
    left_index=max(which(x>alkane$RT))
    right_index=min(which(x<alkane$RT))
    RI=100*(alkane$nC[left_index]+(log10(x)-log10(alkane$RT[left_index]))/(log10(alkane$RT[right_index])-log10(alkane$RT[left_index])))
    return(RI)}
}

trt_tab=read.delim("../../sample_list.txt",header = T,stringsAsFactors = F)
alkane=data.frame(nC=21:30,RT=c(10.727,11.701,12.522,13.220,13.830,14.374,14.910,15.454,16.023,16.624))
file_list=list.files("./",pattern = "CHC") 
dat=lapply(file_list,function(x) read.delim(x,header = T,stringsAsFactors = F)[,-1])
names(dat)=strsplit2(file_list,".txt")[,1]

dat_rename=lapply(dat,function(x){colnames(x)[1]="RT"; return(x)})
dat_rename=lapply(dat_rename,function(x){x[,13]="unknown"; return(x)})
dat_narm=lapply(dat_rename,na.omit)
dat_filter=lapply(dat_narm,function(x) x[as.numeric(x$S.N)>3|is.na(as.numeric(x$S.N)),])
dat_RI=lapply(dat_filter,function(x) {
  tmp=apply(x,1,function(y) {tmp=RI_calculation(y[1],alkane)})
  x$RI=tmp
  return(x)
})
dat_poi=lapply(dat_RI,function(x){
  x[x$RI>=2097&x$RI<=2103,13]="n-C21"
  x[x$RI>=2183.5&x$RI<=2189.5,13]="C22:1?"
  x[x$RI>=2191.5&x$RI<=2197.5,13]="c-VA"
  x[x$RI>=2197&x$RI<=2203,13]="n-C22"
  x[x$RI>=2276&x$RI<=2282,13]="9-T"
  x[x$RI>=2282&x$RI<=2288,13]="7-T"
  x[x$RI>=2291&x$RI<=2297,13]="5-T"
  x[x$RI>=2297&x$RI<=2303,13]="n-C23"
  x[x$RI>=2380&x$RI<=2386,13]="6-C24:1?"
  x[x$RI>=2388&x$RI<=2394,13]="5-C24:1?"
  x[x$RI>=2397&x$RI<=2403,13]="n-C24"
  x[x$RI>=2462&x$RI<=2468,13]="7,11-PD"
  x[x$RI>=2475&x$RI<=2481,13]="9-P"
  x[x$RI>=2483&x$RI<=2489,13]="7-P"
  x[x$RI>=2497&x$RI<=2503,13]="n-C25"
  x[x$RI>=2662&x$RI<=2668,13]="7,11-HD"
  x[x$RI>=2697&x$RI<=2703,13]="n-C27"
  x[x$RI>=2860.3&x$RI<=2866.3,13]="7,11-ND"
  x[x$RI>=2997&x$RI<=3003,13]="IS"
  return(x)
  })
sapply(dat_poi,function(x) unique(x$Name))

quantity_matrix=sapply(dat_poi,function(x) {
  tmp=c()
  for (i in c("n-C21","C22:1?","c-VA","n-C22","9-T","7-T","5-T","n-C23","6-C24:1?","5-C24:1?","n-C24","7,11-PD","9-P","7-P","n-C25","7,11-HD","n-C27","7,11-ND","IS")){
    if(i%in%x$Name) tmp=c(tmp,as.numeric(max(x$Area[x$Name%in%i])))
    else tmp=c(tmp,0)
  }
  return(tmp)
})
rownames(quantity_matrix)=c("n-C21","C22:1?","c-VA","n-C22","9-T","7-T","5-T","n-C23","6-C24:1?","5-C24:1?","n-C24","7,11-PD","9-P","7-P","n-C25","7,11-HD","n-C27","7,11-ND","IS")
conc_matrix=t(quantity_matrix)/quantity_matrix["IS",]*0.1  #0.1mg/ml
relative=function(x) x/sum(x)
CLR=function(x) {
  tmp=c()
  y=na.omit(x)
  y=y[y>0]
  tmp[x>0&!is.na(x)]=log(y/prod(y)^(1/length(y)))
  tmp[x==0|is.na(x)]=NA
  return(tmp)
}
relative_matrix=apply(conc_matrix[,c(2,4:18)],1,relative)
CLR_matrix=apply(relative_matrix,2,CLR)
rownames(CLR_matrix)=colnames(conc_matrix)[c(2,4:18)]
CLR_matrix_no_na=na.omit(CLR_matrix)
####overlay####
setwd("/Volumes/Temp1/shengkai/CHC/")
Bm=read.csv("B1m-1.csv",stringsAsFactors = F)
chromatograph_m=Bm[2:1097,]
baseline_m=Bm[-1:-1098,]
plot(as.numeric(chromatograph_m[,1]),chromatograph_m[,2]-baseline_m[,2],type="l",xlab="Retention time (min)",ylab="Signal intensity",ylim=c(0,400000))
abline(h=max(0.025*(chromatograph_m[,2]-baseline_m[,2])),lty=2)
#points(as.numeric(baseline[,1])/1000/60,2*baseline[,2],type="l",lty=3)

plot(as.numeric(chromatograph_m[,1]),log(chromatograph_m[,2])-log(baseline_m[,2]),type="l",xlab="Retention time (min)",ylab="Log Intensity")
#abline(h=log(max(0.05*(chromatograph_m[,2]-baseline_m[,2]))),lty=2)
abline(h=log(2),lty=3)

Bf=read.csv("B1f-1.csv",stringsAsFactors = F)
chromatograph_f=Bf[2:1097,]
baseline_f=Bf[-1:-1098,]
plot(as.numeric(chromatograph_f[,1]),chromatograph_f[,2]-baseline_f[,2],type="l",xlab="Retention time (min)",ylab="Signal intensity",ylim=c(0,400000))
abline(h=max(0.025*(chromatograph_f[,2]-baseline_f[,2])),lty=2)
#points(as.numeric(baseline[,1])/1000/60,2*baseline[,2],type="l",lty=3)

plot(as.numeric(chromatograph_f[,1]),log(chromatograph_f[,2])-log(baseline_f[,2]),type="l",xlab="Retention time (min)",ylab="Log Intensity")
#abline(h=log(max(0.05*(chromatograph_f[,2]-baseline_f[,2]))),lty=2)
abline(h=log(2),lty=3)

png("overlay.png",height = 8,width = 16,unit="cm",res = 600,pointsize = 6)
plot(as.numeric(chromatograph_m[,1]),(chromatograph_m[,2]-baseline_m[,2])/10^3,type="l",las=1,
     xlab="Retention time (min)",ylab="Signal intensity",xlim=c(11,17),ylim=c(-300,300),col="royalblue",lwd=1.5)
points(as.numeric(chromatograph_f[,1]),-(chromatograph_f[,2]-baseline_f[,2])/10^3,type="l",col="salmon",lwd=1.5)
legend("topright", legend=c("male","female"),col=c("royalblue","salmon"), lty=1,lwd=2,bty="n")
text(dat_poi$CHC066$RT[-c(16,20)],
     (chromatograph_m[,2]-baseline_m[,2])[dat_poi$CHC066[-c(16,20),7]]/10^3,
     1:length(dat_poi$CHC066$RT[-c(16,20)]),pos = 3)
text(12.450,300,5)
#abline(h=max(0.025*(chromatograph_m[,2]-baseline_m[,2])),lty=2,col="royalblue")
#abline(h=-max(0.025*(chromatograph_f[,2]-baseline_f[,2])),lty=2,col="salmon")
dev.off()

####modeling hot-base only####
evo=substr(trt_tab$pop[-c(16,65,112,113,114,115)],1,1)
pop=trt_tab$pop[-c(16,65,112,113,114,115)]
sex=trt_tab$sex[-c(16,65,112,113,114,115)]
conc_matrix_HB=conc_matrix[evo%in%c("B","H"),]
relative_matrix_HB=t(relative_matrix[,evo%in%c("B","H")])
CLR_matrix_HB=t(CLR_matrix[,evo%in%c("B","H")])

avg_relative_matrix_HB=apply(relative_matrix_HB,2,function(x) tapply(x,paste(evo,sex)[evo%in%c("B","H")],mean))
se_relative_matrix_HB=apply(relative_matrix_HB,2,function(x) tapply(x,paste(evo,sex)[evo%in%c("B","H")],function(a) sd(a)/sqrt(length(a))))
logFC_relative_matrix_HB=log2(avg_relative_matrix_HB[c(3,4),]/avg_relative_matrix_HB[c(1,2),])
pheatmap(t(logFC_relative_matrix_HB),cellwidth = 20,cellheight = 20,
        cluster_rows = F,cluster_cols = F,labels_col = c("female","male"),
        color = myColor,breaks = myBreaks,
        filename = "/Volumes/Temp1/shengkai/CHC/Figure1c_CHC.png",res=600)

pca_HB=prcomp(t(na.omit(t(CLR_matrix_HB))))
ve=pca_HB$sdev^2/sum(pca_HB$sdev^2)

png("/Volumes/Temp1/shengkai/CHC/PCA_HB_12_consistent_peaks.png",width = 8.7,height = 8.7,units = "cm",res=600,pointsize = 8)
plot(pca_HB$x,col=ifelse(evo[colnames(CLR_matrix)%in%colnames(na.omit(t(CLR_matrix_HB)))]=="H","salmon","forestgreen"),
     pch=as.numeric(as.factor(pop[colnames(CLR_matrix)%in%colnames(na.omit(t(CLR_matrix_HB)))])),
     xlab=paste0("PC1 (",round(ve[1],4)*100,"%)"),ylab=paste0("PC2 (",round(ve[2],4)*100,"%)"),asp=1)
dev.off()
png("/Volumes/Temp1/shengkai/CHC/PCA_HB_12_consistent_peaks_simple_label.png",width = 8.7,height = 8.7,units = "cm",res=600,pointsize = 8)
plot(pca_HB$x,pch=ifelse(evo[colnames(CLR_matrix)%in%colnames(na.omit(t(CLR_matrix_HB)))]=="H",19,1),
     col=ifelse(sex[colnames(CLR_matrix)%in%colnames(na.omit(t(CLR_matrix_HB)))]=="m","purple","orange"),
     xlab=paste0("PC1 (",round(ve[1],4)*100,"%)"),ylab=paste0("PC2 (",round(ve[2],4)*100,"%)"),asp=1,ylim=c(-1,1))
legend("topleft",legend=c("Anc. M.","Anc. F.", "Evo. M.","Evo. F."),col=c("purple","orange","purple","orange"),
       pch=c(1,1,19,19),bty="n")
dev.off()

#loading plot
png("/Volumes/Temp1/shengkai/CHC/PCA_HB_12_consistent_peaks_loading.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mar=c(5,6,2,1),bg="NA")
barplot(t(pca_HB$rotation[rev(c(2,5,10,1,3,6,4,7,9,8,11,12)),1:2]),beside = T,
        legend.text = c("PC1 (F./M.)","PC2 (Anc./Evo.)"),
        horiz = T,las=1,xlab="loadings",args.legend = list(x="bottomright"),xlim=c(-0.6,0.6),
        names.arg = rev(c("n-C22","n-C23","n-C25","C22:1","7-T","9-T","6-C24:1","5-C24:1","7-P","7,11-PD","7,11-HD","7,11-ND")))
dev.off()

p_evo_HB_aov=apply(CLR_matrix_HB[,],2,function(x) anova(lm(x~sex[evo%in%c("B","H")]*evo[evo%in%c("B","H")]))[1:3,5])
p_pop_HB_aov=apply(CLR_matrix_HB[evo[evo%in%c("B","H")]%in%"H",],2,function(x) anova(lm(x~sex[evo%in%c("H")]*pop[evo%in%c("H")]))[1:3,5])
write.csv(apply(p_evo_HB_aov,1,function(x) p.adjust(x,"BH")),"/Volumes/Temp1/shengkai/CHC/anova_HB_model1.csv")
write.csv(apply(p_pop_HB_aov,1,function(x) p.adjust(x,"BH")),"/Volumes/Temp1/shengkai/CHC/anova_HB_model2.csv")

anova(lm(pca_HB$x[,1]~sex[evo%in%c("B","H")]*evo[evo%in%c("B","H")]))
anova(lm(pca_HB$x[,2]~sex[evo%in%c("B","H")]*evo[evo%in%c("B","H")]))

mean_PC_score=apply(pca_HB$x,2,function(x) tapply(x,group[evo%in%c("B","H")],mean))
SE_PC_score=apply(pca_HB$x,2,function(x) tapply(x,group[evo%in%c("B","H")],function(y) sd(y)/sqrt(length(y))))

png("/Volumes/Temp1/shengkai/CHC/mean_PC_scores_HB.png",height = 15.1,width = 15.1,units = "cm",res=600,pointsize = 8)
par(mfrow=c(2,2),mar=c(4.5,4.5,1.5,1.5))
for (i in 1:4){
  plot(NA,type="n",xlab="population",ylab=paste0("mean PC", i, " score (", round(ve[i],4)*100,"%)"),xlim=c(0.5,2.5),
       ylim=range(mean_PC_score[,i]+1.1*c(-SE_PC_score[,i],SE_PC_score[,i])),xaxt="n",cex.axis=1.25)
  axis(1,at=1:2,labels = c("Anc.","Evo."))
  points(1:2,mean_PC_score[c(1,3),i],type = "b",pch=19,col="salmon")
  points(1:2,mean_PC_score[c(2,4),i],type = "b",pch=19,col="royalblue")
  arrows(c(1,2,1,2),mean_PC_score[c(1,3),i]-SE_PC_score[c(1,3),i],c(1,2,1,2),mean_PC_score[c(1,3),i]+SE_PC_score[c(1,3),i],
         angle=90,length=0.05,code=3,col="salmon")
  arrows(c(1,2,1,2),mean_PC_score[c(2,4),i]-SE_PC_score[c(2,4),i],c(1,2,1,2),mean_PC_score[c(2,4),i]+SE_PC_score[c(2,4),i],
         angle=90,length=0.05,code=3,col="royalblue")
}
dev.off()

png("/Volumes/Temp1/shengkai/CHC/mean_PC_scores_HB_PC1&2.png",height = 8.7,width = 8.7,units = "cm",res=600,pointsize = 8)
par(mfrow=c(1,2),mar=c(4.5,4.5,1.5,1.5))
for (i in 1:2){
  plot(NA,type="n",xlab="population",ylab=paste0("PC", i, " score (", round(ve[i],4)*100,"%)"),xlim=c(0.5,2.5),
       ylim=range(mean_PC_score[,i]+1.1*c(-SE_PC_score[,i],SE_PC_score[,i])),xaxt="n",cex.axis=1)
  axis(1,at=1:2,labels = c("Anc.","Evo."))
  points(1:2,mean_PC_score[c(1,3),i],type = "b",pch=19,col="salmon")
  points(1:2,mean_PC_score[c(2,4),i],type = "b",pch=19,col="royalblue")
  arrows(c(1,2,1,2),mean_PC_score[c(1,3),i]-SE_PC_score[c(1,3),i],c(1,2,1,2),mean_PC_score[c(1,3),i]+SE_PC_score[c(1,3),i],
         angle=90,length=0.05,code=3,col="salmon")
  arrows(c(1,2,1,2),mean_PC_score[c(2,4),i]-SE_PC_score[c(2,4),i],c(1,2,1,2),mean_PC_score[c(2,4),i]+SE_PC_score[c(2,4),i],
         angle=90,length=0.05,code=3,col="royalblue")
}
dev.off()

mean_PC_score[3:4,2]/mean_PC_score[1:2,2]

png("/Volumes/Temp1/shengkai/CHC/boxplot_model1_all.png",height = 17.4,width = 17.4,units = "cm",res=600,pointsize = 8)
par(mfrow=c(4,4))
for(i in c("n-C22","n-C23","n-C24","n-C25","n-C27","C22:1?","9-T","7-T","5-T","6-C24:1?","5-C24:1?","9-P","7-P","7,11-PD","7,11-HD","7,11-ND")){
  boxplot(CLR_matrix_HB[,i]~evo[evo%in%c("B","H")]+sex[evo%in%c("B","H")],las=1,main=i,ylab="CLR-transformed Rel. Conc.")
}
dev.off()
  
png("/Volumes/Temp1/shengkai/CHC/boxplot_model1_alkane.png",height = 4.4,width = 17.4,units = "cm",res=600,pointsize = 8)
par(mfrow=c(1,5))
for(i in c("n-C22","n-C23","n-C24","n-C25","n-C27")){
  boxplot(CLR_matrix_HB[,i]~evo[evo%in%c("B","H")]+sex[evo%in%c("B","H")],las=1,main=i,ylab="CLR-transformed Rel. Conc.")
}
dev.off()

png("/Volumes/Temp1/shengkai/CHC/boxplot_model1_monoene.png",height = 8.7,width = 17.4,units = "cm",res=600,pointsize = 8)
par(mfrow=c(2,4))
for(i in c("C22:1?","9-T","7-T","5-T","6-C24:1?","5-C24:1?","9-P","7-P")){
  boxplot(CLR_matrix_HB[,i]~evo[evo%in%c("B","H")]+sex[evo%in%c("B","H")],las=1,main=i,ylab="CLR-transformed Rel. Conc.")
}
dev.off()

png("/Volumes/Temp1/shengkai/CHC/boxplot_model1_diene.png",height = 4.4,width = 13.05,units = "cm",res=600,pointsize = 8)
par(mfrow=c(1,3))
for(i in c("7,11-PD","7,11-HD","7,11-ND")){
  boxplot(CLR_matrix_HB[,i]~evo[evo%in%c("B","H")]+sex[evo%in%c("B","H")],las=1,main=i,ylab="CLR-transformed Rel. Conc.")
}
dev.off()

png("/Volumes/Temp1/shengkai/CHC/boxplot_model2_all.png",height = 17.4,width = 17.4,units = "cm",res=600,pointsize = 8)
par(mfrow=c(4,4))
for(i in c("n-C22","n-C23","n-C24","n-C25","n-C27","C22:1?","9-T","7-T","5-T","6-C24:1?","5-C24:1?","9-P","7-P","7,11-PD","7,11-HD","7,11-ND")){
  boxplot(CLR_matrix_HB[evo[evo%in%c("B","H")]%in%c("H"),i]~pop[evo%in%c("H")]+sex[evo%in%c("H")],las=1,main=i,ylab="CLR-transformed Rel. Conc.")
}
dev.off()

png("/Volumes/Temp1/shengkai/CHC/boxplot_model2_7-P.png",height = 8.7,width = 17.4,units = "cm",res=600,pointsize = 8)
boxplot(CLR_matrix_HB[evo[evo%in%c("B","H")]%in%c("H"),"7-P"]~sex[evo%in%c("H")]+pop[evo%in%c("H")],las=1,main="7-P",ylab="CLR-transformed Rel. Conc.")
dev.off()

fm_mm=lmer(CLR_matrix_HB[,8]~sex[evo%in%c("B","H")]+evo[evo%in%c("B","H")]+(evo[evo%in%c("B","H")]|pop[evo%in%c("B","H")]))

####mating assay####
setwd("/Volumes/Temp1/shengkai/mate_choice_CGE/")
####experiment: Base v.s. Hot evolved ####
dat_BH=read.csv("./mating_assay_data_HB.csv",header = T,stringsAsFactors = F)
dat_BH$Y_idx=(sqrt(dat_BH[,4]*dat_BH[,5])-sqrt(dat_BH[,6]*dat_BH[,7]))/(sqrt(dat_BH[,4]*dat_BH[,5])+sqrt(dat_BH[,6]*dat_BH[,7]))
dat_BH$pHm=(dat_BH$HH+dat_BH$HB)/rowSums(dat_BH[,4:7])
dat_BH$pHf=(dat_BH$BH+dat_BH$HH)/rowSums(dat_BH[,4:7])

fm=with(dat_BH,lm(Y_idx~pair+coloring))
anova(fm)
tapply(dat_BH$Y_idx,dat_BH$pair,function(x) t.test(x,mu = 0))
avg_y_BH=tapply(dat_BH$Y_idx,dat_BH$pair,mean)
se_y_BH=tapply(dat_BH$Y_idx,dat_BH$pair,function(x) sd(x)/sqrt(length(x)))
t.test(dat_BH$Y_idx,mu=0)

fm_pHm=with(dat_BH,lm(pHm~pair+coloring))
anova(fm_pHm)
fm_pHf=with(dat_BH,lm(pHf~pair+coloring))
anova(fm_pHf)
tapply(dat_BH$pHm,dat_BH$pair,function(x) t.test(x,mu = 0.5))
tapply(dat_BH$pHf,dat_BH$pair,function(x) t.test(x,mu = 0.5))

fisher.test(matrix(colSums(dat_BH[,c(4,6,7,5)]),2,2))
matrix_cont_table=apply(dat_BH[,c(4,6,7,5)],2,function(x) tapply(x,dat_BH$pair,sum))
res_chisq=apply(matrix_cont_table,1,function(x) chisq.test(matrix(x,2,2)))
res_FET=apply(matrix_cont_table,1,function(x) fisher.test(matrix(x,2,2)))

overall_avg_y_BH=mean((sqrt(sapply(res_FET,function(x) x$estimate))-1)/(sqrt(sapply(res_FET,function(x) x$estimate))+1))
overall_se_y_BH=sd((sqrt(sapply(res_FET,function(x) x$estimate))-1)/(sqrt(sapply(res_FET,function(x) x$estimate))+1))/sqrt(3)
t.test(tapply(dat_BH$Y_idx,dat_BH$pair,mean),mu=0)
t.test(overall_avg_y_BH,mu=0)

avg_pHm=rowSums(matrix_cont_table[,3:4])/rowSums(matrix_cont_table)
avg_pHf=rowSums(matrix_cont_table[,c(2,4)])/rowSums(matrix_cont_table)
t.test(avg_pHm,mu=0.5)
t.test(avg_pHf,mu=0.5)


####experiment: cross replicates####
dat_cross_repl=read.csv("./mating_assay_data_cross_repl.csv",header = T,stringsAsFactors = F)
dat_cross_repl$Y_idx=(sqrt(dat_cross_repl[,4]*dat_cross_repl[,5])-sqrt(dat_cross_repl[,6]*dat_cross_repl[,7]))/(sqrt(dat_cross_repl[,4]*dat_cross_repl[,5])+sqrt(dat_cross_repl[,6]*dat_cross_repl[,7]))
fm2=with(dat_cross_repl,lm(Y_idx~pair+coloring))
anova(fm2)
tapply(dat_cross_repl$Y_idx,dat_cross_repl$pair,function(x) t.test(x,mu = 0))
avg_y_CR=tapply(dat_cross_repl$Y_idx,dat_cross_repl$pair,mean)
se_y_CR=tapply(dat_cross_repl$Y_idx,dat_cross_repl$pair,function(x) sd(x)/sqrt(length(x)))
t.test(dat_cross_repl$Y_idx,mu=0)
overall_avg_y_CR=mean(dat_cross_repl$Y_idx)
overall_se_y_CR=sd(dat_cross_repl$Y_idx)/sqrt(length(dat_cross_repl$Y_idx))

fisher.test(matrix(colSums(dat_cross_repl[,c(4,6,7,5)]),2,2))
matrix_cont_table2=apply(dat_cross_repl[,c(4,6,7,5)],2,function(x) tapply(x,dat_cross_repl$pair,sum))
res_chisq2=apply(matrix_cont_table2,1,function(x) chisq.test(matrix(x,2,2)))
res_FET2=apply(matrix_cont_table2,1,function(x) fisher.test(matrix(x,2,2)))
overall_avg_y_CR=mean(dat_cross_repl$Y_idx)
overall_se_y_CR=sd(dat_cross_repl$Y_idx)/sqrt(length(dat_cross_repl$Y_idx))


####plotting####
png("Figure2b.png",width = 8.7,height = 8.7*4/3,units = "cm",res = 600,pointsize = 8)
bp=barplot(c(avg_y_BH,avg_y_CR),names.arg = c("A/H1","A/H3","A/H6","H1/H3","H1/H6","H3/H6"),
           ylim=c(-1,1),xlab="Mating combination",ylab = "Mean Y",
           col=rep(c("white","grey60"),each=3))
abline(h=0)
arrows(bp,c(avg_y_BH,avg_y_CR)-c(se_y_BH,se_y_CR),
       bp,c(avg_y_BH,avg_y_CR)+c(se_y_BH,se_y_CR),code = 3,length = 0.05,angle = 90)
text(1.9,0.8,labels = expression(2.8%*%10^-5))
text(5.5,0.8,labels = 0.66)
dev.off()

png("avg_Y.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 6)
par(mfrow=c(1,2))
bp=barplot(avg_y_BH,names.arg = c("A/H1","A/H3","A/H6"),ylim=c(-1,1),xlab="Mating combination",ylab = "Mean Y")
abline(h=0)
arrows(bp,avg_y_BH-se_y_BH,bp,avg_y_BH+se_y_BH,code = 3,length = 0.05,angle = 90)

bp2=barplot(avg_y_CR,names.arg = c("H1/H3","H1/H6","H3/H6"),ylim=c(-1,1),xlab="Mating combination",ylab = "Mean Y")
abline(h=0)
arrows(bp2,avg_y_CR-se_y_CR,bp,avg_y_CR+se_y_CR,code = 3,length = 0.05,angle = 90)
dev.off()

png("y_visualization.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 6)
par(mfrow=c(1,2))
bp1=barplot((sqrt(sapply(res_FET,function(x) x$estimate))-1)/(sqrt(sapply(res_FET,function(x) x$estimate))+1),
            names.arg = c("A/H1","A/H3","A/H6"),ylim=c(-1,1),xlab="Mating combination",ylab = "Mean Y")
abline(h=0)
text(bp1,(sqrt(sapply(res_FET,function(x) x$estimate))-1)/(sqrt(sapply(res_FET,function(x) x$estimate))+1)+0.1,
     labels = scientific(sapply(res_FET,function(x) x$p.value),digits = 2))
text(1.9,0.8,labels = "2.8e-05")
bp2=barplot((sqrt(sapply(res_FET2,function(x) x$estimate))-1)/(sqrt(sapply(res_FET2,function(x) x$estimate))+1),
            names.arg = c("H1/H3","H1/H6","H3/H6"),ylim=c(-1,1),xlab="Mating combination",ylab = "Mean Y")
abline(h=0)
text(bp2,(sqrt(sapply(res_FET2,function(x) x$estimate))-1)/(sqrt(sapply(res_FET2,function(x) x$estimate))+1)+c(0.1,-.1,-.1),
     labels = scientific(sapply(res_FET2,function(x) x$p.value),digits = 2))
#text(1.9,0.8,labels = "6.6e-01")
dev.off()


png("overall_Y.png",width = 8.7,height = 8.7,units = "cm",res = 600,pointsize = 6)
bp3=barplot(c(overall_avg_y_BH,overall_avg_y_CR),names.arg = c("ancestral/evolved","cross-evolved replicates"),
            ylim=c(-1,1),xlab="Mating combination",ylab = "Mean Y")
abline(h=0)
arrows(bp3,c(overall_avg_y_BH,overall_avg_y_CR)-c(overall_se_y_BH,overall_se_y_CR),
       bp3,c(overall_avg_y_BH,overall_avg_y_CR)+c(overall_se_y_BH,overall_se_y_CR),
       code = 3,length = 0.05,angle = 90)
dev.off()

####progeny assay####
dat=read.csv("/Volumes/Temp1/shengkai/mate_choice_CGE/progeny_assay_all.csv",header = T,stringsAsFactors = F)

dat_new=data.frame(P=rep(dat$P,2),M=rep(dat$M,2),repl=rep(dat$repl,2),transfer=rep(c("A","B"),each=45),
                   m=c(dat$T2_m,dat$T3_m),f=c(dat$T2_f,dat$T3_f))
dat_new$total=dat_new$m+dat_new$f
dat_new$cross[paste0(dat_new$P,dat_new$M)%in%"H1H1"]=1
dat_new$cross[paste0(dat_new$P,dat_new$M)%in%"H3H3"]=2
dat_new$cross[paste0(dat_new$P,dat_new$M)%in%"H6H6"]=3
dat_new$cross[paste0(dat_new$P,dat_new$M)%in%c("H1H3","H3H1")]=4
dat_new$cross[paste0(dat_new$P,dat_new$M)%in%c("H1H6","H6H1")]=5
dat_new$cross[paste0(dat_new$P,dat_new$M)%in%c("H3H6","H6H3")]=6

#dat_new$ratio=log2(dat_new$m/dat_new$f)

fm=with(dat_new,lm(total~paste(P,M)))
anova(fm)
LSD_res = LSD.test(dat_new$total,paste(dat_new$P,dat_new$M),DFerror = 81,MSerror = 57.257,alpha = 0.05,p.adj = "BH")
kru_res=kruskal(dat_new$total,paste(dat_new$P,dat_new$M),0.05,"BH")
fm1=with(dat_new,lm(total~P+M+transfer))
anova(fm1)
# fm2=with(dat_new,lm(abs(ratio)~paste(P,M)))
# anova(fm2)
# kru_res2=kruskal(abs(dat_new$ratio),paste(dat_new$P,dat_new$M),0.05,"BH")
fm3=with(dat_new,lm(total~as.factor(P==M)))
anova(fm3)
kru_res3=kruskal(dat_new$total,as.factor(dat_new$P==dat_new$M),0.05,"BH")

avg_total = (dat_new$total[1:45]+dat_new$total[-(1:45)])/2
fm4=lm(avg_total~dat_new$cross[1:45])
anova(fm4)
kru_res4=kruskal(avg_total,dat_new$cross[1:45],0.05,"none")

png("/Volumes/Temp1/shengkai/CHC/FigureS8.png",height = 8.7,width = 8.7,units = "cm",res=600,pointsize = 8)
par(mar=c(4,5,2,2))
boxplot((dat_new$total[1:45]+dat_new$total[-(1:45)])/2~dat_new$cross[1:45],
        names = NA,ylab = "Number of viable progeny",xlab ="",col = rep(c("white","grey60"),each = 3))
mtext(c("Pop. 1","H1","H3","H6","H1","H1","H3"),at = c(0,1:6),side = 1,line = 1)
mtext(c("Pop. 2","H1","H3","H6","H3","H6","H6"),at = c(0,1:6),side = 1,line = 2.5)
dev.off()

tapply((dat_new$total[1:45]+dat_new$total[-(1:45)])/2,dat_new$cross[1:45],mean)
#compatibility
#1-3: 32.10/(32.20+37.20)*2=0.925
#1-6: 29.55/(32.20+34.10)*2=0.891
#3-6: 33.30/(37.20+34.10)*2=0.934
#distance
#1-3: 1.43
#1-6: 1.58
#3-6: 1.27
plot(x = c(1.43,1.58,1.27),y=c(0.925,0.891,0.934))
cor.test(c(0.925,0.891,0.934),c(1.43,1.58,1.27))
gg_heatmap_data=data.frame(var1=rep(c("1","3","6"),each=3),var2=rep(c("1","3","6"),3),
                           value=c(NA,0.925,0.891,1.43,NA,0.934,1.58,1.27,NA))
gg_label_data=data.frame(var1=c("1","3","6"),var2=c("1","3","6"),
                           value=c(1,3,6))

ggheatmap <- ggplot(na.omit(gg_heatmap_data), aes(var2, var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 1.25, limit = c(0.89,1.61), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+
  coord_fixed()

png("/Volumes/Temp1/shengkai/CHC/FigureS5.png",height = 8.7,width = 8.7,units = "cm",res=600)
ggheatmap + 
  geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
  labs(y="Distance", x = "Compatibility")+
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.text = element_blank())+
  geom_text(aes(x = var2, y = var1, label=value),data = gg_label_data, color = "black", size = 6) 
dev.off()  

####simple comparison####
with(dat_new,t.test(total[P==M],total[P!=M],alternative = "greater"))
with(dat_new,wilcox.test(total[P==M],total[P!=M],alternative = "greater"))
with(dat_new,t.test(abs(ratio[P==M]),abs(ratio[P!=M])))
with(dat_new,wilcox.test(abs(ratio[P==M]),abs(ratio[P!=M])))
with(dat_new,t.test((total[P==M][1:15]+total[P==M][-(1:15)])/2,(total[P!=M][1:30]+total[P!=M][-(1:30)])/2,alternative = "greater"))
with(dat_new,wilcox.test((total[P==M][1:15]+total[P==M][-(1:15)])/2,(total[P!=M][1:30]+total[P!=M][-(1:30)])/2,alternative = "greater"))

png("/Volumes/Temp1/shengkai/CHC/Figure4b.png",height = 8.7,width = 8.7,units = "cm",res=600,pointsize = 8)
par(mar=c(4,5,2,2))
with(dat_new,boxplot(total[P==M],total[P!=M],names=c("Within-replicate","cross-replicate"),
                     ylab="Numbers of viable progenies",col=c("white","grey60")))
text(1.5,55,labels = "p = 0.024")
dev.off()

png("/Volumes/Temp1/shengkai/CHC/Figure4b_new.png",height = 8.7,width = 8.7,units = "cm",res=600,pointsize = 8)
par(mar=c(4,5,2,2))
with(dat_new,boxplot((total[P==M][1:15]+total[P==M][-(1:15)])/2,
                     (total[P!=M][1:30]+total[P!=M][-(1:30)])/2,
                     names=c("Within-replicate","cross-replicate"),
                     ylab="Numbers of viable progeny",col=c("white","grey60")))
text(1.5,49,labels = "p = 0.031")
dev.off()

png("./progeny_assay_v1.png",height = 8.7,width = 15.4,units = "cm",res=600,pointsize = 8)
par(mfrow=c(1,2))
with(dat_new,boxplot(total[P==M],total[P!=M],names=c("matched","mismatched"),main="Viability",
                     ylab="Numbers of viable progenies"))
with(dat_new,boxplot(abs(ratio[P==M]),abs(ratio[P!=M]),names=c("matched","mismatched"),
                     main="Sex allocation",ylab=expression(abs(log[2]("Male"/"Female")))))
dev.off()
