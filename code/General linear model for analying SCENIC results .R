
## Read SCENIC AUC results and correct barcodes ##
MSC=
library(stringr)
AUC=read.csv('F:/AUC_eachCell.csv')
meta=MSC@meta.data

BM1_origin=str_sub(row.names(subset(meta,group=='BM1')),start=5)
BM2_origin=str_sub(row.names(subset(meta,group=='BM2')),start=5)
BM3_origin=str_sub(row.names(subset(meta,group=='BM3')),start=5)
AD1_origin=str_sub(row.names(subset(meta,group=='AD1')),start=5)
AD2_origin=str_sub(row.names(subset(meta,group=='AD2')),start=5)
AD3_origin=str_sub(row.names(subset(meta,group=='AD3')),start=5)
PM1_origin=str_sub(row.names(subset(meta,group=='PM1')),start=5)
PM2_origin=str_sub(row.names(subset(meta,group=='PM2')),start=5)
PM3_origin=str_sub(row.names(subset(meta,group=='PM3')),start=5)
UC1_origin=str_sub(row.names(subset(meta,group=='UC1')),start=5)
UC2_origin=str_sub(row.names(subset(meta,group=='UC2')),start=5)
UC3_origin=str_sub(row.names(subset(meta,group=='UC3')),start=5)

ADbarcode1=intersect(paste(AD1_origin,'.1',sep = ''),colnames(AUC))
ADbarcode2=intersect(paste(AD2_origin,'.2',sep = ''),colnames(AUC))
ADbarcode3=intersect(paste(AD3_origin,'.3',sep = ''),colnames(AUC))
BMbarcode1=intersect(paste(BM1_origin,'.4',sep = ''),colnames(AUC))
BMbarcode2=intersect(paste(BM2_origin,'.5',sep = ''),colnames(AUC))
BMbarcode3=intersect(paste(BM3_origin,'.6',sep = ''),colnames(AUC))
PMbarcode1=intersect(paste(PM1_origin,'.7',sep = ''),colnames(AUC))
PMbarcode2=intersect(paste(PM2_origin,'.8',sep = ''),colnames(AUC))
PMbarcode3=intersect(paste(PM3_origin,'.9',sep = ''),colnames(AUC))
UCbarcode1=intersect(paste(UC1_origin,'.10',sep = ''),colnames(AUC))
UCbarcode2=intersect(paste(UC2_origin,'.11',sep = ''),colnames(AUC))
UCbarcode3=intersect(paste(UC3_origin,'.12',sep = ''),colnames(AUC))

row.names(meta)=c(ADbarcode1,ADbarcode2,ADbarcode3,BMbarcode1,BMbarcode2,BMbarcode3,
                  PMbarcode1,PMbarcode2,PMbarcode3,UCbarcode1,UCbarcode2,UCbarcode3)


##¢ÙCluster 123 vs Other
meta_S123vsOther=read.csv('meta_for_GLM_S123vsOther.csv',row.names = 'barcode')
AUC_S123vsOther=AUC[,as.vector(row.names(meta_S123vsOther))]
AUC_S123vsOther=cbind(t(AUC_S123vsOther),meta_S123vsOther$S123vsOther)
colnames(AUC_S123vsOther)[349]='S123_vs_other'

colnames(AUC_S123vsOther)=gsub("[ (.*)]", "", colnames(AUC_S123vsOther))
AUC_S123vsOther=as.data.frame(AUC_S123vsOther)

compare_t_list=list()
for(i in colnames(AUC_S123vsOther[,c(1:348)])){
  f=as.formula(paste(i,'S123_vs_other',sep='~'))
  model=glm(f,data = AUC_S123vsOther)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:348)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}
write.csv(t(compare_t_mat),'F:/S123_vs_other.csv')
rm(AUC_S123vsOther,meta_S123vsOther,compare_t_list,compare_t_mat)

##¢ÚCluster 4 vs 123
meta_S4vsS123=read.csv('F:/meta_for_GLM_S4vsS123.csv',row.names = 'barcode')
AUC_S4vsS123=AUC[,as.vector(row.names(meta_S4vsS123))]
AUC_S4vsS123=cbind(t(AUC_S4vsS123),meta_S4vsS123$S4vsS123)
colnames(AUC_S4vsS123)[349]='S4_vs_S123'

colnames(AUC_S4vsS123)=gsub("[ (.*)]", "", colnames(AUC_S4vsS123))
AUC_S4vsS123=as.data.frame(AUC_S4vsS123)

compare_t_list=list()
for(i in colnames(AUC_S4vsS123[,c(1:348)])){
  f=as.formula(paste(i,'S4_vs_S123',sep='~'))
  model=glm(f,data = AUC_S4vsS123)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:348)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}
write.csv(t(compare_t_mat),'F:/S4_vs_S123.csv')
rm(AUC_S4vsS123,meta_S4vsS123,compare_t_list,compare_t_mat)

##¢ÛCluster 5 vs 4
meta_S5vsS4=read.csv('F:/meta_for_GLM_S5vsS4.csv',row.names = 'barcode')
AUC_S5vsS4=AUC[,as.vector(row.names(meta_S5vsS4))]
AUC_S5vsS4=cbind(t(AUC_S5vsS4),meta_S5vsS4$S5vsS4)
colnames(AUC_S5vsS4)[349]='S5_vs_S4'

colnames(AUC_S5vsS4)=gsub("[ (.*)]", "", colnames(AUC_S5vsS4))
AUC_S5vsS4=as.data.frame(AUC_S5vsS4)

compare_t_list=list()
for(i in colnames(AUC_S5vsS4[,c(1:348)])){
  f=as.formula(paste(i,'S5_vs_S4',sep='~'))
  model=glm(f,data = AUC_S5vsS4)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:348)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}
write.csv(t(compare_t_mat),'F:/S5_vs_S4.csv')
rm(AUC_S5vsS4,meta_S5vsS4,compare_t_list,compare_t_mat)

##¢ÛCluster 6 vs 4
meta_S6vsS4=read.csv('F:/meta_for_GLM_S6vsS4.csv',row.names = 'barcode')
AUC_S6vsS4=AUC[,as.vector(row.names(meta_S6vsS4))]
AUC_S6vsS4=cbind(t(AUC_S6vsS4),meta_S6vsS4$S6vsS4)
colnames(AUC_S6vsS4)[349]='S6_vs_S4'

colnames(AUC_S6vsS4)=gsub("[ (.*)]", "", colnames(AUC_S6vsS4))
AUC_S6vsS4=as.data.frame(AUC_S6vsS4)

compare_t_list=list()
for(i in colnames(AUC_S6vsS4[,c(1:348)])){
  f=as.formula(paste(i,'S6_vs_S4',sep='~'))
  model=glm(f,data = AUC_S6vsS4)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:348)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}
write.csv(t(compare_t_mat),'F:/S6_vs_S4.csv')
rm(AUC_S6vsS4,meta_S6vsS4,compare_t_list,compare_t_mat)

##¢ÛCluster 7 vs 4
meta_S7vsS4=read.csv('F:/meta_for_GLM_S7vsS4.csv',row.names = 'barcode')
AUC_S7vsS4=AUC[,as.vector(row.names(meta_S7vsS4))]
AUC_S7vsS4=cbind(t(AUC_S7vsS4),meta_S7vsS4$S7vsS4)
colnames(AUC_S7vsS4)[349]='S7_vs_S4'

colnames(AUC_S7vsS4)=gsub("[ (.*)]", "", colnames(AUC_S7vsS4))
AUC_S7vsS4=as.data.frame(AUC_S7vsS4)

compare_t_list=list()
for(i in colnames(AUC_S7vsS4[,c(1:348)])){
  f=as.formula(paste(i,'S7_vs_S4',sep='~'))
  model=glm(f,data = AUC_S7vsS4)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:348)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}
write.csv(t(compare_t_mat),'F:/S7_vs_S4.csv')
rm(AUC_S7vsS4,meta_S7vsS4,compare_t_list,compare_t_mat)








