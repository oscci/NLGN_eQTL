.
#Read in file from twin_discordant analysis
doubledata<-read.csv('doubledata_orig.csv')
colnames(doubledata)[1]<-'twin_id' #this had been wrongly labelled as rand_id, so for clarity redo it.
#Read in file from Nuala
nualadat<-read.csv('Twin1and2_fusion_data.csv')
colnames(nualadat)[2:3]<-c('id','NLGN4Xexp')
#remove DB and underscore
nualadat$id<-gsub('DB','',nualadat$id)
nualadat$id<-gsub('_','',nualadat$id)

#now add gene exp column to doubledata
doubledata$NLGN4Xexp<-NA
for(i in 1:nrow(nualadat)){
  w<-which(doubledata$twin_id==nualadat$id[i])
  if(length(w)>0){
    doubledata$NLGN4Xexp[w]<-nualadat$NLGN4Xexp[i]
  }
}
doubledata$NLGN4Xexp2<-NA
half1<-1:(nrow(doubledata)/2)
half2<-half1+nrow(doubledata)/2
doubledata$NLGN4Xexp2[half1]<-doubledata$NLGN4Xexp[half2]
doubledata$NLGN4Xexp2[half2]<-doubledata$NLGN4Xexp[half1]
plot(doubledata$NLGN4Xexp[half1],doubledata$NLGN4Xexp2[half1],col=doubledata$zygosity)

doubledata$NLGN4diff<-abs(doubledata$NLGN4Xexp-doubledata$NLGN4Xexp2)
doubledata$NLGN4conc<-doubledata$NLGN4diff
w<-which(abs(doubledata$NLGN4diff)>0)
doubledata$NLGN4disc[w]<-1
table(doubledata$NLGN4disc,doubledata$twin)
#shows twin1 or 2 in column and different (0/1) on NLGN4expression in rows.

cor(doubledata$NLGN4Xexp,doubledata$NLGN4Xexp2,use="complete.obs")

#Ah, it turns out we have MZ and DZ pairs! We only have same-sex (because only females!)

MZ<-doubledata[doubledata$zygosity==1,]
DZ<-doubledata[doubledata$zygosity==2,]

plot(MZ$NLGN4Xexp,MZ$NLGN4Xexp2)
plot(DZ$NLGN4Xexp,DZ$NLGN4Xexp2)

#Need to add phenotypes 
genopheno<-read.csv('twin_with_genopheno.csv') #file saved from wrangling data from all.data in analysis Double_hit_results.Rmd

#now add pheno columns
doubledata$nonword_rep_ss<-NA
doubledata$langfactor<-NA
doubledata$global_neurodev<-NA
for(i in 1:nrow(genopheno)){
  w<-which(doubledata$twin_id==genopheno$ID[i])
  if(length(w)>0){
    doubledata$nonword_rep_ss[w]<-genopheno$nonword_rep_ss[i]
    doubledata$langfactor[w]<-genopheno$langfactor[i]
    doubledata$global_neurodev[w]<-genopheno$global_neurodev[i]
  }
}
#now make double entry
doubledata$nonword_rep_ss2<-NA
doubledata$langfactor2<-NA
doubledata$global_neurodev2<-NA

#find column numbers
colrange<-(ncol(doubledata)-2):ncol(doubledata)

doubledata[half1,colrange]<-doubledata[half2,(colrange-3)]
doubledata[half2,colrange]<-doubledata[half1,(colrange-3)]

MZrow<-which(doubledata$zygosity==1)
DZrow<-which(doubledata$zygosity>1)

#Correlations for nonword rep
cor(doubledata$nonword_rep_ss[MZrow],doubledata$nonword_rep_ss2[MZrow],use='complete.obs')
cor(doubledata$nonword_rep_ss[DZrow],doubledata$nonword_rep_ss2[DZrow],use='complete.obs')

cor(doubledata$langfactor[MZrow],doubledata$langfactor2[MZrow],use='complete.obs')
cor(doubledata$langfactor[DZrow],doubledata$langfactor2[DZrow],use='complete.obs')
plot(doubledata$langfactor,doubledata$langfactor2,col=doubledata$zygosity)

cor(doubledata$global_neurodev[MZrow],doubledata$global_neurodev2[MZrow],use='complete.obs')
cor(doubledata$global_neurodev[DZrow],doubledata$global_neurodev2[DZrow],use='complete.obs')
table(doubledata$global_neurodev[MZrow],doubledata$global_neurodev2[MZrow])
table(doubledata$global_neurodev[DZrow],doubledata$global_neurodev2[DZrow])
doubledata$nonword_rep_diff<-abs(doubledata$nonword_rep_ss-doubledata$nonword_rep_ss2)
doubledata$langfactor_diff<-abs(doubledata$langfactor-doubledata$langfactor2)

doubledata$neurodev_diff<-abs(doubledata$global_neurodev-doubledata$global_neurodev2)

plot(doubledata$NLGN4diff[DZrow],doubledata$nonword_rep_diff[DZrow])
plot(doubledata$NLGN4diff[DZrow],doubledata$langfactor_diff[DZrow])

plot(doubledata$NLGN4diff[DZrow],doubledata$neurodev_diff[DZrow])

