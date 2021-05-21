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

doubledata$NLGN4diff<-doubledata$NLGN4Xexp-doubledata$NLGN4Xexp2
w<-which(abs(doubledata$NLGN4diff)>0)
doubledata$NLGN4diff[w]<-1
table(doubledata$NLGN4diff,doubledata$twin)
#shows twin1 or 2 in column and different (0/1) on NLGN4expression in rows.

cor(doubledata$NLGN4Xexp,doubledata$NLGN4Xexp2,use="complete.obs")

#Ah, it turns out we have MZ and DZ pairs! We only have same-sex (because only females!)

MZ<-doubledata[doubledata$zygosity==1,]
DZ<-doubledata[doubledata$zygosity==2,]

plot(MZ$NLGN4Xexp,MZ$NLGN4Xexp2)
plot(DZ$NLGN4Xexp,DZ$NLGN4Xexp2)

