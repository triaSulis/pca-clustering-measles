library(psych)
library(ggplot2)
library(PerformanceAnalytics)
library(dplyr)
library(mvnormtest)
library(cluster) #Clustering
library(factoextra) #Clustering & Viz
library(tidyverse) #Data Manipulation
library(dplyr)
library(corrplot)
library(fpc) #Clustering
library(clValid) #Choose c
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(stats)

#Data Preparation
df=read.csv("C:/Users/Tria/Documents/Download/S2-Stat/Analisis Data C/W12/campak.csv",header=T,sep=";")
rownames(df)=df$Provinsi #kolom kecamatan didrop
df=df[,c(2:13)]#mendrop variable "kesling" karena berpotensi sebagai variabel kategorik.
head(df)
summary(df)
win.graph()
boxplot(df)

#Multivariate Normality
dft=t(df)
mshapiro.test(dft) #tdk normal, karena multivariat dan skala beda2

#Bartlett Test (Homogeneity)
bartlett.test(df) #variansinya sama atau ngga, karena multivariat dan skala beda2

#Correlation/Multicolinearity, PCA bisa --> kalau ada multiko
r=cor(df)
r
corrplot(r)
chart.Correlation(df) #***korelasi signifikan, banyak itu ada multiko, bisa dilanjut
cortest.bartlett(r,n=nrow(df))

#KMO Test
KMO(r) #overall>0.5, stop if <0.5; MSA>0.5, drop if <0.5 (1 by 1 pake backward, sampai memenhi semua, kecuali STBM) ==>bisa diPCA/FA
#Setelah mendrop variabel STBM (overall naik 0.79) maka MSA each item sudah diatas 0.5

#PCA
df.pr=prcomp(df,center=TRUE, scale=TRUE) #discaling variabel dengan mean
df.pr
summary(df.pr)

#Eigen Value #punya eigen vektor dari banyaknya pcn =menentukan dengan eigenvalue>1
round(df.pr$sdev^2,2) 
screeplot(df.pr, type = "l", npcs=7, main="Screeplot")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"), col=c("red"), lty=5, cex=0.6)

#Variance Explained
exp_var=df.pr$sdev^2/sum(df.pr$sdev^2)
exp_var
win.graph()
fviz_eig(df.pr, addlabels=TRUE, ylim=c(0,100))

#Cumulative Variance Plot
cumpro=cumsum(df.pr$sdev^2/sum(df.pr$sdev^2))
cumpro
plot(cumpro[0:7], xlab="PC#", ylab="Amount of explained variance", main="Cumulative Variance Plot")
bline(v = 2, col="blue", lty=5)
abline(h = 0.6653290, col="blue", lty=5)
legend("bottomright", legend=c("Cut-off @PC2"), col=c("blue"), lty=5, cex=0.6)

#PC Eq
round(df.pr$rotation[,1],3)
df.pr$x[,1]

#FA menggunakan metode MLE, penentuan berapa factor
df.fa=factanal(df,1, scores="regression", rotation="varimax") #berapa faktor yang didapat sewaktu pca, rotation untuk lebih mendetailkan setiap variabel masuk mana karena bedanya diatas 0.5.
df.fa$loadings
c.mle=solve(r)%*%df.fa$loadings
df.fa$scores

#FA menggunakan metode PC
ncomp=2
rawLoadings=df.pr$rotation[,1:ncomp] %*% diag(df.pr$sdev, ncomp, ncomp) #unrotated data
rawLoadings
rotatedLoadings=varimax(rawLoadings)$loadings #rotated data
rotatedLoadings
L=as.matrix(rawLoadings)
L2=c(-0.9327449,0.9704237,-0.8557196,-0.9533302,-0.9792645,-0.9795471,-0.9847029,-0.9625463,-0.9616820,-0.9502218,-0.9653107,-0.8693857)
L2=matrix(L2,nrow=12,ncol=1)
L2
C=L%*%solve(t(L)%*%L)
C.rt=L2%*%solve(t(L2)%*%L2)
sc.fa=as.matrix(scale(df))%*%C
sc.fa.rt=as.matrix(scale(df))%*%C.rt

#New df
df=read.csv("C:/Users/Tria/Documents/Download/S2-Stat/Analisis Data C/W12/campak.csv",header=T,sep=";")
df.pc=data.frame(df[,-1],df.pr$x[,1]) #ada tambahan variabel Jumlah Kasus DBD "r" ;rendah, "t":tinggi
df.pc
head(df.pc)
summary(df.pc)
win.graph()
boxplot(df.pc)
head(df.pc)
rpc=cor(df.pc[,1:2])
cortest.bartlett(rpc,n=nrow(df.pc[,1:2])) #tidak ada multiko, tolak H0

#Heat Map
win.graph()
fviz_pca_var(df.pr,col.var="coord",gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),repel=T, axes=c(1,2)) #coord untuk visualnya lingkaran, panjang panah sesuai jari2nya

#Group by PC
fviz_pca_ind(df.pr,geom.ind="point",pointshape=21,pontsize=2,fill.ind=as.factor(df.pc$df.Provinsi),col.ind="black",palette="jco",addEllipses=T,label="var",col.var="black",repel=T)

#Biplot
fviz_pca_biplot(df.pr, repel=TRUE, col.var = "#0f0501", col.ind=as.factor(df.pc$df.Provinsi))
ggbiplot(df.pr, labels=rownames(df), ellipse=TRUE, groups=as.factor(df.pc$df.Provinsi))

#Clustering K-Means
dataclus=df.pc
dataclus
df.pc
win.graph()
fviz_nbclust(dataclus, kmeans, method="wss") #Elbow Method
fviz_nbclust(dataclus, kmeans, method="silhouette") #Silhouette Method
intern=clValid(dataclus, nClust = 2:5,clMethods = c("hierarchical","kmeans","pam"), validation = "internal")
summary(intern)
km_fit = kmeans(dataclus,centers=2)
print(km_fit)
df.clus=data.frame(dataclus,km_fit$cluster) #Adding Cluster to DF
win.graph()
fviz_cluster(km_fit,data=dataclus) #Cluster Plot
table(km_fit$cluster) #Number of members in each clusters
km_fit$centers #Represented object from each clusters
df.clus %>%
  mutate(cluster=km_fit.cluster) %>%
  group_by(cluster) %>%
  summarise_all("mean") #Desc of each clusters

#Model Criterions
score_KM=function(k){
  km=kmeans(dataclus, centers=k)
  ss=silhouette(km$cluster, dist(dataclus))
  mean(ss[,3])
}
k=2
avg_sil=sapply(k,score_KM) #Silhouette Score
avg_sil
mod_cri=function(Data, nc, c) #penghitungan kriteria kebaikan
{
  n = dim(Data)[1]
  p = dim(Data)[2]
  X = Data[,1:(p-1)]
  Group = Data[,p]
  p = dim(X)[2]
  Mean.X = matrix(ncol = p, nrow = (nc+1))
  for (i in 1:nc)
  {
    for (j in 1:p)
    {
      Mean.X[i,j] = mean(X[which(Group==i),j])
      Mean.X[(nc+1),j] = mean(X[,j])
    }
  }
  SST = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      SST[i,j] = (X[i,j] - Mean.X[(nc+1),j])^2
    }
  }
  SST = sum(sum(SST))
  SSE = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      for (k in 1:nc)
      {
        if (Group[i]==k)
        {
          SSE[i,j] = (X[i,j] - Mean.X[k,j])^2
        }
      }
    }
  }
  SSE = sum(sum(SSE))
  Rsq = (SST-SSE)/SST
  icdrate = 1-Rsq
  Pseudof = (Rsq/(c-1))/((icdrate)/(nc-c))
  ssb=SST-SSE
  list(SSW=SSE, SST=SST, SSB=ssb, Rsq=Rsq, icdrate=icdrate, pseudof=Pseudof)
}
mod_cri(df.clus,length(df.clus),2)



