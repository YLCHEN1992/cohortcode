#CSTATISTICS VERSION 202408221453
#version 20240722931
#Genetic analysis annotation by 
#vector genetype as character class and split by comma
comg=function(datap){
N=length(datap)
T=table(datap)/length(datap)
data=unlist(strsplit(datap,","))
D=table(data)/length(data);PIC=0
D2=D^2;for(i in 1:length(D2)){for(j in 1:length(D2)){if(j>i){PIC=PIC+2*D2[i]*D2[j]}}}
PIC=as.numeric(1-sum(D2)-PIC)
n=length(D)
F=matrix(0,n,n)
rownames(F)=names(D)
colnames(F)=names(D)
idn=c()
for(i in 1:n){
for(j in 1:n){
F[i,j]=D[i]*D[j]
idn=c(idn,paste(rownames(F)[i],colnames(F)[j],sep=","))}}
egf=as.vector(F);names(egf)=idn
egf=egf[!duplicated(egf)]
for(i in 1:length(egf)){tem=unlist(strsplit(names(egf)[i],","));if(tem[1]!=tem[2]){egf[i]=2*egf[i]}}
smg=function(gnm){paste(sort(unlist(strsplit(gnm,","))),collapse=",")}
names(egf)=sapply(names(egf),smg)
names(T)=sapply(names(T),smg)
fdata=data.frame(oct=as.vector(T),exp=egf[match(names(T),names(egf))])
df=nrow(fdata)
if(N>29){ka=sum((fdata[,1]-fdata[,2])^2/fdata[,2])}else{ka=sum((abs(fdata[,1]-fdata[,2])-0.5)^2/fdata[,2])}
hwp=1-pchisq(ka,df)
tf=-2*N*((1-D)*log(1-D)/D);tl=-2*N*D*log(D)/(1-D)
TE=egf[match(names(T),names(egf))]
hoE=sum(TE[sapply(names(TE),function(x){tmp=unlist(strsplit(x,","));tmp[1]==tmp[2]})])
heE=sum(TE[!sapply(names(TE),function(x){tmp=unlist(strsplit(x,","));tmp[1]==tmp[2]})])
ho=sum(T[sapply(names(T),function(x){tmp=unlist(strsplit(x,","));tmp[1]==tmp[2]})])
he=sum(T[!sapply(names(T),function(x){tmp=unlist(strsplit(x,","));tmp[1]==tmp[2]})])
fis=(heE-he)/heE;ne=1/ho;ibd=1/(1-fis)
list("samplenumber"=N,"GenetypeFrq"=T,"GeneFrq"=D,"HWEP"=hwp,"PIC"=PIC,
"He"=he,"Ho"=ho,"Ne"=ne,"Tfix"=tf,"Tloss"=tl,"FIS"=fis,"IBD"=1/(1-fis))}
#locus linkage calucation via two genetype vectors
gld=function(x,y){
xf=table(x)/length(x);yf=table(y)/length(y)
z=c(names(xf)[1],names(yf)[1]);y=y[grep(z[1],x)];y=y[grep(z[2],y)]
pab=length(y)/length(x);dd=pab-xf[1]*yf[1]
dmax=max(c(-xf[1]*yf[1],-(1-xf[1])*(1-yf[1])));dmin=min(c(xf[1]*(1-yf[1]),(1-xf[1])*yf[1]))
if(dd>0){ddd=dd/dmin}else if(dd<0){ddd=dd/dmax}else{ddd=dd}
rr=dd^2/(xf[1]*(1-xf[1])*yf[1]*(1-yf[1]))
list("info"=z,"D"=as.numeric(dd),"D'"=as.numeric(ddd),"R^2"=as.numeric(rr))}
#mutiple locus linkage calucation via genetype matrix 
mgld=function(mdata,i=4){apply(mdata,2,function(y){apply(mdata,2,function(x){gld(x,y)[[i]]})})}

#version 20240722931
#Feature analysis annotation by 
#PCA classifying via row
PCA=function(mydata){
mydata=scale(mydata, center = T, scale = T)
mtcars.pca=prcomp(mydata, center = F,scale. = F)
sm.importance=as.data.frame(summary(mtcars.pca)$importance)
pc1=paste("PC1(",format(sm.importance$PC1[2]*100),"%)",sep="")
pc2=paste("PC2(",format(sm.importance$PC2[2]*100),"%)",sep="")
pc3=paste("PC3(",format(sm.importance$PC3[2]*100),"%)",sep="")
layout=as.data.frame(mtcars.pca$x)
list(layout=layout,PC1=pc1,pc2=pc2,pc3=pc3)}

#umap and tsne analysis by packages
fumap=function(data,seed=123,ps=T){
if (!require("umap", quietly = TRUE)){install.packages("umap")}
library(umap);set.seed(seed)
result=umap(data, preserve.seed=ps)
result}
ftsne=function(data){
if (!require("tsne", quietly = TRUE)){install.packages("tsne")}
library(tsne)
result=tsne(data,max_iter=1000,min_cost=0,k=2,initial_dims = 30)
result}

#Simple PCoA
PCoA=function(df,method='bray',nf=2){
if (!require("ade4", quietly = TRUE)){install.packages("ade4")}
if (!require("vegan", quietly = TRUE)){install.packages("vegan")}
library(ade4);library(vegan)
df.dist=vegdist(df,method=method)
pcoa=dudi.pco(df.dist,scannf=F,nf=nf)
layout=as.data.frame(pcoa$li)
pcoa$eig=pcoa$eig/sum(pcoa$eig)
pcoa1=pcoa$eig[1];pcoa2=pcoa$eig[2];pcoa3=pcoa$eig[3]
list(layout=layout,pcoa1=pcoa1,pcoa2=pcoa2,pcoa3=pcoa3)}

#Latent Class Analysis
LCA=function(data,colns,nclass,maxiter=3000){
if (!require("poLCA", quietly = TRUE)){install.packages("poLCA")}
if (!require("reshape2", quietly = TRUE)){install.packages("reshape2")}
library(poLCA)
colns=paste(colns,collapse=",")
formulas=as.formula(paste0("cbind(",colns,")~1"))
M1=poLCA(formula=formulas,data=data,nclass=nclass,graph=F,maxiter=maxiter)
table=reshape2::melt(M1$probs)
list(layout=table,aic=M1$aic,bic=M1$bic,predclass=M1$predclass,posterior=M1$posterior)}

#Latent Profile Analysis
LPAs=function(LPA_data,colns,G,modelNames=NULL){
if (!require("mclust", quietly = TRUE)){install.packages("mclust")}
library(mclust)
LPA_data=LPA_data[,colnames(LPA_data)%in%colns]
Mclust(LPA_data,G=G,modelNames=modelNames)}
LPA=function(LPA_data,colns){
if (!require("mclust", quietly = TRUE)){install.packages("mclust")}
library(mclust)
LPA_data=LPA_data[,colnames(LPA_data)%in%colns]
BIC=mclustBIC(LPA_data);ICL=mclustICL(LPA_data)
PBIC=unlist(strsplit(names(summary(BIC))[1],","))
PICL=unlist(strsplit(names(summary(ICL))[1],","))
LPA_BIC=LPAs(LPA_data,colns,G=as.numeric(PBIC[2]),modelNames=PBIC[1])
LPA_ICL=LPAs(LPA_data,colns,G=as.numeric(PICL[2]),modelNames=PICL[1])
list(BIC=BIC,ICL=ICL,LPA_BIC=LPA_BIC,LPA_ICL=LPA_ICL,
BICpredclass=LPA_BIC$classification,BICposterior=LPA_BIC$z,
ICLpredclass=LPA_ICL$classification,ICLposterior=LPA_ICL$z,
Boostrap=mclustBootstrapLRT(LPA_data,modelName=PBIC[1]))}

#version 202407221622 
#Model analysis annotation by 
#linear regression model
chen_lm=function(data,x,y,weights=NULL){
formulas=as.formula(paste(y,paste(x,collapse="+"),sep="~"))
mod=summary(lm(formula=formulas,weights=weights,data=data))
mod$fstatistic=as.numeric(mod$fstatistic)
list(coef=mod$coefficients,r=mod$r.squared,adj.r=mod$adj.r.squared,F=mod$fstatistic,
pvalues=pf(mod$fstatistic[1],mod$fstatistic[2],mod$fstatistic[3],lower.tail=F))}

#logistic regression model
chen_glm=function(data,x,y,family="gaussian",weights=NULL){
formulas=as.formula(paste(y,paste(x,collapse="+"),sep="~"))
mod=glm(formula=formulas,family=family,data=data);bic=BIC(mod);mod=summary(mod)
mod$fstatistic=as.numeric(mod$fstatistic)
list(coef=mod$coefficients,aic=mod$aic,bic=bic)}

#LASSO regression
chen_lasso=function(data,x,y,seed=106,family=c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
weights=NULL,maxit=100000){
if (!require("glmnet", quietly = TRUE)){install.packages("glmnet")}
library(glmnet)
set.seed(seed)
data=na.omit(data)
x=data[,colnames(data)%in%x]
y=as.vector(data[,colnames(data)==y])
cvfit=cv.glmnet(as.matrix(x),y,family=family,weights=weights,alpha=1,maxit=maxit)
fit=glmnet(as.matrix(x),y,family=family,weights=weights,alpha=1,maxit=maxit,lambda = cvfit$lambda.min)
coefs=coef(fit,s=cvfit$lambda.min)
index=which(coefs!=0)
actcoef=coefs[index]
lassgene=row.names(coefs)[index]
coefs=data.frame(factor=lassgene,coef=actcoef);coefs}

#ridge regression
chen_ridge=function(data,x,y,seed=106,family= c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
weights=NULL,maxit=100000){
if (!require("glmnet", quietly = TRUE)){install.packages("glmnet")}
library(glmnet)
set.seed(seed)
data=na.omit(data)
x=data[,colnames(data)%in%x]
y=as.vector(data[,colnames(data)==y])
cvfit=cv.glmnet(as.matrix(x),y,family=family,weights=weights,alpha=0,maxit=maxit)
fit=glmnet(as.matrix(x),y,family=family,weights=weights,alpha=0,maxit=maxit,lambda = cvfit$lambda.min)
coefs=coef(fit,s=cvfit$lambda.min)
index=which(coefs!=0)
actcoef=coefs[index]
lassgene=row.names(coefs)[index]
coefs=data.frame(factor=lassgene,coef=actcoef);coefs}

#Support Vector Machine Model
chen_svm=function(data,x,y,seed=1071,kernel="radial"){
if (!require("e1071", quietly = TRUE)){install.packages("e1071")}
library(e1071)
set.seed(seed)
formulas=as.formula(paste(y,paste(x,collapse="+"),sep="~"))
svm(formula=formulas,data=data,kernel=kernel)}

#RandomForeset Model
chen_randomforest=function(data,x,y,seed=1071,ntree=500,kernel="radial"){
if (!require("randomForest", quietly = TRUE)){install.packages("randomForest")}
library(randomForest)
set.seed(seed)
formulas=as.formula(paste(y,paste(x,collapse="+"),sep="~"))
mod=randomForest(formula=formulas,data=data,importance=TRUE,proximity=TRUE,ntree=ntree)
import=importance(mod)
list(mod=mod,importance=import)}

#Model evaluate by values
modevalue=function(truevalues,predictions,ROC=F,method="pearson"){
R2=cor(truevalues, predictions,method=method)^2
RP=cor.test(truevalues, predictions,method=method)$p.value
MSE=mean((truevalues - predictions)^2)
RMSE=sqrt(MSE)
MAE=mean(abs(truevalues - predictions))
if(ROC){library(pROC) 
auc=as.numeric(roc(truevalues, predictions)$auc)
list(auc=auc,r2=R2,rp=RP,mse=MSE,rmse=RMSE,mae=MAE)
}else{list(r2=R2,rp=RP,mse=MSE,rmse=RMSE,mae=MAE)}}
#other functions including step() AIC() BIC()

#Survival Analysis
chen_survival=function(sdata,time,status,x){
if (!require("survival", quietly = TRUE)){install.packages("survival")}
library(survival)
st=Surv(sdata[,colnames(sdata)==time],sdata[,colnames(sdata)==status])
formulas=as.formula(paste0("st~",paste(x,collapse="+")))
logrank=survdiff(formulas,data=sdata)
coxmod=summary(coxph(formulas,data=sdata))
plotfit=survfit(formulas,data=sdata)
list(coxmod=coxmod,coxtable=coxmod$coefficients,
logrankf=logrank$chisq,logrankp=logrank$pvalue,plotfit=plotfit)}

#version 202407230850
#baseline analysis annotation by 
#basetable statistics
chen_basetable=function(data,allvars,targetfactor,fvars,addOverall=F,includeNA=F,
testApprox=chisq.test,testExact=fisher.test,
testNormal=oneway.test,testNonNormal=kruskal.test){
if (!require("tableone",quietly=TRUE)){install.packages("tableone")}
library(tableone)
Svytab=CreateTableOne(vars=allvars,strata=targetfactor,data=data,addOverall=addOverall,
factorVars=fvars,includeNA=includeNA,testApprox=testApprox, testExact=testExact,
testNormal=testNormal, testNonNormal=testNonNormal)
Svytab}

#basetable statistics plus
chen_mbasetable=function(data,x=".",y="",z=F,selec=NA,subset=NULL,method=1){
if (!require("compareGroups", quietly = TRUE)){install.packages("compareGroups")}
library(compareGroups)
formulas=as.formula(paste(y,paste(x,collapse="+"),sep="~"))
restab=descrTable(formulas, data=data,selec=selec,subset=subset,method=method)
if(as.logical(z)){classtab=strataTable(restab, z)
list(restab=restab,classtab=classtab)}else{restab}}
#export2word(restab, file='table1.docx')
