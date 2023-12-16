library(glmnet)
library(survival)
library(foreign)
mydata<-read.spss("C:/Users/CSH/Desktop/ATRT训练集.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
head(mydata)   
str(mydata)
mydata$Race<-relevel(mydata$Race,ref = 'White')
mydata$Age<-relevel(mydata$Age,ref = '1-')
mydata$Sex<-relevel(mydata$Sex,ref = 'Male')
mydata$Householdincome<-relevel(mydata$Householdincome,ref = '75000+')
mydata$Grade<-relevel(mydata$Grade,ref = 'III-IV')
mydata$Site<-relevel(mydata$Site,ref = 'Brain')
mydata$Laterality<-relevel(mydata$Laterality,ref = 'Left')
mydata$SEERstage<-relevel(mydata$SEERstage,ref = 'Localized')
mydata$Tumorsize<-relevel(mydata$Tumorsize,ref = '4+')
mydata$S<-relevel(mydata$S,ref = 'Yes')
mydata$C<-relevel(mydata$C,ref = 'Yes')
mydata$R<-relevel(mydata$R,ref = 'Yes')
mydata$Years<-relevel(mydata$Years,ref = '2010+')
mydata$Treatment<-relevel(mydata$Treatment,ref = 'S')
str(mydata)
x<-as.matrix(mydata[,5:17])
y<-Surv(mydata$times,mydata$OS==1)
y<-Surv(mydata$times,mydata$CSS==1)
#Lasso??ģ
lasso<-glmnet(x,y,family = "cox",alpha = 1)
print(lasso)
plot(lasso, xvar = "lambda", label = TRUE)
lasso.coef<-predict(lasso, s=0.1197188, type = "coefficients")
lasso.coef
#??????֤
set.seed(123)
fitCV<-cv.glmnet(x,y,family= "cox",
                 type.measure = "deviance",
                 nfolds = 15)
plot(fitCV)
fitCV$lambda.1se
coef(fitCV,s= "lambda.1se")
fitCV$lambda.min
coef(fitCV,s= "lambda.min")

#cox回归模型
library(foreign)
library(DynNom)
library(glmnet)
library(rms)
library(survival)
library(survminer)
mydata<-read.spss("C:/Users/CSH/Desktop/ATRT验证集.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
head(mydata)
mydata$OS<-ifelse(mydata$OS=="Dead",1,0)
mydata$CSS<-ifelse(mydata$CSS=="Dead",1,0)
mydata$Race<-relevel(mydata$Race,ref = 'White')
mydata$Age<-relevel(mydata$Age,ref = '1-')
mydata$Sex<-relevel(mydata$Sex,ref = 'Male')
mydata$Grade<-relevel(mydata$Grade,ref = 'III-IV')
mydata$SEERstage<-relevel(mydata$SEERstage,ref = 'Localized')
mydata$S<-relevel(mydata$S,ref = 'Yes')
mydata$C<-relevel(mydata$C,ref = 'Yes')
mydata$R<-relevel(mydata$R,ref = 'Yes')
mydata$Householdincome<-relevel(mydata$Householdincome,ref = '75000-')
mydata$Laterality<-relevel(mydata$Laterality,ref = 'Left')
mydata$Site<-relevel(mydata$Site,ref = 'Brain')
mydata$Tumorsize<-relevel(mydata$Tumorsize,ref = 'others')
mydata$Years<-relevel(mydata$Years,ref = '2010+')
dd<-datadist(mydata)
options(datadist = 'dd')
###
mod<-cph(Surv(times,OS==1)~Age+SEERstage+Tumorsize+S+R+C,x=T,y=T,data = mydata,surv = T)
DynNom(mod) 
DNbuilder(mod)
options(encoding = "UTF-8")
install.packages('rsconnect')    
library(rsconnect)




coxm<-cph(Surv(times,CSS==1)~Age+SEERstage+Tumorsize+S+R+C,x=T,y=T,data = mydata,surv = T)
coxm1<-cph(Surv(times,CSS==1)~SEERstage+C,x=T,y=T,data = mydata,surv = T)
coxm2<-cph(Surv(times,OS==1)~Age+Tumorsize+S,x=T,y=T,data = mydata,surv = T)
library(Hmisc)
surv<-Survival(coxm)
rcorrcens(surv(times,CSS)~predict(coxm),data = mydata)
surv<-Survival(coxm1)
rcorrcens(surv(times,CSS)~predict(coxm1),data = mydata)
surv<-Survival(coxm2)
rcorrcens(surv(times,OS)~predict(coxm2),data = mydata)

surv1<-function(x)surv(1*12,lp=x)
surv2<-function(x)surv(1*24,lp=x)
surv3<-function(x)surv(1*36,lp=x)
nom<-nomogram(coxm,fun = list(surv1,surv2,surv3),lp = F,funlabel = c("1-Year OS",'2-Year OS','3-Year OS'),maxscale = 100,fun.at = c('0.95','0.85','0.80','0.70','0.60','0.50','0.4','0.3','0.2','0.1'))
plot((nom),xfrac=.3)
plot(nom,
     xfrac = .35,
     cex.var = 1,
     cex.axis = 0.8,
     tcl = -0.5,
     lmgp = 0.3,
     label.every = 1,
     naxes = 13,
     col.grid = gray(c(0.8,0.95)),
     lplabel = "Liner Predictorlp",
     points.label = 'Points',
     total.points.label = 'Total Points',
     force.label = T)



cal<-calibrate(coxm,cmethod = 'KM',method = 'boot',u=22,m=20,B=1000)
plot(cal,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim=c(0.0,1),ylim=c(0.0,1),
     xlab="Nomogram-predicted Probability of 1-Year OS",
     ylab="Actual 1-Year OS(proportion)",
     col=c(rgb(192,98,83,maxColorValue = 255)))
lines(cal[,c("mean.predicted","KM")],type = "b",lwd=2,col=c(rgb(192,98,83,maxColorValue = 255)),pch=16) 
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue = 255)))


#####################???????????��?##############################
library('nomogramFormula')
ret = nomogramFormula::TotalPoints.rms(rd=mydata,fit=coxm,nom=nom)
ret$'total points'
write.csv(ret,'ATRT??????SEER.csv')
#####################璁＄畻姣忎釜寰楀???##############################

#ROC
library(timeROC)
mydata$model1 = predict(coxm,data = mydata)
mydata$model2 = predict(coxm1,data = mydata)
mydata$model3 = predict(coxm2,data = mydata)
roc1 = timeROC(T=mydata$times,
               delta = mydata$OS,
               marker = mydata$model1,
               cause = 1,
               times = 12,
               ROC =TRUE,
               iid = TRUE
)
roc1$AUC
roc2 = timeROC(T=mydata$times,
               delta = mydata$OS,
               marker = mydata$model2,
               cause = 1,
               times = 12,
               ROC =TRUE,
               iid = TRUE
)
roc2$AUC     
roc3 = timeROC(T=mydata$times,
               delta = mydata$OS,
               marker = mydata$model3,
               cause = 1,
               times = 12,
               ROC =TRUE,
               iid = TRUE
)
roc3$AUC   
pdf('ROC.pdf')
plot(roc1,time=12,col='red',title ='ROC',lw=2)
plot(roc2,time=12,add=T,col='green',lw=2)
plot(roc3,time=12,add=T,col='blue',lw=2)
legend('bottomright',
       c(paste('model1 AUC',round(roc1$AUC,3)[2],sep=':'),
         paste('model2 AUC',round(roc2$AUC,3)[2],sep=':'))
       ,lw=2,col=c('red','green','blue'))

dev.off()
compare(roc1,roc2)




##NRI
library(nricens)
library(survival)
library(foreign)
library(rms)
mydata<-read.spss("C:/Users/CSH/Desktop/???ʰ????ⲿ??֤??.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
mydata$OS<-ifelse(mydata$OS=="Dead",1,0)
names(mydata)
time<-mydata$times
status<-mydata$OS
j1<-as.matrix(subset(mydata,select = c(clinical)))
j2<-as.matrix(subset(mydata,select = c(Age,Site,clinical,Tstage,N,Treatment)))
mod.std<-coxph(Surv(time,status)~.,data.frame(time,status,j1),x=TRUE)
mod.new<-coxph(Surv(time,status)~.,data.frame(time,status,j2),x=TRUE)
p.std=get.risk.coxph(mod.std,t0=60)
p.new=get.risk.coxph(mod.new,t0=60)
nricens(mdl.std = mod.std,mdl.new = mod.new,t0=60,cut = c(0.2,0.4),
        niter = 100,alpha = 0.05,updown = 'category')


##IDIָ??
library(survIDINRI)
t0=60
x<-IDI.INF(mydata[,2:1],p.std,p.new,t0,npert = 100)
IDI.INF.OUT(x)
IDI.INF.GRAPH(x)



##计算临床决策曲线 失败???
library(dcurves)
library(foreign)
library(survival)
library(dplyr)
mydata<-read.spss("C:/Users/CSH/Desktop/下咽.sav",use.value.labels = F,to.data.frame = T)
mydata<-na.omit(mydata)
mydata$Age<-as.factor(mydata$Age)
mydata$Race<-as.factor(mydata$Race)
mydata$Site<-as.factor(mydata$Site)
mydata$Tstage<-as.factor(mydata$Tstage)
mydata$N<-as.factor(mydata$N)
mydata$M<-as.factor(mydata$M)
mydata$Treatment<-as.factor(mydata$Treatment)
f1<-coxph(Surv(Months,Status)~Tstage+N+M,mydata)
f2<-coxph(Surv(Months,Status)~Tstage+N+M+Age+Race+Site+Treatment,mydata)
mydata$pr_failuref136=c(1-(summary(survfit(f1,newdata=mydata),Months=36)$surv))
mydata$pf236=c(1-(summary(survfit(f1,newdata=mydata),Months=36)$surv))
mydata$pf336=c(1-(summary(survfit(f1,newdata=mydata),Months=36)$surv))



##?߼??ع???DCA
library(rms)
library(ggDCA)
library(survival)  
library(foreign)
rm(list = ls())     
mydata<-read.spss("C:/Users/CSH/Desktop/超声.sav")
mydata<-as.data.frame(mydata)
head(mydata)
attach(mydata)
dd<-datadist(mydata)
options(datadist='dd')
model1<-coxph(Surv(Months,Status==1)~Tstage+N+M+Age+Race+Site+Treatment,data=mydata)        
dca1<-dca(model1,
          new.data = NULL,
          Months=60)      


#????????DCA
library(rms)
library(foreign)
library(ggDCA)
library(ggplot2)
library(survival)
library(survminer)
mydata<-read.spss("C:/Users/CSH/Desktop/ATRTѵ��??.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
head(mydata)
mydata$OS<-ifelse(mydata$OS=="Dead",1,0)
mydata$CSS<-ifelse(mydata$CSS=="Dead",1,0)
mydata$Race<-relevel(mydata$Race,ref = 'White')
mydata$Age<-relevel(mydata$Age,ref = '1-')
mydata$Sex<-relevel(mydata$Sex,ref = 'Male')
mydata$Grade<-relevel(mydata$Grade,ref = 'III-IV')
mydata$SEERstage<-relevel(mydata$SEERstage,ref = 'Localized')
mydata$S<-relevel(mydata$S,ref = 'Yes')
mydata$C<-relevel(mydata$C,ref = 'Yes')
mydata$R<-relevel(mydata$R,ref = 'Yes')
mydata$Householdincome<-relevel(mydata$Householdincome,ref = '75000-')
mydata$Laterality<-relevel(mydata$Laterality,ref = 'Left')
mydata$Site<-relevel(mydata$Site,ref = 'Brain')
mydata$Tumorsize<-relevel(mydata$Tumorsize,ref = 'others')
mydata$Years<-relevel(mydata$Years,ref = '2010+')
mydata$scoreOS<-relevel(mydata$scoreOS,ref = 'Low-risk')
mydata$scoreCSS<-relevel(mydata$scoreCSS,ref = 'Low-risk')
dd<-datadist(mydata)
options(datadist = 'dd')
M1<-cph(Surv(times,Status)~Tstage+N,data=mydata)
M2<-cph(Surv(times,Status)~Age+Site+Race+Tstage+N+Treatment,data=mydata)
d<-dca(M1,M2,
       times=c(1))

ggplot(d,linetype=1)  

#KM???????߻???
fit<-survfit(Surv(times,OS)~scoreOS,data = mydata)
plot(fit)
p1<-ggsurvplot(fit)
p1
#??????
ggsurvplot(fit,data = mydata,
           legend.title="Risk group",
           legend.labs=c("Low-risk","High-risk"),
           pval=TRUE,
           pval.method=TRUE,
           risk.table=TRUE,
           tables.height=0.2,
           tables.theme=theme_cleantable("Risk group"),
           palette=c("#E7b800","#2E9FDF"),
           ggtheme=theme_bw())
#??????
ggsurvplot(fit,data = mydata,
           pval=TRUE,
           pval.method=TRUE,
           risk.table=TRUE,
           risk.table.col="RiskscoreOS",
           linetype = "RiskscoreOS",
           tables.theme=theme_cleantable("RiskscoreOS"),
           palette=c("#E7b800","#2E9FDF","#00AFBB","#B53DB2","#FF5044","#B53DB2","#7764F6"),
           ggtheme=theme_bw())

#??????
ggsurvplot(fit,data = mydata,
           pval=TRUE,
           pval.method=TRUE,
           risk.table=TRUE,
           risk.table.col="Age",
           linetype = "Age",
           tables.theme=theme_cleantable("Age"),
           palette=c("#E7b800","#2E9FDF","#00AFBB","#FF5044","#7764F6","#B53DB2"),
           ggtheme=theme_bw())



#?????ࣨ??????50%OS?ĺ??߼?95%?????????䣩
ggsurvplot(fit,data = mydata,
           surv.median.line = "hv",
           pval=TRUE,
           pval.method=TRUE,
           conf.int=TRUE,
           risk.table=TRUE,
           risk.table.col="scoreOS",
           linetype = "scoreOS",
           tables.theme=theme_cleantable("scoreOS"),
           palette=c("#E7b800","#2E9FDF","#00AFBB","#B53DB2","#7764F6","#FF5044"),
           ggtheme=theme_bw())



##ɭ??ͼ##
library(survminer)
library(survival)
library(foreign)
library(DynNom)
library(rms)
mydata<-read.spss("C:/Users/CSH/Desktop/ATRTѵ��??.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
head(mydata)
mydata$OS<-ifelse(mydata$OS=="Dead",1,0)
mydata$CSS<-ifelse(mydata$CSS=="Dead",1,0)
mydata$Race<-relevel(mydata$Race,ref = 'White')
mydata$Age<-relevel(mydata$Age,ref = '1-')
mydata$Sex<-relevel(mydata$Sex,ref = 'Male')
mydata$Householdincome<-relevel(mydata$Householdincome,ref = '75000+')
mydata$Grade<-relevel(mydata$Grade,ref = 'III-IV')
mydata$Site<-relevel(mydata$Site,ref = 'Brain')
mydata$Laterality<-relevel(mydata$Laterality,ref = 'Left')
mydata$SEERstage<-relevel(mydata$SEERstage,ref = 'Localized')
mydata$Tumorsize<-relevel(mydata$Tumorsize,ref = '4+')
mydata$S<-relevel(mydata$S,ref = 'Yes')
mydata$C<-relevel(mydata$C,ref = 'Yes')
mydata$R<-relevel(mydata$R,ref = 'Yes')
mydata$Years<-relevel(mydata$Years,ref = '2010+')
mydata$Treatment<-as.numeric(mydata$Treatment)
str(mydata)
dd<-datadist(mydata)
options(datadist = 'dd')

model <- coxph( Surv(times, CSS) ~ Age + SEERstage + Tumorsize + S + R + C,
                data = mydata )
model
ggforest(model)
ggforest(model,
         main = "Hazard ratio CSS", # ???ñ???
         cpositions = c(0.06, 0.2, 0.35), # ????ǰ???е????Ծ???
         fontsize = 1.0, # ??????????С
         refLabel = "reference",
         noDigits = 2) #???ñ???С????λ??

#DCA????
setwd("C:/Users/CSH/Desktop")
source("stdca.R") #һ??Ҫ??stdca.R????֮ǰ?趨????ʼĿ¼?С?
dev<-read.csv("data.csv")
head(dev)
str(dev)
library(rms)
library(foreign)
library(survival)

Srv=Surv(dev$time,dev$Status)
coxmod=coxph(Srv~Gender+Site+SEERstage+Tumorsize+R+S+Age,data=dev)
dev$one.years.Survival.Probabilitynew = c(1- (summary(survfit(coxmod,newdata=dev),times=12)$surv))  #????1????????
dev$two.years.Survival.Probabilitynew = c(1- (summary(survfit(coxmod,newdata=dev),times=24)$surv))  #????2????????
dev$three.years.Survival.Probabilitynew = c(1- (summary(survfit(coxmod,newdata=dev),times=36)$surv))  #????3????????
write.csv(dev, "devnnew.csv") 
stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=12,predictors="one.years.Survival.Probabilitynew",xstop=0.6,smooth=TRUE)
stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=24,predictors="two.years.Survival.Probabilitynew",xstop=0.6,smooth=TRUE)
stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=36,predictors="three.years.Survival.Probabilitynew",xstop=0.6,smooth=TRUE)

#��ģ?ͱȽ?
coxmod1<-coxph(Srv ~Gender+Site+SEERstage+Tumorsize+R+S+Age,data=dev)
coxmod2<-coxph(Srv ~SEERstage,data=dev)

dev$model1<-c(1- (summary(survfit(coxmod1,newdata=dev),times=36)$surv))

dev$model2<-c(1- (summary(survfit(coxmod2,newdata=dev),times=36)$surv))




stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=36,predictors=c("model1"),smooth=TRUE)
stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=36,predictors=c("model1","model2"),xstop = 0.8,smooth=TRUE)

