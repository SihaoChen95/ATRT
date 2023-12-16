library(glmnet)
library(survival)
library(foreign)
library(ggplot2)
library(randomForestSRC)
mycol<-c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
         "#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
         "#E41A1C","#377EB8","#4DAF4A","#FF7F00","#FFFF33" )
mydata<-read.spss("C:/Users/CSH/Desktop/ATRT训练集（随机森林算法）.sav")
mydata<-as.data.frame(mydata)
head(mydata)   
str(mydata)
rf.model<- rfsrc(Surv(times,OS) ~ ., data =mydata,
                 ntree = 1000, 
                 splitrule="logrank",
                 importance = TRUE,
                 nodesize = 15)
rf.model

plot(get.tree(rf.model,3))
## print results of trained forest
print(rf.model)
## plot results of trained forest
plot(rf.model)

rf.top<-var.select(rf.model)
rf.top

rf.top2<-data.frame(
  Feature=rf.top$topvars,
  vimp=rf.top$varselect[rf.top$topvars,2])

ggplot(aes(x=reorder(Feature,vimp),y=vimp,fill=Feature),data=rf.top2)+
  geom_col()+
  coord_flip()+
  theme_bw()+
  labs(x="")+
  ggtitle("VIP gene in Randomforest")+
  scale_fill_manual(values = mycol)+
  theme(legend.position = "")
