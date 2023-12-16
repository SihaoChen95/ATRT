library(shinyPredict)
library(DynNom)
library(foreign)
library(glmnet)
library(rms)
library(survival)
library(survminer)
library(survival) # 加载包
mydata<-read.csv("C:/Users/CSH/Desktop/ATRT.csv")
View(mydata) # 预览数据集
mydata$status <- as.numeric(mydata$status)
mydata$age <- factor(mydata$age) 
mydata$seerstage <- factor(mydata$seerstage) 
mydata$tumorsize <- factor(mydata$tumorsize) 
mydata$surgery <- factor(mydata$surgery) 
mydata$radiation <- factor(mydata$radiation) 
mydata$chemotherapy <- factor(mydata$chemotherapy) 


tmp.m3 <- coxph(Surv(time , status ) ~ age + seerstage + tumorsize + surgery + radiation + chemotherapy,
                data=mydata, 
                model = FALSE, y=FALSE)

tmp.m4 <- coxph(Surv(time , status ) ~ age + seerstage + tumorsize,
                data = mydata,
                model = FALSE, y=FALSE)

tmp.m5 <- coxph(Surv(time , status ) ~ age + seerstage + tumorsize + surgery + radiation + chemotherapy,
                data=mydata, 
                model = FALSE, y=FALSE)

options(encoding = "UTF-8")
deployApp(appName = "DynNomappforATRTinOS")
install.packages('rsconnect')    
library(rsconnect)
rsconnect::setAccountInfo(name='atrt', token='DC7AF2A7080BD08B35F86782579AF8F5', secret='hY4UBje4UfslaHCbnv5WCz/4RgjUKuUgoVj/3GsG')

shinyPredict(models=list("Model 1"= tmp.m3), 
             data=mydata[, c("time","status","age","seerstage","tumorsize","surgery","radiation","chemotherapy")], # 创建动态列线图的数据
             path = "C:\\Users\\CSH\\Documents", # 动态列线图shiny app文件存放位置
             title="Predicting ATRT CSS mortality") # 列线图的名称

