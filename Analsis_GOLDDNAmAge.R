
setwd("~/Desktop/BioageV2/Result4/")
################################
### https://ngdc.cncb.ac.cn/ewas/datahub/download
load("DNAmAge/Data/age_methylation_v1/age_methylation_v1.RData")
head(age_download[,1:5])
dim(age_download)

nas_stat <- function(x){
  length(x[is.na(x)])
}

nl=apply(age_download, 1, nas_stat)
length(nl[nl<838)]) ### 10% missing rate
age_download=age_download[which(nl<838),]
dim(age_download)

load("DNAmAge/Data/sample_age_methylation_v1/sample_age.RData")
head(sample_age)
dim(sample_age)

a= -9.367827
b= 0.08002944
h=exp(a)*exp(b*(sample_age$age))
summary(h)

cl=which(sample_age$age>=18)
length(cl)

k=1
pl=list()
chunk_size <- 10000
n_rows <- nrow(age_download)
ii <- split(1:n_rows, ceiling(seq_along(1:n_rows) / chunk_size))
ii

for(i in ii){
  tt=age_download[i, cl]
  #tt=apply(tt, 2, as.numeric)
  print(k)
  pl[[k]]=apply(tt,1, function(x) WGCNA::cor(as.numeric(x), sample_age$age[cl],use="pairwise.complete.obs"))
  k=k+1
}

pl1=unlist(pl)
length(pl1[which(abs(pl1)>0.3)]) ## 0.1
pl2=pl1[which(abs(pl1)>0.3)]
head(pl2)
length(pl2)

length(cl)
summary(sample_age$age[cl])
sd(sample_age$age[cl])

###########
#setwd("~/Desktop/BioageV2/Result4/")
#load("DNAmAge/Data/age_methylation_v1/age_methylation_v1.RData")
#head(age_download[,1:5])
#dim(age_download)

load("DNAmAge/Data/sample_age_methylation_v1/sample_age.RData")
head(sample_age)
dim(sample_age)

a= -9.367827
b= 0.08002944
h=exp(a)*exp(b*(sample_age$age))
summary(h)

#cpgs=names(pl2)[2:length(pl2)]
#saveRDS(cpgs,"Resi/Included.Cpgs3.rds")

cpgs=readRDS("Resi/Included.Cpgs3.rds")
length(cpgs)

library(biglasso)
x=age_download[rownames(age_download) %in% cpgs,cl]
dim(x)

x=data.frame(t(x))
#saveRDS(x,"Resi/Included.Cpgs.Beta2.rds")
#saveRDS(x,"Resi/Included.Cpgs.Beta.rds")

#x=readRDS("Resi/Included.Cpgs.Beta.rds")
x=imputeMissings::impute(x)
head(x[,1:5])

X.bm <- as.big.matrix(as.matrix(x))
head(X.bm[,1:5])

y=log(h)[cl]
summary(y)

fit <- biglasso(X.bm, y)
plot(fit)

cvfit <- tryCatch(
  {
    cv.biglasso(X.bm, y, seed = 1234, nfolds = 5, ncores = 8,lambda.min = 1e-3, nlambda = 100)
  },
  error = function(cond) {
    cv.biglasso(X.bm, y, seed = 1234, nfolds = 5, ncores = 8,lambda.min = 1e-3, nlambda = 100)
  }
)

summary(cvfit)
saveRDS(cvfit,"Resi/Lasso.cvfit3.rds")

cvfit=readRDS("Resi/Lasso.cvfit3.rds")
pdf("Resi/Lasso.meth2.pdf",height = 4,width = 6)
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1) ,mgp = c(2.5, 0.5, 0))
plot(cvfit, type = "all",lambda.min = 1e-5)
dev.off()

yhat=predict(cvfit,X.bm)
plot(exp(yhat[,1]),exp(y))

pd=data.frame(yh=yhat[,1],DNAmAge,age=sample_age$age[cl])
l=lm(yh~age,data=pd)
s=summary(l)
s
a1=s$coefficients[1,1]
b1=s$coefficients[2,1]
resi=yhat[,1]-b1*sample_age$age[cl]

plot(sample_age$age[cl],resi)
plot(sample_age$age[cl],yhat[,1])
plot(sample_age$age[cl],exp(yhat[,1]))

DNAmAge=log(exp(yhat[,1])/exp(a1))/b1
plot(DNAmAge,sample_age$age[cl])
cor.test(DNAmAge,sample_age$age[cl])
################ 


y1=resi

k=1
pl=list()
chunk_size <- 10000
n_rows <- nrow(age_download)
ii <- split(1:n_rows, ceiling(seq_along(1:n_rows) / chunk_size))
ii

for(i in ii){
  tt=age_download[i, cl]
  #tt=apply(tt, 2, as.numeric)
  print(k)
  pl[[k]]=apply(tt,1, function(x) WGCNA::cor(as.numeric(x), y1,use="pairwise.complete.obs"))
  k=k+1
}

sapply(pl, length)
pl1=unlist(pl)
length(pl1[which(abs(pl1)>0.1)]) ## 0.1
pl2=pl1[which(abs(pl1)>0.1)]
head(pl2)
length(pl2)

############################
library(biglasso)
cpgs1=names(pl2)[2:length(pl2)]
cpgs1=cpgs1[grep("^cg",cpgs1)]
length(cpgs1)

x1=age_download[rownames(age_download) %in% cpgs1,cl]
dim(x1)

x1=data.frame(t(x1))
x1=imputeMissings::impute(x1)

X.bm1 <- as.big.matrix(as.matrix(x1))
cvfit2 <- tryCatch(
  {
    cv.biglasso(X.bm1, y1, seed = 1234, nfolds = 5, ncores = 8,lambda.min = 1e-3, nlambda = 100)
  },
  error = function(cond) {
    cv.biglasso(X.bm1, y1, seed = 1234, nfolds = 5, ncores = 8,lambda.min = 1e-3, nlambda = 100)
  }
)

summary(cvfit2)
#saveRDS(cvfit2,"Resi/Lasso.cvfit.resi.rds")

cvfit2=readRDS("Resi/Lasso.cvfit.resi.rds")
pdf("Resi/Lasso.meth3.pdf",height = 4,width = 6)
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1) ,mgp = c(2.5, 0.5, 0))
plot(cvfit2, type = "all",lambda.min = 1e-5)
dev.off()

####################
coef=coef(cvfit)[,1]
coef=coef[coef!=0]
coef2=coef(cvfit2)[,1]
coef2=coef2[coef2!=0]

resi1=predict(cvfit2,X.bm1)
GDNAmAge=sample_age$age[cl] + (resi1-a1)/b1
pd$GDNAmAge=GDNAmAge[,1]
pd$id=sample_age$sample_id[cl]
pd$geo=rownames(sample_age)[cl]
pd$sex=sample_age$sex[cl]

plot(sample_age$age[cl],GDNAmAge)
plot(sample_age$age[cl],(resi1-a1)/b1)
cor.test(sample_age$age[cl],as.numeric((resi1-a1)/b1))

coef3=coef2
coef3[1]=(coef3[1]-a1)
coef3=coef3/b1

x2=as.matrix(x1[,match(names(coef3)[2:length(coef3)],colnames(x1))])
x2=apply(x2, 2, as.numeric)
resi2=x2 %*% coef3[2:length(coef3)] + coef3[1]
summary(resi2)

plot(sample_age$age[cl],sample_age$age[cl]+resi2)

x3=as.matrix(X.bm[,match(names(coef)[2:length(coef)],colnames(X.bm))])
x3=apply(x3, 2, as.numeric)
dnamage=log(exp(x3 %*% coef[2:length(coef)] + coef[1])/exp(a1))/b1
summary(dnamage)

plot(sample_age$age[cl],dnamage)

saveRDS(coef,"Resi/DNAmAge.cpgs.rds")
saveRDS(coef3,"Resi/GOLDDNAmAge.cpgs.rds")

write.table(coef,"Resi/DNAmAge.cpgs.txt",col.names = F,sep = "\t")
write.table(coef3,"Resi/GOLDDNAmAge.cpgs.txt",col.names = F,sep="\t")

##################
val.sam=read.table("Resi/sample_cancer.txt")

val.sam=val.sam[which(val.sam$age>=18),]
dim(val.sam)
tx=val.sam1$sample_id[val.sam1$sample_type=="control"]
length(tx)
rep=tx[tx %in% colnames(age_download)]

val.sam=val.sam[!(val.sam$sample_id %in% rep),]
dim(val.sam)
summary(val.sam$age)

val.meth=read.csv("Resi/methylation_matrix_large.csv")
colnames(val.meth)

nl1=apply(val.meth, 1, nas_stat)
summary(nl1)
val.meth=val.meth[which(nl1<92),]
dim(val.meth)
com.sam=val.meth[,1][val.meth[,1] %in% val.sam$sample_id]
length(com.sam)

val.sam1=val.sam[match(com.sam,val.sam$sample_id),]
val.meth=val.meth[match(com.sam,val.meth[[1]]),]
dim(val.sam1)
dim(val.meth)

val.meth=imputeMissings::impute(val.meth)

cpgs.coef3=names(coef3)[2:length(coef3)]
coef31=cpgs.coef3[cpgs.coef3 %in% colnames(val.meth)]
length(coef31)

x2=as.matrix(val.meth[,match(coef31,colnames(val.meth))])
x2=apply(x2, 2, as.numeric)

resi2=x2 %*% coef3[match(coef31,names(coef3))] + coef3[1]
summary(resi2)

plot(val.sam1$age,resi2)
plot(val.sam1$age[val.sam1$sample_type=="control"],
     resi2[val.sam1$sample_type=="control"])

val.sam1$resi=resi2

########################################################
cpgs.coef=names(coef)[2:length(coef)]
coef11=cpgs.coef[cpgs.coef %in% colnames(val.meth)]
length(coef11)

x3=as.matrix(val.meth[,match(coef11,colnames(val.meth))])
x3=apply(x3, 2, as.numeric)
dnamage=log(exp(x3 %*% coef[match(coef11,names(coef))] + coef[1])/exp(a1))/b1
summary(dnamage)


plot(val.sam1$age,dnamage)
plot(val.sam1$age,val.sam1$age+val.sam1$resi)

val.sam1$dnamage=dnamage

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

x1 <- val.sam1$age[val.sam1$sample_type=="control"]
y1 <- (val.sam1$age+val.sam1$resi)[val.sam1$sample_type=="control"]
x2 <- val.sam1$age[val.sam1$sample_type=="control"]
y2 <- dnamage[val.sam1$sample_type=="control"]
x3 <- val.sam1$age[val.sam1$sample_type=="control"]
y3 <- val.sam1$resi[val.sam1$sample_type=="control"]
x4 <- val.sam1$age[val.sam1$sample_type=="control"]
y4 <- (dnamage-val.sam1$age)[val.sam1$sample_type=="control"]

add_cor_text <- function(x, y, pos = "topleft") {
  cor_val <- cor(x, y, use = "complete.obs")
  p_val <- cor.test(x, y)$p.value
  p_text <- ifelse(p_val < 0.001, "p < 0.001", 
                   ifelse(p_val < 0.01, paste("p =", round(p_val, 3)),
                          paste("p =", round(p_val, 2))))
  legend("topleft", pos, legend = bquote(italic(r) ~ "=" ~ .(round(cor_val, 3)) ~ .(p_text)), 
         bty = "n", cex = 1,inset = c(-0.1, 0.02))
}

plot(val.sam1$age[val.sam1$sample_type=="control"],
     (val.sam1$age+val.sam1$resi)[val.sam1$sample_type=="control"],
     xlab = "Chronological Age (years)",
     ylab = "GOLD DNAmAge (years)",
     pch = 16, col = "blue", cex = 0.6)
add_cor_text(x1, y1)

plot(val.sam1$age[val.sam1$sample_type=="control"],
     dnamage[val.sam1$sample_type=="control"],
     xlab = "Chronological Age (years)",
     ylab = "DNAmAge (years)",
     pch = 16, col = "red", cex = 0.6)
add_cor_text(x2, y2)

plot(val.sam1$age[val.sam1$sample_type=="control"],
     val.sam1$resi[val.sam1$sample_type=="control"],
     xlab = "Chronological Age (years)",
     ylab = "Residual (years)",
     pch = 16, col = "purple", cex = 0.6)
add_cor_text(x3, y3)

plot(val.sam1$age[val.sam1$sample_type=="control"],
     (dnamage-val.sam1$age)[val.sam1$sample_type=="control"],
     xlab = "Chronological Age (years)",
     ylab = "Acceleration (years)",
     pch = 16, col = "green", cex = 0.6)
add_cor_text(x4, y4)

par(mfrow = c(1, 1))
#######################


head(pd)
pd1=data.frame(DNAmAge=c(pd$DNAmAge,val.sam1$dnamage),
               age=c(pd$age,val.sam1$age),
               GDNAmAge=c(pd$GDNAmAge,val.sam1$resi+val.sam1$age),
               type=c(rep("Train",nrow(pd)),
                      rep("Test",nrow(val.sam1))))
pd1$type=factor(pd1$type,levels = c("Train","Test"))
library(ggpubr)
library(ggplot2)

p1=ggplot(pd1, aes(age, GDNAmAge)) +
  geom_point(aes(color = type), size = 0.5) +
  stat_cor(aes(group = type, color = type)) +
  theme_classic() + 
  ylab("GOLD DNAmAge") + 
  xlab("Age") +
  scale_y_continuous(expand = c(0.05,0.15)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.2)) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))
  
p2=ggplot(pd1,aes(age,DNAmAge)) +
  geom_point(aes(color = type), size = 0.5) +
  stat_cor(aes(group = type, color = type)) +
  theme_classic() + 
  ylab("DNAmAge") + 
  xlab("Age") +
  scale_y_continuous(expand = c(0.05,0.15)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.2)) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))


####################################
library(survival)
library(survminer)

val.sam2=val.sam1[!is.na(val.sam2$vital_status),]
dim(val.sam2)
table(val.sam2$vital_status)
val.sam2$vital_status=as.numeric(as.factor(val.sam2$vital_status)) -1
table(val.sam2$vital_status)

val.sam2$race[grep("Indian",val.sam2$race)]="Asian"
val.sam2$race[grep("Asian",val.sam2$race)]="Asian"

val.sam2$race[grep("black",val.sam2$race)]="Black"
val.sam2$race[grep("white",val.sam2$race)]="White"
table(val.sam2$race)
summary(val.sam2$overall_survival)
table(val.sam2$vital_status)
summary(val.sam2$age)
sd(val.sam2$age)

cc=coxph(Surv(overall_survival,vital_status)~age+resi+race+sex,data=val.sam2)
summary(cc)

cc1=coxph(Surv(overall_survival,vital_status)~age+dnamage+race+sex,data=val.sam2)
summary(cc1)

val.sam2$gbioage=val.sam2$resi #+ val.sam2$age
val.sam2$Groups=as.factor(cut(val.sam2$gbioage,breaks=quantile(val.sam2$gbioage)))
f <- survfit(Surv(overall_survival,vital_status) ~ Groups,data=val.sam2)
library(ggsci)


p3=ggsurvplot(f,
              data = val.sam2,
              xlim=c(0,5),
              ylim=c(0,0.75),
              break.x.by = 1,
              fun = "event",
              pval = TRUE,
              pval.coord = c(0.1,0.65),
              ylab = "Incidence",
              xlab = "Follow-ups (years)",
              legend.title = "Residual",
              legend.labs = paste("Q",1:4,sep=""),
              legend = c(0.8,0.2),
              censor  = FALSE,
              palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"),
              ggtheme = theme_classic())


val.sam2$dnamaccel=val.sam2$dnamage-val.sam2$age
val.sam2$Groups1=as.factor(cut(val.sam2$dnamaccel,breaks=quantile(val.sam2$dnamaccel)))
f <- survfit(Surv(overall_survival,vital_status) ~ Groups1,data=val.sam2)

p4=ggsurvplot(f,
              data = val.sam2,
              xlim=c(0,5),
              ylim=c(0,0.75),
              break.x.by = 1,
              fun = "event",
              pval = TRUE,
              pval.coord = c(0.1,0.65),
              ylab = "Incidence",
              xlab = "Follow-ups (years)",
              legend.title = "Accel",
              legend.labs = paste("Q",1:4,sep=""),
              legend = c(0.8,0.2),
              censor  = FALSE,
              palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"),
              ggtheme = theme_classic())

p31=p3$plot+theme(legend.key.size  = unit(0.7,"lines")) + guides(color=guide_legend(ncol = 2))
p41=p4$plot+theme(legend.key.size = unit(0.7,"lines"))+ guides(color=guide_legend(ncol = 2))

pdf("Resi/Figure.pdf",height = 5,width = 6)
ggarrange(p1,NULL,p2,
          p31,NULL,p41,
          labels = c("A","","B","C","","D"),
          widths = c(1,0.05,1),
          ncol = 3,nrow = 2)
dev.off()

train.sam=pd[,c("id","age","sex","GDNAmAge","DNAmAge")]
head(train.sam)

val.sam3=val.sam2[,c("sample_id","age","sex","resi","dnamage")]
val.sam3$resi=val.sam3$resi + val.sam3$age
colnames(val.sam3)=c("id","age","sex","GDNAmAge","DNAmAge")
head(val.sam3)

write.csv(train.sam,"Resi/Info.train.samples.csv",row.names = F)
write.csv(val.sam3,"Resi/Info.cancer.samples.csv",row.names = F)

###################
library(methylclockData)
library(methylclock)
setwd("~/Desktop/BioageV2/Result4/")
val.sam3=read.csv("Resi/Info.cancer.samples.csv")
#load("Resi/cancer_methylation_v1.RData")
cancer_download=readr::read_csv("Resi/methylation_matrix_large2.csv")
head(cancer_download[,1:5])
length(which(colnames(cancer_download) %in% val.sam3$id))

#cancer_download[,1]=rownames(cancer_download)
common.cpgs=readRDS("Resi/Common.clocks.cpgs.rds")
head(common.cpgs)
length(common.cpgs)

#cancer_download=cancer_download[rownames(cancer_download) %in% common.cpgs,]
cancer_download=cancer_download[,c(1,which(colnames(cancer_download) %in% val.sam3$id))]
dim(cancer_download)

cancer_download[,2:ncol(cancer_download)]=apply(cancer_download[,2:ncol(cancer_download)],2,as.numeric)
cancer_download[is.na(cancer_download)]

chunk_size <- 1000
n_cols <- ncol(cancer_download)
ii <- split(2:n_cols, ceiling(seq_along(1:n_cols) / chunk_size))
ii
DNAmClocks=list()
k=1
for(i in ii){
  DNAmClocks[[k]]=DNAmAge(cancer_download[,c(1,i)])
  print(k)
  k=k+1
}

DNAmClocks1 <- do.call(rbind, DNAmClocks)
head(DNAmClocks1)

val.sam=read.table("Resi/sample_cancer.txt")
val.sam=val.sam[match(val.sam3$id,val.sam$sample_id),]
val.sam3$time=val.sam$overall_survival
val.sam3$status=val.sam$vital_status
val.sam3$status=as.numeric(as.factor(val.sam3$status)) -1
val.sam3$sex=val.sam$sex
val.sam3$race=val.sam$race
val.sam3$race[grep("Indian",val.sam3$race)]="Asian"
val.sam3$race[grep("Asian",val.sam3$race)]="Asian"
val.sam3$race[grep("black",val.sam3$race)]="Black"
val.sam3$race[grep("white",val.sam3$race)]="White"
val.sam3=val.sam3[match(DNAmClocks1[[1]],val.sam3$id),]
val.sam3=cbind(val.sam3,data.frame(DNAmClocks1))


library(survival)
library(survival)
library(forestplot)

# 存储结果的列表
results <- list()
nclocks=c("GDNAmAge","DNAmAge","Horvath","Hannum","Levine","BNN","EN")
for(ck in nclocks){
  fl <- paste("Surv(time, status) ~ age +", ck, "+ race + sex", sep = "")
  cc1 <- coxph(as.formula(fl), data = val.sam3)
  s <- summary(cc1)
  ck_result <- s$coefficients[ck, ]
  results[[ck]] <- data.frame(
    Clock = ck,
    HR = exp(ck_result["coef"]),
    Lower_CI = exp(ck_result["coef"] - 1.96 * ck_result["se(coef)"]),
    Upper_CI = exp(ck_result["coef"] + 1.96 * ck_result["se(coef)"]),
    P_value = ck_result["Pr(>|z|)"]
  )
}

forest_data <- do.call(rbind, results)
rownames(forest_data) <- NULL
results = forest_data 

library(ggplot2)
results$Clock[1]="GOLD DNAmAge"
results$Clock=factor(results$Clock,levels = rev(results$Clock))

pdf("Resi/HR.pdf",height = 3,width = 4)
ggplot(results, aes(x = HR, y = Clock)) +
  geom_point(size = 3, color = "navy") +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(x = "Hazard Ratio (95% CI)",
       y = "Epigenetic Clock") +
  theme_classic2() + xlim(1,1.06) +
  geom_text(aes(x=1.02,label = sprintf("%.3f (%.3f-%.3f)", 
                                HR, Lower_CI, Upper_CI)),
            hjust = -0.2, size = 3)
dev.off()


library(ggplot2)
library(ggpubr)
library(patchwork)

plot_list <- list()

for(ck in nclocks[3:7]) {
  # 计算相关系数用于标题
  cor_test <- cor.test(val.sam3$age, val.sam3[[ck]], method = "pearson")
  cor_label <- sprintf("r = %.2f, p = %.3f", cor_test$estimate, cor_test$p.value)
  
  p <- ggplot(val.sam3, aes_string(x = "age", y = ck)) +
    geom_point(color = "grey80", size = 0.5, alpha = 0.6) +
    stat_cor(method = "pearson", size = 2.5, 
             label.x.npc = 0.2, label.y.npc = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", size = 0.5) +
    xlab("Chronological age") +
    ylab(ck) +
    ggtitle(ck) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          axis.title = element_text(size = 9))
  
  plot_list[[ck]] <- p
}

wrap_plots(plot_list, ncol = 3) + 
  plot_annotation(#title = "Correlation between Chronological Age and Epigenetic Clocks",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

vars=c(c("id","age","sex","time","status"),nclocks)
vars
val.sam31=val.sam3[,vars]
head(val.sam31)
write.csv(val.sam31,"Resi/Info.cancer.samples.csv",row.names = F)
