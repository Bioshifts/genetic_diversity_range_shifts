#Figures: effect size of |VIS|*|GD|
################################################################################
#required packages
list.of.packages <- c(
    "ggplot2","gridExtra","grid","lattice") 


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

#data
setwd("/home/rbertrand/W/Bioshift/GD_study/boot_analysis/REMLF2/") #change path according to its location on your own computer
t1=read.csv2("coeff_allW.csv",sep=";",dec=".",h=T)
nB=unique(t1$nB)
v1=unique(t1$var)
t2=t1
res=data.frame(var="Intercept",Estimate=t2$Estimate[t2$var=="(Intercept)"],nB,ordV=NA,Param='CE') #balanced and presence/absence data
res=rbind(res,
          data.frame(var="Intercept",Estimate=t2$Estimate[t2$var=="(Intercept)"]+t2$Estimate[t2$var=="ParamLE"],nB,ordV=NA,Param='LE'),
          data.frame(var="Intercept",Estimate=t2$Estimate[t2$var=="(Intercept)"]+t2$Estimate[t2$var=="ParamTE"],nB,ordV=NA,Param='TE'))
res=rbind(res,
          data.frame(var="|VIS|",Estimate=t2$Estimate[t2$var=="scale(vel_abs)"],nB,ordV=4,Param='CE'), 
          data.frame(var="|VIS|",Estimate=t2$Estimate[t2$var=="scale(vel_abs)"]+t2$Estimate[t2$var=="scale(vel_abs):ParamLE"],nB,ordV=7,Param='LE'),
          data.frame(var="|VIS|",Estimate=t2$Estimate[t2$var=="scale(vel_abs)"]+t2$Estimate[t2$var=="scale(vel_abs):ParamTE"],nB,ordV=1,Param='TE'))
res=rbind(res,
          data.frame(var="GD",Estimate=t2$Estimate[t2$var=="scale(GD)"],nB,ordV=5,Param='CE'), 
          data.frame(var="GD",Estimate=t2$Estimate[t2$var=="scale(GD)"]+t2$Estimate[t2$var=="scale(GD):ParamLE"],nB,ordV=8,Param='LE'),
          data.frame(var="GD",Estimate=t2$Estimate[t2$var=="scale(GD)"]+t2$Estimate[t2$var=="scale(GD):ParamTE"],nB,ordV=2,Param='TE'))

res=rbind(res,
          data.frame(var="|VIS|:GD",Estimate=t2$Estimate[t2$var=="scale(vel_abs):scale(GD)"],nB,ordV=6,Param='CE'), 
          data.frame(var="|VIS|:GD",Estimate=t2$Estimate[t2$var=="scale(vel_abs):scale(GD)"]+t2$Estimate[t2$var=="scale(vel_abs):scale(GD):ParamLE"],nB,ordV=9,Param='LE'),
          data.frame(var="|VIS|:GD",Estimate=t2$Estimate[t2$var=="scale(vel_abs):scale(GD)"]+t2$Estimate[t2$var=="scale(vel_abs):scale(GD):ParamTE"],nB,ordV=3,Param='TE'))

res$ID=paste0(res$var,"_",res$Param)
s1=10 #the seed to make reproducible random staff
r2=read.csv2("R2_allW_sel10000.csv",sep=";",dec=".",h=T) #10000 bootstrapped models

res=merge(res,data.frame(nB=r2$nB),by.x="nB",by.y="nB")
res=subset(res,is.na(ordV)==F)
res=res[order(res$ordV,decreasing=F),]
boxplot(Estimate~as.factor(ordV),var,data=res)
rr3=res
#bootstrap test: H0= coefficient equal to 0; H1= coefficient higher than 0
test1=function(x,mu=0){
    return(1-(length(x[x>mu])/length(x)))
}

#bootstrap test: H0= coefficient equal to 0; H1= coefficient lower than 0
test2=function(x,mu=0){
    return(1-(length(x[x<mu])/length(x)))
}

res=data.frame(var=names(tapply(rr3$Estimate,rr3$ID,mean)),moy=tapply(rr3$Estimate,rr3$ID,mean),sd=tapply(rr3$Estimate,rr3$ID,sd),median=tapply(rr3$Estimate,rr3$ID,median),q025=tapply(rr3$Estimate,rr3$ID,quantile,probs=0.025),
               q975=tapply(rr3$Estimate,rr3$ID,quantile,probs=0.975), pv.inf0=tapply(rr3$Estimate,rr3$ID,test2,mu=0),pv.sup0=tapply(rr3$Estimate,rr3$ID,test1,mu=0))
rr3$expEst=exp(rr3$Estimate)
boxplot(expEst~as.factor(ordV),var,data=rr3)
a=0
for(i in unique(res$var)){
    rr3a=subset(rr3,Estimate>=res$q025[res$var==i] & Estimate<=res$q975[res$var==i] & ID==i)
    if(a==0){
        rr4=rr3a
    }else{
        rr4=rbind(rr4,rr3a)
    }
    a=a+1
}

#figure v1
theme_set(
    theme_classic() +
        theme(legend.position = "top")
)

rr3$ordV2=1
rr3$ordV2[rr3$Param=="CE"]=2
rr3$ordV2[rr3$Param=="LE"]=3
rr4$ordV2=1
rr4$ordV2[rr4$Param=="CE"]=2
rr4$ordV2[rr4$Param=="LE"]=3
rr4$ordVX=1
rr4$ordVX[rr4$var=="GD"]=2
rr4$ordVX[rr4$var=="|VIS|:GD"]=3
rr4$ordVX=as.factor(rr4$ordVX)
rr4$ordV2=as.factor(rr4$ordV2)

gg1 <- ggplot(rr4, aes(x = ordV2, y = expEst,fill = ordVX))+
    geom_hline(yintercept=1, linetype = "dashed")

gg1 <- gg1 + geom_violin(trim = T,position = position_dodge(0.8), width=1.5 )+
    stat_summary(fun = "median",
                 data=rr4,
                 geom = "point",
                 color = "black",
                 position = position_dodge(0.8),
                 size=0.65) +
    scale_fill_manual(values = c("blanchedalmond", "steelblue1", "tan1"),
                      labels=c("Climate change velocity","Genetic diversity","Climate change velocity:Genetic diversity"))+
    theme(legend.position = c(0.6, 0.88),legend.title=element_blank(),legend.text = element_text(family="Arial",size = 8,colour="black"),
          legend.key.size = unit(0.5, 'cm'))+
    labs(x="", y="Effect size") +
    theme(axis.text = element_text(family="Arial",size = 9,colour="black")) +
    theme(axis.title = element_text(family="Arial",size = 10,colour="black")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(plot.title = element_text(hjust = 0.5,size=11))

gg1

gg2= gg1 + scale_x_discrete(breaks=c(1:3),labels=c("Trailing","Centroïd","Leading"))+
    scale_y_continuous(breaks=c(0.5,0.75,1,1.25,1.5,1.75),limits=c(0.45,1.85),expand=c(0,0),labels=as.character(c(0.5,0.75,1,1.25,1.5,1.75)))
gg2

png("/home/rbertrand/W/Bioshift/GD_study/boot_analysis/REMLF2/fig2_v1.png",unit="cm",width=9,height=9,res=300)#,width=547,height=360
gg2
dev.off()

#figure v2
############################"
#|VIS| effect
rr4a=subset(rr4,var=="|VIS|")
rr3a=subset(rr3,var=="|VIS|")
gg1 <- ggplot(rr4a, aes(x = ordV2, y = expEst,fill=Param))+
    geom_hline(yintercept=1, linetype = "dashed")+
    geom_vline(xintercept=1:3, linetype = "dashed",col="lightgrey")

gg1 <- gg1 + geom_violin(data=rr4a,aes(x = ordV2, y = expEst,fill = Param),trim = T) +
    #theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"))+
    coord_flip(ylim=c(0.78,1.9))+ 
    stat_summary(fun = "median",
                 data=rr3a,
                 geom = "point",
                 color = "black") +
    #scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    scale_fill_manual(values = c("blanchedalmond", "steelblue1", "tan1"))+
    #geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)
    ggtitle("Climate change velocity")+
    theme(legend.position = "none") +
    labs(x="", y="") +
    theme(axis.text = element_text(size = 11)) +
    theme(axis.title = element_text(size = 11)) +
    #theme(legend.position = "right") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(plot.title = element_text(hjust = 0.5,size=12))

gg1= gg1 + scale_x_discrete(breaks=c(1:3),labels=c("Trailing","Centroïd","Leading")) +
    scale_y_continuous(breaks=c(0.8,1,1.2,1.4,1.6,1.8),limits=c(0,3),expand=c(0,0),labels=as.character(c(0.8,1,1.2,1.4,1.6,1.8)))+
    labs(tag = 'a')
#labs(tag = '(a)')

#gg1=gg1+      geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)
gg1

#GD effect
rr4a=subset(rr4,var=="GD")
rr3a=subset(rr3,var=="GD")
gg2 <- ggplot(rr4a, aes(x = ordV2, y = expEst,fill=Param))+
    geom_hline(yintercept=1, linetype = "dashed")+
    geom_vline(xintercept=1:3, linetype = "dashed",col="lightgrey")

gg2 <- gg2 + geom_violin(data=rr4a,aes(x = ordV2, y = expEst,fill=Param),trim = T) +
    #theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"))+
    coord_flip(ylim=c(0.72,1.28))+ 
    stat_summary(fun = "median",
                 data=rr3a,
                 geom = "point",
                 color = "black") +
    #scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    scale_fill_manual(values = c("blanchedalmond", "steelblue1", "tan1"))+
    ggtitle("Genetic diversity")+
    theme(legend.position = "none") +
    labs(x="", y="Estimate") +
    theme(axis.text = element_text(size = 11)) +
    theme(axis.title = element_text(size = 11)) +
    #theme(legend.position = "right") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line.y=element_blank(),axis.line.x = element_line(colour = "black"))+
    theme(plot.title = element_text(hjust = 0.5,size=12))

gg2 = gg2 +scale_x_discrete(breaks=c(1:3),labels=c("Trailing","Centroïd","Leading"))+
    scale_y_continuous(breaks=seq(0.7,1.3,by=0.1),limits=c(0,3),expand=c(0,0),labels=as.character(c(0.7,0.8,0.9,1,1.1,1.2,1.3)))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(tag = 'b')
#labs(tag = '(a)')
gg2

#|VIS|:GD effect
rr4a=subset(rr4,var=="|VIS|:GD")
rr3a=subset(rr3,var=="|VIS|:GD")
gg3 <- ggplot(rr4a, aes(x = ordV2, y = expEst,fill=Param))+
    geom_hline(yintercept=1, linetype = "dashed")+
    geom_vline(xintercept=1:3, linetype = "dashed",col="lightgrey")

gg3 <- gg3 + geom_violin(data=rr4a,aes(x = ordV2, y = expEst,fill=Param),trim = T) +
    #theme(panel.grid.major.x = element_line(colour = "black", linetype = "dashed"))+
    coord_flip(ylim=c(0.48,1))+ 
    stat_summary(fun = "median",
                 data=rr3a,
                 geom = "point",
                 color = "black") +
    #scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    scale_fill_manual(values = c("blanchedalmond", "steelblue1", "tan1"))+
    ggtitle("Climate change velocity:Genetic diversity")+
    theme(legend.position = "none") +
    labs(x="", y="") +
    theme(axis.text = element_text(size = 11)) +
    theme(axis.title = element_text(size = 11)) +
    #theme(legend.position = "right") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line.y=element_blank(),axis.line.x = element_line(colour = "black"))+
    theme(plot.title = element_text(hjust = 0.5,size=12))

gg3 = gg3 + scale_x_discrete(breaks=c(1:3),labels=c("Trailing","Centroïd","Leading"))+
    scale_y_continuous(breaks=seq(0.5,1,by=0.1),limits=c(0,3),expand=c(0,0),labels=as.character(c(0.5,0.6,0.7,0.8,0.9,1)))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(tag = 'c')
gg3

##############pooling figures together
#Créer un emplacement vide
blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

png("/home/rbertrand/W/Bioshift/GD_study/boot_analysis/REMLF2/fig2_v2.png",unit="cm",width=27,height=11,res=300)#,width=547,height=360
#gs <- lapply(1:4, function(ii) 
#grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
lay <- rbind(c(1,2,3,4))
grid.arrange(gg1,gg2,gg3, 
             ncol=3, nrow=1, widths=c(10.2,9,9), heights=c(11),layout_matrix = lay)
dev.off()

