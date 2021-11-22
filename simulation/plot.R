setwd("E:/GWASmeta/sim-plot")
df1<-read.csv("sim1.csv")
colnames(df1)[3]<-"FEM"
colnames(df1)[4]<-"REM"
colnames(df1)[5]<-"EXH"
colnames(df1)[6]<-"MCMC"
colnames(df1)[7]<-"SSS"
library(ggplot2)

patte <- c("#333366","#CC6633","#669900","#666666","#CC6633","#993300")
my_theme <- theme(panel.background = element_blank(),
                  panel.border=element_rect(fill=NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background=element_blank(),
                  axis.text.x=element_text(colour="black",size=18),#angle = 20,hjust = 1),
                  axis.text.y=element_text(colour="black",size=25),
                  axis.title=element_text(colour="black",size=25,vjust=-1),
                  axis.ticks=element_line(colour="black"),
                  plot.margin=unit(c(1,1,1,1),"line"),
                  legend.title=element_blank(),
                  legend.key = element_rect(fill = "white"),
                  #legend.position = c(0.12,0.98),
                  #legend.justification = c(0.02,0.98),
                  legend.position = "right",
                  legend.text=element_text(size=25),
                  legend.key.size=unit(1, "cm"))
colours <- c("ABF(EXH)"="#996600","ABF(MCMC)"="#CC6633","ABF(SSS)"="#669900",
             "p-value(FEM)"="#666666","p-value(REM)"="#FFCC33")
boost <- 1600
ggplot(df1)+
  geom_point(aes(trueOR,MCMC,colour="ABF(MCMC)"),size=1)+
  geom_line(aes(trueOR,MCMC,colour="ABF(MCMC)"),size=1)+
  geom_point(aes(trueOR,EXH,colour="ABF(EXH)"),size=1)+
  geom_line(aes(trueOR,EXH,colour="ABF(EXH)"),size=1)+
  geom_point(aes(trueOR,SSS,colour="ABF(SSS)"),size=1)+
  geom_line(aes(trueOR,SSS,colour="ABF(SSS)"),size=1)+
  geom_point(aes(trueOR,FEM*boost,colour="p-value(FEM)"),size=1)+
  geom_line(aes(trueOR,FEM*boost,colour="p-value(FEM)"),size=1)+
  geom_point(aes(trueOR,REM*boost,colour="p-value(REM)"),size=1)+
  geom_line(aes(trueOR,REM*boost,colour="p-value(REM)"),size=1)+
  labs(x="true OR",y="log ABF")+
  scale_y_continuous(sec.axis=sec_axis(~./boost,name="p-value"))+
  scale_colour_manual(name="class",values=colours)+
  my_theme+
  #geom_line(aes(trueOR,log(1000),colour="ABF=1000"),size=1.5)+
  scale_x_continuous(breaks=seq(0.6,2.0,0.2))+
  guides(shape=guide_legend(override.aes=list(size=10)))+
  #scale_colour_manual(values=c("#CC6633","#333366","#669900","#666666","#FFCC33"))+
  scale_linetype_manual(values = c("ABF(EXH)"='dotted')) 

df2<-read.csv("sim2.csv",head=TRUE)
dfmean<-data.frame(value=unlist(df2[1,],use.names = FALSE))
dfmean$name <- c(rep(c("exh,indep,sigma=0.5","exh,fixed,sigma=0.5",
                       "exh,sigma=0.5,rho=0.7","exh,sigma=0.5,rho=0.3",
                       "exh,sigma=0.8,rho=0.7"),2),rep(c("mcmc,indep,sigma=0.5",
                                                   "mcmc,fixed,sigma=0.5","mcmc,sigma=0.5,rho=0.7",
                                                   "mcmc,sigma=0.5,rho=0.3","mcmc,sigma=0.8,rho=0.7"),2),
                 rep(c("shotgun,indep,sigma=0.5","shotgun,fixed,sigma=0.5",
                       "shotgun,sigma=0.5,rho=0.7","shotgun,sigma=0.5,rho=0.3",
                       "shotgun,sigma=0.8,rho=0.7"),2))
dfmean$algorithm <- c(rep("EXH",10),rep("MCMC",10),rep("SSS",10))
dfmean$prior <- rep(c("indep","fixed","corr 1","corr 2","corr 3"),3)
dfmeanvalue <- dfmean[c(seq(1,5),seq(11,15),seq(21,25)),]
dfmeantime <- dfmean[c(seq(6,10),seq(16,20),seq(26,30)),]
ggplot(dfmeanvalue,aes(prior,value))+
  geom_bar(stat="identity", aes(fill = algorithm),position="dodge")+
  #coord_flip()+
  ylab("mean logABF")+
  #xlab("method")+
  my_theme+
  scale_fill_manual(values=patte)
ggplot(dfmeantime,aes(prior,value))+
  geom_bar(stat="identity", aes(fill = algorithm),position="dodge")+
  #coord_flip()+
  ylab("mean time")+
  #xlab("method")+
  my_theme+
  scale_fill_manual(values=patte)

df2de <- read.csv("sim2_detail.csv",head=TRUE,row.names= 1)
df2de

library(forcats)
df3 <- read.csv("sim3_sum.csv",row.names= 1)
df3
df3plot <- data.frame(value=unlist(df3[1,seq(1,20,2)],use.names = FALSE),time=unlist(df3[1,seq(2,20,2)],use.names = FALSE))
df3plot$name <- c("max","mcmc,1000","mcmc,5000","mcmc,10000",
                  "mcmc,20000","shotgun,100","shotgun,200",
                  "shotgun,500","shotgun,1000","shotgun,2000")
df3plot$name <- fct_inorder(df3plot$name)
df3plot
df3plot$iter <- c("max","1000m","5000m","10000m","20000m","100s","200s","500s","1000s","2000s")
df3plot$iter <- factor(df3plot$iter,levels= c("max","1000m","5000m","10000m","20000m","100s","200s","500s","1000s","2000s"))
df3plot$algorithm <- c("EXH",rep("MCMC",4),rep("SSS",5))
ggplot(df3plot)+geom_bar(aes(iter,value,fill = algorithm),stat="identity")+
  labs(y="mean logABF",x="iteration")+
  #coord_flip()+
  my_theme+
  scale_fill_manual(values=patte)+
  scale_x_discrete(labels=c("max",1000,5000,10000,20000,100,200,500,1000,2000))

ggplot(df3plot)+geom_bar(aes(iter,time,fill = algorithm),stat="identity")+
  labs(y="time",x="iteration")+
  #coord_flip()+
  my_theme+
  scale_fill_manual(values=patte)+
  scale_x_discrete(labels=c("max",1000,5000,10000,20000,100,200,500,1000,2000))

df3de <- read.csv("sim3_detail.csv",head=TRUE,row.names=1)
mcmc10000 <- unlist(df3de[3,],use.names = FALSE)
plot(mcmc10000,type="l",xlab="Experiment No.",ylab="log ABF",main="mcmc,n. iter=10000",ylim=c(0,40))
abline(h=df3[1,1],col="blue")
shotgun200 <- unlist(df3de[6,],use.names = FALSE)
plot(shotgun200,type="l",xlab="Experiment No.",ylab="log ABF",main="shotgun, n.iter=200",ylim=c(0,40))
abline(h=df3[1,1],col="blue")
# hist(mcmc10000,xlab="log ABF",breaks=15,main="mcmc, n.iter=10000",xlim=c(0,25))
# abline(v=df3[1,1],lwd=2,col="blue")
# hist(shotgun500,xlab="log ABF",breaks=15,main="shotgun, n.iter=500")
# abline(v=df3[1,1],lwd=2,col="blue")
mcmc20000 <- unlist(df3de[4,],use.names = FALSE)
plot(mcmc20000,type="l",xlab="Experiment No.",ylab="log ABF",main="mcmc, n.iter=20000",ylim=c(0,40))
abline(h=df3[1,1],col="blue")
shotgun500 <- unlist(df3de[7,],use.names = FALSE)
plot(shotgun500,type="l",xlab="Experiment No.",ylab="log ABF",main="shotgun, n.iter=500",ylim=c(0,40))
abline(h=df3[1,1],col="blue")
# hist(mcmc20000,xlab="log ABF",breaks=15,main="mcmc, n.iter=20000",xlim=c(0,25))
# abline(v=df3[1,1],lwd=2,col="blue")
# hist(shotgun1000,xlab="log ABF",breaks=15,main="shotgun, n.iter=1000")
# abline(v=df3[1,1],lwd=2,col="blue")
shotgun1000 <- unlist(df3de[8,],use.names = FALSE)
plot(shotgun1000,type="l",xlab="Experiment No.",ylab="log ABF",main="shotgun, n.iter=1000",ylim=c(0,40))
abline(h=df3[1,1],col="blue")

df8 <- read.csv("sim3_sum.csv",row.names= 1)
df8
df8plot <- data.frame(value=unlist(df8[1,seq(1,22,2)],use.names = FALSE),time=unlist(df8[1,seq(2,22,2)],use.names = FALSE))
df8plot$name <- c("max","mcmc,1000","mcmc,5000","mcmc,10000",
                  "mcmc,20000","mcmc,50000","shotgun,100","shotgun,200",
                  "shotgun,500","shotgun,1000","shotgun,2000")
df8plot$name <- fct_inorder(df3plot$name)
df8plot
ggplot(df8plot)+geom_col(aes(name,value))+
  labs(y="mean logABF",x="method")+
  coord_flip()+
  my_theme

ggplot(df8plot)+geom_col(aes(name,time))+
  labs(y="time",x="method")+
  coord_flip()+
  my_theme

df1_1 <- read.csv("sim1-1.csv",row.names=1)
colnames(df1_1)[2]<-"FEM"
colnames(df1_1)[3]<-"REM"
colnames(df1_1)[4]<-"EXH"
colnames(df1_1)[5]<-"MCMC"
colnames(df1_1)[6]<-"SSS"
boost <- 1800

ggplot(df1_1)+
  geom_point(aes(trueOR,MCMC,colour="ABF(MCMC)"),size=1)+
  geom_line(aes(trueOR,MCMC,colour="ABF(MCMC)"),size=1)+
  geom_point(aes(trueOR,EXH,colour="ABF(EXH)"),size=1)+
  geom_line(aes(trueOR,EXH,colour="ABF(EXH)"),size=1)+
  geom_point(aes(trueOR,SSS,colour="ABF(SSS)"),size=1)+
  geom_line(aes(trueOR,SSS,colour="ABF(SSS)"),size=1)+
  geom_point(aes(trueOR,FEM*boost,colour="p-value(FEM)"),size=1)+
  geom_line(aes(trueOR,FEM*boost,colour="p-value(FEM)"),size=1)+
  geom_point(aes(trueOR,REM*boost,colour="p-value(REM)"),size=1)+
  geom_line(aes(trueOR,REM*boost,colour="p-value(REM)"),size=1)+
  labs(x="true OR",y="log ABF")+
  scale_y_continuous(sec.axis=sec_axis(~./boost,name="p-value"))+
  scale_colour_manual(name="class",values=colours)+
  my_theme+
  #geom_line(aes(trueOR,log(1000),colour="ABF=1000"),size=1.5)+
  scale_x_continuous(breaks=seq(0.6,2.0,0.2))+
  guides(shape=guide_legend(override.aes=list(size=10)))
  #scale_colour_manual(values=c("#CC6633","#333366","#669900","#666666","#FFCC33"))+
  #scale_linetype_manual(values = c("ABF(EXH)"='dotted')) 
