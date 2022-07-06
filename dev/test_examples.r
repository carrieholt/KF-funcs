#==============================================================
#Functions to run TMB stock recruitment functions
#Author: Catarina Wor
#Date: 14th of June 2018
#==============================================================


library(ggplot2)


#============================================
devtools::document()
devtools::load_all()
#update online documentation
#pkgdown::build_site()

#equivalent to ctrl + b in Rstudio 
#devtools::build(binary= TRUE)

#devtools::build()


#detach("KFfuncs", unload=TRUE)
#devtools::unload("KFfuncs")

#path.install <-  "C:/Users/worc/Documents/KFfuncs_0.0.0.1000.tar.gz"
#path.install <- "/Users/catarinawor/Documents/timevar/KFfuncs_0.0.0.1000.tar.gz"

system(paste0("Rcmd.exe INSTALL --preclean --no-multiarch --with-keep.source ", path.install))


library(KFfuncs)




#===============================================================
#devtools::load_all("../KF-funcs")
#devtools::document("../KF-funcs")
#devtools::install_local("../KF-funcs")


#===============================================================
x <-Stellako$ETS
y <-log(Stellako$Rec/Stellako$ETS)


initial <- list()
initial$mean.a <- lm(y~x)$coefficients[1]
initial$var.a <- 1
initial$b <- -lm(y~x)$coefficients[2]
initial$ln.sig.e <- log(1)
initial$ln.sig.w <- log(1)
initial$Ts <- 0
initial$EstB <- TRUE

Stel<-kf.rw(initial=initial,x=x,y=y)
names(Stel)
Stel$smoothe.y


steltmbdat<-data.frame(S=Stellako$ETS,
                  R=Stellako$Rec)

Steltmb <- kfTMB(steltmbdat)


Steltmb$tmb_obj$report()
summary(sdreport(Steltmb$tmb_obj))


Steltmbmc <- kfTMBmcmc(data=steltmbdat,iter=3000 )
names(Steltmbmc)

summary(Steltmbmc$fitmcmc)

#use/test package


mydata<-list(S=Stellako$ETS,
  R=Stellako$Rec
  )

resu <- rbTMB(data=mydata, priorratiovar=c(1,1), silent = FALSE, control = TMBcontrol())

resu$sd_report
resu$tmb_obj$report()
sdrep<-summary(sdreport(resu$tmb_obj))
lowatmb<-sdrep[which(rownames(sdrep)=="alpha"),1]-1.96*sdrep[which(rownames(sdrep)=="alpha"),2]
higatmb<-sdrep[which(rownames(sdrep)=="alpha"),1]+1.96*sdrep[which(rownames(sdrep)=="alpha"),2]

rekf <- kfTMB(data=mydata, silent = FALSE, control = TMBcontrol())
names(rekf)

cbind(sqrt(rekf$tmb_obj$report()$smoothevara),
  kfrep[which(rownames(kfrep)=="smoothemeana"),2])

kfrep <- summary(sdreport(rekf$tmb_obj))
kfrep[which(rownames(kfrep)=="smoothemeana"),1]
kfrep[which(rownames(kfrep)=="smoothevara"),1]

lowkftmb<-kfrep[which(rownames(kfrep)=="smoothemeana"),1]-1.96*sqrt(kfrep[which(rownames(kfrep)=="smoothevara"),1])
higkftmb<-kfrep[which(rownames(kfrep)=="smoothemeana"),1]+1.96*sqrt(kfrep[which(rownames(kfrep)=="smoothevara"),1])

lowkftmbapprox<-kfrep[which(rownames(kfrep)=="smoothemeana"),1]-1.96*kfrep[which(rownames(kfrep)=="smoothemeana"),2]
higkftmbapprox<-kfrep[which(rownames(kfrep)=="smoothemeana"),1]+1.96*kfrep[which(rownames(kfrep)=="smoothemeana"),2]



mdf<-data.frame(est_a=c(Stel$smoothe.mean.a,
  sdrep[which(rownames(sdrep)=="alpha"),1],
  kfrep[which(rownames(kfrep)=="smoothemeana"),1]),
  type=rep(c("Kalman-filter","recursive Bayes","Kalman-filter TMB"),each=length(x)),
    ind=rep(1:length(x),3))

mdf$lower<-c(Stel$smoothe.mean.a-1.96*sqrt(Stel$smoothe.var.a),lowatmb,lowkftmb)
mdf$upper<-c(Stel$smoothe.mean.a+1.96*sqrt(Stel$smoothe.var.a),higatmb,higkftmb)

Stel$smoothe.var.a


p <- ggplot(mdf)
p <- p + geom_line(aes(x=ind,y=est_a, col=type),size=1.4)
p <- p + geom_ribbon(aes(x=ind,ymin=lower, ymax=upper, fill=type),alpha=0.2)
p <- p + theme_bw(18)
p

#===============================================
#compare analytical se and the TMB sd approximation 

kfdf<-data.frame(est_a=c(Stel$smoothe.mean.a,
  kfrep[which(rownames(kfrep)=="smoothemeana"),1],
  kfrep[which(rownames(kfrep)=="smoothemeana"),1]),
  type=rep(c("Kalman-filter","Kalman-filter TMB", "Kalman-filter TMB approx"),each=length(x)),
    ind=rep(1:length(x),3))

kfdf$lower<-c(Stel$smoothe.mean.a-1.96*sqrt(Stel$smoothe.var.a),lowkftmb, lowkftmbapprox)
kfdf$upper<-c(Stel$smoothe.mean.a+1.96*sqrt(Stel$smoothe.var.a),higkftmb,higkftmbapprox)


pp <- ggplot(kfdf)
pp <- pp + geom_line(aes(x=ind,y=est_a, col=type),size=1.4)
pp <- pp + geom_ribbon(aes(x=ind,ymin=lower, ymax=upper, fill=type),alpha=0.2)
pp <- pp + theme_bw(18)
pp






library(TMB)
library(KFfuncs)
compile("../src/Rickerkf.cpp", '-O1 -g', DLLFLAGS='', framework = 'TMBad')
dyn.load(dynlib("../src/Rickerkf"))


mydata<-list(S=Stellako$ETS,
  R=Stellako$Rec
  )
tmb_data <- list(
    x = mydata$S,
    y = log(mydata$R/mydata$S)
  )

 tmb_params <- list(
    initmeana   = lm(y~x, data=tmb_data)$coefficients[[1]],
    #loginitvara = log(1),
    b           = lm(y~x, data=tmb_data)$coefficients[[2]],
    logsige     = log(1),
    logsigw     = log(1)
  )

obj <- MakeADFun(tmb_data,tmb_params,DLL="Rickerkf")#,lower = -Inf, upper = Inf)
newtonOption(obj, smartsearch=FALSE)

  opt <- nlminb(obj$par,obj$fn,obj$gr)
  obj$rep()
kfrep <- summary(sdreport(obj))
   
