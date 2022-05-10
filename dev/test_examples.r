#==============================================================
#Functions to run TMB stock recruitment functions
#Author: Catarina Wor
#Date: 14th of June 2018
#==============================================================


library(ggplot2)
devtools::document()

#update online documentationH
#pkgdown::build_site()

#equivalent to ctrl + b in Rstudio 
#devtools::build(binary= TRUE)

devtools::build()


detach("package:KF-funcs", unload=TRUE)
path.install <-  "C:/Users/worc/Documents/KFfuncs_0.0.0.1000.tar.gz"

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
initial$EstB <- "TRUE"

Stel<-kf.rw(initial=initial,x=x,y=y)
names(Stel)




#use/test package


mydata<-list(S=Stellako$ETS,
  R=Stellako$Rec
  )


    obs_logRS = log(Stellako$Rec/Stellako$ETS),

summary(Stellako)


resu <- rbTMB(data=mydata, priorratiovar=c(1,1), silent = FALSE, control = TMBcontrol())

names(resu)
resu$sd_report
resu$tmb_obj$report()
sdrep<-summary(sdreport(resu$tmb_obj))
lowatmb<-sdrep[which(rownames(sdrep)=="alpha"),1]-1.96*sdrep[which(rownames(sdrep)=="alpha"),2]
higatmb<-sdrep[which(rownames(sdrep)=="alpha"),1]+1.96*sdrep[which(rownames(sdrep)=="alpha"),2]


mdf<-data.frame(est_a=c(Stel$smoothe.mean.a,sdrep[which(rownames(sdrep)=="alpha"),1]),
  type=rep(c("Kalman-filter","recursive Bayes"),each=length(x)),
    ind=rep(1:length(x),2))

mdf$lower<-c(Stel$smoothe.mean.a-1.96*Stel$"sig.w",lowatmb)
mdf$upper<-c(Stel$smoothe.mean.a+1.96*Stel$"sig.w",higatmb)



p <- ggplot(mdf)
p <- p + geom_line(aes(x=ind,y=est_a, col=type),size=1.4)
p <- p + geom_ribbon(aes(x=ind,ymin=lower, ymax=upper, fill=type),alpha=0.2)
p <- p + theme_bw(18)
p










