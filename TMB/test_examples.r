#==============================================================
#Functions to run TMB stock recruitment functions
#Author: Catarina Wor
#Date: 14th of June 2018
#==============================================================


library(TMB)
#library(ggplot)

library(tmbstan)
devtools::load_all()

#===============================================================
#Stellako
x <-Stellako$ETS
y <-log(Stellako$Rec/Stellako$ETS)




mydata<-list(obs_logR=log(Stellako$Rec),
  obs_S=Stellako$ETS,
  prbeta1=1.0,
  prbeta2=1.0
  )

compile("Rickerkf_ratiovar.cpp", "-O1 -g", DLLFLAGS="")

dyn.load("Rickerkf_ratiovar")


parameters <- list(
  alphao=lm(y~x)$coefficients[1],
  logbeta = log(-lm(y~x)$coefficients[2]),
  rho=.2,
  logvarphi= 0.1,
  alpha=rep(0.9,length(mydata$obs_logR))  )

obj<-MakeADFun(mydata,parameters,random="alpha",DLL="Rickerkf_ratiovar")
newtonOption(obj, smartsearch=FALSE)

opt<-nlminb(obj$par,obj$fn,obj$gr)
rep<-obj$report()


#MCMC
fitmcmc1 <- tmbstan(obj, chains=3,
              iter=1000000, init="random",
              lower=c(0.1,-1.0,0.0,-3.0),
              upper=c(5.0,8.,1.0,5.0),
               control = list(adapt_delta = 0.98))

    mc <- extract(fitmcmc1, pars=names(obj$par),
              inc_warmup=TRUE, permuted=FALSE)
    


    fit_summary <- summary(fitmcmc1)

    



#=====================================
#Kalman filter

initial <- list()
initial$mean.a <- lm(y~x)$coefficients[1]
initial$var.a <- 1
initial$b <- -lm(y~x)$coefficients[2]
initial$ln.sig.e <- log(1)
initial$ln.sig.w <- log(1)
initial$Ts <- 1
initial$EstB <- "True"

Stel<-kf.rw(initial=initial,x=x,y=y)

mdf<-data.frame(est_a=c(Stel$smoothe.mean.a,rep$alpha,fit_summary$summary[5:65,6]),
  type=rep(c("Kalman-filter","Recursive Bayes MLE","Recursive Bayes"),each=length(x)),
    ind=rep(1:length(x),3))

mdf$lower<-c(rep(NA,length(x)*2),fit_summary$summary[5:65,4])
mdf$upper<-c(rep(NA,length(x)*2),fit_summary$summary[5:65,8])



p <- ggplot2::ggplot(mdf)
p <- p + ggplot2::geom_line(aes(x=ind,y=est_a, col=type),size=1.4)
p <- p + ggplot2::geom_ribbon(aes(x=ind,ymin=lower, ymax=upper, fill=type),alpha=0.2)
p <- p + theme_bw(18)
p



#===========================================
#Harrison

har<-read.csv("../data/Harrison_Apr18.csv")


hx <-har$S
hy <-log(har$R/har$S)


mydatah<-list(obs_logR=log(har$R),
  obs_S=har$S,
  prbeta1=1.0,
  prbeta2=1.0
  )

parametersh <- list(
  alphao=lm(hy~hx)$coefficients[1],
  logbeta = log(-lm(hy~hx)$coefficients[2]),
  rho=.2,
  logvarphi= 0.1,
  alpha=rep(0.9,length(mydata$obs_logR))  )


objh<-MakeADFun(mydatah,parametersh,random="alpha",DLL="Rickerkf_ratiovar")
newtonOption(objh, smartsearch=FALSE)

opth<-nlminb(objh$par,objh$fn,objh$gr)
reph<-objh$report()



#=====================================
#Kalman filter

initialh <- list()
initialh$mean.a <- lm(hy~hx)$coefficients[1]
initialh$var.a <- 1
initialh$b <- -lm(hy~hx)$coefficients[2]
initialh$ln.sig.e <- log(1)
initialh$ln.sig.w <- log(1)
initialh$Ts <- 1
initialh$EstB <- "True"

HR<-kf.rw(initial=initialh,x=hx,y=hy)






#===================================================
#old functions these won't work with this example
runTMB<-function(A,comps=TRUE){

	####
	#This function compiles and runs a standard TMB model
	###
  setwd(A$DIR)

  cppfile<-paste(A$dll,".cpp",sep="")

  if(comps==TRUE){
    compile(cppfile,libtmb=FALSE, "-O1 -g", DLLFLAGS="")
    dyn.load(dynlib(A$dll))
  }
  
  obj<-MakeADFun(A$dat,A$params,random=A$rndm,DLL=A$dll)
  newtonOption(obj, smartsearch=FALSE)
  
  opt<-nlminb(obj$par,obj$fn,obj$gr)
  rep<-obj$report()

  return(obj)
}




posteriorsdf<-function(B){

	####
	#This function runs TMB stan and produces clean posteriors in data-frame format
	###

    fitmcmc1 <- tmbstan(B$obj, chains=B$nchain,
              iter=B$iter, init="random",
              lower=B$lowbd, upper=B$hibd,
               control = list(adapt_delta = 0.98))

    mc <- extract(fitmcmc1, pars=names(B$obj$par),
              inc_warmup=TRUE, permuted=FALSE)
    

    fit_summary <- summary(fitmcmc1)

    posterior <- as.array(fitmcmc1)

    mainrun <- melt(posterior)

    poslist <- list(fit=fitmcmc1,
        fit_summary=fit_summary,
        posteriors=mainrun,
        mcmcobj=mc
      ) 

    return(poslist) 
}



plot_posteriors<-function(df,salvar=FALSE,DIR="",nome=""){

	####
	#This function plots posterior distributions by chain
	###

  pm<-ggplot(df)
  pm<-pm+geom_density(aes(x=value, color=chains))
  pm<-pm+facet_wrap(~parameters, scales="free")
  print(pm)

  if(salvar){
    setwd(DIR)
    ggsave(nome, plot=pm, width=10,height=7)

  }

}



results_table<-function(D){

	####
	#This function produces tex tables with SR function results
	####
  if(sum(!is.na(D$MCMC))){
    tab<-data.frame(Parameter=D$param_names,
                      MLE=D$MLE,
                      Median=c(apply(D$MCMC,2,function(x) quantile(x, .5)),D$other[2,]),
                      Lower=c(apply(D$MCMC,2,function(x) quantile(x, .025)),D$other[1,]),
                      Upper=c(apply(D$MCMC,2,function(x) quantile(x, .975)),D$other[3,]))
  }else{
    tab<-data.frame(Parameter=D$param_names,
                      MLE=D$MLE,
                      Median=c(D$other[2,]),
                      Lower=c(D$other[1,]),
                      Upper=c(D$other[3,]))
  }

    setwd(D$DIR)
    tabtmp<-xtable(tab, digits=D$digits,caption = D$caption,label=D$labs)
    digits(tabtmp)<-D$digits

    print(tabtmp,sanitize.text.function = function(x) {x},
      include.rownames = FALSE, 
  file=D$filename,caption.placement = "top")

}



model_pred_plot<-function(M, salvar=FALSE,DIR="",filename=""){

    sumfit<-apply(M$predBayes,2,function(x) quantile(x, probs=c(0.025,.5,0.975)))
    fitdf<-as.data.frame(t(sumfit))
    names(fitdf)<-c("lower","estimate","upper")
    fitdf<-cbind(fitdf,M$orig_data)
    fitdf$type<-"Bayesian"


    fitdf1<-fitdf
    fitdf1$estimate<-M$predFreq
    fitdf1$lower<-NA
    fitdf1$upper<-NA
    fitdf1$type<-"MLE"

  fitdf<-fitdf[order(SR$S_adj),]
  fitdf1<-fitdf1[order(SR$S_adj),]

  fitdf<-rbind(fitdf1,fitdf)

    p<-ggplot(fitdf)
    p<-p+geom_line(aes(x=S_adj,y=estimate, col=type),size=1.5)
    p<-p+geom_ribbon(aes(x=S_adj,ymin=lower,ymax=upper, fill=type),alpha=0.4)
    p <- p + geom_text(aes(x=S_adj,y=R,label=BroodYear ),hjust=0, vjust=0)
    p <- p + theme_bw(16)
    p <- p + ylab("Recruits") + xlab("Spawners")
    print(p)

    if(salvar){
      setwd(DIR)
      ggsave(filename, plot=p, width=10,height=7)
    }

}







