#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //Kalman filter translated to TMB from Carrie holt R-code
  DATA_VECTOR(x);   // observed log recruitment
  DATA_VECTOR(y);    // observed  Spawner
  
 
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value
  //rho         -> Proportion of total variance associated with obs error.
  //varphi      -> Total precision
  //alpha       -> Time-varying alpha

  PARAMETER(initmeana);
  PARAMETER(loginitvara);
  PARAMETER(b);
  PARAMETER(logsige);
  PARAMETER(logsigw);

  Type sige = exp(logsige);
  Type sigw = exp(logsigw);
  Type initvara = exp(loginitvara);

  int Tmax = x.size();

  vector<Type> priormeana(Tmax),  priorvara(Tmax), yhat(Tmax),
  f(Tmax), v(Tmax), postmeana(Tmax), postvara(Tmax), filtery(Tmax),
   smoothemeana(Tmax), smoothevara(Tmax), smoothey(Tmax);


  vector<Type> pstar(Tmax-1);
  
  Type ans = Type(0.0);

  
  
  //first year
  priormeana(0) = initmeana;
  priorvara(0) = initvara;
  yhat(0) = priormeana(0) + b * x(0);
  v(0) = y(0) - yhat(0);
  f(0) = priorvara(0) + sige*sige;
     
  postmeana(0) = priormeana(0) + (priorvara(0) * (v(0)/f(0)));
  postvara(0) = priorvara(0) - (priorvara(0)*priorvara(0)/f(0));
  filtery(0) = postmeana(0) + b * x(0);
  ans +=  (log(f(0)) + (v(0)*v(0)/f(0)))/Type(2.);



  for(int t=1; t<Tmax; t++){
     
    priormeana(t) = postmeana(t-1);
    priorvara(t) = postvara(t-1) + sigw*sigw;
    
    //Step 2: Generate predicted value for y(t) given y(t-1) and error
    yhat(t) = priormeana(t) + b * x(t);
    v(t) = y(t) - yhat(t);
    f(t) = priorvara(t) + sige*sige;
    
    // Step 3: Generate posterior distribution for intercept (a):
    postmeana(t) = priormeana(t) + (priorvara(t) * (v(t)/f(t)));
    postvara(t) = priorvara(t) - (priorvara(t)*priorvara(t)/f(t));
    filtery(t) = postmeana(t) + b * x(t);
    ans += (log(f(t)) + (v(t)*v(t)/f(t)))/Type(2.);

  }



  //smoothing part

  
  smoothemeana(Tmax-1) = postmeana(Tmax-1);
  smoothevara(Tmax-1) = postvara(Tmax-1);
  smoothey(Tmax-1) = smoothemeana(Tmax-1) + b * x(Tmax-1);
  

  for(int i=Tmax-2; i>=0; --i){
      
      pstar(i) = postvara(i)/priorvara(i + 1);
      smoothemeana(i) = postmeana(i) + pstar(i) * (smoothemeana(i + 1) - priormeana(i + 1));
      smoothevara(i) = postvara(i) + pstar(i)*pstar(i) * (smoothevara(i + 1) - priorvara(i + 1));    
      smoothey(i) = smoothemeana(i) + b * x(i);
  }

  
   REPORT(ans)
   REPORT(Tmax)
   REPORT(priormeana)
   REPORT(priorvara)
   REPORT(yhat)
   REPORT(f)
   REPORT(v)
   REPORT(postmeana)
   REPORT(postvara)
   REPORT(filtery)
   REPORT(pstar)
   ADREPORT(smoothemeana)
   REPORT(smoothevara)
   REPORT(smoothey)
   REPORT(initmeana)
   REPORT(initvara)
   REPORT(b)
   REPORT(sige)
   REPORT(sigw)
    return ans;

}

