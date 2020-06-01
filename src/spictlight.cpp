
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll = 0;

  DATA_SCALAR(dt);
  DATA_VECTOR(logobsC);
  DATA_VECTOR(logobsI);
  DATA_VECTOR(ii);
  DATA_VECTOR(ic);
  DATA_VECTOR(nc);
  DATA_VECTOR(priorlogn);

  PARAMETER(logsdb);
  PARAMETER(logsdi);
  PARAMETER(logsdf);
  PARAMETER(logsdc);
  PARAMETER(logn);
  PARAMETER(logm);
  PARAMETER(logK);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logF);
  PARAMETER(logq);

  int ns = logF.size();
  int nobsI = logobsI.size();
  int nobsC = logobsC.size();
  Type m = exp(logm);
  Type n = exp(logn);
  Type gamma = pow(n, n/(n-1)) / (n - 1);
  Type K = exp(logK);
  Type sdf = exp(logsdf);
  Type sdb = exp(logsdb);
  Type sdb2 = pow(sdb,2);
  Type sdc = exp(logsdc);
  Type sdi = exp(logsdi);
  vector<Type> F = exp(logF);
  vector<Type> B = exp(logB);

  // Prior on logn
  if(priorlogn(0) != 1){
    nll -= dnorm(logn, priorlogn(1), priorlogn(2), 1);
  }

  // Fishing mortality as RW
  for(int i=1; i<ns; i++){
    nll -= dnorm(logF(i), logF(i-1), sqrt(dt) * sdf, true);
  }

  // Biomass
  vector<Type> logBpred(ns);
  for(int i=0; i<(ns-1); i++){
    logBpred(i+1) = log(B(i)) + (gamma*m/K - gamma*m/K*pow(B(i)/K, n-1.0) - F(i) - 0.5*sdb2)*dt;
    nll -= dnorm(logBpred(i+1), logB(i+1), sqrt(dt) * sdb, true);
  }

  // Observations
  // catches
  int ind;
  vector<Type> Cpredsub(ns);
  for(int i=0; i<ns; i++){
    Cpredsub(i) =  F(i) * B(i) * dt;
  }
  vector<Type> logCpred(nobsC);
  vector<Type> Cpred(nobsC);
  for(int i=0; i<nobsC; i++){
    Cpred(i) = 0.0;
  }
  for(int i=0; i<nobsC; i++){
    for(int j=0; j<nc(i); j++){
      ind = CppAD::Integer(ic(i)-1) + j;
      Cpred(i) += Cpredsub(ind);
    }
    logCpred(i) = log(Cpred(i));
  }
  for(int i=0; i<nobsC; i++){
    nll -= dnorm(logCpred(i), logobsC(i), sdc, true);
  }

  // index
  vector<Type> logIpred(nobsI);
  for(int i=0; i<nobsI; i++){
    ind = CppAD::Integer(ii(i)-1);
    logIpred(i) = logq + logB(ind);
    nll -= dnorm(logobsI(i), logIpred(i), sdi, true);
  }

  // reporting
  ADREPORT(Cpredsub);

  return nll;
}
