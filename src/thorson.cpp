// Pella-Tomlinson surplus production model reparamterised according to Fletcher
// ala Thorson et al. (2012) "Spawning biomass reference points for exploited
// marine fishes, incorporating taxonomic and body size information"

#include <TMB.hpp>

templat<class Type>
Type objective_function<Type>::operator() ()
{
  // data
  DATA_VECTOR(Cobs); // catch
  DATA_VECTOR(SBobs); // biomass

  // parameters
  PARAMETER(logPhi);
  PARAMETER(logSB0);
  // PARAMETER(logSB1SB0);
  PARAMETER(logY);
  PARAMETER_VECTOR(logSB);
  PARAMETER(logSigma_p);
  PARAMETER(logSigma_m);

  // quantities
  Type sigma_p = exp(logSigma_p);
  Type sigma_m = exp(logSigma_m);
  Type phi = exp(logPhi);
  Type gamma = pow(phi, phi/(phi-1.0)) / (phi - 1.0);
  Type SB0 = exp(logSB0);
  //  Type SB1SB0 = exp(logSB1SB0);
  Type y = exp(logY);
  vector<Type> SB = exp(logSB);
  //Type y = MSY / SB0;
  int nt = SB.size();

  // surplus production
  //  SB(0) = SB1SB0 * SB0;
  //  logSB(0) = log(SB(0));
  vector<Type> S(nt);
  for(int t=0; t<nt; t++){
        // if(phi != 1.0){
      // Pella and Tomlinson model
      S(t) = gamma * y * SB(t) - gamma * y * SB0 * pow(SB(t)/SB0, phi);
    // }else{
    //   // Fox model if phi == 1
    //   S(t) = exp(1) * y * SB(t) * log(SB(t)/SB0);
    // }
    std::cout << "S(" << t << "): " << S(t) << std::endl;
  }

  // initialise nll
  Type nll = 0;

  // process error
  vector<Type> SBpred(nt);
  for(int t=0; t<(nt-1); t++){
    SBpred(t+1) = SB(t) + S(t) - Cobs(t);
    // process error
    if(!isNA(SB(t+1))) nll -= dnorm(log(SBpred(t+1)), logSB(t+1),
                                        sigma_p, true);
    std::cout << "SBpred(" << t << "): " << SBpred(t) << std::endl;
  }

  // observation error (SBobs(t) = SB(t) * exp(tau(t)))
  for(int t=1; t<nt; t++){
    if(!isNA(SBobs(t))) nll -= dnorm(log(SBobs(t)), logSB(t), sigma_m, true);
  }

  // ADreports
  ADREPORT(sigma_p);
  ADREPORT(sigma_m);
  ADREPORT(S);

  // return nll
  return nll;
}
