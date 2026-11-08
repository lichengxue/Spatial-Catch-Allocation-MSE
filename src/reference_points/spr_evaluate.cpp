#include <TMB.hpp>
#include "legacy_wham/all.hpp"
template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(SPR_weights); DATA_ARRAY(log_M); DATA_ARRAY(FAA);
  DATA_IVECTOR(spawn_seasons); DATA_IVECTOR(spawn_regions); DATA_IVECTOR(fleet_regions);
  DATA_IMATRIX(fleet_seasons); DATA_VECTOR(fracyr_seasons);
  DATA_IARRAY(can_move); DATA_IARRAY(must_move); DATA_IVECTOR(mig_type);
  DATA_ARRAY(trans_mu_base); DATA_MATRIX(L); DATA_INTEGER(which_F_age);
  DATA_ARRAY(waa_ssb); DATA_ARRAY(waa_catch); DATA_ARRAY(mature);
  DATA_SCALAR(percentSPR); DATA_ARRAY(NAA); DATA_MATRIX(fracyr_SSB);
  DATA_IVECTOR(years); DATA_IVECTOR(Fbar_ages); DATA_VECTOR(R_XSPR);
  DATA_INTEGER(small_dim); DATA_INTEGER(SPR_weight_type); DATA_INTEGER(bias_correct);
  DATA_ARRAY(marg_NAA_sigma);
  PARAMETER(logF);
  // With one iteration the original WHAM routine evaluates SPR at the supplied
  // F, without taking any Newton steps. R uses a bracketed root search instead.
  vector< array<Type> > res = get_SPR_res(SPR_weights,log_M,FAA,spawn_seasons,spawn_regions,
    fleet_regions,fleet_seasons,fracyr_seasons,can_move,must_move,mig_type,trans_mu_base,L,
    which_F_age,waa_ssb,waa_catch,mature,percentSPR,NAA,fracyr_SSB,exp(logF),
    years,years,years,years,years,years,years,Fbar_ages,R_XSPR,small_dim,
    SPR_weight_type,bias_correct,marg_NAA_sigma,0,1);
  int ns=spawn_regions.size();
  Type ratio=exp(res(3)(0,ns)-res(4)(0,ns));
  vector<Type> Fbar=exp(vector<Type>(res(17).matrix().row(0)));
  vector<Type> SSBref=exp(vector<Type>(res(1).matrix().row(0)));
  REPORT(ratio); REPORT(Fbar); REPORT(SSBref);
  return ratio-Type(.4);
}
