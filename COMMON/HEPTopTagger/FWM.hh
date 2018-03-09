#ifndef __FWM_HH__
#define __FWM_HH__

#include <math.h>
#include <cmath>

#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_legendre.h"

#include "./HEPTopTagger.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Boost.hh"

using namespace std;
using namespace fastjet;  

class FWM {
public:
  FWM();
  FWM(vector<fastjet::PseudoJet> jets);
  FWM(HEPTopTagger::HEPTopTagger htt, int selection);
  
  double U(unsigned order);
  double Pt(unsigned order);
  double Pt(unsigned order, fastjet::PseudoJet ref_pj);
 
private:
  double cos_Omega(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2);
  inline double ATan2(double x, double y);
  double Theta(fastjet::PseudoJet j);
  double legendre(int l, double x) {return gsl_sf_legendre_Pl(l,x);};

  vector<fastjet::PseudoJet> _jets;
  };

//--------------------------------------------------------------------
#endif // __FWM_HH__
