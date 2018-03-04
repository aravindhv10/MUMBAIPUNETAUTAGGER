// Example for a low_pt working point
#ifndef __LOWPT_HH__
#define __LOWPT_HH__

#include <math.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"

#include "./HEPTopTagger.hh"
#include "./FWM.hh"

class LowPt {
 public:
  LowPt();
  
  bool is_tagged(HEPTopTagger::HEPTopTagger htt);
};
#endif
