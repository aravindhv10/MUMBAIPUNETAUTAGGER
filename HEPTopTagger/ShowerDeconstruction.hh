#include "AnalysisParameters.h"

const Parameter listOfParameters[] = {
    CREATE_PARAMETER(R),
    CREATE_PARAMETER(ptmin_jet),
    CREATE_PARAMETER(n_of_jets),
    CREATE_PARAMETER(rapmax_jet),
    CREATE_PARAMETER(eta_cell),
    CREATE_PARAMETER(phi_cell),
    CREATE_PARAMETER(cell_ptcut),
    CREATE_PARAMETER(ptcut_microjets),
    CREATE_PARAMETER(masscut_microjets),
    CREATE_PARAMETER(lambda_mu_ext),
    CREATE_PARAMETER(lambda_theta_ext),
    CREATE_PARAMETER(cut_n_sub),
    CREATE_PARAMETER(fractionISR_of_hardInt_PT),
    CREATE_PARAMETER(noISR),
    CREATE_PARAMETER(topmass),
    CREATE_PARAMETER(topwindow),
    CREATE_PARAMETER(wmass),
    CREATE_PARAMETER(wwindow),
    CREATE_PARAMETER(Higgsmass),
    CREATE_PARAMETER(delta_Higgsmass),
    CREATE_PARAMETER(Rsmall),
    CREATE_PARAMETER(m_min_jet),
    CREATE_PARAMETER(br_wqq),
    CREATE_PARAMETER(useBtag),
    CREATE_PARAMETER(tagprob),
    CREATE_PARAMETER(fakeprob)
  };

AnalysisParameters::AnalysisParameters()
  : Deconstruction::Parameters() {
  init();
}

AnalysisParameters::AnalysisParameters(const std::string &input)
  : Deconstruction::Parameters() {
  init();
  read(input);
}

void AnalysisParameters::init() {
  for (int k = 0; k < sizeof(listOfParameters)/sizeof(Parameter); ++k) {
    insert(listOfParameters[k].m_key, 0);
  }
}

#include "BackgroundModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "Helper.h"

#include "Parameters.h"

using namespace Deconstruction;

Deconstruction::BackgroundModel::BackgroundModel(Parameters &param, Flavour::id flav)
  : Deconstruction::Model::Model(param), m_flav(flav) {
}

Deconstruction::BackgroundModel::~BackgroundModel() {
}

//double Deconstruction::BackgroundModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sum) {
double Deconstruction::BackgroundModel::weight(const StoredJet &jets, fastjet::PseudoJet &sum) {
  // Background model ---> calls make_splitting with gluon

  // the Background model consists simply of a FSRshower
  // keep in mind that the H for the generation of the 
  // fat jet is applied in Deconstruct::deconstruct()

  // check for theta cut if there is no grandmother:
  if (sum.perp() <= sum.m())
    return 0;

  double w = start_splitting(jets, m_flav, Shower::QCD);
  m_calc.clear();

  return w;
}

double Deconstruction::BackgroundModel::hamiltonian(double pTsum) {
  double kH2 = square(pTsum);
  double Hfj_sig = Cte::Npdf_signal*std::pow(Cte::pTmin2/kH2, Cte::Npdf_signal)/kH2;
  return Hfj_sig;
}

#include "BModel.h"

#include <iostream>
#include <vector>
#include <fastjet/PseudoJet.hh>
#include "Helper.h"

#include "Parameters.h"

using namespace Deconstruction;
using namespace std;

Deconstruction::BModel::BModel(Parameters &param, Flavour::id flavourb, Flavour::id flavourt)
  : Deconstruction::Model::Model(param), m_flavourb(flavourb), m_flavourt(flavourt) {
}

Deconstruction::BModel::~BModel() {
}

//double Deconstruction::BModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sum) {
double Deconstruction::BModel::weight(const StoredJet &jets, fastjet::PseudoJet &sum) {
  // top is color connected partner von der ersten emission vom b
  // evtl sollte top bzw. bottom flavor mitgegeben werden - sind ja korrelliert.

  StoredJet topset(jets.store(), m_topset);

#ifdef DEBUGCODE
  printoutput(jets,"input in btopmodel", DEBUG);
  LOG(DEBUG) << "input flavor: " << m_flavourt << endl;
  LOG(DEBUG) << "top: " << endl;
  printoutput(topset, DEBUG);
#endif

  // check for theta cut if there is no grandmother:

  double bglusplitweight = 0;

  StoredJet empty(jets.store());

  // weight is the same for qquark and antiqquark...
  if ( (m_flavourt == Flavour::t) && (m_flavourb == Flavour::b) ) {
    bglusplitweight = make_splitting(jets,empty,topset,topset,Flavour::b,Flavour::t,Flavour::noflav,Flavour::t,Shower::b);  // b radiates right (top right col partner)
  } else if ( (m_flavourt == Flavour::t) && (m_flavourb == Flavour::bbar) ) {
    bglusplitweight = make_splitting(jets,topset,empty,topset,Flavour::b,Flavour::t,Flavour::t,Flavour::noflav,Shower::b); // antib radiates left (top left col.partner)
  } else
    LOG(ERROR) << "ERROR in btopmodel: either tflavor or bflavor wrong...";

#ifdef DEBUGCODE
  LOG(DEBUG) << "bglusplitweight: " << bglusplitweight << endl;
#endif

  return bglusplitweight;
}

//void Deconstruction::BModel::setTopset(const std::vector<fastjet::PseudoJet> &topset) {
void Deconstruction::BModel::setTopset(const StoredJet &topset) {
  m_topset = topset.getList();
}

double Deconstruction::BModel::hamiltonian(double pTsum) {
  return 1.0;
}

#include "Deconstruct.h"
#include "Parameters.h"
#include "AnalysisParameters.h"
#include "Model.h"

#include "Exception.h"

#include "JetInfo.h"

#include "Helper.h"

#include <algorithm>

#include <fastjet/PseudoJet.hh>
#include "StoredJet.h"

#include <map>

using namespace Deconstruction;
using namespace std;
using namespace fastjet;

Deconstruction::Deconstruct::Deconstruct(Parameters &param, Model &signal, Model &background, Model &isr)
  : m_param(param), m_signal(signal), m_background(background), m_isr(isr) {
  std::cout << std::endl;
  std::cout << "Shower Deconstruction" << std::endl;
  std::cout << "=======================" << std::endl;
  std::cout << "If you use this code, please cite the following papers:" << std::endl;
  std::cout << "D. E. Soper and M. Spannowsky, ``Finding physics signals with shower deconstruction,''  Phys. Rev. D 84 (2011) 074002 [arXiv:1102.3480 [hep-ph]]." << std::endl;
  std::cout << "D. E. Soper and M. Spannowsky, ``Finding top quarks with shower deconstruction,'' Phys. Rev. D 87 (2013) 5,  054012 [arXiv:1211.3140 [hep-ph]]." << std::endl;
  std::cout << "In case of questions, please contact Danilo.Enoque.Ferreira.De.Lima@cern.ch or Michael.Spannowsky@cern.ch" << std::endl;
  std::cout << std::endl;
}

Deconstruction::Deconstruct::~Deconstruct() {
}

Model &Deconstruction::Deconstruct::signal() {
  return m_signal;
}

Model &Deconstruction::Deconstruct::background() {
  return m_background;
}

Model &Deconstruction::Deconstruct::isr() {
  return m_isr;
}

Parameters &Deconstruction::Deconstruct::param() {
  return m_param;
}

double Deconstruction::Deconstruct::deconstruct(std::vector<fastjet::PseudoJet> &input,
                              double &wSignal, double &wBackground) {

  m_signalWeight.clear();

  Model::initlevel = -1;// to agree with dave's numbering scheme
  Model::binitlevel = -1;
  Model::tinitlevel = -1;
  Model::Winitlevel = -1;

  Model::level=Model::initlevel;
  Model::btopshower_level=Model::binitlevel;
  Model::tshower_level=Model::tinitlevel;
  Model::Wshower_level=Model::Winitlevel;

  // create an ordering for the input microjets
  // and calculate the fatjet momentum
  fastjet::PseudoJet fatjet(0,0,0,0);
  for (unsigned int i = 0; i < input.size(); ++i) {
    fatjet += input[i];
    // the new UserInfoBase type will be owned by PseudoJet, which deletes it when it goes out of scope
    input[i].set_user_info(new JetInfo(i, input[i].user_index()));
#ifdef DEBUGCODE
    LOG(DEBUG) << "mom " << std::endl;
    printjet(input[i]);
#endif
  }
  // no sort necessary so far ...
  //std::sort(input.begin(), input.end(), lessThanIndex);

  // according to eq. 3.1: Qsquare = pT_fatjet^2 + m_fatjet^2
  double Qsquare = square(fatjet.perp()) + square(fatjet.m());

  // store the input vector
  Storage store;
  store.set(input);

  wBackground = 0;
  wSignal = 0;

  // this embeds the powerset function to calculate ISR/SB combination
  // it will take much less time and memory to do this calculation inline, instead of creating
  // and storing all combinations (they are uncorrelated, so the following can be used)

  // This calculates the set of all combination of elements in v
  // This algorithm takes O(n * 2^n)

  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with v.size() bits with any number of bits activated
  // unsigned long long is 8 bytes long in most computers = 8*8 bits = 64 bits
  // this function works well as long as the maximum number of elements in the set is 64
  // (for more elements, a cleverer way can be thought of using more than one unsigned long long)

  unsigned int nElements = (unsigned int) input.size();
  unsigned int iElements = (unsigned int) param()[p_cut_n_sub];
  if (nElements < iElements)
    iElements = nElements;
  iElements = nElements;
  if (iElements > sizeof(unsigned long long)*8) {
    throw NEW_EXCEPTION("Number of microjets [maxElements] exceeds size of unsigned long long ([maxsize] bits)\
                         in this computer.")
                         .setParam("maxsize", to_string(sizeof(unsigned long long)*8))
                         .setParam("maxElements", to_string(iElements));
  }
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  // we generate the isrWeights and the signal weights first
  // then we sum the signal weight and check whether we should move on
  std::vector<double> isrWeight(nSets, 0);

  // the bits position in index i give the positions of the vector to be used
  for (unsigned long long i = 0; i < nSets; ++i) {
    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    // creating a vector in result itself and then changing it saves copy time

    StoredJet isrJets(store);
    StoredJet hardJets(store);

    fastjet::PseudoJet isrSum(0,0,0,0);
    fastjet::PseudoJet hardSum(0,0,0,0);

    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        isrJets.push_back(k);
        //copyJI(input[k], isrJets.back());
        isrSum += input[k];
      } else { // if a jet is not in the ISR set, it is in the hard proccess
        hardJets.push_back(k);
        //copyJI(input[k], hardJets.back());
        hardSum += input[k];
      }
    }

    // no sort needed so far ...
    //std::sort(isrJets.begin(), isrJets.end(), lessThanIndex);
    //std::sort(hardJets.begin(), hardJets.end(), lessThanIndex);
  
    // now go over this combination of isr/sb jets and calculate weights for ISR, signal and background
    // first for ISR:
    // fractionISR_of_hardInt_PT: value which determines the fraction of ISR_pT vs HardInt_pT
    if (square(isrSum.perp()) >= Qsquare*square(param()[p_fractionISR_of_hardInt_PT]) ) {
      continue; // demand that the ISR p_T is less than fractionISR_of_hardInt_PT times the virtuality
                // of the hard proccess
    }

#ifdef DEBUGCODE
    LOG(DEBUG) << "accepted ISR vs Signal branch" << endl;
    LOG(DEBUG) << "ISR_PT**2 : " << square(isrSum.perp()) << endl;
    LOG(DEBUG) << "Qsquare : " << Qsquare << endl;
    LOG(DEBUG) << "Qsquare*pow(fractionISR_of_hardInt_PT,2) : " << Qsquare*square(param()[p_fractionISR_of_hardInt_PT]) << endl;
    LOG(DEBUG) << square(isrSum.perp()) << " < " << Qsquare*square(param()[p_fractionISR_of_hardInt_PT]) << endl;

    LOG(DEBUG) << "ISRbranch" << endl; printoutput(isrJets);
    LOG(DEBUG) << "SBbranch" << endl; printoutput(hardJets);

    LOG(DEBUG) << "pT_ISR^2: " << square(isrSum.perp()) << "  <  " << Qsquare/4 << endl;
    LOG(DEBUG) << "input ISR: " << endl;
    printoutput(isrJets);
#endif

    double thisISRWeight = 1.0;
    if (param()[p_noISR] < 0.1) {
      m_isr.setQsquare(Qsquare);
      thisISRWeight = m_isr.weight(isrJets, hardSum)*m_isr.hamiltonian();
    }
    isrWeight[i] = thisISRWeight;

#ifdef DEBUGCODE
    LOG(DEBUG) << "signal hypothesis: " << endl;
    LOG(DEBUG) << "input hard interaction" << endl;
    printoutput(hardJets);
#endif

    m_signal.setQsquare(Qsquare);
    double thisSignalWeight = m_signal.weight(hardJets, hardSum);
    double thisSignalH = m_signal.hamiltonian(hardSum.perp());
    wSignal += thisSignalH*thisISRWeight*thisSignalWeight;

    m_signalWeight.insert(std::pair<double, std::vector<fastjet::PseudoJet> >(thisSignalH*thisISRWeight*thisSignalWeight, (std::vector<fastjet::PseudoJet>) (hardJets)));

#ifdef DEBUGCODE
    LOG(DEBUG) << "input Hard" << endl;
    printoutput(hardJets,DEBUG);
    LOG(DEBUG) << "isrweight: " << isrWeight[i] << endl;
    LOG(DEBUG) << "vH = topweight: " << thisSignalWeight << endl;
    LOG(DEBUG) << "topweight*isrweight: " << isrWeight[i]*thisSignalWeight << endl;
    LOG(DEBUG) << "Hamiltonian signal: " << thisSignalH << endl;
    LOG(DEBUG) << "Result of this history (topweight*Hamiltonina= h* vH): " << thisSignalH*thisSignalWeight << endl;
    LOG(DEBUG) << "Result of this history including ISR (topweight*Hamiltonina*isrweight): " << thisSignalH*thisSignalWeight*isrWeight[i] << endl;

    LOG(DEBUG) << "------------------------------------------------ " << endl;
    LOG(DEBUG) << endl;

    LOG(DEBUG) << endl;
#endif


  // don't bother calculating the background weight if the signal weight is negative
  // need to create two loops for this: is it worth it if the sample is signal-enriched?
  // commenting it out for now: it should be faster if the sample has mostly signal
  }

  if (!(wSignal > 0)) {
    wBackground = -1;
    return wSignal/wBackground;
  }

  // background weights
  for (unsigned long long i = 0; i < nSets; ++i) {

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    // creating a vector in result itself and then changing it saves copy time

    StoredJet hardJets(store);
    StoredJet isrJets(store);

    fastjet::PseudoJet isrSum(0,0,0,0);
    fastjet::PseudoJet hardSum(0,0,0,0);

    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        isrJets.push_back(k);
        //copyJI(input[k], isrJets.back());
        isrSum += input[k];
      } else { // if a jet is not in the ISR set, it is in the hard proccess
        hardJets.push_back(k);
        //copyJI(input[k], hardJets.back());
        hardSum += input[k];
      }
    }

    // no sort needed so far ...
    //std::sort(isrJets.begin(), isrJets.end(), lessThanIndex);
    //std::sort(hardJets.begin(), hardJets.end(), lessThanIndex);
  
    // now go over this combination of isr/sb jets and calculate weights for ISR, signal and background
    // first for ISR:
    // fractionISR_of_hardInt_PT: value which determines the fraction of ISR_pT vs HardInt_pT
    if (square(isrSum.perp()) >= Qsquare*square(param()[p_fractionISR_of_hardInt_PT]) ) {
      continue; // demand that the ISR p_T is less than fractionISR_of_hardInt_PT times the virtuality
                // of the hard proccess
    }

#ifdef DEBUGCODE
    LOG(DEBUG) << "background hypothesis: " << endl;
    LOG(DEBUG) << "input hard interaction" << endl; printoutput(hardJets);
#endif


    double thisBackgroundWeight = m_background.weight(hardJets, hardSum);
    double thisBackgroundH = m_background.hamiltonian(hardSum.perp());
    wBackground += thisBackgroundH*isrWeight[i]*thisBackgroundWeight;

#ifdef DEBUGCODE
    LOG(DEBUG) << "input ISR" << endl; printoutput(isrJets);
    LOG(DEBUG) << "isrweight: " <<  isrWeight[i]  << endl;
    LOG(DEBUG) << "backgroundweight: " <<  thisBackgroundWeight << endl;
    LOG(DEBUG) << "Hamiltonian backg: " << thisBackgroundH << endl;
    LOG(DEBUG) << "background*Hamiltonina: " << thisBackgroundH*thisBackgroundWeight << endl;
    LOG(DEBUG) << "------------------------------------------------ " << endl;
    LOG(DEBUG) << endl;
    LOG(DEBUG) << "isrweight*background*Hamiltonian: " << thisBackgroundH*thisBackgroundWeight*isrWeight[i]  << endl;
#endif
  }

  double finalWeight = 0;
  if (wBackground > 0)
    finalWeight = wSignal/wBackground;

  return finalWeight;
}

const std::multimap<double, std::vector<fastjet::PseudoJet> > &Deconstruction::Deconstruct::signalWeight() const {
  return m_signalWeight;
}

#include "Exception.h"

#include <string>
#include <map>
#include <iostream>

Deconstruction::Exception::Exception(const std::string &msg, int lineno, const std::string &sourceFile)
  : m_msg(msg), m_lineno(lineno), m_sourceFile(sourceFile) {
}

Deconstruction::Exception &Deconstruction::Exception::setParam(const std::string &key, const std::string &value) {
  m_param[key] = value;
  return *this;
}

std::string Deconstruction::Exception::what() const {
  std::string reason;
  reason += "Deconstruction::Exception in ";
  reason += m_sourceFile;
  reason += ", line ";
  reason += m_lineno;
  reason += ":\n";
  reason += fixedMsg();

  return reason;
}

std::string Deconstruction::Exception::fixedMsg() const {
  std::string msg = m_msg;

  for (std::map<std::string, std::string>::const_iterator it = m_param.begin(); it != m_param.end(); ++it) {
    std::string lookFor = "[";
    lookFor += it->first;
    lookFor += "]";
    findAndReplace(msg, lookFor, it->second);
  }
  return msg;
}

void Deconstruction::Exception::findAndReplace(std::string msg, const std::string &key, const std::string &value) const {
  size_t pos = 0;
  while ((pos = msg.find(key, pos)) != std::string::npos) {
     msg.replace(pos, msg.length(), value);
     pos += key.length();
  }
}

#include "HBBModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "Helper.h"

#include "Parameters.h"
#include "AnalysisParameters.h"

using namespace Deconstruction;

Deconstruction::HBBModel::HBBModel(Parameters &param)
  : Deconstruction::Model::Model(param) {
}

Deconstruction::HBBModel::~HBBModel() {
}

double Deconstruction::HBBModel::weight(const StoredJet &jets, fastjet::PseudoJet &sum) {

  double signalweight = 0;

  // there must be 2 microjets for the Higgs decay
  if (jets.size() < 2) return signalweight;

  if (std::fabs(jets.sum().m() - param()[p_Higgsmass]) > param()[p_delta_Higgsmass]) {
    return signalweight;
  } else { //assign weight H*exp(-S):
    signalweight = 4.0*square(M_PI)/param()[p_Higgsmass]/param()[p_delta_Higgsmass];
    signalweight /= pow(pow(jets.sum().m(), 2) - pow(param()[p_Higgsmass], 2), 2) + pow(param()[p_Higgsmass], 2)*pow(param()[p_delta_Higgsmass], 2);
  }

  double w = 0;

  unsigned int iElements = (unsigned int) jets.size();
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements
  StoredJet empty(jets.store());
  for (unsigned long long i = 1; i < nSets; ++i) {
    StoredJet ubranch(jets.store());
    StoredJet dbranch(jets.store());
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        ubranch.push_back(jets.getIdx(k));
      } else {
        dbranch.push_back(jets.getIdx(k));
      }
    }

    double valueleft = 0;
    double valueright = 0;
    m_calc.clear();
    valueleft = make_splitting(ubranch, empty, dbranch, jets, Flavour::b, Flavour::noflav, Flavour::noflav, Flavour::noflav, Shower::b);
    m_calc.clear();
    valueright = make_splitting(dbranch, ubranch, empty, jets, Flavour::bbar, Flavour::noflav, Flavour::noflav, Flavour::noflav, Shower::b);

    w += valueleft*valueright;
  }

  return w*signalweight;
}

double Deconstruction::HBBModel::hamiltonian(double pTsum) {
  double kH2 = square(pTsum);
  double Higgsmass2 = square(param()[p_Higgsmass]);

  double Hfj_sig = Cte::Npdf_signal*std::pow((Cte::pTmin2 + Higgsmass2)/(kH2 + Higgsmass2),Cte::Npdf_signal)/(kH2 + Higgsmass2);

  return Hfj_sig;
}

#include "Helper.h"

#include <vector>
#include <cmath>
#include "Exception.h"
#include <fastjet/PseudoJet.hh>

#include "JetInfo.h"

#include <iostream>
#include <iomanip>
#include "stdlib.h"

#include <sstream>

#include "StoredJet.h"

using namespace Deconstruction;
using std::cout;
using std::endl;
using std::max;
using std::setw;
using std::vector;

const double Deconstruction::Cte::pi          = 3.141592653589793;
const double Deconstruction::Cte::TR          = 0.5;
const double Deconstruction::Cte::CA          = 3.0;
const double Deconstruction::Cte::CF          = (square(Cte::CA) - 1.0)/2.0/Cte::CA;
const double Deconstruction::Cte::nf          = 5;
const double Deconstruction::Cte::b0          = (33.0 - 2.0*Cte::nf)/12.0/Cte::pi;

// default is 0.5
const double Deconstruction::Cte::cnp         = 1.0;
// default is 0.5
const double Deconstruction::Cte::kappa_np2   = 4.0;
// default is 1
const double Deconstruction::Cte::kappa_p2    = 4.0;
// default is 2
const double Deconstruction::Cte::nnp         = 1.5;
// default is 2
const double Deconstruction::Cte::cr          = 2.0;
// default is 4
const double Deconstruction::Cte::nr          = 1.0;

// minpT for top quark
const double Deconstruction::Cte::pTmin2      = square(400.0);
// new default is 1.4
const double Deconstruction::Cte::Npdf_backg  = 2.0;
const double Deconstruction::Cte::Npdf_signal = 2.0;

// sudakov top factor a
const double Deconstruction::Cte::alphattg    = 1.0;
// sudakov top factor b
const double Deconstruction::Cte::betattg     = 1.0;

const double Deconstruction::Cte::smalldouble = 1e-7;

//std::string to_string(unsigned int i) {
//  std::stringstream ss;
//  ss << i;
//  return ss.str();
//}

void Deconstruction::dummy() {
  std::vector<fastjet::PseudoJet>::iterator it;
  std::vector<std::vector<fastjet::PseudoJet> > v;
  std::vector<fastjet::PseudoJet> x;
  power_set_orig<fastjet::PseudoJet,  std::vector<fastjet::PseudoJet>::iterator >(it, it, v, 0, 0);
  find_subset_orig(it, it, 1, x, v);
}

template <class T , class Iter>
void  Deconstruction::find_subset_orig( Iter first , Iter last , int n , vector<T>& foo, vector<vector<T> > & result ) {

  if (n == 0) {
    result.push_back(foo);        
  } else {
    for ( Iter iter = first ; iter != last ; ++iter ) {
      foo.push_back( *iter ) ;
      ++iter ;
      find_subset_orig( iter , last , n-1 , foo, result ) ;
      --iter ;
      foo.pop_back() ;
    }
  }
}

template <class T , class Iter>
void  Deconstruction::power_set_orig( Iter first , Iter last, vector < vector < T > > & result, int start, int end) {
  vector<T>  sets ;
  for ( int i = start ; i <= end ; ++i ) {
    sets.resize(0) ;
    find_subset_orig( first , last , i , sets, result ) ;
  }
}


void Deconstruction::copyJI(const fastjet::PseudoJet &from, fastjet::PseudoJet &to) {
  to.set_user_info(new JetInfo(from.user_info<JetInfo>().i(), from.user_info<JetInfo>().user_index()));
}

/// does the actual work for printing out a jet
void Deconstruction::printoutput(const std::vector<fastjet::PseudoJet> &leftbranch) {
  LOG(DEBUG) << " { ";
  for(unsigned ii=0; ii<leftbranch.size(); ii++) {
    LOG(DEBUG) << leftbranch[ii].user_info<JetInfo>().i() << " ";
  }
  LOG(DEBUG) << " } " << endl;
}

void Deconstruction::printoutput(const std::vector<fastjet::PseudoJet> &leftbranch, MsgLevel x) {
  LOG(x) << " { ";
  for(unsigned ii=0; ii<leftbranch.size(); ii++) {
    LOG(x) << leftbranch[ii].user_info<JetInfo>().i() << " ";
  }
  LOG(x) << " } " << endl;
}
void Deconstruction::printoutput(const std::vector<fastjet::PseudoJet> &leftbranch, const std::string &s) {
  LOG(DEBUG) << s<< endl;
  LOG(DEBUG) << " { ";
  for(unsigned ii=0; ii<leftbranch.size(); ii++) {
    LOG(DEBUG) << leftbranch[ii].user_info<JetInfo>().i() << " ";
  }
  LOG(DEBUG) << " } " << endl;
}

void Deconstruction::printoutput(const std::vector<fastjet::PseudoJet> &leftbranch, const std::string &s, MsgLevel x) {
  if (x == FORCE) {
    std::cout << s<< endl;
    std::cout << " { ";
    for(unsigned ii=0; ii<leftbranch.size(); ii++) {
      std::cout << leftbranch[ii].user_info<JetInfo>().i() << " ";
    }
    std::cout << " } " << endl;
    return;
  }
  LOG(x) << s<< endl;
  LOG(x) << " { ";
  for(unsigned ii=0; ii<leftbranch.size(); ii++) {
    LOG(x) << leftbranch[ii].user_info<JetInfo>().i() << " ";
  }
  LOG(x) << " } " << endl;
}

void Deconstruction::printoutput(const StoredJet &leftbranch) {
  std::vector<fastjet::PseudoJet> lj = (std::vector<fastjet::PseudoJet>) leftbranch;
  printoutput(lj);
}

void Deconstruction::printoutput(const StoredJet &leftbranch, MsgLevel x) {
  std::vector<fastjet::PseudoJet> lj = (std::vector<fastjet::PseudoJet>) leftbranch;
  printoutput(lj, x);
}
void Deconstruction::printoutput(const StoredJet &leftbranch, const std::string &s) {
  std::vector<fastjet::PseudoJet> lj = (std::vector<fastjet::PseudoJet>) leftbranch;
  printoutput(lj, s);
}

void Deconstruction::printoutput(const StoredJet &leftbranch, const std::string &s, MsgLevel x) {
  std::vector<fastjet::PseudoJet> lj = (std::vector<fastjet::PseudoJet>) leftbranch;
  printoutput(lj, s, x);
}

void Deconstruction::printjet (const fastjet::PseudoJet &jet) {
  LOG(DEBUG) << "E, px, py, pz = "
       << " " << setw(10) << jet.e()
       << " " << setw(10) << jet.px()
       << " " << setw(10) << jet.py()
       << " " << setw(10) << jet.pz() << endl;
}

void Deconstruction::printjet (const fastjet::PseudoJet &jet, const std::string &s, Deconstruction::MsgLevel x) {
  if (x == FORCE) {
    std::cout << s <<endl;
    std::cout << "M, Pt, y, phi = "
       << " " << setw(10) << jet.m()
       << " " << setw(10) << jet.perp()
       << " " << setw(10) << jet.rap()
       << " " << setw(10) << jet.phi() << endl;
    return;
  }
  LOG(x) << s <<endl;
  LOG(x) << "M, Pt, y, phi = "
       << " " << setw(10) << jet.m()
       << " " << setw(10) << jet.perp()
       << " " << setw(10) << jet.rap()
       << " " << setw(10) << jet.phi() << endl;
}

// This calculates the set of all combination of elements in v
// This algorithm takes O(n * 2^n)
template <class T>
void Deconstruction::powerset(const std::vector<T> &v,
                              std::vector<std::vector<T> > &result,
                              const unsigned int minElements,
                              const unsigned int maxElements) {
  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with v.size() bits with any number of bits activated
  // unsigned long long is 8 bytes long in most computers = 8*8 bits = 64 bits
  // this function works well as long as the maximum number of elements in the set is 64
  // (for more elements, a cleverer way can be thought of using more than one unsigned long long)
  unsigned int nElements = (unsigned int) v.size();
  if (maxElements > sizeof(unsigned long long)*8) {
    throw NEW_EXCEPTION("Using powerset function with maxElements = [maxElements] and number of elements \
                         in vector = [nElements]. The maximum size of sets to be used in the powerset \
                         function is [maxsize], given by the size in bits of the unsigned long long in \
                         this computer.")
                         .setParam("maxsize", to_string(sizeof(unsigned long long)*8))
                         .setParam("maxElements", to_string(maxElements))
                         .setParam("nElements", to_string(nElements));
  }
  unsigned int iElements = std::min(maxElements, nElements);
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  // the bits position in index i give the positions of the vector to be used
  for (unsigned long long i = 0; i <= nSets; ++i) {
    // calculate if the sum of bits in i is >= minElements
    if (minElements >= 1) {
      // ok, we really need to test this
      unsigned int sumBits = numberOfSetBits(i);
      if (sumBits < minElements)
        continue;
    }
    // } else { } // no minimum number of elements ... ignore the previous check

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    // creating a vector in result itself and then changing it saves copy time
    result.push_back(std::vector<T>());
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0)
        result.back().push_back(v[k]);
    }
  }
}

unsigned int Deconstruction::numberOfSetBits(unsigned long long i) {
  // use popcount algorithm to count number of bits set in i
  // some processors have an instruction that does that using a lookup table
  // if this is not available gcc has a very efficient implementation of this
  // This requires gcc >= 3.4
  return __builtin_popcount(i);
}

double Deconstruction::square(const double &x) {
  return x*x;
}

fastjet::PseudoJet Deconstruction::sum(const std::vector<fastjet::PseudoJet> &v) {
  fastjet::PseudoJet result(0,0,0,0);
  for (int i = 0; i < v.size(); ++i) {
    result += v[i];
  }
  return result;
}

double Deconstruction::square_p1minusp2(const fastjet::PseudoJet & p1, const fastjet::PseudoJet & p2) {
  return ((p1.e()-p2.e())*(p1.e()-p2.e()) - (p1.px()-p2.px())*(p1.px()-p2.px()) - (p1.py()-p2.py())*(p1.py()-p2.py()) - (p1.pz()-p2.pz())*(p1.pz()-p2.pz()));
}

double Deconstruction::scalprod(const fastjet::PseudoJet & p1, const fastjet::PseudoJet & p2) {
  return (p1.e()*p2.e() - p1.px()*p2.px() - p1.py()*p2.py() - p1.pz() * p2.pz());
}

double Deconstruction::alphas(const double scale2) {
  double mZ2= square(91.1876);
  double alphasMZ = 0.118;
  
  double als = alphasMZ/(1.0 + alphasMZ*Cte::b0*std::log(scale2/mZ2));
  return als;
}

double dot(const fastjet::PseudoJet &x, const fastjet::PseudoJet &y) {
  return x.e()*y.e() - x.px()*y.px() - x.py()*y.py() - x.pz()*y.pz();
}

#include "ISRModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "Helper.h"

#include "Parameters.h"

#include <memory>
#include <algorithm>

#include "partition_deconstruction.h"

using namespace Deconstruction;
using namespace std;
using namespace fastjet;

Deconstruction::ISRModel::ISRModel(Parameters &param)
  : Deconstruction::Model::Model(param) {
}

Deconstruction::ISRModel::~ISRModel() {
}

//double Deconstruction::ISRModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::ISRModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {

  double Qsquare = m_Qsquare;

  double isr_sudakov = 1.0; // simplified version -- cancels anyway between signal and background!
  // 0) calculate total Sudakov, which tells how likely it is
  //    that vector<PseudoJet*> input ended up in the cone.
  //    This Sudakov cancels in S/B...

  if (jets.size() == 0)
    return isr_sudakov;

  vector<int> v = jets.getList();

  int iisize((int)jets.size());

  vector< vector< vector<int> > > partonsetmix;
  partonsetmix.clear();
  try
    {
      partition::iterator it(iisize);
      while(true)
        {
          // cout << it << " : " << it.subsets() << " : ";
          auto_ptr<vector<vector<int> > > part = it[v];
          partonsetmix.push_back(*part);

          ++it;
        }
    }
  catch (overflow_error&)
    {
      //return(0.0);
    }


  // now we have all combinations for the ISR which enter the FSR shower.

  // man hat: H*FSR_shower fuer jede moegliche ISR-Abstrahlung-- 
  //          der Sudakov ist nach unserer neuen definition global 
  //          gleich und cancelt zwischen Signal und Untergrund.
  //          ISR_sudakov damit = 1;

  double ISRweight = 0.0; // muss summiert werden
  for(unsigned ii=0; ii<partonsetmix.size(); ii++)
    {
      if(partonsetmix[ii].size()==0) continue; //faengt addition von ISRtmp=1 ab
      double ISRtmp(1.0); // muss multipliziert werden
      for(unsigned jj=0; jj<partonsetmix[ii].size(); jj++)
        {
          // all weights for each vector in partonsetmix is summed up

          // weight for each ISR radiation:
          // H*Sudakov_tot*FSRshower
          StoredJet jets_ij(jets.store(), partonsetmix[ii][jj]);
          PseudoJet tmp_J = jets_ij.sum();

          // theta cut in ISR hamiltonian:
          // has to jump out of the loop with ISRtmp=0.0
          // because if in {{1} {2}} the cut is failed by {2} only, the whole
          // shower history is rejected....
          if(tmp_J.perp() <= tmp_J.m() )
            {
              ISRtmp=0.0;
              break;
            }


          // ISR/UE hamiltonian, non-pert part und pert part
          double hamiltonian(8.0*M_PI*Cte::CA*alphas(square(tmp_J.perp())+Cte::kappa_p2)
                             /(square(tmp_J.perp()) + Cte::kappa_p2)
                             /pow(1.0+Cte::cr*tmp_J.perp()/sqrt(Qsquare),Cte::nr)
                             + 16.0*M_PI*Cte::cnp*pow(Cte::kappa_np2,Cte::nnp-1)/pow(square(tmp_J.perp()) + Cte::kappa_np2,Cte::nnp));

          //ISRtmp *= hamiltonian * start_splitting(partonsetmix[ii][jj], Flavour::g, Shower::QCD);
          ISRtmp *= hamiltonian * start_splitting(jets_ij, Flavour::g, Shower::QCD);
        }

      ISRweight += ISRtmp;
    }

  return ISRweight*isr_sudakov;

}

double Deconstruction::ISRModel::hamiltonian(double pTsum) {
  return 1.0;
}

#include "JetInfo.h"

Deconstruction::JetInfo::JetInfo()
  : fastjet::PseudoJet::UserInfoBase(), m_i(0), m_user_index(0) {
}

Deconstruction::JetInfo::JetInfo(unsigned int i, int user_index)
  : fastjet::PseudoJet::UserInfoBase(), m_i(i), m_user_index(user_index) {
}

Deconstruction::JetInfo::~JetInfo() {
}

unsigned int &Deconstruction::JetInfo::i() {
  return m_i;
}

int &Deconstruction::JetInfo::user_index() {
  return m_user_index;
}

const unsigned int &Deconstruction::JetInfo::i() const {
  return m_i;
}

const int &Deconstruction::JetInfo::user_index() const {
  return m_user_index;
}

bool Deconstruction::lessThanIndex(const fastjet::PseudoJet &i, const fastjet::PseudoJet &j) {
  return (i.user_info<Deconstruction::JetInfo>().i() < j.user_info<Deconstruction::JetInfo>().i());
}

#include "Message.h"
#include <streambuf>

Deconstruction::nullbuf Deconstruction::MsgInterface::_nullbuf = Deconstruction::nullbuf();

Deconstruction::nullstream Deconstruction::MsgInterface::nullStream = Deconstruction::nullstream();
Deconstruction::Message Deconstruction::MsgInterface::msg = Deconstruction::Message();

Deconstruction::nullbuf::nullbuf()
  : std::streambuf() {
}

Deconstruction::nullbuf::nullbuf(const Deconstruction::nullbuf &c)
  : std::streambuf() {
}

Deconstruction::nullstream::nullstream()
: std::ostream(&Deconstruction::MsgInterface::_nullbuf) {
}

Deconstruction::nullstream::nullstream(const Deconstruction::nullstream &c)
: std::ostream(&Deconstruction::MsgInterface::_nullbuf) {
}

std::ostream &Deconstruction::nullstream::m(const Deconstruction::MsgLevel x) {
  return *this;
}

Deconstruction::Message::Message()
  : std::ostream(std::cout.rdbuf()) {
  m_level = INFO;
}

Deconstruction::Message::Message(const Deconstruction::Message &c)
  : std::ostream(std::cout.rdbuf()) {
  m_level = c.m_level;
}

Deconstruction::Message::~Message() {
}

void Deconstruction::Message::setLevel(const Deconstruction::MsgLevel level) {
  m_level = level;
}

const Deconstruction::MsgLevel &Deconstruction::Message::level() const {
  return m_level;
}

std::ostream &Deconstruction::Message::m(const Deconstruction::MsgLevel level) {
  if (level >= m_level)
    return *this;
  return Deconstruction::MsgInterface::nullStream;
}

#include "Model.h"

#include "Helper.h"

#include "Parameters.h"
#include "AnalysisParameters.h"

#include "TopModel.h"

#include "Exception.h"

#include <fastjet/PseudoJet.hh>

#include <string>
#include <iostream>

#include <iomanip>

#include <cmath>

#include "JetInfo.h"

#include "StoredJet.h"
#include "StoredCalculations.h"

using namespace Deconstruction;

using std::endl;
using std::max;
using std::fabs;
using std::min;
using std::vector;
using namespace fastjet;

int Model::level = -1;
int Model::btopshower_level= -1;
int Model::tshower_level = -1;
int Model::Wshower_level = -1;

int Model::initlevel = -1;
int Model::binitlevel = -1;
int Model::tinitlevel = -1;
int Model::Winitlevel = -1;

StoredCalculations Model::m_calc = StoredCalculations();

Deconstruction::Model::Model(Parameters &p)
  : m_param(p) {

  nrbreitw = 2.0;

  wwidth = std::fabs(square(param()[p_wmass]+param()[p_wwindow]) - square(param()[p_wmass]))/nrbreitw/param()[p_wmass];
  delta_wmass = nrbreitw*wwidth;

  topwidth = std::fabs(square(param()[p_topmass]+param()[p_topwindow]) - square(param()[p_topmass]))/nrbreitw/param()[p_topmass];
  delta_topmass = nrbreitw*topwidth;

  dR2_outside = square(param()[p_R]);

  lambda_mu = param()[p_lambda_mu_ext];
  topmass = param()[p_topmass];

  m_useBtag = param()[p_useBtag];
  m_tagprob = param()[p_tagprob];
  m_fakeprob = param()[p_fakeprob];
}

void Deconstruction::Model::setQsquare(double q) {
  m_Qsquare = q;
}
 
Deconstruction::Model::~Model() {
}

Parameters &Deconstruction::Model::param() const {
  return m_param;
}

double Deconstruction::Model::start_splitting(const StoredJet &input, const Flavour::id flavour, const Shower::id shower) {
  m_calc.clear();

  StoredJet empty(input.store());
  return make_splitting(input, empty, empty, empty, flavour, Flavour::noflav, Flavour::noflav, Flavour::noflav, shower);
}

double Deconstruction::Model::make_splitting(const StoredJet &input, const StoredJet &leftcolpartner, const StoredJet &rightcolpartner, const StoredJet &grandmother, const Flavour::id flavor, const Flavour::id granflavor, const Flavour::id leftcolflav, const Flavour::id rightcolflav, const Shower::id shower) {

  // check if the calculation was already done
  std::vector<int> a = input.getList();
  std::sort(a.begin(), a.end());
  std::vector<int> b = leftcolpartner.getList();
  std::sort(b.begin(), b.end());
  std::vector<int> c = rightcolpartner.getList();
  std::sort(c.begin(), c.end());
  std::vector<int> d = grandmother.getList();
  std::sort(d.begin(), d.end());
  StoredKey sk;
  for (int i = 0; i < a.size(); ++i) sk.push_back(a[i]);
  sk.push_back(-1);
  for (int i = 0; i < b.size(); ++i) sk.push_back(b[i]);
  sk.push_back(-1);
  for (int i = 0; i < c.size(); ++i) sk.push_back(c[i]);
  sk.push_back(-1);
  for (int i = 0; i < d.size(); ++i) sk.push_back(d[i]);
  sk.push_back(-1);
  sk.push_back(flavor);
  sk.push_back(granflavor);
  sk.push_back(leftcolflav);
  sk.push_back(rightcolflav);
  sk.push_back(shower);
  double w = 0;
  if (m_calc.check(sk, w)) {
    return w;
  }

  ++Model::level;
  if (shower==Shower::W) ++Model::Wshower_level;
  if (shower==Shower::b) ++Model::btopshower_level;
  if (shower==Shower::t) ++Model::tshower_level;

#ifdef DEBUGCODE
  LOG(DEBUG) << "input make_splitting - level: " << Model::level<<  endl;
  if(shower==Shower::W) LOG(DEBUG) << "Wshower_level: " << Model::Wshower_level << endl;
  if(shower==Shower::b) LOG(DEBUG) << "Model::btopshower_level: " << Model::btopshower_level << endl;
  if(shower==Shower::t) LOG(DEBUG) << "tshower_level: " << Model::tshower_level << endl;

  printoutput(input,"input");
  printoutput(leftcolpartner,"leftcolpartner");
  printoutput(rightcolpartner,"rightcolpartner");
  printoutput(grandmother,"grandmother");
  LOG(DEBUG) << "flavor: " << flavor << endl;
  LOG(DEBUG) << "granflavor: " << granflavor << endl;
  LOG(DEBUG) << "leftcolflavor: " << leftcolflav << endl;
  LOG(DEBUG) << "rightcolflavor: " << rightcolflav << endl;
  LOG(DEBUG) << "shower: " << shower << endl;

#endif

  /// top has to decay into at least 3 particles /////////
  if ((flavor==Flavour::t || flavor==Flavour::tbar) && (input.size() < 3) ) {
    --Model::level;
    if (shower==Shower::W) --Model::Wshower_level;
    if (shower==Shower::b) --Model::btopshower_level;
    if (shower==Shower::t) --Model::tshower_level;
    m_calc.store(sk, 0);
    return 0;
  }
  //////////////////////////////////////////////////////////////////////////////
  /////////// if last particle of one branch //////////////////////////////////
  if (input.size() < (unsigned) 2) {
    if (input.size() == 0) {
      --Model::level;
      if (shower==Shower::W) --Model::Wshower_level;
      if (shower==Shower::b) --Model::btopshower_level;
      if (shower==Shower::t) --Model::tshower_level;
      m_calc.store(sk, 0);
      return 0;
    }

    PseudoJet tmptot(input.sum());
    PseudoJet tmpgrandm(grandmother.sum());
    PseudoJet tmpmrightcol(rightcolpartner.sum());
    PseudoJet tmpmleftcol(leftcolpartner.sum());

    double weight_lastparton(0.0);

    if (flavor==Flavour::g)
      weight_lastparton = sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor, granflavor, leftcolflav,rightcolflav,shower);
    else if ((flavor == Flavour::q) || (flavor == Flavour::qbar) || (flavor == Flavour::b) || (flavor == Flavour::bbar))
      weight_lastparton = sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower);
    else
      LOG(ERROR) << "error in sudakov make_splitting: wrong flavor" << endl;

    // for radiation off b-quark after top decay we need to know 
    // the flavor of the color connected partner (here top)

    LOG(DEBUG) << "sudakov value last parton: " << weight_lastparton << endl;

    double btaggerweight(1.0);
    if (m_useBtag) {
      double isNotB = 0;
      if ( (flavor != Flavour::b) && (flavor != Flavour::bbar)) isNotB = 1.0;
      if (tmptot.user_index() == 0)   btaggerweight = 1;     // not tried to be b-tagged
      else if (tmptot.user_index() < 0) btaggerweight = 1 - (m_tagprob*(1-isNotB) + m_fakeprob*isNotB );
      else if (tmptot.user_index() > 0) btaggerweight = m_tagprob*(1-isNotB) + m_fakeprob*isNotB;
    }
    // construct here sudakov if input is last parton (microjet)
	  
    LOG(DEBUG) << "exp(-sudakov)*btag last parton: " << weight_lastparton*btaggerweight << endl;

    --Model::level;
    if(shower==Shower::W) --Model::Wshower_level;
    if(shower==Shower::b) --Model::btopshower_level;
    if(shower==Shower::t) --Model::tshower_level;
 
    m_calc.store(sk, weight_lastparton*btaggerweight);
    return(weight_lastparton*btaggerweight); // this has to be changed to the sudakov
  }
  
  // return value
  double valuetotal(0.0);
  
  /////////////////// the top quark is special: it decays always before hadronization 
  /////////////////// thus we do not allow the top to radiate down to the hadronization scale
  /////////////////// in every step it has to be checked if the top decays and the weight has
  /////////////////// to be added to the t -> t g splitting weight.

  double topdecayweight(0.0);
  double antitopdecayweight(0.0);

  if(flavor == Flavour::t && (input.size()>2) ) {
    PseudoJet tmptot(input.sum());
      
    if (std::fabs(square(tmptot.m())-square(topmass))<topmass*delta_topmass) {
      PseudoJet tmpgrandm(grandmother.sum());
      PseudoJet tmpmrightcol(rightcolpartner.sum());
      PseudoJet tmpmleftcol(leftcolpartner.sum());
	  
      topdecayweight = sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,Flavour::t, granflavor,leftcolflav,rightcolflav,shower);                             // multiply sudakov for top before decay
      double sudtop=topdecayweight;
      double topdecaymodel = 0;
      TopModel tm(param(), flavor);
      topdecaymodel=tm.weight(input, tmptot);
      topdecayweight *= topdecaymodel;

#ifdef DEBUGCODE
      LOG(DEBUG) << "input in top decay model: " << endl;
      printoutput(input);
      LOG(DEBUG) << "results: " << endl;
      LOG(DEBUG) << "exp(-St) top vor decay: " << sudtop << endl;
      LOG(DEBUG) << "Htopdecaymodel: " << topdecaymodel << endl;
      LOG(DEBUG) << "Htopdecaymodel*exp(-St): " << topdecayweight << endl;
#endif

      valuetotal += topdecayweight;
    }
  }


  if(flavor == Flavour::tbar && (input.size()>2)) {
    PseudoJet tmptot(input.sum());
      
    if (std::fabs(square(tmptot.m())-square(topmass))<topmass*delta_topmass) {
      PseudoJet tmpgrandm(grandmother.sum());
      PseudoJet tmpmrightcol(rightcolpartner.sum());
      PseudoJet tmpmleftcol(leftcolpartner.sum());
	  
      antitopdecayweight = sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,Flavour::tbar, granflavor,leftcolflav,rightcolflav,shower);

      // multiply sudakov for antitop before decay 
      
      TopModel tm(param(), flavor);
      antitopdecayweight *= tm.weight(input, tmptot);
      valuetotal += antitopdecayweight;
    }
  }
  ///////////// until here decays /////////////////////////////

  StoredJet empty(input.store());

  /////////// now start splittings ////////////////////////////
  /////// split into left and right branches /////////////////////////

  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with input.size() bits with any number of bits activated
  unsigned int iElements = (unsigned int) input.size();
  // iElements cannot be greater than the maximum,
  // but this has already been tested before, when generating the powerset for ISR/SB
  // there is no point in testing it again, since these sets must be smaller or equal in size
  // than the subsets of the powerset of all microjets
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  // the bits position in index i give the positions of the vector to be used
  // we are skipping the empty set combination here
  // that's why i starts at one and ends at nSets-1
  for (unsigned long long i = 1; i < nSets; ++i) {
    //std::vector<fastjet::PseudoJet> leftbranch;
    //std::vector<fastjet::PseudoJet> rightbranch;

    StoredJet leftbranch(input.store());
    StoredJet rightbranch(input.store());

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        leftbranch.push_back(input.getIdx(k));
      } else {
        rightbranch.push_back(input.getIdx(k));
      }
    }

    PseudoJet tmptot(input.sum());
    PseudoJet tmpL(leftbranch.sum());
    PseudoJet tmpR(rightbranch.sum());
    PseudoJet tmpgrandm(grandmother.sum());
    PseudoJet tmpmrightcol(rightcolpartner.sum());
    PseudoJet tmpmleftcol(leftcolpartner.sum());
      
#ifdef DEBUGCODE
    LOG(DEBUG) << "split left right: " << endl;
    LOG(DEBUG) << "left: " << endl;
    printoutput(leftbranch);
    LOG(DEBUG) << "right: " << endl;
    printoutput(rightbranch);
#endif

    // next make_splittings have to be called with input (mother) as argument for grandmother...
    double valueleft(0.0);
    double valueright(0.0);
    double valuegbbleft(0.0);
    double valuegbbright(0.0);
    double valuegqqleft(0.0);
    double valuegqqright(0.0);
    
    double splitfactorgqq(0.0),splitfactorggg(0.0),splitfactorgbb(0.0);
    double splitfactorbbg(0.0),splitfactorqqg(0.0),splitfactorantibbg(0.0),splitfactorantiqqg(0.0);
    double splitfactorttg(0.0),splitfactorantittg(0.0);
    
    double checksudtot(0.0);
    double checksplittot(0.0);

    double gbb(0.0),ggg(0.0),gqq(0.0),qqg(0.0),antiqqg(0.0),bbg(0.0),antibbg(0.0);
    double gtotal(0.0), ttg(0.0),antittg(0.0);
     
    // IMPORTANT
    // if Flavour::g splits to quarks, the anti-quark goes left and the quark right
    // if Higgs splits to quarks, the anti-quark goes right and the quark left
    
    if(flavor == Flavour::g) // Flavour::g splits
    {
      //msp: ATTENTION! for massive b I might need to include new check for splitting (including massive b to calculate virtuality)

      //check if splitting allowed
	  
      if( ((square(tmpL.m())/tmpL.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5)) ||
          ((square(tmpR.m())/tmpR.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5)) )
        continue;

      // calc sudakov for this "propagator line"
      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor, granflavor,leftcolflav,rightcolflav,shower));

      /// g-> g g splitting
#ifdef DEBUGCODE
      LOG(DEBUG) << "do g->g g splitting " << endl;
#endif

      splitfactorggg = HGluonGG(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

#ifdef DEBUGCODE
      LOG(DEBUG) << "for g->gg splitting: " << endl;
      LOG(DEBUG) << "Hgqq: " << splitfactorggg << endl;
      LOG(DEBUG) << "sudakov: " << sudtot << endl;
      LOG(DEBUG) << "Hggg*sudakov: " << sudtot*splitfactorggg << endl;
#endif

      valueleft = make_splitting(leftbranch,leftcolpartner, rightbranch,input,Flavour::g,flavor,leftcolflav,Flavour::g,shower);
      valueright = make_splitting(rightbranch,leftbranch, rightcolpartner,input,Flavour::g,flavor,Flavour::g,rightcolflav,shower);
	  
      ggg = valueleft*valueright*splitfactorggg;

      // g->q q splitting
#ifdef DEBUGCODE
      LOG(DEBUG) << "do g->q q splitting " << endl;	       
#endif

      // TO speed up the code HGluonQQ can be replaced by HGluonBB and simply multiplied by nf-1
      // HGluonQQ=HGluonBB*(nf-1)
      splitfactorgqq = (Cte::nf-1)*HGluonQQ(tmptot, tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor, granflavor, leftcolflav , rightcolflav, shower);

#ifdef DEBUGCODE	 
      LOG(DEBUG) << "for g->qq splitting: " << endl;
      LOG(DEBUG) << "Hgqq: " << splitfactorgqq/Cte::nf << endl;
      LOG(DEBUG) << "nf*Hgqq: " << splitfactorgqq << endl;
      LOG(DEBUG) << "sudakov: " << sudtot << endl;
      LOG(DEBUG) << "nf*Hgqq*sudakov: " << sudtot*splitfactorgqq << endl;
#endif

      valuegqqleft = make_splitting(leftbranch,leftcolpartner, empty,input,Flavour::qbar,flavor,leftcolflav,Flavour::noflav,shower);
      valuegqqright = make_splitting(rightbranch, empty, rightcolpartner,input,Flavour::q,flavor,Flavour::noflav,rightcolflav,shower);
	  
      gqq = valuegqqleft*valuegqqright*splitfactorgqq;
	
      // g->b b splitting
      LOG(DEBUG) << "do g->b b splitting " << endl;
      splitfactorgbb = HGluonBB(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor, granflavor, leftcolflav , rightcolflav, shower); 
      LOG(DEBUG) << "H g->bb : " << splitfactorgbb << endl;
      LOG(DEBUG) << "H exp(-S) : " << sudtot*splitfactorgbb << endl;

      valuegbbleft = make_splitting(leftbranch, leftcolpartner, empty, input,Flavour::bbar,flavor,leftcolflav,Flavour::noflav,shower);
      valuegbbright = make_splitting(rightbranch, empty, rightcolpartner,input,Flavour::b,flavor,Flavour::noflav,rightcolflav,shower);
	      	      
      gbb = valuegbbleft*valuegbbright*splitfactorgbb;

      gtotal = (gbb+gqq+ggg)*sudtot;

#ifdef DEBUGCODE
      LOG(DEBUG) << "flavor: " << flavor << endl;
      printoutput(input, "total branch");
      printoutput(leftbranch,"left branch");
      printoutput(rightbranch,"right branch");
      LOG(DEBUG) << "in Flavour::gsplitting: " << endl;
      LOG(DEBUG) << "ggg: valueleft - " << valueleft << endl;
      LOG(DEBUG) << "ggg: valueright - " << valueright << endl;
      LOG(DEBUG) << "sudakov factor for this propagator line: " << sudtot << endl;
      LOG(DEBUG) << "ggg: Hggg - " << splitfactorggg << endl;
      LOG(DEBUG) << "exp(-S) : " << sudtot << endl;
      LOG(DEBUG) << "Hggg*exp(-s): " << splitfactorggg*sudtot << endl;
      LOG(DEBUG) << "valueleft*valueright*H " << ggg << endl;

      LOG(DEBUG) << "gqq: valuegqqleft - " << valuegqqleft << endl;
      LOG(DEBUG) << "gqq: valuegqqright - " << valuegqqright << endl;
      LOG(DEBUG) << "gqq: Hgqq - " << splitfactorgqq << endl;
      LOG(DEBUG) << "exp(-S) : " << sudtot << endl;
      LOG(DEBUG) << "Hgqq*exp(-S) - " << splitfactorgqq*sudtot << endl;
      LOG(DEBUG) << "valueleft*valueright*H " << gqq << endl;

      LOG(DEBUG) << "Hgbb: " << splitfactorgbb << endl;
      LOG(DEBUG) << "Hgqq: " << splitfactorgqq << endl;
      LOG(DEBUG) << "Hggg: " << splitfactorggg << endl;
      LOG(DEBUG) << "Hggg*exp(-S): " << splitfactorggg*sudtot << endl;
      LOG(DEBUG) << "Hgbb*exp(-S): " << splitfactorgbb*sudtot << endl;
      LOG(DEBUG) << "Hgqq*exp(-S): " << splitfactorgqq*sudtot << endl;
      LOG(DEBUG) << "(Hggg+nf*Hgqq)*exp(-S): " << (splitfactorgqq+splitfactorgbb+splitfactorggg)*sudtot << endl;
      LOG(DEBUG) << "gbb : " << gbb << endl;
      LOG(DEBUG) << "gqq : " << gqq << endl;
      LOG(DEBUG) << "ggg : " << ggg << endl;
      LOG(DEBUG) << "sudakov: " << sudtot << endl;
      LOG(DEBUG) << "gbb*exp(-S) : " << gbb*sudtot << endl;
      LOG(DEBUG) << "gqq*exp(-S) : " << gqq*sudtot << endl;
      LOG(DEBUG) << "ggg*exp(-S) : " << ggg*sudtot << endl;
      LOG(DEBUG) << "(gqq+ggg)*exp(-S): " << (gqq+ggg)*sudtot << endl;
      LOG(DEBUG) << "for gtotal: " << gtotal << endl;
#endif

    } else if(flavor == Flavour::t) {
      // top can only radiate if for the decay 3 subjets remain
      if (input.size() < 3)
        continue;

      // definition of virtuality $\mu_J^2 = (p_A + p_B)^2 - m_J^2$
      double muJ2(square(tmptot.m())-square(topmass));
      double muh2(std::fabs(square(tmpL.m())-square(topmass))); //top

      //msp changed in last version check!! if(muh2 < topmass*topwidth) muh2 = topmass*topwidth;
      // muh2 = topmass*topwidth;
      double mus2(square(tmpR.m())); //Flavour::g
	  
      //check splitting with massive daughters
      //check if splitting allowed

      if (((std::fabs(muh2)/tmpL.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5)) ||
          ((std::fabs(mus2)/tmpR.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5))) 
	 continue;

      // new cut:
      /*
      if(granflavor == noflav) {
        if((2*lambda_mu*(muJ2/tmptot.perp())) > (2*(muJ2+square(topmass))/tmptot.perp())) continue;
      } else if(granflavor == Flavour::t) {
        if((2*lambda_mu*(muJ2/tmptot.perp())) > (square(tmpgrandm.m())-square(topmass))/tmpgrandm.perp()) continue;
      } else
        LOG(DEBUG) << "ERROR grandmother no Flavour::t flavor although jet is Flavour::t" << endl;
      */
      // hier muss man immer beide gewichte addieren: splitting t -> t + g und top decay

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));

      valueleft = make_splitting(leftbranch, leftcolpartner, rightbranch,input,Flavour::t,flavor,Flavour::noflav,Flavour::g,shower);
      valueright = make_splitting(rightbranch,leftbranch, rightcolpartner,input,Flavour::g,flavor,Flavour::t,rightcolflav,shower);

      splitfactorttg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      ttg = valueleft*valueright*splitfactorttg*sudtot;
    } else if(flavor == Flavour::tbar) {
      // antitop can only radiate if for the decay 3 subjets remain 
      //if (input.size() < 3)
      if (input.size() < 3)
        continue;

      // definition of virtuality $\mu_J^2 = (p_A + p_B)^2 - m_J^2$
      double muJ2(square(tmptot.m())-square(topmass));
      double muh2(std::fabs(square(tmpR.m())-square(topmass))); //top
      // if(muh2 < topmass*topwidth) muh2 = topmass*topwidth;
      // muh2 = topmass*topwidth;
      double mus2(square(tmpL.m())); //Flavour::g

      //check splitting with massive daughters
      //check if splitting allowed
      if (((std::fabs(muh2)/tmpR.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5)) ||
         ((std::fabs(mus2)/tmpL.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5)) )
        continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));

      valueleft = make_splitting(leftbranch, leftcolpartner, rightbranch,input,Flavour::g,flavor,leftcolflav,Flavour::tbar,shower);
      valueright = make_splitting(rightbranch,leftbranch, empty,input,Flavour::tbar,flavor,Flavour::g,Flavour::noflav,shower);
	  
      splitfactorantittg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      antittg = valueleft*valueright*splitfactorantittg*sudtot;
    } else if(flavor == Flavour::bbar) // antib splits always antib (right) Flavour::g (left)
    {
      // definition of virtuality $\mu_J^2 = (p_A + p_B)^2 - m_J^2$
      // however take b massless 
      double muJ2(square(tmptot.m()));
      double muh2(square(tmpR.m())); //antib
      double mus2(square(tmpL.m())); //Flavour::g
	   
      if (((std::fabs(muh2)/tmpR.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5)) ||
          ((std::fabs(mus2)/tmpL.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5)) )
	      continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));

      // left splitting inherits leftcolpartner from mother
      // and the right branch is the right color partner
      valueleft = make_splitting(leftbranch, leftcolpartner, rightbranch,input,Flavour::g,flavor,leftcolflav,Flavour::bbar,shower);
      // right splitting inherits rightcolpartner from mother
      // and the left branch is the left color partner
      valueright = make_splitting(rightbranch, leftbranch, empty,input,Flavour::bbar,flavor,Flavour::g,Flavour::noflav,shower);

      splitfactorantibbg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      antibbg = valueleft*valueright*splitfactorantibbg*sudtot;
    } else if(flavor == Flavour::b) // b splits always b (right) Flavour::g (left)
    {
      // definition of virtuality $\mu_J^2 = (p_A + p_B)^2 - m_J^2$
      double muJ2(square(tmptot.m()));
      double muh2(square(tmpL.m())); //b
      double mus2(square(tmpR.m())); //Flavour::g
	  
      //check splitting with massive daughters
      //check if splitting allowed
      if (((std::fabs(muh2)/tmpL.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5)) ||
          ((std::fabs(mus2)/tmpR.perp()) >= lambda_mu*(muJ2/tmptot.perp()*0.5)) )
        continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));
	  
      valueleft = make_splitting(leftbranch, empty, rightbranch,input,Flavour::b,flavor,Flavour::noflav,Flavour::g,shower);
      valueright = make_splitting(rightbranch,leftbranch, rightcolpartner,input,Flavour::g,flavor,Flavour::b,rightcolflav,shower);

      splitfactorbbg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      bbg = valueleft*valueright*splitfactorbbg*sudtot;
		
#ifdef DEBUGCODE
      LOG(DEBUG) << " flavor: " << flavor << endl;
      LOG(DEBUG) << " sudakov: " << sudtot << endl;
      LOG(DEBUG) << " leftbranch: " << endl;
      printoutput(leftbranch);
      LOG(DEBUG) << " rightbranch: " << endl;
      printoutput(rightbranch);
      LOG(DEBUG) << "splitfactorbbg: " << splitfactorbbg << endl;
      LOG(DEBUG) << "splitfactorbbg*sudtot: " << splitfactorbbg*sudtot << endl;
      LOG(DEBUG) << "bbg: " << bbg << endl;
#endif
    }  else if(flavor == Flavour::qbar) // antiq splits antiq (right) Flavour::g(left)
    {
      if (((square(tmpL.m())/tmpL.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5)) ||
          ((square(tmpR.m())/tmpR.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5)) )
        continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));

      splitfactorantiqqg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

      LOG(DEBUG) << "H_qbar->qbarg : " << splitfactorantiqqg << endl;
      LOG(DEBUG) << "exp(-S) : " << sudtot << endl;
      LOG(DEBUG) << "H exp(-S) : " << sudtot*splitfactorantiqqg << endl;

      // left splitting inherits leftcolpartner from mother
      // and the right branch is the right color partner
      valueleft = make_splitting(leftbranch, leftcolpartner, rightbranch,input,Flavour::g,flavor,leftcolflav,Flavour::qbar,shower);
      // right splitting inherits rightcolpartner from mother
      // and the left branch is the left color partner
      valueright = make_splitting(rightbranch, leftbranch, empty,input,Flavour::qbar,flavor,Flavour::g,Flavour::noflav,shower);

      antiqqg=valueleft*valueright*splitfactorantiqqg*sudtot;
    } else if(flavor == Flavour::q) // b splits (which is left branch)
    {

      if ( ((square(tmpL.m())/tmpL.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5)) ||
         ((square(tmpR.m())/tmpR.perp()) >= lambda_mu*(square(tmptot.m())/tmptot.perp()*0.5)) )
        continue;

      double sudtot(sudakov(tmptot,tmpmleftcol,tmpmrightcol,tmpgrandm,flavor,granflavor, leftcolflav, rightcolflav,shower));
      splitfactorqqg = HQuark(tmptot,tmpmleftcol,tmpmrightcol,tmpL,tmpR,tmpgrandm,flavor,granflavor,leftcolflav,rightcolflav,shower);

#ifdef DEBUGCODE
      LOG(DEBUG) << "H_qbar->qbarg : " << splitfactorqqg << endl;
      LOG(DEBUG) << "exp(-S) : " << sudtot << endl;
      LOG(DEBUG) << "H exp(-S) : " << sudtot*splitfactorqqg << endl;
#endif

      valueleft=make_splitting(leftbranch, empty, rightbranch,input,Flavour::q,flavor,Flavour::noflav,Flavour::g,shower);
      valueright=make_splitting(rightbranch,leftbranch, rightcolpartner,input,Flavour::g,flavor,Flavour::q,rightcolflav,shower);

      qqg=valueleft*valueright*splitfactorqqg*sudtot;

#ifdef DEBUGCODE
      LOG(DEBUG) << " flavor: " << flavor << endl;
      LOG(DEBUG) << " sudakov: " << sudtot << endl;
      LOG(DEBUG) << " leftbranch: " << endl;
      printoutput(leftbranch);
      LOG(DEBUG) << " rightbranch: " << endl;
      printoutput(rightbranch);
      LOG(DEBUG) << "splitfactorqqg: " << splitfactorqqg << endl;
      LOG(DEBUG) << "splitfactorqqg*sudtot: " << splitfactorqqg*sudtot << endl;
      LOG(DEBUG) << "qqg: " << qqg << endl;
#endif

    } else
      LOG(ERROR) << "in make_splittin: WRONG FLAVOR: " << flavor << " -- serious error" << endl;

    valuetotal += gtotal+bbg+antibbg+qqg+antiqqg+ttg+antittg;

#ifdef DEBUGCODE      
    LOG(DEBUG) << "in make_splitting value total: " << valuetotal << endl;
    LOG(DEBUG) << "flavor: " << flavor << endl;
    LOG(DEBUG) << "gtotal: " << gtotal << endl;
    LOG(DEBUG) << "bbg: " << bbg << endl;
    LOG(DEBUG) << "antibbg: " << antibbg << endl;
    LOG(DEBUG) << "qqg: " << qqg << endl;
    LOG(DEBUG) << "antiqqg: " << antiqqg << endl;
    LOG(DEBUG) << "ttg: " << ttg << endl;
    LOG(DEBUG) << "antittg: " << antittg << endl;
    LOG(DEBUG) << "topdecayweight: " << topdecayweight << endl;
    LOG(DEBUG) << "antitopdecayweight: " << antitopdecayweight << endl;
      
    //if((flavor == Flavour::t || flavor == Flavour::tbar) && (level == 0)) 
    LOG(DEBUG) << "level: " << level << endl;
    LOG(DEBUG) << "splitting combination nr. " << i << " of " << nSets<< endl;
    LOG(DEBUG) << "flavor: " << flavor << endl;
    LOG(DEBUG) << "left branch: " << endl;
    printoutput(leftbranch);
    LOG(DEBUG) << "right branch: " << endl;
    printoutput(rightbranch);
    LOG(DEBUG) << "checksudtot: " << checksudtot << endl;
    LOG(DEBUG) << "checksplittot: " << checksplittot << endl;
    LOG(DEBUG) << "topdecayweight: " << topdecayweight << endl;
    LOG(DEBUG) << "antitopdecayweight: " << antitopdecayweight << endl;
    LOG(DEBUG) << "ttg: " << ttg << endl;

    LOG(DEBUG) << "antittg: " << antittg << endl;
    LOG(DEBUG) << "gtotal: " << gtotal << endl;
    LOG(DEBUG) << "bbg: " << bbg << endl;
    LOG(DEBUG) << "antibbg: " << antibbg << endl;
    LOG(DEBUG) << "qqg: " << qqg << endl;
    LOG(DEBUG) << "antiqqg: " << antiqqg << endl;
    LOG(DEBUG) << "valuetotal: " << valuetotal << endl;

#endif

  } // end loop over left branches
  
  // if((flavor == Flavour::t || flavor == Flavour::tbar) && (level == 0))
#ifdef DEBUGCODE
  LOG(DEBUG) << "level: " << level << endl;
  printoutput(input);
  LOG(DEBUG) << "tshower: " << Shower::t << endl;
  LOG(DEBUG) << "topdecay: " << topdecayweight << endl;
  LOG(DEBUG) << "radiation off top: " << valuetotal - topdecayweight << endl;
  LOG(DEBUG) << "topdecay + radiation: " << valuetotal << endl << endl;
  LOG(DEBUG) << "____________" << endl;
  LOG(DEBUG) << "____________" << endl;
#endif

  --Model::level;
  if(shower==Shower::W) --Model::Wshower_level;
  if(shower==Shower::b) --Model::btopshower_level;
  if(shower==Shower::t) --Model::tshower_level;

  m_calc.store(sk, valuetotal);

  return valuetotal;
}

double Deconstruction::Model::HQuark(const fastjet::PseudoJet &tmptot,
                              const fastjet::PseudoJet & tmpleftcol,
                              const fastjet::PseudoJet & tmprightcol,
                              const fastjet::PseudoJet & tmpL,
                              const fastjet::PseudoJet & tmpR,
                              const fastjet::PseudoJet & tmpgrandm,
                              const Flavour::id flavour,
                              const Flavour::id granflavour,
                              const Flavour::id tmpleftcolflav,
                              const Flavour::id tmprightcolflav,
                              const Shower::id shower) {
  double Hamiltonian = 0.0;

  double anglefac = 1.0;  // should be one here because if there is no col. connected partner we take it 1;

  //note for fermions, left=hard, right=soft
  //     antifemrions, left=soft, right=hard

  double gfunc = 0.0;

  double mu2 = square(tmptot.m());

  if (flavour == Flavour::t) {
    if (shower != Shower::t)
      throw NEW_EXCEPTION("flavor = Flavour::t but not tshower, cant be...");
    if (tmprightcol.e() > Cte::smalldouble) {
      anglefac = ((tmpL.squared_distance(tmpR) + square(param()[p_topmass])/square(tmpL.perp()))*(tmpL.squared_distance(tmprightcol) + square(param()[p_topmass])/square(tmpL.perp()))
                  - square(param()[p_topmass])/square(tmpL.perp())*tmpR.squared_distance(tmprightcol));
      anglefac *= 1.0/((tmpL.squared_distance(tmpR) + square(param()[p_topmass])/square(tmpL.perp()))*(tmpR.squared_distance(tmprightcol) + tmpR.squared_distance(tmpL) + square(param()[p_topmass])/square(tmpL.perp())));
    } else {
      anglefac = tmpL.squared_distance(tmpR) / (tmpL.squared_distance(tmpR)+(square(param()[p_topmass])/square(tmpL.perp())));
    }

    mu2 -= square(param()[p_topmass]); // top virtuality is p^2-mt^2

    if (mu2 < 0.0)
      throw NEW_EXCEPTION("ERROR: MSP top in HQuark, mu2 cannot be smaller 0 ");

    Hamiltonian = 8.0*M_PI*Cte::CF*alphas(mu2)/mu2*tmptot.perp()/tmpR.perp()*(1.0+square(tmpL.perp()/tmptot.perp()))*anglefac;
  } else if (flavour == Flavour::tbar) {
    if (shower != Shower::t)
      throw NEW_EXCEPTION("flavour = Flavour::tbar but not tshower, cant be...");

    if (tmpleftcol.e() > Cte::smalldouble) {
      anglefac = ((tmpR.squared_distance(tmpL) + square(param()[p_topmass])/square(tmpR.perp()))*(tmpR.squared_distance(tmpleftcol) + square(param()[p_topmass])/square(tmpR.perp()))
                 - square(param()[p_topmass])/square(tmpR.perp())*tmpL.squared_distance(tmpleftcol));
      anglefac *= 1.0/((tmpR.squared_distance(tmpL) + square(param()[p_topmass])/square(tmpR.perp()))*(tmpL.squared_distance(tmpleftcol) + tmpL.squared_distance(tmpR) + square(param()[p_topmass])/square(tmpR.perp())));
    } else {
      anglefac = tmpL.squared_distance(tmpR) / (tmpL.squared_distance(tmpR)+(square(param()[p_topmass])/square(tmpR.perp())));
    }

    mu2 -= square(param()[p_topmass]); // top virtuality is p^2-mt^2
    if (mu2 < 0.0)
      throw NEW_EXCEPTION("ERROR: MSP antitop in HQuark, mu2 cannot be smaller 0 ");
    Hamiltonian = 8.0*M_PI*Cte::CF*alphas(mu2)/mu2*tmptot.perp()/tmpL.perp()*(1.0+square(tmpR.perp()/tmptot.perp()))*anglefac;

  } else if (flavour == Flavour::q || flavour == Flavour::b) //assume massless Flavour::bs, thus Flavour::q = Flavour::b for kinematic
  {

    if ((tmprightcolflav == Flavour::t) || (tmprightcolflav == Flavour::tbar)) {
      if (tmprightcol.e() < Cte::smalldouble)
        throw NEW_EXCEPTION("HQuark: color conn. right flavor is top but no energy?? cant be ... ");

      double numerator = (tmprightcol.squared_distance(tmpR)+square(param()[p_topmass])/square(tmprightcol.perp()))*
                         (tmprightcol.squared_distance(tmpL)+square(param()[p_topmass])/square(tmprightcol.perp()))
                         - (square(param()[p_topmass])/square(tmprightcol.perp()))*tmpR.squared_distance(tmpL);
      double denom = 0.0;
      if (shower == Shower::t) {
        denom = (tmpR.squared_distance(tmprightcol) + square(param()[p_topmass])/square(tmprightcol.perp()))*(tmpR.squared_distance(tmpL)+tmpR.squared_distance(tmprightcol) + square(param()[p_topmass])/square(tmprightcol.perp()));
      } else if (shower == Shower::b) {
        denom = square(tmprightcol.squared_distance(tmpR)+square(param()[p_topmass])/square(tmprightcol.perp()));
      } else
        throw NEW_EXCEPTION("ERROR in HQuark for quark: if top color connected has to be tshower or btopshower!! ");

      anglefac = numerator/denom;
    } else {
      if (tmprightcolflav == Flavour::noflav) {
        anglefac = 1.0;
      } else {
        anglefac = tmprightcol.squared_distance(tmpL)/(tmpR.squared_distance(tmpL) + tmpR.squared_distance(tmprightcol));
      }
    }

    Hamiltonian = 8.0*M_PI*Cte::CF*alphas(mu2)/mu2*tmptot.perp()/tmpR.perp()*(1.0+square(tmpL.perp()/tmptot.perp()))*anglefac;
  } else if (flavour == Flavour::qbar || flavour == Flavour::bbar) {
    //note for fermions, left=hard, right=soft
    //     antifemrions, left=soft, right=hard

    if ((tmpleftcolflav == Flavour::t) || (tmpleftcolflav == Flavour::tbar)) {
      if (tmpleftcol.e() < Cte::smalldouble)
        throw NEW_EXCEPTION("ERROR in HQuark: color conn. left flavor is top but no energy?? cant be ... ");

      double numerator = (tmpleftcol.squared_distance(tmpL)+square(param()[p_topmass])/square(tmpleftcol.perp()))*
                         (tmpleftcol.squared_distance(tmpR)+square(param()[p_topmass])/square(tmpleftcol.perp()))
                       - ((square(param()[p_topmass])/square(tmpleftcol.perp()))*tmpL.squared_distance(tmpR));

      double denom = 0;

      if (shower == Shower::t) {
        denom = ((tmpL.squared_distance(tmpleftcol) + square(param()[p_topmass])/square(tmpleftcol.perp()))*(tmpR.squared_distance(tmpL)+tmpL.squared_distance(tmpleftcol) + square(param()[p_topmass])/square(tmpleftcol.perp())));
      } else if (shower == Shower::b) {
        denom = square(tmpleftcol.squared_distance(tmpL)+square(param()[p_topmass])/square(tmpleftcol.perp()));
      } else
        throw NEW_EXCEPTION("ERROR in HQuark for antiquark: if antitop color connected has to be tshower or btopshower!! ");

      anglefac=numerator/denom;
    } else {
      if (tmpleftcolflav == Flavour::noflav) {
        anglefac = 1.0;
      } else {
        anglefac = tmpleftcol.squared_distance(tmpR)/(tmpR.squared_distance(tmpL) + tmpL.squared_distance(tmpleftcol));
      }
    }

    Hamiltonian = 8.0*M_PI*Cte::CF*alphas(mu2)/mu2*tmptot.perp()/tmpL.perp()*(1.0+square(tmpR.perp()/tmptot.perp()))*anglefac;
  }

  if (Hamiltonian < 0.0)
    throw NEW_EXCEPTION("Error - HQuark cannot be smaller 0.");

  return Hamiltonian;
}

double Deconstruction::Model::HGluonQQ(const fastjet::PseudoJet & tmptot,
                                const fastjet::PseudoJet & tmpmleftcol,
                                const fastjet::PseudoJet & tmpmrightcol,
                                const fastjet::PseudoJet & tmpL,
                                const fastjet::PseudoJet & tmpR,
                                const fastjet::PseudoJet & tmpgrandm,
                                const Flavour::id flavour,
                                const Flavour::id granflavour,
                                const Flavour::id tmpleftcolflav,
                                const Flavour::id tmprightcolflav,
                                const Shower::id shower) {
  double Hgqq = 8.0*M_PI*Cte::TR*alphas(square(tmptot.m()))*(square(tmpL.perp()) + square(tmpR.perp()))/square(tmptot.perp())/square(tmptot.m());

  return Hgqq;
}

double Deconstruction::Model::HGluonBB(const fastjet::PseudoJet & tmptot, 
				 const fastjet::PseudoJet & tmpmleftcol, 
				 const fastjet::PseudoJet & tmpmrightcol,
				 const fastjet::PseudoJet & tmpL,
				 const fastjet::PseudoJet & tmpR, 
				 const fastjet::PseudoJet & tmpgrandm,
				 const Flavour::id flavor,
				 const Flavour::id granflavor,
				 const Flavour::id tmpleftcolflav,
				 const Flavour::id tmprightcolflav,
				 const Shower::id shower)
{
  
  double mu2(square(tmptot.m())); // should be (q_b+q_bbar)^2 with q_b^2=q_bbar^2=mb^2 if bmass is considered
  //////// Hgbb //////////////////////////////
  double Hgbb(8.0*M_PI*Cte::TR*alphas(mu2)*(square(tmpL.perp()) + square(tmpR.perp()))/square(tmptot.perp())/mu2);

  if(Hgbb < 0.0)
    LOG(ERROR) << "MSP: Error - Hgbb cannot be smaller 0" << endl;

  return Hgbb;
}


double Deconstruction::Model::HGluonGG(const fastjet::PseudoJet & tmptot, 
				    const fastjet::PseudoJet & tmpmleftcol, 
				    const fastjet::PseudoJet & tmpmrightcol,
				    const fastjet::PseudoJet & tmpL,
				    const fastjet::PseudoJet & tmpR, 
				    const fastjet::PseudoJet & tmpgrandm,
				    const Flavour::id flavor,
				    const Flavour::id granflavor,
				    const Flavour::id tmpleftcolflav,
				    const Flavour::id tmprightcolflav,
				    const Shower::id shower)
{

  double Hggg(0.0);
  
  
  //////////////   Hggg  /////////////////////////////////

  double anglefac(1.0);

  fastjet::PseudoJet soft;
  fastjet::PseudoJet hard;
  fastjet::PseudoJet colcon;
  Flavour::id colconflav;

  if(tmpL.perp() < tmpR.perp())
    {
      soft = tmpL;
      hard = tmpR;
      colcon = tmpmleftcol;
      colconflav = tmpleftcolflav;
    }
  else
    {
      soft = tmpR;
      hard = tmpL;
      colcon = tmpmrightcol;
      colconflav = tmprightcolflav;
    }
    

  if(colconflav == Flavour::noflav)
    {
      anglefac = 1.0;
    }
  else
    {
    
      if(shower==Shower::t)
	{
	  if((colconflav == Flavour::t) || (colconflav == Flavour::tbar))
	    {
	      anglefac = ((soft.squared_distance(colcon)+square(param()[p_topmass])/square(colcon.perp()))*(hard.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp())) - (square(param()[p_topmass])/square(colcon.perp()))*soft.squared_distance(hard));
	      
	      anglefac *= 1/((soft.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp())) * (soft.squared_distance(hard) + soft.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp())));	  
	    }
	  else
	    {
	      anglefac = hard.squared_distance(colcon)/(soft.squared_distance(hard)+soft.squared_distance(colcon));
	    }
	  
	}
      else if(shower == Shower::b)
	{
	  if((colconflav == Flavour::t) || (colconflav == Flavour::tbar))
	    {

	      anglefac = ((soft.squared_distance(colcon)+square(param()[p_topmass])/square(colcon.perp()))*(hard.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp())) - (square(param()[p_topmass])/square(colcon.perp()))*soft.squared_distance(hard));
	      
	      anglefac *= 1/square(soft.squared_distance(colcon) + square(param()[p_topmass])/square(colcon.perp()));
	    }
	  else
	    {
	      anglefac = hard.squared_distance(colcon)/(soft.squared_distance(hard)+soft.squared_distance(colcon));
	    }
	}
      else
	{
	  if((colconflav == Flavour::t) || (colconflav == Flavour::tbar)) LOG(ERROR) << "colcon = top in wrong shower : " << shower << endl;
	  
	  anglefac = hard.squared_distance(colcon)/(soft.squared_distance(hard)+soft.squared_distance(colcon));
	}

    } // end if no color partner else
  
  Hggg = 8.0*M_PI*Cte::CA*alphas(square(tmptot.m()))*square(tmptot.perp())/(square(tmptot.m())*soft.perp()*hard.perp()) * square(1.0 - soft.perp()*hard.perp()/square(tmptot.perp())) * anglefac;

  if(Hggg < 0.0) LOG(ERROR) << "MSP: Error - Hggg cannot be small 0" << endl;
  
  return Hggg;
}  




double Deconstruction::Model::sudakovGluon(const fastjet::PseudoJet & tmptot, 
					const fastjet::PseudoJet & tmpmleftcol, 
					const fastjet::PseudoJet & tmpmrightcol,
					const fastjet::PseudoJet & tmpgrandm,
					const Flavour::id flavor,
					const Flavour::id granflavor,
					const Flavour::id tmpleftcolflav,
					const Flavour::id tmprightcolflav,
					const Shower::id shower)
{
  
  double mothermass2(square(tmptot.m()));  
  double kappaK(0.0);
  double mu2max(0.0);
  //following line checks for existence of grandmother
  if(tmpgrandm.e() < Cte::smalldouble)
    {
      kappaK = 2.0*square(tmptot.perp())/tmptot.perp();
    }
  else 
    {
      if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	{
	  kappaK = (square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp();
	}
      else
	{
	  kappaK = square(tmpgrandm.m())/tmpgrandm.perp();
	}
    }
  

  double Sggg(0.0), Sgqq(0.0);
  double deltaSggg(0.0);

  double colconAmass(0.0);
  double colconBmass(0.0);


  //////// Sggg //////////////////////////////////////////////////////////////////////
  double theta_kA(0.0);
  if(tmpmleftcol.e()<Cte::smalldouble)
    {
      theta_kA = sqrt(dR2_outside);
    }
  else
    {
      if(tmpleftcolflav == Flavour::t || tmpleftcolflav == Flavour::tbar) colconAmass = param()[p_topmass]; 
      theta_kA=sqrt(tmptot.squared_distance(tmpmleftcol) + square(colconAmass)/square(tmpmleftcol.perp()));
    }
  
  double theta_kB(0.0);
  if(tmpmrightcol.e()<Cte::smalldouble)
    {
      theta_kB = sqrt(dR2_outside);
    }
  else 
    {
      if(tmprightcolflav == Flavour::t || tmprightcolflav == Flavour::tbar) colconBmass = param()[p_topmass];
      theta_kB=sqrt(tmptot.squared_distance(tmpmrightcol) + square(colconBmass)/square(tmpmrightcol.perp()));
    }


  //// BparaL, BparaR, AparaL, AparaR

  double BparaL(0.0), BparaR(0.0), AparaL(0.0), AparaR(0.0);
  
  
  double fbargg((log(2.0)-11.0/12.0)*Cte::CA);
  double f0gg(Cte::CA);
  

  // sudakov for massive and massless color connected parnters the same except B-> B' and A->A'
  // the primed ones are for massive partners.


  //MSP: where do we need kappa?


     double SILmin(0.0);
     double SILmax(0.0);
     double SIRmin(0.0);
     double SIRmax(0.0);

     // left radiation:
     //mu2star= AparaL*square(tmptot.perp())*exp(BparaL+fbargg/f0gg);
     

     // right radiation:


  
  // if grandmass2 =0 ->error
 
  Sggg = Cte::CA/M_PI/square(Cte::b0)*(log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*kappaK/2.0))*(1.0/alphas(theta_kA*theta_kB*square(tmptot.perp()))-11.0*Cte::b0/12.0) + 1.0/alphas(std::fabs(mothermass2)) - 1.0/alphas(tmptot.perp()*kappaK/2.0));
  
  
  if(Sggg<0) Sggg =0;
  
  if(shower == Shower::b)
    {
      if( tmpleftcolflav == Flavour::t || tmpleftcolflav == Flavour::tbar)
	{
	  deltaSggg = Cte::CA/2.0/M_PI/Cte::b0 *(log(square(theta_kA)/(square(param()[p_topmass])/square(tmpmleftcol.perp()))) + 1.0)*log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*kappaK/2.0));
	}
      else if(tmprightcolflav == Flavour::t || tmprightcolflav == Flavour::tbar)
	{
	  deltaSggg = Cte::CA/2.0/M_PI/Cte::b0 *(log(square(theta_kB)/(square(param()[p_topmass])/square(tmpmrightcol.perp()))) + 1.0)*log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*kappaK/2.0));
	}
      
      /*
      if(deltaSggg<0.0)
	{
	  LOG(DEBUG) << "ERROR: MSP deltaSggg<0 cannot be! "<< endl;
	}
      */

      if(deltaSggg < 0.0) deltaSggg = 0;

      //      if(deltaSggg < 0.0) LOG(DEBUG) << "MSP: ERROR deltaSggg < 0. Cannot be!" << endl;

      Sggg = Sggg + deltaSggg ;
    }

  ///////// Sgqq //////////////////////////////////////////////////////////////////////
  Sgqq = Cte::TR/3.0/M_PI/Cte::b0*log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*kappaK/2.0));

  // flavor factor g-> all 5 quarks
  Sgqq = Cte::nf*Sgqq;  // Sgqq = Sgbb

  double Sg = Sggg + Sgqq;

  return(1.0);
  //  return(exp(-Sg));

}



double Deconstruction::Model::sudakovQuark(const fastjet::PseudoJet & tmptot, 
					const fastjet::PseudoJet & tmpmleftcol, 
					const fastjet::PseudoJet & tmpmrightcol,
					const fastjet::PseudoJet & tmpgrandm,
					const Flavour::id flavor,
					const Flavour::id granflavor,
					const Flavour::id tmpleftcolflav,
					const Flavour::id tmprightcolflav,
					const Shower::id shower)
{

  double kappaK(0);
  double virt2(square(tmptot.m()));
  if(flavor == Flavour::t || flavor == Flavour::tbar) virt2 -= square(param()[p_topmass]);

  //following line checks for existence of grandmother

  // use fix starting scale from top decay in W shower:
  if(shower == Shower::W)
    {
      if(Model::Wshower_level == Model::Winitlevel+1)
	{
	  if((tmpgrandm.e() < Cte::smalldouble)) LOG(ERROR) << "MSP: ERROR! Grandm shouldnt be empty in W shower" << endl;
	  kappaK = 2.0*square(param()[p_wmass])/tmpgrandm.perp(); // grandperp should be perp of W
	}
      else
	{
	  kappaK=square(tmpgrandm.m())/tmpgrandm.perp();
	}
    }


  // use fix starting scale from top decay in b shower:  
  if(shower == Shower::b)
    {
      if(Model::btopshower_level == Model::binitlevel+1)
	{
	  if((tmpgrandm.e() < Cte::smalldouble)) LOG(ERROR) << "MSP: ERROR! Grandm shouldnt be empty in btop shower" << endl;
	  kappaK = 2.0*(square(param()[p_topmass])-square(param()[p_wmass]))/tmpgrandm.perp(); //tmpgrand is top here      
	}
      else
	{
	  kappaK=square(tmpgrandm.m())/tmpgrandm.perp();
	}
    }

  if(shower==Shower::QCD) {      
    if(tmpgrandm.e() < Cte::smalldouble) kappaK=2.0*square(tmptot.m())/tmptot.perp();
    else kappaK=square(tmpgrandm.m())/tmpgrandm.perp();}


  // check if topshower at the beginning:
  if(shower==Shower::t)
    { 
      
      if(Model::tshower_level == Model::tinitlevel+1)
	{
	  kappaK = 2.0*(square(tmptot.perp())+square(param()[p_topmass]))/tmptot.perp();
	}
      else
	{
	  // MSP Following needs to be changed... minus param()[p_topmass]^2 is wrong!!
	  if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	    {
	      kappaK = (square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp();
	    }
	  else if(granflavor == Flavour::noflav)
	    { 
	      LOG(ERROR) << "in sudakov Quark: find no flav for GrandM ERROR" << endl;
	    }
	  else kappaK = square(tmpgrandm.m())/tmpgrandm.perp();
	}

    }
  
  double Sqqg(0.0);


  if(flavor == Flavour::t)
    {
      if(shower != Shower::t) LOG(ERROR) << "flavor Flavour::t but shower not tshower, cant be..." << endl;
      double drtotcol2(dR2_outside);
      if(tmpmrightcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else drtotcol2=tmptot.squared_distance(tmpmrightcol);
      
     
      // msp - check here the mass in the alphas term
      Sqqg = Cte::CF/M_PI/Cte::b0*log(alphas(std::fabs(virt2))/alphas(tmptot.perp()*kappaK/2.0))*log((Cte::alphattg*drtotcol2*square(tmptot.perp())+2.0*Cte::betattg*square(param()[p_topmass]))/square(param()[p_topmass]));
      
      
      if(Sqqg <0) Sqqg=0;

      double Sdecay(0.0);
      
      if(std::fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
	{
	  if(kappaK/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	    {
	      Sdecay = log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/topwidth));
	    }
	  else Sdecay = log(atan(delta_topmass/topwidth));
	  
	  Sdecay -= log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth));
	}

            
      if(Sdecay < 0)
	{
	  LOG(ERROR) << "MSP: ERROR Sdecay for top cannot be < 0 " << endl;
	  LOG(DEBUG) << "splitting sudakov: " << Sqqg << endl;
	  LOG(DEBUG) << "kappaK: " << kappaK << endl;
	  LOG(DEBUG) << "topmass: " << param()[p_topmass] << endl;
	  LOG(DEBUG) << "topwidth: " << topwidth << endl;
	  LOG(DEBUG) << "delta_topmass: " << delta_topmass << endl;
	  LOG(DEBUG) << "jetmass: " << tmptot.m() << endl;
	  LOG(DEBUG) << "delta_top: " << delta_topmass << endl;
	  LOG(DEBUG) << "pT jet: " << tmptot.perp() << endl;
	  LOG(DEBUG) << "log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/delta_topmass)): " << log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/delta_topmass)) << endl;
	  LOG(DEBUG) << "log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth)) " << log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth)) << endl;
	}
	

    
      Sqqg += Sdecay;

    }
  else if(flavor == Flavour::tbar)
    {
      
      if(shower != Shower::t) LOG(ERROR) << "flavor Flavour::tbar but shower not tshower, cant be..." << endl;
      double drtotcol2(dR2_outside);
      if(tmpmleftcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else drtotcol2=tmptot.squared_distance(tmpmleftcol);

      // msp - check here the mass in the alphas term
      Sqqg = Cte::CF/M_PI/Cte::b0*log(alphas(std::fabs(virt2))/alphas(tmptot.perp()*kappaK/2.0))*log((drtotcol2*square(tmptot.perp())+2*square(param()[p_topmass]))/square(param()[p_topmass]));


      if(Sqqg <0) Sqqg=0;



      double Sdecay(0.0);
      
      if(std::fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
	{
	  if(kappaK/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	    {
	      Sdecay = log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/delta_topmass));
	    }
	  else Sdecay = log(atan(delta_topmass/topwidth));
	  
	  Sdecay -= log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth));
	}
      
      if(Sdecay < 0.0) LOG(ERROR) << "MSP: ERROR Sdecay cannot be < 0 " << endl;
      
      Sqqg += Sdecay;


    }
  else if(flavor == Flavour::q || flavor == Flavour::b)
    {
      
      double drtotcol2(dR2_outside);
      if(tmpmrightcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else 
	{
	  if(tmprightcolflav == Flavour::t ||tmprightcolflav == Flavour::tbar)
	    {
	      drtotcol2=tmptot.squared_distance(tmpmrightcol) + square(param()[p_topmass])/square(tmpmrightcol.perp());
	    }
	  else
	    {
	      drtotcol2=tmptot.squared_distance(tmpmrightcol);
	    }
	}

      //LOG(DEBUG) << "drtotcol2 : " << drtotcol2 << endl;
      Sqqg = Cte::CF/M_PI/square(Cte::b0)*(log(alphas(std::fabs(virt2))/alphas(tmptot.perp()*kappaK/2.0))*
			       (1.0/alphas(drtotcol2*square(tmptot.perp())) - 3.0*Cte::b0/4.0) + 1.0/alphas(std::fabs(virt2))-1.0/alphas(tmptot.perp()*kappaK/2.0));

      if(Sqqg <0) Sqqg=0;

      if((shower == Shower::b) && (tmprightcolflav == Flavour::t ||tmprightcolflav == Flavour::tbar))
	{
	  double deltaSqqg = Cte::CF/M_PI/Cte::b0*(log((tmpmrightcol.squared_distance(tmptot)+square(param()[p_topmass])/square(tmpmrightcol.perp()))/(square(param()[p_topmass])/square(tmpmrightcol.perp()))) + 1.0);
	  deltaSqqg *= log(alphas(std::fabs(virt2))/(alphas(tmptot.perp()*kappaK/2.0)));

	  /*
	  if(deltaSqqg < 0.0)
	    {
	      LOG(DEBUG) << "MSP: ERROR deltaSqqg < 0: " << deltaSqqg << " Cannot be!" << endl;
	      LOG(DEBUG) << "flavor: "  << flavor  << endl;
	      LOG(DEBUG) << "granflavor: " << granflavor << endl;
	      LOG(DEBUG) << "tmpmleftcolflav: " << tmpleftcolflav << endl;
	      LOG(DEBUG) << "tmpmrightcolflav: " << tmprightcolflav << endl;
	      LOG(DEBUG) << "shower: " << shower << endl;
	      LOG(DEBUG) << "tmptot: " << endl;
	      printjet(tmptot);
	      LOG(DEBUG) << "tmpmleftcol: " << endl;
	      printjet(tmpmleftcol);
	      LOG(DEBUG) << "tmprightcol: " << endl;
	      printjet(tmpmrightcol);
	      LOG(DEBUG) << "kappaK : "<< kappaK << endl;
	    }
	  */
	  if(deltaSqqg < 0.0) deltaSqqg = 0.0;

	  Sqqg = Sqqg + deltaSqqg;
	}

    }
  else if(flavor == Flavour::qbar || flavor == Flavour::bbar)
    {

      double drtotcol2(dR2_outside);
      if(tmpmleftcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else 
	{
	  if(tmpleftcolflav == Flavour::t ||tmpleftcolflav == Flavour::tbar)
	    {
	      drtotcol2=tmptot.squared_distance(tmpmleftcol) + square(param()[p_topmass])/square(tmpmleftcol.perp());
	    }
	  else
	    {
	      drtotcol2=tmptot.squared_distance(tmpmleftcol);
	    }
	}

      //LOG(DEBUG) << "drtotcol2 : " << drtotcol2 << endl;
      Sqqg = Cte::CF/M_PI/square(Cte::b0)*(log(alphas(std::fabs(virt2))/alphas(tmptot.perp()*kappaK/2.0))*
			       (1.0/alphas(drtotcol2*square(tmptot.perp())) - 3.0*Cte::b0/4.0) + 1.0/alphas(std::fabs(virt2))-1.0/alphas(tmptot.perp()*kappaK/2.0));

      if(Sqqg <0) Sqqg=0;


      if((shower == Shower::b) && (tmpleftcolflav == Flavour::t ||tmpleftcolflav == Flavour::tbar))
	{
	  double deltaSqqg = Cte::CF/M_PI/Cte::b0*(log((tmpmleftcol.squared_distance(tmptot)+square(param()[p_topmass])/square(tmpmleftcol.perp()))/(square(param()[p_topmass])/square(tmpmleftcol.perp()))) + 1.0);
	  deltaSqqg *= log(alphas(square(tmptot.m()))/(alphas(tmptot.perp()*kappaK/2.0)));
	  
	  if(deltaSqqg< 0.0) deltaSqqg = 0.0;
	  Sqqg = Sqqg + deltaSqqg;
	}

    }



  //  return(exp(-Sqqg));
  return(1.0);

}




double Deconstruction::Model::sudakovTopEND(const fastjet::PseudoJet & tmptot, 
				      const fastjet::PseudoJet & tmpmleftcol, 
				      const fastjet::PseudoJet & tmpmrightcol,
				      const fastjet::PseudoJet & tmpgrandm,
				      const Flavour::id flavor,
				      const Flavour::id granflavor,
				      const Flavour::id tmpleftcolflav,
				      const Flavour::id tmprightcolflav,
				      const Shower::id shower)
{

  double grandmass2(square(tmpgrandm.m())-square(param()[p_topmass]));
  double grandperp(tmpgrandm.perp());
  //  double mothermass2(topmass*topwidth);
  double mothermass2(square(tmptot.m()));

  //following line checks for existence of grandmother
  if(tmpgrandm.e() < Cte::smalldouble)
    {
      grandmass2 = 2.0*(square(tmptot.perp())+square(param()[p_topmass]));
      grandperp = tmptot.perp();
    }

  double Sqqg(0.0);

  double drtotcol2(dR2_outside);
  if(flavor == Flavour::t)
    {
      if(tmpmrightcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else drtotcol2=tmptot.squared_distance(tmpmrightcol);
    }
  else if(flavor == Flavour::tbar)
    {
      if(tmpmleftcol.e()<Cte::smalldouble)
	{
	  drtotcol2 = dR2_outside;
	}
      else drtotcol2=tmptot.squared_distance(tmpmleftcol);
     
    }
  else LOG(ERROR) << "in sudakovTopEND: wrong flavor " << endl;

  Sqqg = Cte::CF/M_PI/Cte::b0*log(alphas(std::fabs(mothermass2))/alphas(tmptot.perp()*grandmass2/2.0/grandperp)) * log((drtotcol2*square(tmptot.perp())+2.0*square(param()[p_topmass]))/square(param()[p_topmass]));
  
  if(Sqqg <0) Sqqg=0;

  double Sdecay(0.0);
  
  if(std::fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
    {
      if(grandmass2/grandperp/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	{
	  Sdecay = log(atan(grandmass2/grandperp/2.0*tmptot.perp()/param()[p_topmass]/delta_topmass));
	}
      else Sdecay = log(atan(delta_topmass/topwidth));
      
      Sdecay -= log(atan(std::fabs(square(tmptot.m())-square(param()[p_topmass]))/param()[p_topmass]/topwidth));
    }

  if(Sdecay < 0.0) LOG(ERROR) << "MSP: Error - Sdecay cannot be smaller 0" << endl;

  //  return(exp(-Sqqg+Sdecay));
  return(1.0);
}




double Deconstruction::Model::sudakov(const fastjet::PseudoJet & tmptot, 
				   const fastjet::PseudoJet & tmpmleftcol, 
				   const fastjet::PseudoJet & tmpmrightcol,
				   const fastjet::PseudoJet & tmpgrandm,
				   const Flavour::id flavor,
				   const Flavour::id granflavor,
				   const Flavour::id tmpleftcolflav,
				   const Flavour::id tmprightcolflav,
				   const Shower::id shower)
{

  // calculate virtuality of J
  double mu(tmptot.m());
  if(flavor == Flavour::t || flavor == Flavour::tbar) {
    mu = std::sqrt(std::fabs(square(tmptot.m())-square(param()[p_topmass])));
  }

  // calculate mumax (maximum scale), depends if first splitting after decay,... (many cases)
  double mumax(0.0);
    
  if(shower == Shower::W)
    {
      if(Model::Wshower_level == Model::Winitlevel+1)
	{
	  if((tmpgrandm.e() < Cte::smalldouble)) LOG(ERROR) << "MSP: ERROR! Grandm shouldnt be empty in W shower" << endl;
	  mumax = sqrt(tmptot.perp()*square(param()[p_wmass])/tmpgrandm.perp()); // grandperp should be perp of W
	}
      else
	{
	  mumax=sqrt(tmptot.perp()*square(tmpgrandm.m())/(2.0*tmpgrandm.perp()));
	}
    }
  else if(shower == Shower::b)
    {
      if(Model::btopshower_level == Model::binitlevel+1)
	{
	  if((tmpgrandm.e() < Cte::smalldouble)) LOG(ERROR) << "MSP: ERROR! Grandm shouldnt be empty in btop shower" << endl;
	  mumax = sqrt(tmptot.perp()*(square(param()[p_topmass])-square(param()[p_wmass]))/tmpgrandm.perp()); //tmpgrand is top here      
	}
      else
	{
	  mumax=sqrt(tmptot.perp()*square(tmpgrandm.m())/(2.0*tmpgrandm.perp()));
	}
    }
  else if(shower==Shower::QCD) 
    {      
      if(tmpgrandm.e() < Cte::smalldouble)
	{
	  mumax=tmptot.perp();
	}
      else
	{
	  mumax=sqrt(tmptot.perp()/2.0*square(tmpgrandm.m())/tmpgrandm.perp());
	}
    }
  else if(shower==Shower::t)
    { 
      
      if(Model::tshower_level == Model::tinitlevel+1)
	{
	  mumax = sqrt(tmptot.perp()*(square(tmptot.perp())+square(param()[p_topmass]))/tmptot.perp());
	}
      else
	{
	  if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	    {
	      mumax = sqrt(tmptot.perp()/2.0*(square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp());
	    }
	  else if(granflavor == Flavour::noflav)
	    { 
	      LOG(ERROR) << "in sudakov Quark: find no flav for GrandM ERROR" << endl;
	    }
	  else mumax = sqrt(tmptot.perp()/2.0*square(tmpgrandm.m())/tmpgrandm.perp());
	}

    }
  else LOG(ERROR) << "in sudakov(): Wrong shower" << shower << endl;
  
  double f0(0.0);
  double fbar(0.0);
  
  double SIL(0.0);
  double SIR(0.0);
  double SILmumax(0.0);
  double SIRmumax(0.0);
  double Sgqq(0.0);

  double stotal(0.0);

  double stopdecay(0.0);
  double stoprad(0.0);

  if(flavor == Flavour::tbar)
    {
      // antitop can radiate to the left
      // antitop can decay
      double theta2kJL(0.0);
 
      if(tmpmleftcol.e()<Cte::smalldouble)
	{
	  theta2kJL = dR2_outside;
	}
      else theta2kJL=tmptot.squared_distance(tmpmleftcol);
      
      double rhotildeq(square(param()[p_topmass])/square(tmptot.perp()));
      double x0(2.0-1.0*theta2kJL/(theta2kJL+rhotildeq) + 0.1 *log(rhotildeq/(theta2kJL+rhotildeq)));
      double x1(-3.0-1.0*theta2kJL/(theta2kJL+rhotildeq) + 0.4 *log(rhotildeq/(theta2kJL+rhotildeq)));

      double Iphiz(2.0*Cte::CF*((theta2kJL+2.0*rhotildeq)/(theta2kJL+rhotildeq) * log(2.0 + theta2kJL/rhotildeq) - 1.0));

      double tscale0((theta2kJL*square(tmptot.perp())+square(param()[p_topmass]))*exp(x0));
      double tscale1((theta2kJL*square(tmptot.perp())+square(param()[p_topmass]))*exp(x1));

      if(tscale0 <= square(mu))
	{
	  //	  LOG(DEBUG) << "1" << endl;
	  SIL = 0.0;
	}
      else if((tscale1 < square(mu)) &&
	      (tscale0 >= square(mu)))
	{
	  //	  LOG(DEBUG) << "2" << endl;
	  SIL = 1.0/alphas(square(mu)) - 1.0/(alphas(tscale0));
	  SIL += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mu))/(alphas(tscale0)));
	  SIL *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;
	}
      else if(tscale1 >= square(mu))
	{
	  //	  LOG(DEBUG) << "3" << endl;
	  SIL = 1.0/alphas(tscale1) - 1.0/(alphas(tscale0));
	  SIL += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(tscale1)/(alphas(tscale0)));
	  SIL *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;

	  SIL += 1/(2.0*M_PI*Cte::b0)*Iphiz*log(alphas(square(mu))/alphas(tscale1));
	}
      else LOG(ERROR) << "error in sudakov for antitop: Scales dont match" << endl;

      ////// now max part:

      if(tscale0 <= square(mumax))
	{
	  //	  LOG(DEBUG) << "1 max" << endl;
	  SILmumax = 0.0;
	}
      else if((tscale1 < square(mumax)) &&
	      (tscale0 >= square(mumax)))
	{
	  
	  //	  LOG(DEBUG) << "2 max" << endl;

	  SILmumax = 1.0/alphas(square(mumax)) - 1.0/(alphas(tscale0));
	  SILmumax += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mumax))/(alphas(tscale0)));
	  SILmumax *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;

	}
      else if(tscale1 >= square(mumax))
	{
	  //	  LOG(DEBUG) << "3 max" << endl;

	  SILmumax = 1.0/alphas(square(mumax)) - 1.0/(alphas(tscale0));
	  SILmumax += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mumax))/(alphas(tscale0)));
	  SILmumax *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;
	  SILmumax += 1/(2.0*M_PI*Cte::b0)*Iphiz*log(alphas(square(mumax))/alphas(tscale1));
	  
	}
      else LOG(ERROR) << "error in sudakov for antitop: Scales dont match" << endl;

      stoprad = SIL - SILmumax;
	  
      /// antitop decay
      double grandmass2(square(tmpgrandm.m())-square(param()[p_topmass]));
      double grandperp(tmpgrandm.perp());
      //  double mothermass2(topmass*topwidth);
      double mothermass2(square(tmptot.m()));
      

      
      double kappaK(0.0);
      
      if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	{
	  kappaK = (square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp();
	}
      else kappaK = square(tmpgrandm.m())/tmpgrandm.perp();


      if(fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
	{
	  if(kappaK/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	    {
	     stopdecay = log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/topwidth));
	    }
	  else stopdecay = log(atan(delta_topmass/topwidth));
	  
	  stopdecay -= std::log(std::atan(std::fabs(square(tmptot.m())-(param()[p_topmass]*param()[p_topmass]))/param()[p_topmass]/topwidth));
	}
      

      stotal = max(stoprad,0.0) + stopdecay;

    }
  else if(flavor == Flavour::t)
    {
       
      // top can radiate to the right
      // top can decay
      double theta2kJR(0.0);
 
      if(tmpmrightcol.e()<Cte::smalldouble)
	{
	  theta2kJR = dR2_outside;
	}
      else theta2kJR=tmptot.squared_distance(tmpmrightcol);
      
      
      double rhotildeq(square(param()[p_topmass])/square(tmptot.perp()));
      double x0(2.0-1.0*theta2kJR/(theta2kJR+rhotildeq) + 0.1 *log(rhotildeq/(theta2kJR+rhotildeq)));
      double x1(-3.0-1.0*theta2kJR/(theta2kJR+rhotildeq) + 0.4 *log(rhotildeq/(theta2kJR+rhotildeq)));

      double Iphiz(2.0*Cte::CF*((theta2kJR+2.0*rhotildeq)/(theta2kJR+rhotildeq) * log(2.0 + theta2kJR/rhotildeq) - 1.0));
  
      double tscale0((theta2kJR*square(tmptot.perp())+square(param()[p_topmass]))*exp(x0));
      double tscale1((theta2kJR*square(tmptot.perp())+square(param()[p_topmass]))*exp(x1));

      if( tscale0 <= square(mu))
	{
	  SIR=0.0;
	}      
      else if((tscale1 < square(mu)) &&
	      (tscale0 >= square(mu)) )
	{
	  SIR = 1.0/alphas(square(mu)) - 1.0/(alphas(tscale0));
	  SIR += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mu))/(alphas(tscale0)));
	  SIR *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;

	}
      else if(tscale1 >= square(mu))
	{
	  SIR = 1.0/alphas(tscale1) - 1.0/(alphas(tscale0));
	  SIR += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(tscale1)/(alphas(tscale0)));
	  SIR *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;
	  
	  SIR += 1/(2.0*M_PI*Cte::b0)*Iphiz*log(alphas(square(mu))/alphas(tscale1));

	}
      else LOG(ERROR) << "ERROR in sudakovtop: scales do not match " << endl;
      

      if( tscale0 <= square(mumax))
	{
	  SIRmumax=0.0;
	}      
      else if((tscale1 < square(mumax)) &&
	      (tscale0 >= square(mumax)) )
	{
	  SIRmumax = 1.0/alphas(square(mumax)) - 1.0/(alphas(tscale0));
	  SIRmumax += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(square(mumax))/(alphas(tscale0)));
	  SIRmumax *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;

	}
      else if(tscale1 >= square(mumax))
	{
	  SIRmumax = 1.0/alphas(tscale1) - 1.0/(alphas(tscale0));
	  SIRmumax += (1/alphas(tscale0) + x0*Cte::b0)*log(alphas(tscale1)/(alphas(tscale0)));
	  SIRmumax *= 1/2.0/M_PI/(x0-x1)/square(Cte::b0)*Iphiz;
	  
	  SIRmumax += 1/(2.0*M_PI*Cte::b0)*Iphiz*log(alphas(square(mumax))/alphas(tscale1));
	}
      else LOG(ERROR) << "ERROR in sudakovtop: scales do not match " << endl;

      stoprad = SIR - SIRmumax;

      /// top decay:

      double grandmass2(square(tmpgrandm.m())-square(param()[p_topmass]));
      double grandperp(tmpgrandm.perp());
      //  double mothermass2(topmass*topwidth);
      double mothermass2(square(tmptot.m()));
      
      //following line checks for existence of grandmother

      double kappaK(0.0);
      
      if(granflavor == Flavour::t || granflavor == Flavour::tbar)
	{
	  kappaK = (square(tmpgrandm.m())-square(param()[p_topmass]))/tmpgrandm.perp();
	}
      else kappaK = square(tmpgrandm.m())/tmpgrandm.perp();


      if(fabs(square(tmptot.m())-square(param()[p_topmass]))<param()[p_topmass]*delta_topmass)
	{
	  if(kappaK/2.0*tmptot.perp() < param()[p_topmass]*delta_topmass)
	    {
	      stopdecay = log(atan(kappaK/2.0*tmptot.perp()/param()[p_topmass]/topwidth));
	    }
	  else stopdecay = std::log(std::atan(delta_topmass/topwidth));
	  
	  stopdecay -= std::log(std::atan(std::fabs(square(tmptot.m())-(param()[p_topmass]*param()[p_topmass]))/param()[p_topmass]/topwidth));
	}

      stotal = max(stoprad,0.0) + stopdecay;

    }
  else if(flavor == Flavour::q || flavor == Flavour::qbar || flavor == Flavour::b || flavor == Flavour::bbar || flavor == Flavour::g)
    {
      
      if(flavor == Flavour::g)
	{
	  f0 = Cte::CA;
	  fbar = (log(2.0)-11.0/12.0)*Cte::CA;
	}
      else if(flavor == Flavour::q || flavor == Flavour::qbar || flavor == Flavour::b || flavor == Flavour::bbar)
	{
	  f0 = 2.0*Cte::CF;
	  fbar = -3.0/2.0*Cte::CF;
	}
      else LOG(ERROR) << "in sudakov: Wrong flavor: " << flavor << endl;

      
      double rhoL(0.0);
      double theta2kJL(0.0);
      double rhoR(0.0);
      double theta2kJR(0.0);

      double AparaR(0.0);
      double AparaL(0.0);
      double BparaR(0.0);
      double BparaL(0.0);
      /// left color connected partner
      if(tmpleftcolflav == Flavour::t || tmpleftcolflav == Flavour::tbar) // massive
	{
	  rhoL = square(param()[p_topmass])/square(tmpmleftcol.perp());

	  if(tmpmleftcol.e()<Cte::smalldouble)
	    {
	      theta2kJL = dR2_outside;
	    }
	  else theta2kJL=tmptot.squared_distance(tmpmleftcol);

	  // q-t dipole (btopshower) or not
	  if(shower == Shower::b) // a q-t dipole
	    {
	      AparaL=square(rhoL+theta2kJL)/rhoL;
	      BparaL=-1;
	    }
	  else
	    {
	      AparaL=square(rhoL+theta2kJL)/(2*rhoL+theta2kJL);
	      BparaL=rhoL/(rhoL+theta2kJL)*log(rhoL/(2*rhoL+theta2kJL));
	    }
	}
      else ///// massless
	{
	  rhoL = 0.0;
	  
	  if(tmpmleftcol.e()<Cte::smalldouble)
	    {
	      theta2kJL = dR2_outside;
	    }
	  else theta2kJL=tmptot.squared_distance(tmpmleftcol);
	  
	  AparaL=theta2kJL;
	  BparaL=0.0;
	  
	}

      // right color connected partner
      if(tmprightcolflav == Flavour::t || tmprightcolflav == Flavour::tbar) // massive
	{
	  rhoR = square(param()[p_topmass])/square(tmpmrightcol.perp());

	  
	  if(tmpmrightcol.e()<Cte::smalldouble)
	    {
	      theta2kJR = dR2_outside;
	    }
	  else theta2kJR=tmptot.squared_distance(tmpmrightcol);

	  // q-t dipole (btopshower) or not
	  if(shower == Shower::b) // a q-t dipole
	    {
	      AparaR=square(rhoR+theta2kJR)/rhoR;
	      BparaR=-1;
	    }
	  else
	    {
	      AparaR=square(rhoR+theta2kJR)/(2*rhoR+theta2kJR);
	      BparaR=rhoR/(rhoR+theta2kJR)*log(rhoR/(2*rhoR+theta2kJR));
	    }
	}
      else  /// massless
	{
	  rhoR = 0.0;

	  if(tmpmrightcol.e()<Cte::smalldouble)
	    {
	      theta2kJR = dR2_outside;
	    }
	  else theta2kJR=tmptot.squared_distance(tmpmrightcol);
	  
	  AparaR=theta2kJR;
	  BparaR=0.0;
	}

      double mustarL(AparaL*square(tmptot.perp())*exp(BparaL+fbar/f0));
      double mustarR(AparaR*square(tmptot.perp())*exp(BparaR+fbar/f0));

      if((flavor == Flavour::g || flavor == Flavour::qbar || flavor == Flavour::bbar) && (mustarL > square(mu)))
	{	  
	  SIL = f0/2.0/M_PI/square(Cte::b0)/alphas(AparaL*square(tmptot.perp()))+(fbar+BparaL*f0)/2.0/M_PI/Cte::b0;
	  SIL *= log(alphas(square(mu))/alphas(mustarL));
	  SIL -= f0/2.0/M_PI/square(Cte::b0)*(1/alphas(mustarL)-1/alphas(square(mu))); 
	}
	 
      if((flavor == Flavour::g || flavor == Flavour::q || flavor == Flavour::b) && (mustarR > square(mu)))
	{
	  SIR = f0/2.0/M_PI/square(Cte::b0)/alphas(AparaR*square(tmptot.perp()))+(fbar+BparaR*f0)/2.0/M_PI/Cte::b0;
	  SIR *= log(alphas(square(mu))/alphas(mustarR));
	  SIR -= f0/2.0/M_PI/square(Cte::b0)*(1/alphas(mustarR)-1/alphas(square(mu)));
	}

      if((flavor == Flavour::g || flavor == Flavour::qbar || flavor == Flavour::bbar)  && (mustarL > square(mumax)))
	{
	  SILmumax = f0/2.0/M_PI/square(Cte::b0)/alphas(AparaL*square(tmptot.perp()))+(fbar+BparaL*f0)/2.0/M_PI/Cte::b0;
	  SILmumax *= log(alphas(square(mumax))/alphas(mustarL));
	  SILmumax -= f0/2.0/M_PI/square(Cte::b0)*(1/alphas(mustarL)-1/alphas(square(mumax))); 
	}
	 
      if((flavor == Flavour::g || flavor == Flavour::q || flavor == Flavour::b)  && (mustarR > square(mumax)))
	{
	  SIRmumax = f0/2.0/M_PI/square(Cte::b0)/alphas(AparaR*square(tmptot.perp()))+(fbar+BparaR*f0)/2.0/M_PI/Cte::b0;
	  SIRmumax *= log(alphas(square(mumax))/alphas(mustarR));
	  SIRmumax -= f0/2.0/M_PI/square(Cte::b0)*(1/alphas(mustarR)-1/alphas(square(mumax)));
	}

      if(flavor == Flavour::g)
	{
	  
	  Sgqq = Cte::TR/3.0/M_PI/Cte::b0*log(alphas(square(mu))/alphas(square(mumax)));
	  
	  // flavor factor g-> all 5 quarks
	  Sgqq = Cte::nf*Sgqq;  // Sgqq = Sgbb

	  stotal = max(SIL-SILmumax + SIR-SIRmumax,0.0) + Sgqq;
	}

      if(flavor == Flavour::q || flavor == Flavour::b) stotal = max(SIR-SIRmumax,0.0);

      if(flavor == Flavour::qbar || flavor == Flavour::bbar) stotal = max(SIL-SILmumax,0.0);
    }
  else 
    {
      LOG(ERROR) << "In sudakov: Wrong flavor!!" << endl;
    } 

  // Question to Dave:
  // Eq. 47 and Eq. 82 are both supposed to be good for massive color connected partner. The only difference is that Eq. 82 is valid in t-shower and Eq. 47 is valid in the b-shower??

  return(exp(-stotal));
}


double Deconstruction::Model::M2(std::vector<fastjet::PseudoJet> &transformed) {
  return 1.0;
}

#include "Parameters.h"

#include <string>
#include <fstream>
#include <sstream>

#include "Exception.h"

#include <algorithm>

Deconstruction::Parameters::Parameters() {
}

Deconstruction::Parameters::Parameters(const std::string &input) {
  read(input);
}

Deconstruction::Parameters::~Parameters() {
}

int Deconstruction::Parameters::read(const std::string &input) {
  int nReadKeys = 0;

  std::ifstream infile(input.c_str());
  if (!infile.is_open()) {
    throw NEW_EXCEPTION("Failed to open [file] with parameters.").setParam("file", input);
  }
  while (infile.good()) {
    std::stringstream ss;
    infile.get(*ss.rdbuf(), '\n');
    if (infile.get() == EOF)
      break;

    std::string key;
    double value;
    ss >> key >> value;
    (*this)[key] = value;
  }
}

int Deconstruction::Parameters::insert(const std::string &key, const double value) {
  m_param.push_back(value);
  m_keys.push_back(key);
  return m_param.size() - 1;
}

double &Deconstruction::Parameters::operator[](const std::string &key) {
  std::vector<std::string>::iterator it = std::find(m_keys.begin(), m_keys.end(), key);
  if (it == m_keys.end()) {
    return (*this)[insert(key, 0)];
  }
  int index = (int) (it - m_keys.begin());
  return (*this)[index];
}

double &Deconstruction::Parameters::operator[](const int key) {
  return m_param[key];
}

#include "ParseUtils.h"

#include <getopt.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <vector>
#include <algorithm>

void buildOptions(struct extendedOption *extOpt, struct option *&opt, std::string &shortOpts) {
  int numberOfOpts = 0;
  if (extOpt == 0) {
    opt = 0;
    return;
  }
  struct extendedOption *copyOfExtOpt = extOpt;
  while (copyOfExtOpt->name != 0) {
    ++copyOfExtOpt;
    ++numberOfOpts;
  }
  ++numberOfOpts; // for the zero-entry
  std::vector<int> check;
  opt = new option[numberOfOpts];
  std::stringstream ss;
  ss.str("");
  for (int i = 0; i < numberOfOpts; ++i) {
    if (std::find(check.begin(), check.end(), extOpt[i].val) != check.end()) {
      std::cout << "buildOptions: WARNING -- Repeated option identifier (identifier \'" << (char) extOpt[i].val << "\' [" << (int) extOpt[i].val << "])." << std::endl; 
    }

    if (extOpt[i].has_arg == no_argument) { // set a flag
    } else if (extOpt[i].has_arg == required_argument) { // have arguments -- pass them on shortOpts
      ss << (char) extOpt[i].val << ":";
      check.push_back(extOpt[i].val);
    }
    opt[i].name = extOpt[i].name;
    opt[i].has_arg = extOpt[i].has_arg;
    opt[i].flag = extOpt[i].flag;
    opt[i].val = extOpt[i].val;
  }
  shortOpts = ss.str();
}

void dumpHelp(const std::string &nameOfProgram, struct extendedOption *extOpt, const std::string &description) {
  std::cout << description << std::endl;
  std::cout << "Usage: " << nameOfProgram << " [options]" << std::endl;
  std::cout << std::endl;
  std::cout << "where [options] can be a combination of:" << std::endl;
  while (extOpt->name != 0) {
    if (extOpt->has_arg == no_argument) {
      std::cout << "  --" << std::left << std::setw(38) << extOpt->name << "    " << extOpt->description << std::endl;
    } else if (extOpt->has_arg == required_argument) {
      std::cout << "  --" << std::left << std::setw(25) << extOpt->name << "|-" << (char) extOpt->val << "   (arg)      " << extOpt->description << "  (default: ";
      if (extOpt->type == extendedOption::eOTFloat) {
        float *f = (float *) extOpt->pointerToValue;
        std::cout << *f << " [float]";
      } else if (extOpt->type == extendedOption::eOTString) {
        std::string *s = (std::string *) extOpt->pointerToValue;
        std::cout << "\"" << *s << "\" [string]";
      } else if (extOpt->type == extendedOption::eOTInt) {
        int *k = (int *) extOpt->pointerToValue;
        std::cout << *k << " [int]";
      } else {
        std::cout << "unknown";
      }
      std::cout << ")" << std::endl;
    }
    ++extOpt;
  }
  std::cout << std::endl;
  std::cout << "Choose wisely :)" << std::endl;
  std::cout << std::endl;
}

void dumpOptions(struct extendedOption *extOpt) {
  while (extOpt->name != 0) {
    std::cout << std::left << std::setw(25) << extOpt->name << " = ";
    if (extOpt->type == extendedOption::eOTFloat) {
      float *f = (float *) extOpt->pointerToValue;
      std::cout << *f << " [float]";
    } else if (extOpt->type == extendedOption::eOTString) {
      std::string *s = (std::string *) extOpt->pointerToValue;
      std::cout << "\"" << *s << "\" [string]";
    } else if (extOpt->type == extendedOption::eOTInt) {
      int *k = (int *) extOpt->pointerToValue;
      std::cout << *k << " [int]";
    } else {
      std::cout << "unknown";
    }
    std::cout << std::endl;
    ++extOpt;
  }
  std::cout << std::endl;
}

bool parseArguments(int argc, char **argv, struct extendedOption *extOpt) {
  struct option *long_options;
  std::string shortOpts = "";
  buildOptions(extOpt, long_options, shortOpts);

  int c;
  while (true) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, shortOpts.c_str(), long_options, &option_index);
    if (c == -1) break; // end of parsing arguments

    // treat unrecognised option
    if (c == '?') return false;

    if (c == 0) { // flag setting option -- they are set automatically
    } else {
      struct extendedOption *eO = extOpt;
      while (eO->name != 0) { // look for its identifier in the array
        if (c == eO->val) { // if I found the identifier
          //std::cout << "Found ID " << (char) c << " for " << eO->name << " with contents \"" << optarg << "\", type == " << (int) eO->type << std::endl;
          // convert to the appropriate type and copy to the contents of the pointer storing it
          if (eO->type == extendedOption::eOTFloat) {
            std::stringstream ss;
            float *f = (float *) eO->pointerToValue;
            ss << optarg;
            ss >> *f;
          } else if (eO->type == extendedOption::eOTString) {
            std::stringstream ss;
            std::string s;
            ss << optarg;
            ss >> s;
            std::string *sp = (std::string *) eO->pointerToValue;
            *sp = s;
          } else if (eO->type == extendedOption::eOTInt) {
            std::stringstream ss;
            int *k = (int *) eO->pointerToValue;
            ss << optarg;
            ss >> *k;
          }
          break;
        } // if I have found the identifier
        ++eO;
      } // end of loop for finding the identifier of this option in the extended options array
    } // if it is NOT a flag setting option
  } // until they are all parsed

  delete [] long_options;

  return true;
}

#include "Storage.h"
#include "Exception.h"

Deconstruction::Storage::Storage() {
}

Deconstruction::Storage::~Storage() {
}

std::vector<fastjet::PseudoJet> Deconstruction::Storage::get(const std::vector<int> &list) const {
  std::vector<fastjet::PseudoJet> out;
  for (int i = 0; i < list.size(); ++i) {
    if (list[i] >= m_input.size())
      throw NEW_EXCEPTION("Trying to access element out of bound in Storage.");

    out.push_back(m_input[list[i]]);
  }
  return out;
}

fastjet::PseudoJet Deconstruction::Storage::sum(const std::vector<int> &list) const {
  fastjet::PseudoJet out(0,0,0,0);
  for (int i = 0; i < list.size(); ++i) {
    if (list[i] >= m_input.size())
      throw NEW_EXCEPTION("Trying to access element out of bound in Storage.");
    out += m_input[list[i]];
  }
  return out;
}

int Deconstruction::Storage::size() const {
  return m_input.size();
}

void Deconstruction::Storage::set(const std::vector<fastjet::PseudoJet> &input) {
  m_input = input;
}

std::vector<fastjet::PseudoJet> &Deconstruction::Storage::get() {
  return m_input;
}

const std::vector<fastjet::PseudoJet> &Deconstruction::Storage::get() const {
  return m_input;
}

fastjet::PseudoJet &Deconstruction::Storage::operator[](int i) {
  return m_input[i];
}

const fastjet::PseudoJet &Deconstruction::Storage::operator[](int i) const {
  return m_input[i];
}


#include "StoredCalculations.h"
#include "StoredJet.h"

#include <algorithm>
#include <cmath>

using namespace Deconstruction;

Deconstruction::StoredCalculations::StoredCalculations() {
}
Deconstruction::StoredCalculations::~StoredCalculations() {
}

bool Deconstruction::StoredCalculations::check(StoredKey k, double &w) {
  std::map<StoredKey, double>::iterator it = m_table.find(k);
  if (it != m_table.end()) {
    w = it->second;
    return true;
  }
  return false;
}

void Deconstruction::StoredCalculations::store(StoredKey k, double w) {
  m_table.insert(std::pair<StoredKey, double>(k, w));
}

void Deconstruction::StoredCalculations::clear() {
  m_table.clear();
}

#include "StoredJet.h"

#include "Exception.h"

Deconstruction::StoredJet::StoredJet()
  :m_store(0) {
}

Deconstruction::StoredJet::StoredJet(const Storage &store)
  :m_store(&store) {
}

Deconstruction::StoredJet::StoredJet(const Storage *store)
  :m_store(store) {
}

Deconstruction::StoredJet::StoredJet(const Storage &store, std::vector<int> list)
  :m_store(&store), m_list(list) {
}

Deconstruction::StoredJet::StoredJet(const StoredJet &s)
  :m_store(s.m_store), m_list(s.m_list) {
}

Deconstruction::StoredJet::~StoredJet() {
}

std::vector<int> Deconstruction::StoredJet::getList() const {
  return m_list;
}

const Deconstruction::Storage &Deconstruction::StoredJet::store() const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  return *m_store;
}

int Deconstruction::StoredJet::getIdx(int i) const {
  return m_list[i];
}

const fastjet::PseudoJet &Deconstruction::StoredJet::operator[](int i) const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  int j = m_list[i];
  return (*m_store)[j];
}

const fastjet::PseudoJet &Deconstruction::StoredJet::at(int i) const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  int j = m_list.at(i);
  return (*m_store)[j];
}

int Deconstruction::StoredJet::size() const {
  return m_list.size();
}

fastjet::PseudoJet Deconstruction::StoredJet::sum() const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  return m_store->sum(m_list);
}

void Deconstruction::StoredJet::push_back(int j) {
  m_list.push_back(j);
}

Deconstruction::StoredJet::operator std::vector<fastjet::PseudoJet>() const {
  if (!m_store) throw NEW_EXCEPTION("No Storage available!");
  return m_store->get(m_list);
}

#include "TopGluonModel.h"

#include <iostream>

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "Helper.h"

#include "Parameters.h"
#include "AnalysisParameters.h"

using namespace Deconstruction;
using namespace fastjet;
using namespace std;

using std::endl;

Deconstruction::TopGluonModel::TopGluonModel(Parameters &param)
  : Deconstruction::Model::Model(param) {
}

Deconstruction::TopGluonModel::~TopGluonModel() {
}

//double Deconstruction::TopGluonModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::TopGluonModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {

  double result = 0;

  if (sumJets.perp() <= sumJets.m())
    return 0;

#ifdef DEBUGCODE
  LOG(DEBUG) << "input of TopGluonModel: " << endl;
  printoutput(jets, DEBUG);
#endif

  //vector<PseudoJet> firstleftcol,firstrightcol,grandmother;
  //firstleftcol.clear();
  //firstrightcol.clear();
  //grandmother.clear();
  //result = make_splitting(jets, firstleftcol, firstrightcol, grandmother, Flavour::t, Flavour::noflav, Flavour::noflav, Flavour::noflav, Shower::t);

  result = start_splitting(jets, Flavour::t, Shower::t);
  m_calc.clear();

#ifdef DEBUGCODE
  LOG(DEBUG) << "input of TopGluonModel: " << endl;
  printoutput(jets, DEBUG);
  LOG(DEBUG) << "GluTopSplitweight: " << result << endl;
#endif

  return result;

}

double Deconstruction::TopGluonModel::hamiltonian(double pTsum) {
  double kH2 = square(pTsum);
  double resmass2 = square(param()[p_topmass]);
  double Hfj_sig = Cte::Npdf_signal*std::pow((Cte::pTmin2 + resmass2)/(kH2 + resmass2), Cte::Npdf_signal)/(kH2 + resmass2);
  return Hfj_sig;
}



#include "TopModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "Helper.h"
#include "JetInfo.h"

#include "Parameters.h"
#include "AnalysisParameters.h"

#include "Message.h"
#include "BModel.h"
#include "WModel.h"

#include <iostream>

#include <vector>

using namespace Deconstruction;
using namespace fastjet;
using namespace std;

Deconstruction::TopModel::TopModel(Parameters &param, Flavour::id flavour)
  : Deconstruction::Model::Model(param), m_flavour(flavour) {
}

Deconstruction::TopModel::~TopModel() {
}

//double Deconstruction::TopModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::TopModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {
  double result = 0;

#ifdef DEBUGCODE
  LOG(DEBUG) << "In TopModel: " << endl;
  printoutput(jets,"input in TopModel", DEBUG);
#endif

  /////////////////////////////////////////////////////////////
  //////// Make Top decay:   ////////////////////////////////

  // checks if input is already smaller than top window
  // msp: include mass check after cross check with dave
  if ( (jets.size()<3) /*|| (offTop.m() < (topmass-delta_topmass)) */) {
#ifdef DEBUGCODE
    LOG(DEBUG) << "nr microjets < 3 or too light -> no top decay" << endl;
#endif
    return result;
  }

#ifdef DEBUGCODE
  LOG(DEBUG) << "passed lower topmass cut: " << sumJets.m() << endl;
#endif

  /////// split into left and right branches /////////////////////////

  /*
  std::vector<std::vector<fastjet::PseudoJet> > bbranches;

  int start=1;
  int end=jets.size()-1;   // msp: check if input.size()-1 or not..
  power_set_orig<fastjet::PseudoJet,std::vector<fastjet::PseudoJet>::iterator >(jets.begin(), jets.end(), bbranches, start, end);

  for(unsigned i=0; i<bbranches.size(); i++) {
    vector<PseudoJet> Wbranch(jets.size()-bbranches[i].size());
    vector<PseudoJet> bbranch(bbranches[i]);

    sort(jets.begin(), jets.end(),  lessThanIndex);
    sort(bbranch.begin(), bbranch.end(),  lessThanIndex);
    set_difference(jets.begin(), jets.end(), bbranch.begin(), bbranch.end(), Wbranch.begin(), lessThanIndex);
  */

  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with p->mom().size() bits with any number of bits activated
  unsigned int iElements = (unsigned int) jets.size();
  // iElements cannot be greater than the maximum,
  // but this has already been tested before, when generating the powerset for ISR/SB
  // there is no point in testing it again, since these sets must be smaller or equal in size
  // than the subsets of the powerset of all microjets
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  // the bits position in index i give the positions of the vector to be used
  // we are skipping the empty set combination here
  // that's why i starts at one and ends at nSets-1
  for (unsigned long long i = 1; i < nSets; ++i) {

    // reject combinations with less than two partons for the W-boson earlier
    if (numberOfSetBits(i) < 2)
      continue;

    //std::vector<fastjet::PseudoJet> Wbranch;
    //std::vector<fastjet::PseudoJet> bbranch;
    StoredJet Wbranch(jets.store());
    StoredJet bbranch(jets.store());

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        Wbranch.push_back(jets.getIdx(k));
      } else {
        bbranch.push_back(jets.getIdx(k));
      }
    }

    PseudoJet topjet = sumJets;
    PseudoJet Wjet = Wbranch.sum();
    PseudoJet bjet = bbranch.sum();

    //msp: should be checked before branches are splitted
    if (std::fabs(square(topjet.m()) - square(param()[p_topmass])) >= param()[p_topmass]*delta_topmass)
      continue;

    if (std::fabs(square(Wjet.m()) - square(param()[p_wmass])) >= param()[p_wmass]*delta_wmass) continue;

    // already done previously
    //if (Wbranch.size() < 2)
    //  continue;

    // first decay factor for t -> bW

    double decayfactortop = 8.0*square(M_PI)*square(param()[p_topmass])*param()[p_topmass]*topwidth/std::atan(std::fabs(square(topjet.m()) - square(param()[p_topmass]))/(param()[p_topmass]*topwidth));

    decayfactortop /= (square(param()[p_topmass])-square(param()[p_wmass]))*(square(square(topjet.m())-square(param()[p_topmass]))+square(param()[p_topmass])*square(topwidth));

    double bweight = 0.0;
    double Wweight = 0.0;

#ifdef DEBUGCODE
    LOG(DEBUG) << "input b model: " << endl;
    printoutput(bbranch, DEBUG);

    LOG(DEBUG) << "input W model: " << endl;
    printoutput(Wbranch, DEBUG);
#endif

    if (m_flavour == Flavour::t) {
      BModel bm(param(), Flavour::b,Flavour::t);
      bm.setTopset(jets);
      bweight = bm.weight(bbranch, bjet); 
      WModel wm(param(), Flavour::Wp);
      wm.setTopset(jets);
      Wweight = wm.weight(Wbranch, Wjet);
    } else if (m_flavour == Flavour::tbar) {
      BModel bm(param(), Flavour::bbar,Flavour::tbar);
      bm.setTopset(jets);
      bweight = bm.weight(bbranch, bjet); 
      WModel wm(param(), Flavour::Wm);
      wm.setTopset(jets);
      Wweight = wm.weight(Wbranch, Wjet);
    } else {
      LOG(ERROR) << "TopModel: WRONG FLAVOR namely: " << m_flavour << endl;
    }

#ifdef DEBUGCODE
    LOG(DEBUG) << "tflavor: " << m_flavour << endl;
    LOG(DEBUG) << "decayfactortop: " << decayfactortop << endl;
    LOG(DEBUG) << "bweight: " << bweight << endl;
    LOG(DEBUG) << "W shower*decayfactorW: " << Wweight << endl;
    LOG(DEBUG) << "topweight: " << decayfactortop*Wweight*bweight << endl ;

    LOG(DEBUG) << "topweight: " << decayfactortop*Wweight*bweight << endl;
    LOG(DEBUG) << "topweight sum: " << result << endl << endl;
#endif

    result += decayfactortop*Wweight*bweight;
  }

#ifdef DEBUGCODE
  LOG(DEBUG) << "final signal weight in topmodel: " << result << endl << endl;
#endif

  return result;
}

double Deconstruction::TopModel::hamiltonian(double pTsum) {
  return 1.0;
}

#include "WModel.h"

#include <vector>
#include <fastjet/PseudoJet.hh>
#include "Helper.h"
#include "AnalysisParameters.h"

#include "Parameters.h"

#include "JetInfo.h"
#include <iostream>
#include <algorithm>

using namespace Deconstruction;
using namespace std;
using namespace fastjet;

Deconstruction::WModel::WModel(Parameters &param, Flavour::id flavour, bool includeAngleFactor)
  : Deconstruction::Model::Model(param), m_flavour(flavour), m_includeAngleFactor(includeAngleFactor) {
}

Deconstruction::WModel::~WModel() {
}

//void Deconstruction::WModel::setTopset(const std::vector<fastjet::PseudoJet> &topset) {
void Deconstruction::WModel::setTopset(const StoredJet &topset) {
  m_topset = topset.getList();
}

//double Deconstruction::WModel::weight(const std::vector<fastjet::PseudoJet> &jets, fastjet::PseudoJet &sumJets) {
double Deconstruction::WModel::weight(const StoredJet &jets, fastjet::PseudoJet &sumJets) {

  double signalweight = 0;

  StoredJet topset(jets.store(), m_topset);

  PseudoJet Wjet = sumJets;
  PseudoJet tjet = topset.sum();


  // possible combinations for the set in [first, last] are given by permutations of the binary
  // number with p->mom().size() bits with any number of bits activated
  unsigned int iElements = (unsigned int) jets.size();
  // iElements cannot be greater than the maximum,
  // but this has already been tested before, when generating the powerset for ISR/SB
  // there is no point in testing it again, since these sets must be smaller or equal in size
  // than the subsets of the powerset of all microjets
  unsigned long long nSets = (unsigned long long) (0x1 << iElements); // 2^iElements

  StoredJet empty(jets.store());

  // the bits position in index i give the positions of the vector to be used
  // we are skipping the empty set combination here
  // that's why i starts at one and ends at nSets-1
  for (unsigned long long i = 1; i < nSets; ++i) {

    //std::vector<fastjet::PseudoJet> ubranch;
    //std::vector<fastjet::PseudoJet> dbranch;
    StoredJet ubranch(jets.store());
    StoredJet dbranch(jets.store());

    // check the bits set in i and include the elements corresponding to
    // that bit position in result
    for (unsigned int k = 0; k < iElements; ++k) {
      if ( (i & (0x1 << k)) != 0) {
        ubranch.push_back(jets.getIdx(k));
      } else {
        dbranch.push_back(jets.getIdx(k));
      }
    }

    // quark goes left, antiquark goes right
    // if top hypothesis -> W+ > dbar u
    // grandmother is set to be empty because the scale where new shower starts is set by hand

    double valueleft = 0;
    double valueright = 0;

    if (m_flavour == Flavour::Wp) {
      valueleft = make_splitting(ubranch, empty, dbranch, jets, Flavour::q, Flavour::noflav, Flavour::noflav, Flavour::qbar, Shower::W);
      valueright = make_splitting(dbranch, ubranch, empty, jets, Flavour::qbar, Flavour::noflav, Flavour::q, Flavour::noflav,Shower::W);
    } else if (m_flavour == Flavour::Wm) {
      valueleft = make_splitting(dbranch, empty, ubranch, jets, Flavour::q, Flavour::noflav, Flavour::noflav, Flavour::qbar, Shower::W);
      valueright = make_splitting(ubranch, dbranch, empty, jets, Flavour::qbar, Flavour::noflav,Flavour::q, Flavour::noflav, Shower::W);
    } else LOG(ERROR) << "ERROR: wrong flavor in Wmodel: " << m_flavour << endl;

    PseudoJet djet = dbranch.sum();

    double anglefacW = 1;
    if (m_includeAngleFactor) {
      anglefacW = 12.0*scalprod(djet,tjet)*square_p1minusp2(tjet,djet);
      anglefacW /= (square(param()[p_topmass]) - square(param()[p_wmass])) * (square(param()[p_topmass]) + 2*square(param()[p_wmass]));
    }

    double Hwsplit = 8*square(M_PI)*param()[p_wmass]*wwidth/std::atan(std::fabs(square(Wjet.m())-square(param()[p_wmass]))/(param()[p_wmass]*wwidth));
    Hwsplit /= square(square(Wjet.m()) - square(param()[p_wmass])) + square(param()[p_wmass])*square(wwidth);

    ///////// sudakov
    double expSW = 0;
    double Wsud = 0;

    if(std::fabs(square(Wjet.m()) - square(param()[p_wmass])) < param()[p_wmass]*delta_wmass) {
      Wsud=(std::log(std::atan(delta_wmass/wwidth)) - std::log(std::atan(std::fabs(square(Wjet.m())-square(param()[p_wmass]))/(param()[p_wmass]*wwidth))));
      expSW = exp(-Wsud);
    } else
      expSW = 1;

    double decayfactorW = Hwsplit *anglefacW* expSW;

    signalweight += decayfactorW*valueleft*valueright;

#ifdef DEBUGCODE
    LOG(DEBUG) << "top input in Wmodel" << endl;
    printoutput(topset, DEBUG);
    LOG(DEBUG) << "input in Wmodel" << endl;
    printoutput(jets, DEBUG);
    LOG(DEBUG) << "ubranch" << endl;
    printoutput(ubranch, DEBUG);
    LOG(DEBUG) << "dbranch" << endl;
    printoutput(dbranch, DEBUG);
    LOG(DEBUG) << "input flavor: " << m_flavour << endl;
    LOG(DEBUG) << "anglefacW: " << anglefacW << endl;
    LOG(DEBUG) << "Hwsplit: " << Hwsplit << endl;
    LOG(DEBUG) << "Wsud: " << Wsud << endl;
    LOG(DEBUG) << "decayfactorW: " << decayfactorW << endl;
    LOG(DEBUG) << "value left: " << valueleft << endl;

    LOG(DEBUG) << "value right: " << valueright << endl;
    LOG(DEBUG) << "W shower: decayfactorW*valueleft*valueright :  " << decayfactorW*valueleft*valueright << endl;
#endif

  }

#ifdef DEBUGCODE
  LOG(DEBUG) << "W showre total signalweight: " << signalweight << endl;
#endif

  return param()[p_br_wqq]*signalweight;
}

double Deconstruction::WModel::hamiltonian(double pTsum) {
  return 1.0;
}

