#include "AxesFinder.hh"
#include "FWM.hh"
#include "HEPTopTagger.hh"
#include "LowPt.hh"
#include "MeasureFunction.hh"
#include "Njettiness.hh"
#include "NjettinessDefinition.hh"
#include "NjettinessPlugin.hh"
#include "Nsubjettiness.hh"
#include "QHTT.hh"
#include "WinnerTakeAllRecombiner.hh"
//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: AxesFinder.cc 670 2014-06-06 01:24:42Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Functions for minimization.
//
///////

// Given starting axes, update to find better axes by using Kmeans clustering around the old axes
template <int N>
std::vector<LightLikeAxis> AxesFinderFromOnePassMinimization::UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes,
                                  const std::vector <fastjet::PseudoJet> & inputJets) const {
   assert(old_axes.size() == N);

   // some storage, declared static to save allocation/re-allocation costs
   static LightLikeAxis new_axes[N];
   static fastjet::PseudoJet new_jets[N];
   for (int n = 0; n < N; ++n) {
      new_axes[n].reset(0.0,0.0,0.0,0.0);
      new_jets[n].reset_momentum(0.0,0.0,0.0,0.0);
   }

   double precision = _precision;

   /////////////// Assignment Step //////////////////////////////////////////////////////////
   std::vector<int> assignment_index(inputJets.size());
   int k_assign = -1;

   for (unsigned i = 0; i < inputJets.size(); i++){
      double smallestDist = std::numeric_limits<double>::max();  //large number
      for (int k = 0; k < N; k++) {
         double thisDist = old_axes[k].DistanceSq(inputJets[i]);
         if (thisDist < smallestDist) {
            smallestDist = thisDist;
            k_assign = k;
         }
      }
      if (smallestDist > sq(_Rcutoff)) {k_assign = -1;}
      assignment_index[i] = k_assign;
   }

   //////////////// Update Step /////////////////////////////////////////////////////////////
   double distPhi, old_dist;
   for (unsigned i = 0; i < inputJets.size(); i++) {
      int old_jet_i = assignment_index[i];
      if (old_jet_i == -1) {continue;}

      const fastjet::PseudoJet& inputJet_i = inputJets[i];
      LightLikeAxis& new_axis_i = new_axes[old_jet_i];
      double inputPhi_i = inputJet_i.phi();
      double inputRap_i = inputJet_i.rap();

      // optimize pow() call
      // add noise (the precision term) to make sure we don't divide by zero
      if (_beta == 1.0) {
         double DR = std::sqrt(sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i));
         old_dist = 1.0/DR;
      } else if (_beta == 2.0) {
         old_dist = 1.0;
      } else if (_beta == 0.0) {
         double DRSq = sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i);
         old_dist = 1.0/DRSq;
      } else {
         old_dist = sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i);
         old_dist = std::pow(old_dist, (0.5*_beta-1.0));
      }

      // TODO:  Put some of these addition functions into light-like axes
      // rapidity sum
      new_axis_i.set_rap(new_axis_i.rap() + inputJet_i.perp() * inputRap_i * old_dist);
      // phi sum
      distPhi = inputPhi_i - old_axes[old_jet_i].phi();
      if (fabs(distPhi) <= M_PI){
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * inputPhi_i * old_dist );
      } else if (distPhi > M_PI) {
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * (-2*M_PI + inputPhi_i) * old_dist );
      } else if (distPhi < -M_PI) {
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * (+2*M_PI + inputPhi_i) * old_dist );
      }
      // weights sum
      new_axis_i.set_weight( new_axis_i.weight() + inputJet_i.perp() * old_dist );
      // momentum magnitude sum
      new_jets[old_jet_i] += inputJet_i;
   }
   // normalize sums
   for (int k = 0; k < N; k++) {
      if (new_axes[k].weight() == 0) {
         // no particles were closest to this axis!  Return to old axis instead of (0,0,0,0)
         new_axes[k] = old_axes[k];
      } else {
         new_axes[k].set_rap( new_axes[k].rap() / new_axes[k].weight() );
         new_axes[k].set_phi( new_axes[k].phi() / new_axes[k].weight() );
         new_axes[k].set_phi( std::fmod(new_axes[k].phi() + 2*M_PI, 2*M_PI) );
         new_axes[k].set_mom( std::sqrt(new_jets[k].modp2()) );
      }
   }
   std::vector<LightLikeAxis> new_axes_vec(N);
   for (unsigned k = 0; k < N; ++k) new_axes_vec[k] = new_axes[k];
   return new_axes_vec;
}

// Given N starting axes, this function updates all axes to find N better axes.
// (This is just a wrapper for the templated version above.)
std::vector<LightLikeAxis> AxesFinderFromOnePassMinimization::UpdateAxes(const std::vector <LightLikeAxis> & old_axes,
                                      const std::vector <fastjet::PseudoJet> & inputJets) const {
   int N = old_axes.size();
   switch (N) {
      case 1: return UpdateAxesFast<1>(old_axes, inputJets);
      case 2: return UpdateAxesFast<2>(old_axes, inputJets);
      case 3: return UpdateAxesFast<3>(old_axes, inputJets);
      case 4: return UpdateAxesFast<4>(old_axes, inputJets);
      case 5: return UpdateAxesFast<5>(old_axes, inputJets);
      case 6: return UpdateAxesFast<6>(old_axes, inputJets);
      case 7: return UpdateAxesFast<7>(old_axes, inputJets);
      case 8: return UpdateAxesFast<8>(old_axes, inputJets);
      case 9: return UpdateAxesFast<9>(old_axes, inputJets);
      case 10: return UpdateAxesFast<10>(old_axes, inputJets);
      case 11: return UpdateAxesFast<11>(old_axes, inputJets);
      case 12: return UpdateAxesFast<12>(old_axes, inputJets);
      case 13: return UpdateAxesFast<13>(old_axes, inputJets);
      case 14: return UpdateAxesFast<14>(old_axes, inputJets);
      case 15: return UpdateAxesFast<15>(old_axes, inputJets);
      case 16: return UpdateAxesFast<16>(old_axes, inputJets);
      case 17: return UpdateAxesFast<17>(old_axes, inputJets);
      case 18: return UpdateAxesFast<18>(old_axes, inputJets);
      case 19: return UpdateAxesFast<19>(old_axes, inputJets);
      case 20: return UpdateAxesFast<20>(old_axes, inputJets);
      default: std::cout << "N-jettiness is hard-coded to only allow up to 20 jets!" << std::endl;
         return std::vector<LightLikeAxis>();
   }

}

// uses minimization of N-jettiness to continually update axes until convergence.
// The function returns the axes found at the (local) minimum
std::vector<fastjet::PseudoJet> AxesFinderFromOnePassMinimization::getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& seedAxes) const {

   // convert from PseudoJets to LightLikeAxes
   std::vector< LightLikeAxis > old_axes(n_jets, LightLikeAxis(0,0,0,0));
   for (int k = 0; k < n_jets; k++) {
      old_axes[k].set_rap( seedAxes[k].rap() );
      old_axes[k].set_phi( seedAxes[k].phi() );
   }

   // Find new axes by iterating (only one pass here)
   std::vector< LightLikeAxis > new_axes(n_jets, LightLikeAxis(0,0,0,0));
   double cmp = std::numeric_limits<double>::max();  //large number
   int h = 0;
   while (cmp > _precision && h < _halt) { // Keep updating axes until near-convergence or too many update steps
      cmp = 0.0;
      h++;
      new_axes = UpdateAxes(old_axes, inputJets); // Update axes
      for (int k = 0; k < n_jets; k++) {
         cmp += old_axes[k].Distance(new_axes[k]);
      }
      cmp = cmp / ((double) n_jets);
      old_axes = new_axes;
   }

   // Convert from internal LightLikeAxes to PseudoJet
   std::vector<fastjet::PseudoJet> outputAxes;
   for (int k = 0; k < n_jets; k++) {
      fastjet::PseudoJet temp = old_axes[k].ConvertToPseudoJet();
      outputAxes.push_back(temp);
   }

   // this is used to debug the minimization routine to make sure that it works.
   bool do_debug = false;
   if (do_debug) {
      // get this information to make sure that minimization is working properly
      TauComponents seed_tau_components = _measureFunction.result(inputJets, seedAxes);
      double seed_tau = seed_tau_components.tau();
      TauComponents tau_components = _measureFunction.result(inputJets, outputAxes);
      double outputTau = tau_components.tau();
      assert(outputTau <= seed_tau);
   }

   return outputAxes;
}

PseudoJet AxesFinderFromKmeansMinimization::jiggle(const PseudoJet& axis) const {
   double phi_noise = ((double)rand()/(double)RAND_MAX) * _noise_range * 2.0 - _noise_range;
   double rap_noise = ((double)rand()/(double)RAND_MAX) * _noise_range * 2.0 - _noise_range;

   double new_phi = axis.phi() + phi_noise;
   if (new_phi >= 2.0*M_PI) new_phi -= 2.0*M_PI;
   if (new_phi <= -2.0*M_PI) new_phi += 2.0*M_PI;

   PseudoJet newAxis(0,0,0,0);
   newAxis.reset_PtYPhiM(axis.perp(),axis.rap() + rap_noise,new_phi);
   return newAxis;
}


// Repeatedly calls the one pass finder to try to find global minimum
std::vector<fastjet::PseudoJet> AxesFinderFromKmeansMinimization::getAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& seedAxes) const {

   // first iteration
	std::vector<fastjet::PseudoJet> bestAxes = _onePassFinder.getAxes(n_jets, inputJets, seedAxes);

   double bestTau = (_measureFunction.result(inputJets,bestAxes)).tau();

   for (int l = 1; l < _n_iterations; l++) { // Do minimization procedure multiple times (l = 1 to start since first iteration is done already)

      // Add noise to current best axes
      std::vector< PseudoJet > noiseAxes(n_jets, PseudoJet(0,0,0,0));
      for (int k = 0; k < n_jets; k++) {
         noiseAxes[k] = jiggle(bestAxes[k]);
      }

      std::vector<fastjet::PseudoJet> testAxes = _onePassFinder.getAxes(n_jets, inputJets, noiseAxes);
      double testTau = (_measureFunction.result(inputJets,testAxes)).tau();

      if (testTau < bestTau) {
         bestTau = testTau;
         bestAxes = testAxes;
      }
   }

   return bestAxes;
}

// Uses minimization of the geometric distance in order to find the minimum axes.
// It continually updates until it reaches convergence or it reaches the maximum number of attempts.
// This is essentially the same as a stable cone finder.
std::vector<fastjet::PseudoJet> AxesFinderFromGeometricMinimization::getAxes(int /*n_jets*/, const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& currentAxes) const {

   std::vector<fastjet::PseudoJet> seedAxes = currentAxes;
   double seedTau = _function.tau(particles, seedAxes);

   for (int i = 0; i < _nAttempts; i++) {

      std::vector<fastjet::PseudoJet> newAxes(seedAxes.size(),fastjet::PseudoJet(0,0,0,0));

      // find closest axis and assign to that
      for (unsigned int i = 0; i < particles.size(); i++) {

         // start from unclustered beam measure
         int minJ = -1;
         double minDist = _function.beam_distance_squared(particles[i]);

         // which axis am I closest to?
         for (unsigned int j = 0; j < seedAxes.size(); j++) {
            double tempDist = _function.jet_distance_squared(particles[i],seedAxes[j]);
            if (tempDist < minDist) {
               minDist = tempDist;
               minJ = j;
            }
         }

         // if not unclustered, then cluster
         if (minJ != -1) newAxes[minJ] += particles[i];
      }

      // calculate tau on new axes
      seedAxes = newAxes;
      double tempTau = _function.tau(particles, newAxes);

      // close enough to stop?
      if (fabs(tempTau - seedTau) < _accuracy) break;
      seedTau = tempTau;
   }

   return seedAxes;
}

// Go from internal LightLikeAxis to PseudoJet
fastjet::PseudoJet LightLikeAxis::ConvertToPseudoJet() {
    double px, py, pz, E;
    E = _mom;
    pz = (std::exp(2.0*_rap) - 1.0) / (std::exp(2.0*_rap) + 1.0) * E;
    px = std::cos(_phi) * std::sqrt( std::pow(E,2) - std::pow(pz,2) );
    py = std::sin(_phi) * std::sqrt( std::pow(E,2) - std::pow(pz,2) );
    return fastjet::PseudoJet(px,py,pz,E);
}

} //namespace contrib

FASTJET_END_NAMESPACE

double perp(fastjet::PseudoJet v_pj, fastjet::PseudoJet ref_pj) {
  double pt = 0.;
  valarray<double> v = v_pj.four_mom();
  valarray<double> ref = ref_pj.four_mom();
  ref[3] = 0.;
  double mag2 = ref[0]*ref[0] + ref[1]*ref[1] + ref[2]*ref[2];
  double v_ref =  v[0]*ref[0] + v[1]*ref[1] + v[2]*ref[2];
  valarray<double> v_perp= v - (v_ref / mag2) * ref;
  pt = sqrt(v_perp[0]*v_perp[0] + v_perp[1]*v_perp[1] + v_perp[2]*v_perp[2]);
  return pt;
}

FWM::FWM() {};

FWM::FWM(vector<fastjet::PseudoJet> jets) : _jets(jets) {}

FWM::FWM(HEPTopTagger::HEPTopTagger htt, int selection) {
  PseudoJet top = htt.t();
  Unboost rf(top);
  PseudoJet a(-top.px(),-top.py(),-top.pz(),0.);

  vector<PseudoJet> jets;
  if(selection / 1000 == 1)  {
    jets.push_back(a);
  }
  if( (selection%1000)/100 == 1 ) {
    jets.push_back(rf(htt.b()));
  }
  if( (selection%100)/10 == 1 ) {
    jets.push_back(rf(htt.W1()));
  }
  if( (selection%10) == 1 ) {
    jets.push_back(rf(htt.W2()));
  }
  _jets=jets;
}

inline double FWM::ATan2(double y, double x){
  if (x != 0) return  atan2(y, x);
  if (y == 0) return  0;
  if (y >  0) return  M_PI/2;
  else        return -M_PI/2;
}

double FWM::Theta(fastjet::PseudoJet j) {
  return j.px() == 0.0 && j.py() == 0.0 && j.pz() == 0.0 ? 0.0 : ATan2(j.perp(),j.pz());
}

double FWM::cos_Omega(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2) {
  double cos_omega = cos(Theta(jet1)) * cos(Theta(jet2))
    + sin(Theta(jet1)) * sin(Theta(jet2)) * cos(jet1.phi_std() - jet2.phi_std());
  return cos_omega;
}

double FWM::U(unsigned order) {
  double H = 0.;
  double norm = (_jets.size() * _jets.size());
  for(unsigned ii = 0; ii < _jets.size(); ii++){
    for(unsigned jj = 0; jj < _jets.size(); jj++){
      double W = 1.;
      double cos_O;
      if(ii==jj) {
	cos_O=1.0;
      } else {
	cos_O=cos_Omega(_jets[ii], _jets[jj]);
      }
      H += W * legendre(order,cos_O);
    }
  }
  if (norm > 0.) H /= norm;
  return H;
}

double FWM::Pt(unsigned order) {
  fastjet::PseudoJet zaxis(0., 0., 1., 0.);
  return FWM::Pt(order, zaxis);
}

double FWM::Pt(unsigned order, fastjet::PseudoJet ref_pj) {
  double H = 0.;
  double norm = 0.;
  for(unsigned ii = 0; ii < _jets.size(); ii++){
    norm += perp(_jets[ii], ref_pj)*perp(_jets[ii], ref_pj);
    for(unsigned jj = 0; jj < _jets.size(); jj++){
      double W = perp(_jets[ii], ref_pj)*perp(_jets[jj], ref_pj);
      double cos_O;
      if(ii==jj) {
	cos_O=1.0;
      } else {
	cos_O=cos_Omega(_jets[ii], _jets[jj]);
      }
      H += W * legendre(order,cos_O);
    }
  }
  if (norm > 0.) H /= norm;
  return H;
}
namespace HEPTopTagger{
//optimal_R fit
double R_opt_calc_funct(double pt_filt) {
  return 327./pt_filt;
}

bool HEPTopTagger_fixed_R::_first_time = true;

void HEPTopTagger_fixed_R::print_banner() {
  if (!_first_time) {return;}
  _first_time = false;

  std::cout << "#--------------------------------------------------------------------------\n";
  std::cout << "#                           HEPTopTagger 2.0                               \n";
  std::cout << "#                                                                          \n";
  std::cout << "# Please cite JHEP 1506 (2015) 203 [arXiv:1503.05921 [hep-ph]]             \n";
  std::cout << "# and JHEP 1010 (2010) 078 [arXiv:1006.2833 [hep-ph]].                     \n";
  std::cout << "# This code is provided without warranty.                                  \n";
  std::cout << "#--------------------------------------------------------------------------\n";
}

//pt wrt a reference vector
double HEPTopTagger_fixed_R::perp(const fastjet::PseudoJet & vec, const fastjet::PseudoJet & ref) {
  double ref_ref = ref.px() * ref.px() + ref.py() * ref.py() + ref.pz() * ref.pz();
  double vec_ref = vec.px() * ref.px() + vec.py() * ref.py() + vec.pz() * ref.pz();
  double per_per = vec.px() * vec.px() + vec.py() * vec.py() + vec.pz() * vec.pz();
  if (ref_ref > 0.)
    per_per -= vec_ref * vec_ref / ref_ref;
  if (per_per < 0.)
    per_per = 0.;
  return sqrt(per_per);
}

//modified Jade distance
double HEPTopTagger_fixed_R::djademod (const fastjet::PseudoJet& subjet_i, const fastjet::PseudoJet& subjet_j, const fastjet::PseudoJet& ref) {
  double dj = -1.0;
  double delta_phi = subjet_i.delta_phi_to(subjet_j);
  double delta_eta = subjet_i.eta() - subjet_j.eta();
  double delta_R = sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
  dj = perp(subjet_i, ref) * perp(subjet_j, ref) * pow(delta_R, 4.);
  return dj;
}

//minimal |(m_ij / m_123) / (m_w/ m_t) - 1|
double HEPTopTagger_fixed_R::f_rec() {
  double m12 = (_top_subs[0] + _top_subs[1]).m();
  double m13 = (_top_subs[0] + _top_subs[2]).m();
  double m23 = (_top_subs[1] + _top_subs[2]).m();
  double m123 = (_top_subs[0] + _top_subs[1] + _top_subs[2]).m();

  double fw12 = fabs( (m12/m123) / (_mwmass/_mtmass) - 1);
  double fw13 = fabs( (m13/m123) / (_mwmass/_mtmass) - 1);
  double fw23 = fabs( (m23/m123) / (_mwmass/_mtmass) - 1);

  return std::min(fw12, std::min(fw13, fw23));
}

//find hard substructures
void HEPTopTagger_fixed_R::FindHardSubst(const PseudoJet & this_jet, std::vector<fastjet::PseudoJet> & t_parts) {
  PseudoJet parent1(0, 0, 0, 0), parent2(0, 0, 0, 0);
  if (this_jet.m() < _max_subjet_mass || !this_jet.validated_cs()->has_parents(this_jet, parent1, parent2)) {
    t_parts.push_back(this_jet);
  } else {
    if (parent1.m() < parent2.m())
      std::swap(parent1, parent2);
    FindHardSubst(parent1, t_parts);
    if (parent1.m() < _mass_drop_threshold * this_jet.m())
      FindHardSubst(parent2, t_parts);
  }
}

//store subjets as vector<PseudoJet> with [0]->b [1]->W-jet 1 [2]->W-jet 2
void HEPTopTagger_fixed_R::store_topsubjets(const std::vector<PseudoJet>& top_subs) {
  _top_subjets.resize(0);
  double m12 = (top_subs[0] + top_subs[1]).m();
  double m13 = (top_subs[0] + top_subs[2]).m();
  double m23 = (top_subs[1] + top_subs[2]).m();
  double dm12 = fabs(m12 - _mwmass);
  double dm13 = fabs(m13 - _mwmass);
  double dm23 = fabs(m23 - _mwmass);

  if (dm23 <= dm12 && dm23 <= dm13) {
    _top_subjets.push_back(top_subs[0]);
    _top_subjets.push_back(top_subs[1]);
    _top_subjets.push_back(top_subs[2]);
  } else if (dm13 <= dm12 && dm13 < dm23) {
    _top_subjets.push_back(top_subs[1]);
    _top_subjets.push_back(top_subs[0]);
    _top_subjets.push_back(top_subs[2]);
  } else if (dm12 < dm23 && dm12 < dm13) {
    _top_subjets.push_back(top_subs[2]);
    _top_subjets.push_back(top_subs[0]);
    _top_subjets.push_back(top_subs[1]);
  }
  _W = _top_subjets[1] + _top_subjets[2];
  return;
}

//check mass plane cuts
bool HEPTopTagger_fixed_R::check_mass_criteria(const std::vector<PseudoJet> & top_subs) const {
  bool is_passed = false;
  double m12 = (top_subs[0] + top_subs[1]).m();
  double m13 = (top_subs[0] + top_subs[2]).m();
  double m23 = (top_subs[1] + top_subs[2]).m();
  double m123 = (top_subs[0] + top_subs[1] + top_subs[2]).m();
  if (
      (atan(m13/m12) > _m13cutmin && _m13cutmax > atan(m13/m12)
       && (m23/m123 > _rmin && _rmax > m23/m123))
      ||
      (((m23/m123) * (m23/m123) < 1 - _rmin * _rmin* (1 + (m13/m12) * (m13/m12)))
       &&
       ((m23/m123) * (m23/m123) > 1 - _rmax * _rmax * (1 + (m13/m12) * (m13/m12)))
       &&
       (m23/m123 > _m23cut))
      ||
      (((m23/m123) * (m23/m123) < 1 - _rmin * _rmin * (1 + (m12/m13) * (m12/m13)))
       &&
       ((m23/m123) * (m23/m123) > 1 - _rmax * _rmax * (1 + (m12/m13) * (m12/m13)))
       &&
       (m23/m123 > _m23cut))
      ) {
    is_passed = true;
  }
  return is_passed;
}

double HEPTopTagger_fixed_R::nsub(fastjet::PseudoJet jet, int order, fastjet::contrib::Njettiness::AxesMode axes, double beta, double R0) {
  fastjet::contrib::Nsubjettiness nsub(order, axes, beta, R0);
  return nsub.result(jet);
}

HEPTopTagger_fixed_R::HEPTopTagger_fixed_R() : _do_qjets(0),
					       _mass_drop_threshold(0.8), _max_subjet_mass(30.),
					       _mode(Mode(0)), _mtmass(172.3), _mwmass(80.4), _mtmin(150.), _mtmax(200.), _rmin(0.85*80.4/172.3), _rmax(1.15*80.4/172.3),
					       _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), _minpt_tag(200.),
					       _nfilt(5), _Rfilt(0.3), _jet_algorithm_filter(fastjet::cambridge_algorithm), _minpt_subjet(0.),
					       _jet_algorithm_recluster(fastjet::cambridge_algorithm),
					       _zcut(0.1), _rcut_factor(0.5),
					       _q_zcut(0.1), _q_dcut_fctr(0.5), _q_exp_min(0.), _q_exp_max(0.), _q_rigidity(0.1), _q_truncation_fctr(0.0),
					       _debug(false)
{
  _djsum = 0.;
  _delta_top = 1000000000000.0;
  _pruned_mass = 0.;
  _unfiltered_mass = 0.;
  _top_candidate.reset(0., 0., 0., 0.);
  _parts_size = 0;
  _is_maybe_top = _is_masscut_passed = _is_ptmincut_passed = false;
  _top_subs.clear();
  _top_subjets.clear();
  _top_hadrons.clear();
  _top_parts.clear();
  _qweight= -1.;

}

HEPTopTagger_fixed_R::HEPTopTagger_fixed_R(const fastjet::PseudoJet jet) : _do_qjets(0),
									   _jet(jet), _initial_jet(jet),
									   _mass_drop_threshold(0.8), _max_subjet_mass(30.),
									   _mode(Mode(0)), _mtmass(172.3), _mwmass(80.4),  _mtmin(150.), _mtmax(200.), _rmin(0.85*80.4/172.3), _rmax(1.15*80.4/172.3),
									   _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), _minpt_tag(200.),
									   _nfilt(5), _Rfilt(0.3), _jet_algorithm_filter(fastjet::cambridge_algorithm), _minpt_subjet(0.),
									   _jet_algorithm_recluster(fastjet::cambridge_algorithm),
									   _zcut(0.1), _rcut_factor(0.5),
									   _q_zcut(0.1), _q_dcut_fctr(0.5), _q_exp_min(0.), _q_exp_max(0.), _q_rigidity(0.1), _q_truncation_fctr(0.0),
									   _fat(jet),
									   _debug(false)
{}

HEPTopTagger_fixed_R::HEPTopTagger_fixed_R(const fastjet::PseudoJet jet,
					   double mtmass, double mwmass) : _do_qjets(0),
									   _jet(jet), _initial_jet(jet),
									   _mass_drop_threshold(0.8), _max_subjet_mass(30.),
									   _mode(Mode(0)), _mtmass(mtmass), _mwmass(mwmass), _rmin(0.85*80.4/172.3), _rmax(1.15*80.4/172.3),
									   _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), _minpt_tag(200.),
									   _nfilt(5), _Rfilt(0.3), _jet_algorithm_filter(fastjet::cambridge_algorithm), _minpt_subjet(0.),
									   _jet_algorithm_recluster(fastjet::cambridge_algorithm),
									   _zcut(0.1), _rcut_factor(0.5),
									   _q_zcut(0.1), _q_dcut_fctr(0.5), _q_exp_min(0.), _q_exp_max(0.), _q_rigidity(0.1), _q_truncation_fctr(0.0),
									   _fat(jet),
									   _debug(false)
{}

void HEPTopTagger_fixed_R::run() {
  print_banner();

  if ((_mode != EARLY_MASSRATIO_SORT_MASS)
      && (_mode != LATE_MASSRATIO_SORT_MASS)
      && (_mode != EARLY_MASSRATIO_SORT_MODDJADE)
      && (_mode != LATE_MASSRATIO_SORT_MODDJADE)
      && (_mode != TWO_STEP_FILTER) ) {
    std::cout << "ERROR: UNKNOWN MODE" << std::endl;
    return;
  }

  //Qjets
  QjetsPlugin _qjet_plugin(_q_zcut, _q_dcut_fctr, _q_exp_min, _q_exp_max, _q_rigidity, _q_truncation_fctr);
  _qjet_def = fastjet::JetDefinition(&_qjet_plugin);
  _qweight=-1;
  vector<fastjet::PseudoJet> _q_constits;
  ClusterSequence* _qjet_seq;
  PseudoJet _qjet;
  if (_do_qjets){
    _q_constits = _initial_jet.associated_cluster_sequence()->constituents(_initial_jet);
    _qjet_seq = new ClusterSequence(_q_constits, _qjet_def);
    _qjet = sorted_by_pt(_qjet_seq->inclusive_jets())[0];
    _qjet_seq->delete_self_when_unused();
    const QjetsBaseExtras* ext =
      dynamic_cast<const QjetsBaseExtras*>(_qjet_seq->extras());
    _qweight=ext->weight();
    _jet = _qjet;
    _fat = _qjet;
  }

  //initialization
  _djsum = 0.;
  _delta_top = 1000000000000.0;
  _pruned_mass = 0.;
  _unfiltered_mass = 0.;
  _top_candidate.reset(0., 0., 0., 0.);
  _parts_size = 0;
  _is_maybe_top = _is_masscut_passed = _is_ptmincut_passed = false;
  _top_subs.clear();
  _top_subjets.clear();
  _top_hadrons.clear();
  _top_parts.clear();

  //find hard substructures
  FindHardSubst(_jet, _top_parts);

  if (_top_parts.size() < 3) {
    if (_debug) {std::cout << "< 3 hard substructures " << std::endl;}
    return; //such events are not interesting
  }

  // Sort subjets-after-unclustering by pT.
  // Necessary so that two-step-filtering can use the leading-three.
  _top_parts=sorted_by_pt(_top_parts);

  // loop over triples
  _top_parts = sorted_by_pt(_top_parts);
  for (unsigned rr = 0; rr < _top_parts.size(); rr++) {
    for (unsigned ll = rr + 1; ll < _top_parts.size(); ll++) {
      for (unsigned kk = ll + 1; kk < _top_parts.size(); kk++) {

	// two-step filtering
	// This means that we only look at the triplet formed by the
	// three leading-in-pT subjets-after-unclustering.
	if((_mode==TWO_STEP_FILTER) && rr>0)
	  continue;
	if((_mode==TWO_STEP_FILTER) && ll>1)
	  continue;
	if((_mode==TWO_STEP_FILTER) && kk>2)
	  continue;

      	//pick triple
	PseudoJet triple = join(_top_parts[rr], _top_parts[ll], _top_parts[kk]);

	//filtering
	double filt_top_R
	  = std::min(_Rfilt, 0.5*sqrt(std::min(_top_parts[kk].squared_distance(_top_parts[ll]),
					       std::min(_top_parts[rr].squared_distance(_top_parts[ll]),
							_top_parts[kk].squared_distance(_top_parts[rr])))));
	JetDefinition filtering_def(_jet_algorithm_filter, filt_top_R);
	fastjet::Filter filter(filtering_def, fastjet::SelectorNHardest(_nfilt) * fastjet::SelectorPtMin(_minpt_subjet));
	PseudoJet topcandidate = filter(triple);

	//mass window cut
  	if (topcandidate.m() < _mtmin || _mtmax < topcandidate.m()) continue;

	// Sanity cut: can't recluster less than 3 objects into three subjets
	if (topcandidate.pieces().size() < 3)
	  continue;

	// Recluster to 3 subjets and apply mass plane cuts
	JetDefinition reclustering(_jet_algorithm_recluster, 3.14);
	ClusterSequence *  cs_top_sub = new ClusterSequence(topcandidate.pieces(), reclustering);
        std::vector <PseudoJet> top_subs = sorted_by_pt(cs_top_sub->exclusive_jets(3));
	cs_top_sub->delete_self_when_unused();

	// Require the third subjet to be above the pT threshold
	if (top_subs[2].perp() < _minpt_subjet)
	  continue;

	// Modes with early 2d-massplane cuts
	if (_mode == EARLY_MASSRATIO_SORT_MASS      && !check_mass_criteria(top_subs)) {continue;}
	if (_mode == EARLY_MASSRATIO_SORT_MODDJADE  && !check_mass_criteria(top_subs)) {continue;}

	//is this candidate better than the other? -> update
	double deltatop = fabs(topcandidate.m() - _mtmass);
	double djsum = djademod(top_subs[0], top_subs[1], topcandidate)
	  + djademod(top_subs[0], top_subs[2], topcandidate)
	  + djademod(top_subs[1], top_subs[2], topcandidate);
	bool better = false;

	// Modes 0 and 1 sort by top mass
	if ( (_mode == EARLY_MASSRATIO_SORT_MASS)
	     || (_mode == LATE_MASSRATIO_SORT_MASS)) {
	  if (deltatop < _delta_top)
	    better = true;
	}
	// Modes 2 and 3 sort by modified jade distance
	else if ( (_mode == EARLY_MASSRATIO_SORT_MODDJADE)
		  || (_mode == LATE_MASSRATIO_SORT_MODDJADE)) {
	  if (djsum > _djsum)
	    better = true;
	}
	// Mode 4 is the two-step filtering. No sorting necessary as
	// we just look at the triplet of highest pT objects after
	// unclustering
	else if (_mode == TWO_STEP_FILTER) {
	  better = true;
	}
	else {
	  std::cout << "ERROR: UNKNOWN MODE (IN DISTANCE MEASURE SELECTION)" << std::endl;
	  return;
	}

	if (better) {
	  _djsum = djsum;
	  _delta_top = deltatop;
	  _is_maybe_top = true;
	  _top_candidate = topcandidate;
	  _top_subs = top_subs;
	  store_topsubjets(top_subs);
	  _top_hadrons = topcandidate.constituents();
	  // Pruning
	  double _Rprun = _initial_jet.validated_cluster_sequence()->jet_def().R();
	  JetDefinition jet_def_prune(fastjet::cambridge_algorithm, _Rprun);
	  fastjet::Pruner pruner(jet_def_prune, _zcut, _rcut_factor);
	  PseudoJet prunedjet = pruner(triple);
	  _pruned_mass = prunedjet.m();
	  _unfiltered_mass = triple.m();

	  //are all criteria fulfilled?
	  _is_masscut_passed = false;
	  if (check_mass_criteria(top_subs)) {
	    _is_masscut_passed = true;
	  }
	  _is_ptmincut_passed = false;
	  if (_top_candidate.pt() > _minpt_tag) {
	    _is_ptmincut_passed = true;
	  }
	}//end better
      }//end kk
    }//end ll
  }//end rr

  return;
}

void HEPTopTagger_fixed_R::get_info() const {
  std::cout << "#--------------------------------------------------------------------------\n";
  std::cout << "#                          HEPTopTagger Result" << std::endl;
  std::cout << "#" << std::endl;
  std::cout << "# is top candidate: " << _is_maybe_top << std::endl;
  std::cout << "# mass plane cuts passed: " << _is_masscut_passed << std::endl;
  std::cout << "# top candidate mass: " << _top_candidate.m() << std::endl;
  std::cout << "# top candidate (pt, eta, phi): ("
	    << _top_candidate.perp() << ", "
	    << _top_candidate.eta() << ", "
	    << _top_candidate.phi_std() << ")" << std::endl;
  std::cout << "# top hadrons: " << _top_hadrons.size() << std::endl;
  std::cout << "# hard substructures: " << _parts_size << std::endl;
  std::cout << "# |m - mtop| : " << _delta_top << std::endl;
  std::cout << "# djsum : " << _djsum << std::endl;
  std::cout << "# is consistency cut passed: " << _is_ptmincut_passed << std::endl;
  std::cout << "#--------------------------------------------------------------------------\n";
  return;
}

void HEPTopTagger_fixed_R::get_setting() const {
  std::cout << "#--------------------------------------------------------------------------\n";
  std::cout << "#                         HEPTopTagger Settings" << std::endl;
  std::cout << "#" << std::endl;
  std::cout << "# mode: " << _mode << " (0 = EARLY_MASSRATIO_SORT_MASS) " << std::endl;
  std::cout << "#        "         << " (1 = LATE_MASSRATIO_SORT_MASS)  " << std::endl;
  std::cout << "#        "         << " (2 = EARLY_MASSRATIO_SORT_MODDJADE)  " << std::endl;
  std::cout << "#        "         << " (3 = LATE_MASSRATIO_SORT_MODDJADE)  " << std::endl;
  std::cout << "#        "         << " (4 = TWO_STEP_FILTER)  " << std::endl;
  std::cout << "# top mass: " << _mtmass << "    ";
  std::cout << "W mass: " << _mwmass << std::endl;
  std::cout << "# top mass window: [" << _mtmin << ", " << _mtmax << "]" << std::endl;
  std::cout << "# W mass ratio: [" << _rmin << ", " << _rmax << "] (["
	    <<_rmin*_mtmass/_mwmass<< "%, "<< _rmax*_mtmass/_mwmass << "%])"<< std::endl;
  std::cout << "# mass plane cuts: (m23cut, m13min, m13max) = ("
	    << _m23cut << ", " << _m13cutmin << ", " << _m13cutmax << ")" << std::endl;
  std::cout << "# mass_drop_threshold: " << _mass_drop_threshold << "    ";
  std::cout << "max_subjet_mass: " << _max_subjet_mass << std::endl;
  std::cout << "# R_filt: " << _Rfilt << "    ";
  std::cout << "n_filt: " << _nfilt << std::endl;
  std::cout << "# minimal subjet pt: " << _minpt_subjet << std::endl;
  std::cout << "# minimal reconstructed pt: " << _minpt_tag << std::endl;
  std::cout << "# internal jet algorithms (0 = kt, 1 = C/A, 2 = anti-kt): " << std::endl;
  std::cout << "#   filtering: "<< _jet_algorithm_filter << std::endl;
  std::cout << "#   reclustering: "<< _jet_algorithm_recluster << std::endl;
  std::cout << "#--------------------------------------------------------------------------\n";

  return;
}

//uncluster a fat jet to subjets of given cone size
void HEPTopTagger::UnclusterFatjets(const vector<fastjet::PseudoJet> & big_fatjets,
				    vector<fastjet::PseudoJet> & small_fatjets,
				    const ClusterSequence & cseq,
				    const double small_radius) {
  for (unsigned i=0; i < big_fatjets.size(); i++) {
    PseudoJet this_jet = big_fatjets[i];
    PseudoJet parent1(0, 0, 0, 0), parent2(0, 0, 0, 0);
    bool test = cseq.has_parents(this_jet, parent1, parent2);
    double dR = 100;

    if(test) dR = sqrt(parent1.squared_distance(parent2));

    if (!test || dR<small_radius) {
      small_fatjets.push_back(this_jet);
    } else {
      vector<fastjet::PseudoJet> parents;
      parents.push_back(parent1);
      parents.push_back(parent2);
      UnclusterFatjets(parents, small_fatjets, cseq, small_radius);
    }
  }
}

HEPTopTagger::HEPTopTagger() : _do_optimalR(1), _do_qjets(0),
			       _mass_drop_threshold(0.8), _max_subjet_mass(30.),
			       _mode(Mode(0)), _mtmass(172.3), _mwmass(80.4), _mtmin(150.), _mtmax(200.), _rmin(0.85*80.4/172.3), _rmax(1.15*80.4/172.3),
			       _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), _minpt_tag(200.),
			       _nfilt(5), _Rfilt(0.3), _jet_algorithm_filter(fastjet::cambridge_algorithm), _minpt_subjet(0.),
			       _jet_algorithm_recluster(fastjet::cambridge_algorithm),
			       _zcut(0.1), _rcut_factor(0.5),
			       _max_fatjet_R(1.8), _min_fatjet_R(0.5), _step_R(0.1), _optimalR_threshold(0.2),
			       _R_filt_optimalR_calc(0.2), _N_filt_optimalR_calc(10), _r_min_exp_function(&R_opt_calc_funct),
			       _optimalR_mmin(150.), _optimalR_mmax(200.), _optimalR_fw(0.175), _R_opt_diff(0.3),
			       _R_filt_optimalR_pass(0.2), _N_filt_optimalR_pass(5), _R_filt_optimalR_fail(0.3), _N_filt_optimalR_fail(3),
			       _q_zcut(0.1), _q_dcut_fctr(0.5), _q_exp_min(0.), _q_exp_max(0.), _q_rigidity(0.1), _q_truncation_fctr(0.0),
			       _debug(false)
{}

HEPTopTagger::HEPTopTagger(const fastjet::PseudoJet & jet
			   ) : _do_optimalR(1), _do_qjets(0),
			       _jet(jet), _initial_jet(jet),
			       _mass_drop_threshold(0.8), _max_subjet_mass(30.),
			       _mode(Mode(0)), _mtmass(172.3), _mwmass(80.4), _mtmin(150.), _mtmax(200.), _rmin(0.85*80.4/172.3), _rmax(1.15*80.4/172.3),
			       _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), _minpt_tag(200.),
			       _nfilt(5), _Rfilt(0.3), _jet_algorithm_filter(fastjet::cambridge_algorithm), _minpt_subjet(0.),
			       _jet_algorithm_recluster(fastjet::cambridge_algorithm),
			       _zcut(0.1), _rcut_factor(0.5),
			       _max_fatjet_R(jet.validated_cluster_sequence()->jet_def().R()), _min_fatjet_R(0.5), _step_R(0.1), _optimalR_threshold(0.2),
			       _R_filt_optimalR_calc(0.2), _N_filt_optimalR_calc(10), _r_min_exp_function(&R_opt_calc_funct),
			       _optimalR_mmin(150.), _optimalR_mmax(200.), _optimalR_fw(0.175), _R_opt_diff(0.3),
			       _R_filt_optimalR_pass(0.2), _N_filt_optimalR_pass(5), _R_filt_optimalR_fail(0.3), _N_filt_optimalR_fail(3),
			       _q_zcut(0.1), _q_dcut_fctr(0.5), _q_exp_min(0.), _q_exp_max(0.), _q_rigidity(0.1), _q_truncation_fctr(0.0),
			       _fat(jet),
			       _debug(false)
{}

HEPTopTagger::HEPTopTagger(const fastjet::PseudoJet & jet,
			   double mtmass, double mwmass
			   ) : _do_optimalR(1), _do_qjets(0),
			       _jet(jet), _initial_jet(jet),
			       _mass_drop_threshold(0.8), _max_subjet_mass(30.),
			       _mode(Mode(0)), _mtmass(mtmass), _mwmass(mwmass), _mtmin(150.), _mtmax(200.), _rmin(0.85*80.4/172.3), _rmax(1.15*80.4/172.3),
			       _m23cut(0.35), _m13cutmin(0.2), _m13cutmax(1.3), _minpt_tag(200.),
			       _nfilt(5), _Rfilt(0.3), _jet_algorithm_filter(fastjet::cambridge_algorithm), _minpt_subjet(0.),
			       _jet_algorithm_recluster(fastjet::cambridge_algorithm),
			       _zcut(0.1), _rcut_factor(0.5),
			       _max_fatjet_R(jet.validated_cluster_sequence()->jet_def().R()), _min_fatjet_R(0.5), _step_R(0.1), _optimalR_threshold(0.2),
			       _R_filt_optimalR_calc(0.2), _N_filt_optimalR_calc(10), _r_min_exp_function(&R_opt_calc_funct),
			       _optimalR_mmin(150.), _optimalR_mmax(200.), _optimalR_fw(0.175), _R_opt_diff(0.3),
			       _R_filt_optimalR_pass(0.2), _N_filt_optimalR_pass(5), _R_filt_optimalR_fail(0.3), _N_filt_optimalR_fail(3),
			       _q_zcut(0.1), _q_dcut_fctr(0.5), _q_exp_min(0.), _q_exp_max(0.), _q_rigidity(0.1), _q_truncation_fctr(0.0),
			       _fat(jet),
			       _debug(false)
{}

void HEPTopTagger::run() {
  //cout << "--- new Tagger run ---" << endl;

  QjetsPlugin _qjet_plugin(_q_zcut, _q_dcut_fctr, _q_exp_min, _q_exp_max, _q_rigidity, _q_truncation_fctr);
  int maxR = int(_max_fatjet_R * 10);
  int minR = int(_min_fatjet_R * 10);
  int stepR = int(_step_R * 10);
  _qweight=-1;

  if (!_do_optimalR) {
    HEPTopTagger_fixed_R htt(_jet);
    htt.set_mass_drop_threshold(_mass_drop_threshold);
    htt.set_max_subjet_mass(_max_subjet_mass);
    htt.set_filtering_n(_nfilt);
    htt.set_filtering_R(_Rfilt);
    htt.set_filtering_minpt_subjet(_minpt_subjet);
    htt.set_filtering_jetalgorithm(_jet_algorithm_filter);
    htt.set_reclustering_jetalgorithm(_jet_algorithm_recluster);
    htt.set_mode(_mode );
    htt.set_mt(_mtmass);
    htt.set_mw(_mwmass);
    htt.set_top_mass_range(_mtmin, _mtmax);
    htt.set_mass_ratio_range(_rmin, _rmax);
    htt.set_mass_ratio_cut(_m23cut, _m13cutmin, _m13cutmax);
    htt.set_top_minpt(_minpt_tag);
    htt.set_pruning_zcut(_zcut);
    htt.set_pruning_rcut_factor(_rcut_factor);
    htt.set_debug(_debug);
    htt.set_qjets(_q_zcut, _q_dcut_fctr, _q_exp_min, _q_exp_max, _q_rigidity, _q_truncation_fctr);
    htt.run();

    _HEPTopTagger[maxR] = htt;
    _Ropt = maxR;
    _qweight = htt.q_weight();
    _HEPTopTagger_opt = _HEPTopTagger[_Ropt];
  } else {
    _qjet_def = fastjet::JetDefinition(&_qjet_plugin);
    vector<fastjet::PseudoJet> _q_constits;
    ClusterSequence* _qjet_seq;
    PseudoJet _qjet;
    const ClusterSequence* _seq;
    _seq = _initial_jet.validated_cluster_sequence();
    if (_do_qjets){
      _q_constits = _initial_jet.associated_cluster_sequence()->constituents(_initial_jet);
      _qjet_seq = new ClusterSequence(_q_constits, _qjet_def);
      _qjet = sorted_by_pt(_qjet_seq->inclusive_jets())[0];
      _qjet_seq->delete_self_when_unused();
      const QjetsBaseExtras* ext =
	dynamic_cast<const QjetsBaseExtras*>(_qjet_seq->extras());
      _qweight=ext->weight();
      _jet = _qjet;
      _seq = _qjet_seq;
      _fat = _qjet;
    }

    // Do optimalR procedure
    vector<fastjet::PseudoJet> big_fatjets;
    vector<fastjet::PseudoJet> small_fatjets;

    big_fatjets.push_back(_jet);
    _Ropt = 0;

    for (int R = maxR; R >= minR; R -= stepR) {
      UnclusterFatjets(big_fatjets, small_fatjets, *_seq, R / 10.);

      if (_debug) {cout << "R = " << R << " -> n_small_fatjets = " << small_fatjets.size() << endl;}

      _n_small_fatjets[R] = small_fatjets.size();

      // We are sorting by pt - so start with a negative dummy
      double dummy = -99999;

      for (unsigned i = 0; i < small_fatjets.size(); i++) {
	HEPTopTagger_fixed_R htt(small_fatjets[i]);
	htt.set_mass_drop_threshold(_mass_drop_threshold);
	htt.set_max_subjet_mass(_max_subjet_mass);
	htt.set_filtering_n(_nfilt);
	htt.set_filtering_R(_Rfilt);
	htt.set_filtering_minpt_subjet(_minpt_subjet);
	htt.set_filtering_jetalgorithm(_jet_algorithm_filter);
	htt.set_reclustering_jetalgorithm(_jet_algorithm_recluster);
	htt.set_mode(_mode );
	htt.set_mt(_mtmass);
	htt.set_mw(_mwmass);
	htt.set_top_mass_range(_mtmin, _mtmax);
	htt.set_mass_ratio_range(_rmin, _rmax);
	htt.set_mass_ratio_cut(_m23cut, _m13cutmin, _m13cutmax);
	htt.set_top_minpt(_minpt_tag);
	htt.set_pruning_zcut(_zcut);
	htt.set_pruning_rcut_factor(_rcut_factor);
	htt.set_debug(_debug);
	htt.set_qjets(_q_zcut, _q_dcut_fctr, _q_exp_min, _q_exp_max, _q_rigidity, _q_truncation_fctr);

	htt.run();

	if (htt.t().perp() > dummy) {
	  dummy = htt.t().perp();
	  _HEPTopTagger[R] = htt;
	}
      } //End of loop over small_fatjets

      // Only check if we have not found Ropt yet
      if (_Ropt == 0 && R < maxR) {
	// If the new mass is OUTSIDE the window ..
	if (_HEPTopTagger[R].t().m() < (1-_optimalR_threshold)*_HEPTopTagger[maxR].t().m())
	  // .. set _Ropt to the previous mass
	  _Ropt = R + stepR;
      }

      big_fatjets = small_fatjets;
      small_fatjets.clear();
    }//End of loop over R

    // if we did not find Ropt in the loop, pick the last value
    if (_Ropt == 0 && _HEPTopTagger[maxR].t().m() > 0)
      _Ropt = minR;

    //for the case that there is no tag at all (< 3 hard substructures)
    if (_Ropt == 0 && _HEPTopTagger[maxR].t().m() == 0)
      _Ropt = maxR;

    _HEPTopTagger_opt = _HEPTopTagger[_Ropt];

    Filter filter_optimalR_calc(_R_filt_optimalR_calc, SelectorNHardest(_N_filt_optimalR_calc));
    _R_opt_calc = _r_min_exp_function(filter_optimalR_calc(_fat).pt());

    Filter filter_optimalR_pass(_R_filt_optimalR_pass, SelectorNHardest(_N_filt_optimalR_pass));
    Filter filter_optimalR_fail(_R_filt_optimalR_fail, SelectorNHardest(_N_filt_optimalR_fail));
    if(optimalR_type() == 1) {
      _filt_fat = filter_optimalR_pass(_fat);
    } else {
      _filt_fat = filter_optimalR_fail(_fat);
    }
  }
}

//optimal_R type
int HEPTopTagger::optimalR_type() {
  if(_HEPTopTagger_opt.t().m() < _optimalR_mmin || _HEPTopTagger_opt.t().m() > _optimalR_mmax) {
    return 0;
  }
  if(_HEPTopTagger_opt.f_rec() > _optimalR_fw) {
    return 0;
  }
  if(_Ropt/10. - _R_opt_calc > _R_opt_diff) {
    return 0;
  }
  return 1;
}

double HEPTopTagger::nsub_unfiltered(int order, fastjet::contrib::Njettiness::AxesMode axes, double beta, double R0) {
  fastjet::contrib::Nsubjettiness nsub(order, axes, beta, R0);
  return nsub.result(_fat);
}

double HEPTopTagger::nsub_filtered(int order, fastjet::contrib::Njettiness::AxesMode axes, double beta, double R0) {
  fastjet::contrib::Nsubjettiness nsub(order, axes, beta, R0);
  return nsub.result(_filt_fat);
}

HEPTopTagger::~HEPTopTagger(){}
}
//Example for a low_pt working point

LowPt::LowPt() {};

bool LowPt::is_tagged(HEPTopTagger::HEPTopTagger htt) {
  if (htt.is_tagged()) return true;
  if (!htt.is_masscut_passed()) return false;

  double pt = htt.t().pt();
  double m_rec = htt.t().m();
  double m_ratio = (htt.W().m()/htt.t().m()) / (80.4/172.3);
  double m12 = (htt.j1() + htt.j2()).m();
  double m13 = (htt.j1() + htt.j3()).m();
  double m23 = (htt.j2() + htt.j3()).m();
  double atan_13_12 = atan(m12/m13);
  double m23_m123 = m23 / m_rec;

  double FWM_W1W2_U1 = FWM(htt,11).U(1);
  double FWM_pW1W2_U1 = FWM(htt,1011).U(2);
  double FWM_pbW1W2_U1 = FWM(htt,1111).U(2);
  double FWM_pbW2_U1 = FWM(htt,1101).U(1);

  if (pt < 150. || pt > 200.) return false;
  if (m_rec < 108. || m_rec > 282.) return false;
  if (m_ratio < 0.717 || m_ratio > 1.556) return false;
  if (atan_13_12 < 0.441 || atan_13_12 > 0.889) return false;
  if (m23_m123 < 0.412 || m23_m123 > 0.758) return false;
  if(FWM_W1W2_U1 < 0.048 || FWM_W1W2_U1 > 0.373) return false;
  if(FWM_pW1W2_U1 < 0.019 || FWM_pW1W2_U1 > 0.524) return false;
  if(FWM_pbW1W2_U1 < 0.044 || FWM_pbW1W2_U1 > 0.276) return false;
  if(FWM_pbW2_U1 < 0.145 || FWM_pbW2_U1 > 0.445) return false;

  return true;
};
//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: MeasureFunction.cc 670 2014-06-06 01:24:42Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------



FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Measure Function
//
///////

// Return all of the necessary TauComponents for specific input particles and axes
TauComponents MeasureFunction::result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const {

   // first find partition
   // this sets jetPartitionStorage and beamPartitionStorage
   PseudoJet beamPartitionStorage;
   std::vector<fastjet::PseudoJet> jetPartitionStorage = get_partition(particles,axes,&beamPartitionStorage);

   // then return result calculated from partition
   return result_from_partition(jetPartitionStorage,axes,&beamPartitionStorage);
}

std::vector<fastjet::PseudoJet> MeasureFunction::get_partition(const std::vector<fastjet::PseudoJet>& particles,
                                                               const std::vector<fastjet::PseudoJet>& axes,
                                                               PseudoJet * beamPartitionStorage) const {

   std::vector<std::vector<PseudoJet> > jetPartition(axes.size());
   std::vector<PseudoJet> beamPartition;

   // Figures out the partiting of the input particles into the various jet pieces
   // Based on which axis the parition is closest to
   for (unsigned i = 0; i < particles.size(); i++) {

      // find minimum distance; start with beam (-1) for reference
      int j_min = -1;
      double minRsq;
      if (_has_beam) minRsq = beam_distance_squared(particles[i]);
      else minRsq = std::numeric_limits<double>::max(); // make it large value

      // check to see which axis the particle is closest to
      for (unsigned j = 0; j < axes.size(); j++) {
         double tempRsq = jet_distance_squared(particles[i],axes[j]); // delta R distance
         if (tempRsq < minRsq) {
            minRsq = tempRsq;
            j_min = j;
         }
      }

      if (j_min == -1) {
         if (_has_beam) beamPartition.push_back(particles[i]);
         else assert(_has_beam);  // this should never happen.
      } else {
         jetPartition[j_min].push_back(particles[i]);
      }
   }

   // Store beam partition
   if (beamPartitionStorage) {
      *beamPartitionStorage = join(beamPartition);
   }

   // Store jet partitions
   std::vector<PseudoJet> jetPartitionStorage(axes.size(),PseudoJet(0,0,0,0));
   for (unsigned j = 0; j < axes.size(); j++) {
      jetPartitionStorage[j] = join(jetPartition[j]);
   }

   return jetPartitionStorage;
}

// does partition, but only stores index of PseudoJets
std::vector<std::list<int> > MeasureFunction::get_partition_list(const std::vector<fastjet::PseudoJet>& particles,
                                                                 const std::vector<fastjet::PseudoJet>& axes) const {

   std::vector<std::list<int> > jetPartition(axes.size());

   // Figures out the partiting of the input particles into the various jet pieces
   // Based on which axis the parition is closest to
   for (unsigned i = 0; i < particles.size(); i++) {

      // find minimum distance; start with beam (-1) for reference
      int j_min = -1;
      double minRsq;
      if (_has_beam) minRsq = beam_distance_squared(particles[i]);
      else minRsq = std::numeric_limits<double>::max(); // make it large value

      // check to see which axis the particle is closest to
      for (unsigned j = 0; j < axes.size(); j++) {
         double tempRsq = jet_distance_squared(particles[i],axes[j]); // delta R distance
         if (tempRsq < minRsq) {
            minRsq = tempRsq;
            j_min = j;
         }
      }

      if (j_min == -1) {
         assert(_has_beam); // consistency check
      } else {
         jetPartition[j_min].push_back(i);
      }
   }

   return jetPartition;
}


// Uses existing partition and calculates result
// TODO:  Can we cache this for speed up when doing area subtraction?
TauComponents MeasureFunction::result_from_partition(const std::vector<fastjet::PseudoJet>& jet_partition,
                                                     const std::vector<fastjet::PseudoJet>& axes,
                                                     PseudoJet * beamPartitionStorage) const {

   std::vector<double> jetPieces(axes.size(), 0.0);
   double beamPiece = 0.0;

   double tauDen = 0.0;
   if (!_has_denominator) tauDen = 1.0;  // if no denominator, then 1.0 for no normalization factor

   // first find jet pieces
   for (unsigned j = 0; j < axes.size(); j++) {
      std::vector<PseudoJet> thisPartition = jet_partition[j].constituents();
      for (unsigned i = 0; i < thisPartition.size(); i++) {
         jetPieces[j] += jet_numerator(thisPartition[i],axes[j]); //numerator jet piece
         if (_has_denominator) tauDen += denominator(thisPartition[i]); // denominator
      }
   }

   // then find beam piece
   if (_has_beam) {
      assert(beamPartitionStorage); // make sure I have beam information
      std::vector<PseudoJet> beamPartition = beamPartitionStorage->constituents();

      for (unsigned i = 0; i < beamPartition.size(); i++) {
         beamPiece += beam_numerator(beamPartition[i]); //numerator beam piece
         if (_has_denominator) tauDen += denominator(beamPartition[i]); // denominator
      }
   }
   return TauComponents(jetPieces, beamPiece, tauDen, _has_denominator, _has_beam);
}





} //namespace contrib

FASTJET_END_NAMESPACE
//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: Njettiness.cc 677 2014-06-12 18:56:46Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {


///////
//
// Main Njettiness Class
//
///////

Njettiness::Njettiness(const AxesDefinition & axes_def, const MeasureDefinition & measure_def)
: _axes_def(axes_def.create()), _measure_def(measure_def.create()) {
   setMeasureFunctionAndAxesFinder();  // call helper function to do the hard work
}

Njettiness::Njettiness(AxesMode axes_mode, const MeasureDefinition & measure_def)
: _axes_def(createAxesDef(axes_mode)), _measure_def(measure_def.create()) {
   setMeasureFunctionAndAxesFinder();  // call helper function to do the hard work
}

// Convert from MeasureMode enum to MeasureDefinition
// This returns a pointer that will be claimed by a SharedPtr
MeasureDefinition* Njettiness::createMeasureDef(MeasureMode measure_mode, int num_para, double para1, double para2, double para3) const {

   // definition of maximum Rcutoff for non-cutoff measures, changed later by other measures
   double Rcutoff = std::numeric_limits<double>::max();  //large number
   // Most (but not all) measures have some kind of beta value
   double beta = std::numeric_limits<double>::quiet_NaN();
   // The normalized measures have an R0 value.
   double R0 = std::numeric_limits<double>::quiet_NaN();

   // Find the MeasureFunction and set the parameters.
   switch (measure_mode) {
      case normalized_measure:
         beta = para1;
         R0 = para2;
         if(num_para == 2) {
            return new NormalizedMeasure(beta,R0);
         } else {
            throw Error("normalized_measure needs 2 parameters (beta and R0)");
         }
         break;
      case unnormalized_measure:
         beta = para1;
         if(num_para == 1) {
            return new UnnormalizedMeasure(beta);
         } else {
            throw Error("unnormalized_measure needs 1 parameter (beta)");
         }
         break;
      case geometric_measure:
         beta = para1;
         if (num_para == 1) {
            return new GeometricMeasure(beta);
         } else {
            throw Error("geometric_measure needs 1 parameter (beta)");
         }
         break;
      case normalized_cutoff_measure:
         beta = para1;
         R0 = para2;
         Rcutoff = para3; //Rcutoff parameter is 3rd parameter in normalized_cutoff_measure
         if (num_para == 3) {
            return new NormalizedCutoffMeasure(beta,R0,Rcutoff);
         } else {
            throw Error("normalized_cutoff_measure has 3 parameters (beta, R0, Rcutoff)");
         }
         break;
      case unnormalized_cutoff_measure:
         beta = para1;
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in normalized_cutoff_measure
         if (num_para == 2) {
            return new UnnormalizedCutoffMeasure(beta,Rcutoff);
         } else {
            throw Error("unnormalized_cutoff_measure has 2 parameters (beta, Rcutoff)");
         }
         break;
      case geometric_cutoff_measure:
         beta = para1;
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in geometric_cutoff_measure
         if(num_para == 2) {
           return new GeometricCutoffMeasure(beta,Rcutoff);
         } else {
            throw Error("geometric_cutoff_measure has 2 parameters (beta, Rcutoff)");
         }
         break;
      default:
         assert(false);
         break;
   }
   return NULL;
}

// Convert from AxesMode enum to AxesDefinition
// This returns a pointer that will be claimed by a SharedPtr
AxesDefinition* Njettiness::createAxesDef(Njettiness::AxesMode axes_mode) const {

   switch (axes_mode) {
      case wta_kt_axes:
         return new WTA_KT_Axes();
      case wta_ca_axes:
         return new WTA_CA_Axes();
      case kt_axes:
         return new KT_Axes();
      case ca_axes:
         return new CA_Axes();
      case antikt_0p2_axes:
         return new AntiKT_Axes(0.2);
      case onepass_wta_kt_axes:
         return new OnePass_WTA_KT_Axes();
      case onepass_wta_ca_axes:
         return new OnePass_WTA_CA_Axes();
      case onepass_kt_axes:
         return new OnePass_KT_Axes();
      case onepass_ca_axes:
         return new OnePass_CA_Axes();
      case onepass_antikt_0p2_axes:
         return new OnePass_AntiKT_Axes(0.2);
      case onepass_manual_axes:
         return new OnePass_Manual_Axes();
      case min_axes:
         return new MultiPass_Axes(100);
      case manual_axes:
         return new Manual_Axes();
      default:
         assert(false);
         return NULL;
   }
}


// Parsing needed for constructor to set AxesFinder and MeasureFunction
// All of the parameter handling is here, and checking that number of parameters is correct.
void Njettiness::setMeasureFunctionAndAxesFinder() {
   // Get the correct MeasureFunction and AxesFinders
   _measureFunction.reset(_measure_def->createMeasureFunction());
   _startingAxesFinder.reset(_axes_def->createStartingAxesFinder(*_measure_def));
   _finishingAxesFinder.reset(_axes_def->createFinishingAxesFinder(*_measure_def));
}

// setAxes for Manual mode
void Njettiness::setAxes(const std::vector<fastjet::PseudoJet> & myAxes) {
   if (_axes_def->supportsManualAxes()) {
      _currentAxes = myAxes;
   } else {
      throw Error("You can only use setAxes for manual AxesDefinitions");
   }
}

// Calculates and returns all TauComponents that user would want.
// This information is stored in _current_tau_components for later access as well.
TauComponents Njettiness::getTauComponents(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) const {
   if (inputJets.size() <= n_jets) {  //if not enough particles, return zero
      _currentAxes = inputJets;
      _currentAxes.resize(n_jets,fastjet::PseudoJet(0.0,0.0,0.0,0.0));
      _current_tau_components = TauComponents();
      _seedAxes = _currentAxes;
      _currentJets = _currentAxes;
      _currentBeam = PseudoJet(0.0,0.0,0.0,0.0);
   } else {

      _seedAxes = _startingAxesFinder->getAxes(n_jets,inputJets,_currentAxes); //sets starting point for minimization
      if (_finishingAxesFinder) {
         _currentAxes = _finishingAxesFinder->getAxes(n_jets,inputJets,_seedAxes);
      } else {
         _currentAxes = _seedAxes;
      }

      // Find partition and store information
      // (jet information in _currentJets, beam in _currentBeam)
      _currentJets = _measureFunction->get_partition(inputJets,_currentAxes,&_currentBeam);

      // Find tau value and store information
      _current_tau_components = _measureFunction->result_from_partition(_currentJets, _currentAxes,&_currentBeam);  // sets current Tau Values
   }
   return _current_tau_components;
}


// Partition a list of particles according to which N-jettiness axis they are closest to.
// Return a vector of length _currentAxes.size() (which should be N).
// Each vector element is a list of ints corresponding to the indices in
// particles of the particles belonging to that jet.
std::vector<std::list<int> > Njettiness::getPartitionList(const std::vector<fastjet::PseudoJet> & particles) const {
   // core code is in MeasureFunction
   return _measureFunction->get_partition_list(particles,_currentAxes);
}


} // namespace contrib

FASTJET_END_NAMESPACE
//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: NjettinessDefinition.cc 704 2014-07-07 14:30:43Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

   std::string NormalizedMeasure::description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Normalized Measure (beta = " << _beta << ", R0 = " << _R0 << ")";
      return stream.str();
   };

   std::string UnnormalizedMeasure::description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Unnormalized Measure (beta = " << _beta << ", in GeV)";
      return stream.str();
   };

   std::string GeometricMeasure::description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Geometric Measure (beta = " << _beta << ", in GeV)";
      return stream.str();
   };

   std::string NormalizedCutoffMeasure::description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Normalized Cutoff Measure (beta = " << _beta << ", R0 = " << _R0 << ", Rcut = " << _Rcutoff << ")";
      return stream.str();
   };

   std::string UnnormalizedCutoffMeasure::description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Unnormalized Cutoff Measure (beta = " << _beta << ", Rcut = " << _Rcutoff << ", in GeV)";
      return stream.str();
   };

   std::string GeometricCutoffMeasure::description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Geometric Cutoff Measure (beta = " << _beta << ", Rcut = " << _Rcutoff << ", in GeV)";
      return stream.str();
   };


} // namespace contrib

FASTJET_END_NAMESPACE
//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: NjettinessPlugin.cc 663 2014-06-03 21:26:41Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{



std::string NjettinessPlugin::description() const {return "N-jettiness jet finder";}


// Clusters the particles according to the Njettiness jet algorithm
// Apologies for the complication with this code, but we need to make
// a fake jet clustering tree.  The partitioning is done by getPartitionList
void NjettinessPlugin::run_clustering(ClusterSequence& cs) const
{
   std::vector<fastjet::PseudoJet> particles = cs.jets();

   // HACK: remove area information from particles (in case this is called by
   // a ClusterSequenceArea.  Will be fixed in a future FastJet release)
   for (unsigned i = 0; i < particles.size(); i++) {
      particles[i].set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>());
   }


   _njettinessFinder.getTau(_N, particles);

   std::vector<std::list<int> > partition = _njettinessFinder.getPartitionList(particles);

   std::vector<fastjet::PseudoJet> jet_indices_for_extras;

   // output clusterings for each jet
   for (size_t i0 = 0; i0 < partition.size(); ++i0) {
      size_t i = partition.size() - 1 - i0; // reversed order of reading to match axes order
      std::list<int>& indices = partition[i];
      if (indices.size() == 0) continue;
      while (indices.size() > 1) {
         int merge_i = indices.back(); indices.pop_back();
         int merge_j = indices.back(); indices.pop_back();
         int newIndex;
         double fakeDij = -1.0;

         cs.plugin_record_ij_recombination(merge_i, merge_j, fakeDij, newIndex);

         indices.push_back(newIndex);
      }
      double fakeDib = -1.0;

      int finalJet = indices.back();
      cs.plugin_record_iB_recombination(finalJet, fakeDib);
      jet_indices_for_extras.push_back(cs.jets()[finalJet]);  // Get the four vector for the final jets to compare later.
   }

   //HACK:  Re-reverse order of reading to match CS order
   reverse(jet_indices_for_extras.begin(),jet_indices_for_extras.end());

   NjettinessExtras * extras = new NjettinessExtras(_njettinessFinder.currentTauComponents(),jet_indices_for_extras,_njettinessFinder.currentAxes());
   cs.plugin_associate_extras(std::auto_ptr<ClusterSequence::Extras>(extras));

}


} // namespace contrib

FASTJET_END_NAMESPACE
//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: Nsubjettiness.cc 597 2014-04-16 23:07:55Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

//result returns tau_N with normalization dependent on what is specified in constructor
double Nsubjettiness::result(const PseudoJet& jet) const {
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   return _njettinessFinder.getTau(_N, particles);
}

TauComponents Nsubjettiness::component_result(const PseudoJet& jet) const {
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   return _njettinessFinder.getTauComponents(_N, particles);
}

//ratio result uses Nsubjettiness result to find the ratio tau_N/tau_M, where N and M are specified by user
double NsubjettinessRatio::result(const PseudoJet& jet) const {
   double numerator = _nsub_numerator.result(jet);
   double denominator = _nsub_denominator.result(jet);
   return numerator/denominator;
}

} // namespace contrib

FASTJET_END_NAMESPACE

QHTT::QHTT() : _niter(100),   _q_zcut(0.1), _q_dcut_fctr(0.5), _q_exp_min(0.), _q_exp_max(0.), _q_rigidity(0.1), _q_truncation_fctr(0.0)
{};
void QHTT::run(HEPTopTagger::HEPTopTagger htt) {
  htt.get_setting();
  _weight_q1 = -1.;
  _weight_q2 = -1.;
  _m_sum = 0.;
  _m2_sum = 0.;
  _eps_q = 0.;
  _qtags = 0;
  HEPTopTagger::HEPTopTagger _htt_q = htt;
  _htt_q.set_qjets(_q_zcut, _q_dcut_fctr, _q_exp_min, _q_exp_max, _q_rigidity, _q_truncation_fctr);
  _htt_q.do_qjets(true);
  for (int iq = 0; iq < _niter; iq++) {
    _htt_q.run();
    if (_htt_q.is_tagged()) {
      _qtags++;
      _m_sum += _htt_q.t().m();
      _m2_sum += _htt_q.t().m() * _htt_q.t().m();
      if (_htt_q.q_weight() > _weight_q1) {
	_weight_q2 = _weight_q1; _htt_q2 = _htt_q1;
	_weight_q1=_htt_q.q_weight(); _htt_q1 = _htt_q;
      } else if (_htt_q.q_weight() > _weight_q2) {
	_weight_q2=_htt_q.q_weight(); _htt_q2 = _htt_q;
      }
    }
  }
  _eps_q = float(_qtags)/float(_niter);
};
//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: WinnerTakeAllRecombiner.cc 597 2014-04-16 23:07:55Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

std::string WinnerTakeAllRecombiner::description() const {
   return "Winner Take All scheme recombination";
}

// recombine pa and pb by creating pab with energy of the sum of particle energies in the direction of the harder particle
// updated recombiner to use more general form of a metric equal to E*(pT/E)^(alpha), which reduces to pT*cosh(rap)^(1-alpha)
// alpha is specified by the user. The default is alpha = 1, which is the typical behavior. alpha = 2 provides a metric which more
// favors central jets
void WinnerTakeAllRecombiner::recombine(const fastjet::PseudoJet & pa, const fastjet::PseudoJet & pb, fastjet::PseudoJet & pab) const {
   double a_pt = pa.perp(), b_pt = pb.perp(), a_rap = pa.rap(), b_rap = pb.rap();

   // special case of alpha = 1, everything is just pt (made separate so that pow function isn't called)
   if (_alpha == 1.0) {
      if (a_pt >= b_pt) {
         pab.reset_PtYPhiM(a_pt + b_pt, a_rap, pa.phi());
      }
      else if (b_pt > a_pt) {
         pab.reset_PtYPhiM(a_pt + b_pt, b_rap, pb.phi());
      }
   }

   // every other case uses additional cosh(rap) term
   else {
      double a_metric = a_pt*pow(cosh(a_rap), 1.0-_alpha);
      double b_metric = b_pt*pow(cosh(b_rap), 1.0-_alpha);
      if (a_metric >= b_metric) {
   	  double new_pt = a_pt + b_pt*pow(cosh(b_rap)/cosh(a_rap), 1.0-_alpha);
   	  pab.reset_PtYPhiM(new_pt, a_rap, pa.phi());
      }
      if (b_metric > a_metric) {
   	  double new_pt = b_pt + a_pt*pow(cosh(a_rap)/cosh(b_rap), 1.0-_alpha);
   	  pab.reset_PtYPhiM(new_pt, b_rap, pb.phi());
      }
   }
}

} //namespace contrib

FASTJET_END_NAMESPACE
