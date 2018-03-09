#ifndef __QHTT__HH__
#define __QHTT__HH__

#include "HEPTopTagger.hh"

class QHTT {
public:
  QHTT();
  void set_iterations(int niter) {_niter = niter;};
  void set_q_zcut(double zcut) {_q_zcut = zcut;}
  void set_q_dcut_fctr(double dcut_fctr) {_q_dcut_fctr = dcut_fctr;}
  void set_q_exp(double a, double b) {_q_exp_min = a; _q_exp_max = b;}
  void set_q_rigidity(double rigidity) {_q_rigidity = rigidity;}
  void set_q_truncation_fctr(double truncation_fctr) {_q_truncation_fctr = truncation_fctr;}
  
  void run(HEPTopTagger::HEPTopTagger htt);
  HEPTopTagger::HEPTopTagger leading() {return _htt_q1;};
  HEPTopTagger::HEPTopTagger subleading() {return _htt_q2;};
  double weight_leading() {return _weight_q1;};
  double weight_subleading() {return _weight_q2;};
  double eps_q() {return _eps_q;};
  double m_mean() {return _qtags > 0 ? _m_sum/float(_qtags) : 0. ;}
  double m2_mean() {return _qtags > 0 ? _m2_sum/float(_qtags) : 0.;}
  
private:
  int _niter;
  double _q_zcut, _q_dcut_fctr, _q_exp_min, _q_exp_max, _q_rigidity, _q_truncation_fctr;
  int _qtags;
  double _weight_q1, _weight_q2;
  HEPTopTagger::HEPTopTagger _htt_q, _htt_q1, _htt_q2;
  double _m_sum, _m2_sum;
  double _eps_q;
  };
#endif
