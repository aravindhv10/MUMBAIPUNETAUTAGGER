#ifndef NewHEPHeaders2_HH
#define NewHEPHeaders2_HH
#include "CPPFileIO2.hh"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JHTopTagger.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TCanvas.h"
namespace NewHEPHeaders {
    const char junk_address = '0' ;
    const bool DEBUG = false;
    template <typename T> void set_junked ( T * & inptr ) {inptr=(T*)(&junk_address);}
    template <typename T> bool is_junked  ( T * & inptr ) {return (inptr==(T*)(&junk_address));}
    inline size_t shifter  (size_t in) {return (1<<in);}
    inline bool   checkbit (size_t inbits, size_t checkbits) {return ((inbits&checkbits)==checkbits);}
    namespace CONSTANTS   {
        const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
        const double PI2 = PI * 2.0;
    }
    namespace PID         {
        const long DOWN = 1;
        const long UP = 2;
        const long STRANGE = 3;
        const long CHARM = 4;
        const long BOTTOM = 5;
        const long TOP = 6;
        const long GAMMA = 22;
        const long Z = 23;
        const long W = 24;
        const long H0 = 25;
        const long H1 = 35;
        const long A0 = 36;
        const long CHARGED_HIGGS = 37;
        const long GLUON = 21;
        const long ELECTRON_NU = 12;
        const long MUON_NU = 14;
        const long TAU_NU = 16;
        const long ELECTRON = 11;
        const long MUON = 13;
        const long TAU = 15;
        const long CHI10 = 1000022;
        const long CHI20 = 1000023;
    }
    namespace MASS        {
        const double TOP = 173.340;
        const double W = 80.385;
        const double Z = 91.190;
        const double H0 = 125.000;
    }
    namespace DECAY_WIDTH {
        const double W = 2.085 ;
    }
    inline bool detectable (long pid) {
        pid = CPPFileIO::mymod (pid);
        if (pid == PID::ELECTRON_NU) {
            return false;
        }
        if (pid == PID::MUON_NU) {
            return false;
        }
        if (pid == PID::TAU_NU) {
            return false;
        }
        if (pid == PID::CHI10) {
            return false;
        }
        if (pid == PID::CHI20) {
            return false;
        }
        return true;
    }
    inline bool islepton (long pid) {
        pid = CPPFileIO::mymod (pid);
        if (pid == PID::ELECTRON) {
            return true;
        }
        if (pid == PID::MUON) {
            return true;
        }
        if (pid == PID::TAU) {
            return true;
        }
        return false;
    }
    inline bool isblike (long pid) {
        pid=CPPFileIO::mymod(pid);
        if(pid<100) {
            if(pid==PID::BOTTOM) {return true;}
            else {return false;}
        } else {
            while(pid>0) {
                long tmp = pid%10;
                if(tmp==PID::BOTTOM) {return true;}
                pid=pid/10;
            }
            return false ;
        }
    }
    namespace VECTORS {
        const size_t NOTHING = 0;
        const size_t BTAG = shifter(0);
        const size_t TAUTAG = shifter(1);
        const size_t BMESONTAG = shifter(2);
        double VECTOR_EQUALITY_LIMIT = 50.0;
        template < typename TR >                class plane2vector   {
        private:
        public:
            TR x[2];
            void construct (TR _x = 0, TR _y = 0) {
                x[0] = _x;
                x[1] = _y;
            }
            plane2vector (TR _x = 0, TR _y = 0) {
                construct (_x, _y);
            }
            plane2vector (const plane2vector < TR > &c) {
                construct (c.x[0], c.x[1]);
            }
            ~plane2vector () {
            }
            inline TR pt2 () {
                return (x[0] * x[0]) + (x[1] * x[1]);
            }
            inline double pt () {
                return (sqrt (pt2 ()));
            }
            inline TR safenorm2 () {
                TR mag = pt2 ();

                if (CPPFileIO::mymod (mag) < 0.0000000001) {
                    mag = CPPFileIO::mysign (mag) * 0.0000000001;
                }
                return mag;
            }
            inline double phi () {
                double ret = acos (x[0] / sqrt (safenorm2 ()));

                if (x[1] < 0) {
                    ret = CONSTANTS::PI2 - ret;
                }
                return ret;
            }
            inline double dphi (plane2vector < TR > b) {
                double ret = CPPFileIO::mymod (b.phi () - phi ());

                if (ret > CONSTANTS::PI) {
                    ret = CONSTANTS::PI2 - ret;
                }
                return ret;
            }
            inline plane2vector < TR > operator + (plane2vector < TR > b) {
                return plane2vector < TR > (x[0] + b.x[0], x[1] + b.x[1]);
            }
            inline plane2vector < TR > operator - (plane2vector < TR > b) {
                return plane2vector < TR > (x[0] - b.x[0], x[1] - b.x[1]);
            }
            inline TR operator * (plane2vector < TR > b) {
                return (x[0] * b.x[0]) + (x[1] * b.x[1]);
            }
            inline plane2vector < TR > operator * (TR b) {
                return plane2vector < TR > (x[0] * b, x[1] * b);
            }
            inline plane2vector < TR > operator / (TR b) {
                return plane2vector < TR > (x[0] / b, x[1] / b);
            }
            inline double operator  () (plane2vector < TR > b) {
                return dphi (b);
            }
            inline plane2vector < TR > flip () {
                return plane2vector < TR > (-x[0], -x[1]);
            }
            inline plane2vector < TR > dir () {
                plane2vector < TR > ret (x[0], x[1]);
                double mag = sqrt (ret.safenorm2 ());

                ret = ret / mag;
                return ret;
            }
            inline bool operator > (plane2vector < TR > b) {
                return pt2 () > b.pt2 ();
            }
            inline bool operator < (plane2vector < TR > b) {
                return pt2 () < b.pt2 ();
            }
            inline ssize_t operator >> (CPPFileIO::FileFD & f) {
                return f.multiwrite2file (*this);
            }
            inline ssize_t operator << (CPPFileIO::FileFD & f) {
                return f.multiread2file (*this);
            }
            inline TR & operator [] (size_t i) {
                return x[i];
            }
            inline void clearthis () {
                x[0] = 0;
                x[1] = 0;
            }
            inline bool operator == (plane2vector < TR > b) {
                plane2vector < TR > tmp (x[0], x[1]);
                tmp = tmp - b;
                TR diff = tmp.pt2 ();

                diff = CPPFileIO::mymod (diff);
                return diff < VECTOR_EQUALITY_LIMIT;
            }
        };
        template < typename TR >                class euclid3vector  {
        private:
        public:
            plane2vector < TR > xy;
            TR z;

            euclid3vector (TR _x = 0, TR _y = 0, TR _z = 0):xy (_x, _y) {
                z = _z;
            }
            euclid3vector (plane2vector < TR > a, TR _z = 0):xy (a) {
                z = _z;
            }
            euclid3vector (const euclid3vector < TR > &a):xy (a.xy) {
                z = a.z;
            }
            inline double phi () {
                return xy.phi ();
            }
            inline TR pt2 () {
                return xy.pt2 ();
            }
            inline double pt () {
                return xy.pt ();
            }
            inline euclid3vector < TR > operator + (euclid3vector < TR > a) {
                return euclid3vector < TR > (xy + a.xy, z + a.z);
            }
            inline euclid3vector < TR > operator - (euclid3vector < TR > a) {
                return euclid3vector < TR > (xy - a.xy, z - a.z);
            }
            inline TR operator * (euclid3vector < TR > a) {
                return (xy * a.xy) + (z * a.z);
            }
            inline euclid3vector < TR > operator * (TR a) {
                return euclid3vector < TR > (xy * a, z * a);
            }
            inline euclid3vector < TR > operator / (TR a) {
                return euclid3vector < TR > (xy / a, z / a);
            }
            inline TR p2 () {
                return xy.pt2 () + (z * z);
            }
            inline double p () {
                return sqrt (p2 ());
            }
            inline TR & operator [] (size_t ret) {
                if (ret > 1) {
                    return z;
                } else {
                    return xy[ret];
                }
            }
            inline double eta () {
                double tmp_p = p ();

                return 0.5 * log ((tmp_p + z) / (tmp_p - z));
            }
            inline double meta () {
                return CPPFileIO::mymod (eta ());
            }
            inline double cone2 (euclid3vector < TR > b) {
                double tphi = xy.dphi (b.xy);

                tphi = tphi * tphi;
                double teta = eta () - b.eta ();

                teta = teta * teta;
                double ret = teta + tphi;

                return ret;
            }
            inline double cone (euclid3vector < TR > b) {
                return sqrt (cone2 (b));
            }
            inline double dphi (euclid3vector < TR > b) {
                double tphi = xy.dphi (b.xy);

                return tphi;
            }
            inline double operator  () (euclid3vector < TR > b) {
                return cone (b);
            }
            inline TR safenorm2 () {
                TR mag = xy.pt2 () + (z * z);

                if (CPPFileIO::mymod (mag) < 0.0000000001) {
                    mag = CPPFileIO::mysign (mag) * 0.0000000001;
                }
                return mag;
            }
            inline euclid3vector < TR > flip () {
                return euclid3vector < TR > (xy.flip (), -z);
            }
            inline euclid3vector < TR > trans () {
                return euclid3vector < TR > (xy, 0);
            }
            inline euclid3vector < TR > dir () {
                euclid3vector < TR > ret (*this);
                double mag = sqrt (ret.safenorm2 ());

                ret = ret / mag;
                return ret;
            }
            inline bool operator > (euclid3vector < TR > b) {
                return pt2 () > b.pt2 ();
            }
            inline bool operator < (euclid3vector < TR > b) {
                return pt2 () < b.pt2 ();
            }
            inline ssize_t operator >> (CPPFileIO::FileFD & f) {
                return f.multiwrite2file (*this);
            }
            inline ssize_t operator << (CPPFileIO::FileFD & f) {
                return f.multiread2file (*this);
            }
            inline void clearthis () {
                xy = plane2vector < TR > (0, 0);
                z = 0;
            }
            inline bool operator == (euclid3vector < TR > b) {
                euclid3vector < TR > tmp = (*this) - b;
                TR diff = tmp.p2 ();

                diff = CPPFileIO::mymod (diff);
                return diff < VECTOR_EQUALITY_LIMIT;
            }
        };
        template < typename TR >                class lorentz4vector {
        private:
        public:
            euclid3vector < TR > xyz;
            TR t;

            lorentz4vector (TR _x = 0, TR _y = 0, TR _z = 0, TR _t = 0):xyz (_x, _y, _z) {
                t = _t;
            }
            lorentz4vector (euclid3vector < TR > a, TR _t = -1) :xyz (a) {
                if(_t<0) {t=a.p();} else {t=_t;}
            }
            lorentz4vector (plane2vector  < TR > a, TR _z = 0 , TR _t = -1) :xyz (a) {
                if (_t<0) {t=a.pt();} else {t=_t;}
                xyz.z = _z ;
            }
            lorentz4vector (const lorentz4vector < TR > &a):xyz (a.xyz) {
                t = a.t;
            }
            lorentz4vector (const fastjet::PseudoJet & injet) {
                t      = injet.e  () ;
                xyz[2] = injet.pz () ;
                xyz[1] = injet.py () ;
                xyz[0] = injet.px () ;
            }
            inline TR & operator [] (size_t ref) {
                if (ref > 2) {
                    return t;
                } else {
                    return xyz[ref];
                }
            }
            inline TR pt2 () {
                return xyz.pt2 ();
            }
            inline double pt () {
                return xyz.pt ();
            }
            inline TR p2 () {
                return xyz.p2 ();
            }
            inline double p () {
                return xyz.p ();
            }
            inline double phi () {
                return xyz.phi ();
            }
            inline TR m2 () {
                return CPPFileIO::mymod ((t * t) - p2 ());
            }
            inline TR n2 () {
                return CPPFileIO::mymod ((t * t) + p2 ());
            }
            inline double eta () {
                return (0.5 * log ((t + xyz.z) / (t - xyz.z)));
            }
            inline double meta () {
                return CPPFileIO::mymod (eta ());
            }
            inline double peta () {
                return xyz.eta ();
            }
            inline double pmeta () {
                return xyz.meta ();
            }
            inline double m () {
                return sqrt (m2 ());
            }
            inline double n () {
                return sqrt (n2 ());
            }
            inline double dphi (lorentz4vector < TR > b) {
                return xyz.dphi (b.xyz);
            }
            inline lorentz4vector < TR > operator + (lorentz4vector < TR > b) {
                return lorentz4vector < TR > (xyz + b.xyz, t + b.t);
            }
            inline lorentz4vector < TR > operator - (lorentz4vector < TR > b) {
                return lorentz4vector < TR > (xyz - b.xyz, t - b.t);
            }
            inline TR operator * (lorentz4vector < TR > b) {
                return (t * b.t) - (xyz * b.xyz);
            }
            inline lorentz4vector < TR > operator * (TR b) {
                return lorentz4vector < TR > (xyz * b, t * b);
            }
            inline lorentz4vector < TR > operator / (TR b) {
                return lorentz4vector < TR > (xyz / b, t / b);
            }
            inline void operator = ( euclid3vector < TR > other ) {
                xyz = other     ;
                t   = other.p() ;
            }
            inline double pcone2 (lorentz4vector < TR > b) {
                double tphi = xyz.dphi (b.xyz);

                tphi = tphi * tphi;
                double teta = peta () - b.peta ();

                teta = teta * teta;
                double ret = teta + tphi;

                return ret;
            }
            inline double pcone (lorentz4vector < TR > b) {
                return sqrt (pcone2 (b));
            }
            inline double cone2 (lorentz4vector < TR > b) {
                double tphi = xyz.dphi (b.xyz);

                tphi = tphi * tphi;
                double teta = eta () - b.eta ();

                teta = teta * teta;
                double ret = teta + tphi;

                return ret;
            }
            inline double cone (lorentz4vector < TR > b) {
                return sqrt (cone2 (b));
            }
            inline double operator  () (lorentz4vector < TR > b) {
                return cone (b);
            }
            inline lorentz4vector < TR > flip () {
                return lorentz4vector < TR > (xyz.flip (), t);
            }
            inline lorentz4vector < TR > trans () {
                return lorentz4vector < TR > (xyz.trans (), t);
            }
            inline lorentz4vector < TR > dir () {
                return lorentz4vector < TR > (xyz.dir (), t);
            }
            inline bool operator > (lorentz4vector < TR > b) {
                return pt2 () > b.pt2 ();
            }
            inline bool operator < (lorentz4vector < TR > b) {
                return pt2 () < b.pt2 ();
            }
            inline ssize_t operator >> (CPPFileIO::FileFD & f) {
                return f.multiwrite2file (*this);
            }
            inline ssize_t operator << (CPPFileIO::FileFD & f) {
                return f.multiread2file (*this);
            }
            inline bool cleared () {
                return (t < 0);
            }
            inline bool pass () {
                return (t > 0);
            }
            inline void clearthis () {
                t = -1;
                xyz = euclid3vector < TR > (0, 0, 0);
            }
            inline bool operator == (lorentz4vector < TR > b) {
                lorentz4vector < TR > tmp = (*this) - b;
                TR diff = tmp.n2 ();
                diff = CPPFileIO::mymod (diff);
                return diff < VECTOR_EQUALITY_LIMIT;
            }
            inline double gamma () {
                return (double) t / m ();
            }
            inline TR gamma2 () {
                return t * t / m2 ();
            }
            inline double beta () {
                double g = (double) gamma2 ();
                g = 1.0 / g;
                g = 1.0 - g;
                return sqrt (g);
            }
            inline euclid3vector < TR > Velocity () {
                return xyz.dir () * beta ();
            }
            inline lorentz4vector < TR > LorentzBoost (euclid3vector < TR > booster) {
                euclid3vector < TR > parallel = booster * ((xyz * booster) / booster.p2 ());
                euclid3vector < TR > perpendicular = xyz - parallel;
                double gm = booster.p2 ();

                gm = 1.0 - gm;
                gm = (1.0 / gm);
                gm = sqrt (gm);
                lorentz4vector < TR > ret;
                ret.t = gm * (t - (parallel * booster));
                parallel = (parallel - (booster * t)) * (TR) gm;
                ret.xyz = parallel + perpendicular;
                return ret;
            }
            inline lorentz4vector < TR > LorentzBoostGamma (euclid3vector < TR > booster) {
                double gm2   = (double) booster.p2 ()        ;
                double gm    = (double) sqrt       (gm2)     ;
                double beta2 = (double) 1.0    -   (1.0/gm2) ;
                double beta  = (double) sqrt       (beta2)   ;
                euclid3vector < TR > dir = booster / (TR)gm   ;
                euclid3vector < TR > vel = dir     * (TR)beta ;
                euclid3vector < TR > parallel      = dir * (dir*xyz) ;
                euclid3vector < TR > perpendicular = xyz - parallel  ;
                lorentz4vector < TR > ret ; /* Evaluate the vector: */ {
                    ret.t = ( t - (parallel*vel) ) * (TR)gm ;
                    parallel = ( parallel - (vel*t) ) * (TR)gm;
                    ret.xyz = parallel + perpendicular;
                }
                return ret;
            }
            inline fastjet::PseudoJet getpseudojet () {
                return fastjet::PseudoJet (xyz[0], xyz[1], xyz[2], t);
            }
        };
        template < typename TRF, typename TRI > class ParticleNode   {
        private:
            lorentz4vector <TRF> momentum;
            TRI d1, d2, pid;
        public:
            ~ParticleNode () {}
            ParticleNode () {
                d1 = -1;
                d2 = -1;
                pid = 0;
                momentum.clearthis ();
            }
            ParticleNode (TRF _x, TRF _y, TRF _z, TRF _t, TRI _d1, TRI _d2, TRI _pid):momentum (_x, _y, _z, _t) {
                d1 = _d1;
                d2 = _d2;
                pid = _pid;
            }
            ParticleNode (const Pythia8::Particle & part) {
                momentum = lorentz4vector <TRF> (part.px (), part.py (), part.pz (), part.e ());
                pid = part.id ();
                if (part.isFinal ()) {
                    d1 = -1;
                    d2 = -1;
                } else {
                    d1 = part.daughter1 ();
                    d2 = part.daughter2 ();
                }
            }
            inline int id () {
                return pid;
            }
            inline bool isFinal () {
                return (d1 == -1) && (d2 == -1);
            }
            inline int daughter1 () {
                return d1;
            }
            inline int daughter2 () {
                return d2;
            }
            inline float px () {
                return momentum[0];
            }
            inline float py () {
                return momentum[1];
            }
            inline float pz () {
                return momentum[2];
            }
            inline float e () {
                return momentum[3];
            }
            inline float pt () {
                return momentum.pt ();
            }
            inline float eta () {
                return momentum.eta ();
            }
            inline float modeta () {
                return CPPFileIO::mymod (eta());
            }
            inline bool isDetectable () {
                bool ret = (pt () > 0.5) && (modeta () < 6.0);
                ret = ret && detectable (pid);
                return ret;
            }
            inline bool IsGood () {
                return isFinal () && isDetectable ();
            }
            inline bool IsLepton () {
                return islepton (pid);
            }
            inline bool IsBLike () {
                return isblike(pid);
            }
            inline bool IsBMeson () {
                long tmppid = CPPFileIO::mymod(pid) ;
                return ((tmppid>100)&&isblike(tmppid));
            }
            inline bool IsBQuakr () {
                long tmppid = CPPFileIO::mymod(pid) ;
                return (tmppid==PID::BOTTOM);
            }
            inline TRF operator [] (size_t i) {
                return momentum[i] ;
            }
            inline lorentz4vector <TRF> & getvec () {
                return momentum;
            }
            inline lorentz4vector <TRF> & operator () () {
                return momentum;
            }
            inline double operator  () (ParticleNode b) {
                return momentum (b.momentum);
            }
            inline double operator  () (lorentz4vector <TRF> b) {
                return momentum (b);
            }
            inline double pcone (ParticleNode b) {
                return momentum.pcone (b.momentum);
            }
            inline double pcone (lorentz4vector <TRF> b) {
                return momentum.pcone (b);
            }
            inline fastjet::PseudoJet getpseudojet () {
                return momentum.getpseudojet ();
            }
        };
        template < typename TR >                class JetContainer : public lorentz4vector <TR> {
        private:
            inline lorentz4vector <TR> & CstGet (size_t i) {return vectors[0][constituents[i]];}
            std::vector <int> constituents;
            std::vector <lorentz4vector <TR>> *vectors;
        public:
            size_t TAG;
            fastjet::PseudoJet * injet;
            bool tau_tag () {
                if (checkbit(TAG,TAUTAG)) {return true;}
                else {return false;}
            }
            bool bot_tag () {
                if (checkbit(TAG,BTAG)) {return true;}
                else {return false;}
            }
            bool bms_tag () {
                if (checkbit(TAG,BMESONTAG)) {return true;}
                else {return false;}
            }
            JetContainer (std::vector <lorentz4vector<TR>> & _vectors, fastjet::PseudoJet & _injet) {
                TAG = NOTHING;
                /* Get some values: */  {
                    injet = &_injet;
                    vectors = &_vectors;
                    this[0][0] = _injet.px ();
                    this[0][1] = _injet.py ();
                    this[0][2] = _injet.pz ();
                    this[0][3] = _injet.e ();
                }
                /* Get the constituents: */  {
                    constituents.clear ();
                    std::vector < fastjet::PseudoJet > _constituents = _injet.constituents ();
                    for (size_t j = 0; j < _constituents.size (); j++) {
                        constituents.push_back (_constituents[j].user_index ());
                    }
                }
            }
            ~JetContainer () {}
            inline lorentz4vector <TR> & operator  () (size_t i) {return CstGet (i);}
            inline size_t operator  () () {return constituents.size ();}
        };
        class                                         sortorder      {
        private:
        public:
            bool ret;
            sortorder (bool _ret = true) {ret = _ret;}
            ~sortorder () {}
            template <typename T> inline bool operator  () (VECTORS::lorentz4vector <T> & a, VECTORS::lorentz4vector <T> & b) {
                if (a>b) {return ret;}
                else {return !ret;}
            }
            template <typename T> inline bool operator  () (VECTORS::lorentz4vector <T> * a, VECTORS::lorentz4vector <T> * b) {
                if ((*a)>(*b)) {return ret;}
                else {return !ret;}
            }
        };
        template <typename T> void cleanup_vectors ( std::vector <lorentz4vector <T>> & thevectors) {
            std::vector <VECTORS::lorentz4vector <T>> tmp_thevectors;
            for (size_t i = 0; i < thevectors.size (); i++)
                if (!thevectors[i].cleared ()) {
                    tmp_thevectors.push_back (thevectors[i]);
                }
                if (tmp_thevectors.size () > 0) {
                    thevectors.resize (tmp_thevectors.size ());
                    size_t n = thevectors.size () * sizeof (VECTORS::lorentz4vector <T>);
                    memcpy ((void *) &(thevectors[0]), (const void *) &(tmp_thevectors[0]), (size_t) n);
                } else { thevectors.clear (); }
        }
        template <typename TR> class RestTopDecays {
        private:

            TR          ST2 , CT2 , SP2 , CP2 ;
            TR          ST3 , CT3 , SP3 , CP3 ;
            euclid3vector  <TR> Velocity      ;
            lorentz4vector <TR> Daughters [3] ;
            lorentz4vector <TR> InVectors [3] ;

            inline bool ClearAll        () {
                Angle            = -10000.0 ;
                SpinAngle        = -10000.0 ;
                Matrix           =      0.0 ;
                MatrixSpin       =      0.0 ;
                Weight           =      0.0 ;
                WeightSpin       =      0.0 ;
                PhaseSpaceWeight =      0.0 ;
                passed           =    false ;
                return               passed ;
            }
            inline void SetPreliminary  () {
                MT  = MASS::TOP                      ;
                TOP = lorentz4vector <TR> (0,0,0,MT) ;
            }
            inline void SetDirections   () {
                ST2 = sin (T2) ; CT2 = cos (T2) ;
                ST3 = sin (T3) ; CT3 = cos (T3) ;
                SP2 = sin (P2) ; CP2 = cos (P2) ;
                SP3 = sin (P3) ; CP3 = cos (P3) ;
                dir[0] = euclid3vector <TR> (ST2*CP2,ST2*SP2,CT2) ;
                dir[1] = euclid3vector <TR> (ST3*CP3,ST3*SP3,CT3) ;
                Angle  = dir[0] * dir[1]  ;
            }
            inline bool EvalVectorSlave () {
                E3     = ((2.0*E2)-MT) / ( 2.0 * (((E2/MT)*(1.0-Angle))-1.0) ) ;
                passed = ( E3 >= 0 ) ;
                if (!passed) { return ClearAll () ; }
                Daughters [0] = dir[0] * E2                       ;
                Daughters [1] = dir[1] * E3                       ;
                Daughters [2] = TOP - Daughters[0] - Daughters[1] ;
                return passed ;
            }
            inline bool EvalVectors     () {
                SetPreliminary () ;
                passed = ( E2 <= (MT/2.0) ) && ( E2 >= 0 ) ;
                if (!passed) { return ClearAll () ; }
                SetDirections () ;
                return EvalVectorSlave () ;
            }
            inline bool PhaseSpace      () {
                if (!passed) { return ClearAll () ; }
                PhaseSpaceWeight = ( ((1.0-Angle)*E2) - MT ) ;
                PhaseSpaceWeight = 8 * PhaseSpaceWeight * PhaseSpaceWeight ;
                PhaseSpaceWeight = PhaseSpaceWeight / E2*MT*((2*E2)-MT)*ST2*ST3 ;
                return passed ;
            }
            inline bool EvalMatrix      () {
                if (!passed) { return ClearAll () ; }
                double tmp     = MASS::W              / MT ;
                double complex = DECAY_WIDTH::W * tmp / MT ;
                double real    = 1 - (tmp*tmp) - (2*E2/MT) ;
                Matrix         = (E3*MT*MT) * (MT-(2*E3)) / ( (real*real) + (complex*complex) ) ;
                return passed;
            }
            inline bool FinalWeight     () {
                if (!passed) { return ClearAll () ; }
                Weight = Matrix * PhaseSpaceWeight ;
                return passed ;
            }
            inline bool EvalChain       () {
                if (PhaseSpace()) if (EvalMatrix()) if (FinalWeight()) {return true;}
                return ClearAll () ;
            }
            inline bool EvalSpinMatrix ( euclid3vector  <TR> SpinDir ) {
                if (!passed) { return ClearAll () ; }
                MatrixSpin = ( (dir[1]*SpinDir) + 1 ) * Matrix ;
                WeightSpin = MatrixSpin * PhaseSpaceWeight ;
                return passed ;
            }
            inline bool permute        (size_t i, size_t j) {
                size_t k     = 3 - i - j    ;
                Daughters[0] = InVectors[i] ;
                Daughters[1] = InVectors[j] ;
                Daughters[2] = InVectors[k] ;
                dir[0] = Daughters[0].xyz.dir() ;
                dir[1] = Daughters[1].xyz.dir() ;
                Angle  = dir[0] * dir[1] ;
                ST2 = dir[0].pt() ;
                ST3 = dir[1].pt() ;
                E2  = Daughters[0][3] ;
                E3  = Daughters[1][3] ;
                passed = true  ;
                return EvalChain () ;
            }
            inline bool ChangeE2       (double _E2) {
                E2 = _E2 ;
                passed = ( E2 <= (MT/2.0) ) && ( E2 >= 0 ) ;
                if (!passed) { return ClearAll () ; }
                if (EvalVectorSlave()) {return EvalChain ();}
                return ClearAll () ;
            }

            inline bool Init (lorentz4vector <TR> a, lorentz4vector <TR> b, lorentz4vector <TR> c) {
                InVectors[0] = a     ;
                InVectors[1] = b     ;
                InVectors[2] = c     ;
                TOP          = InVectors[2] + InVectors[1] + InVectors[0] ;
                Velocity     = TOP.xyz.dir() ;
                Velocity     = Velocity*(TOP[3]/TOP.m()) ;
                TOP          = TOP.LorentzBoostGamma          (Velocity) ;
                InVectors[0] = InVectors[0].LorentzBoostGamma (Velocity) ;
                InVectors[1] = InVectors[1].LorentzBoostGamma (Velocity) ;
                InVectors[2] = InVectors[2].LorentzBoostGamma (Velocity) ;
                return         permute (0,1)                             ;
            }
            inline void Init (lorentz4vector <TR> *_Daughters) { Init (_Daughters[0],_Daughters[1],_Daughters[2]) ; }
            inline bool Init (TR _E2, TR _T2, TR _P2, TR _T3, TR _P3 ) {
                E2=_E2; T2=_T2; P2=_P2; T3=_T3; P3=_P3;
                if (EvalVectors()) if (PhaseSpace()) if (EvalMatrix()) if (FinalWeight()) {return true;}
                return false ;
            }

        public:

            euclid3vector  <TR> dir        [3]  ;
            euclid3vector  <TR> SpinVector      ;
            lorentz4vector <TR> TOP             ;

            TR E2 , T2 , P2 , T3 , P3 , E3 , MT ;
            TR Angle  , SpinAngle  ;
            TR Matrix , MatrixSpin ;
            TR Weight , WeightSpin ;
            TR PhaseSpaceWeight    ;
            bool passed            ;

            inline TR                  operator () ( TR _E2                                 ) {
                if (!passed) {return 0;}
                ChangeE2(_E2); return Weight;
            }
            inline TR                  operator () ( size_t i , size_t j                    ) {
                if (permute(i,j)) {return Weight ;}
                else {return 0;}
            }
            inline bool                operator () ( lorentz4vector <TR> *_Daughters        ) { return Init ( _Daughters          ) ; }
            inline bool                operator () ( TR _E2, TR _T2, TR _P2, TR _T3, TR _P3 ) { return Init ( _E2,_T2,_P2,_T3,_P3 ) ; }
            inline TR                  operator () ( euclid3vector  <TR> SpinDir            ) {
                if (!passed) { return 0 ; }
                return EvalSpinMatrix (SpinDir) ;
            }
            inline TR                  operator () ( TR TS , TR PS                          ) {
                if (!passed) { return 0 ; }
                TR STS = sin (TS) ; TR CTS = cos (TS) ;
                TR SPS = sin (PS) ; TR CPS = cos (PS) ;
                euclid3vector <TR>     SpinDir  (STS*CPS,STS*SPS,CTS) ;
                return EvalSpinMatrix (SpinDir)                       ;
            }
            inline lorentz4vector <TR> operator [] ( size_t i                               ) { return Daughters [i] ; }

            RestTopDecays(){}
            RestTopDecays ( lorentz4vector <TR> *_Daughters        ) { Init ( _Daughters          ) ; }
            RestTopDecays ( TR _E2, TR _T2, TR _P2, TR _T3, TR _P3 ) { Init ( _E2,_T2,_P2,_T3,_P3 ) ; }
            ~RestTopDecays(){}
        };
    }
    typedef VECTORS::plane2vector   < float       > vector2         ;
    typedef VECTORS::euclid3vector  < float       > vector3         ;
    typedef VECTORS::lorentz4vector < float       > vector4         ;
    typedef VECTORS::RestTopDecays  < float       > topdecayelement ;
    typedef VECTORS::ParticleNode   < float , int > pythia_node     ;
    typedef VECTORS::JetContainer   < float       > JetContainer    ;
    typedef std::vector < vector2 > vector2s ;
    typedef std::vector < vector3 > vector3s ;
    typedef std::vector < vector4 > vector4s ;
    typedef std::vector < JetContainer > JetContainers ;
    typedef std::vector < pythia_node  > pythia_nodes  ;
    typedef std::vector < fastjet::PseudoJet > pseudojets ;
    inline vector4 pseudojet2vector4 (const fastjet::PseudoJet & injet) {
        return vector4 (injet.px(),injet.py(),injet.pz(),injet.e()) ;
    }
    namespace MainEventData {
        class EventData {
        private:
            pythia_nodes store ;
            float event_weight ;
            float event_sigma  ;
            fastjet::ClusterSequence * clust_seq_nrmjet ;
            fastjet::JetDefinition jet_def_nrmjet ;
            inline void clear_jet_cluster () {
                if(!is_junked(clust_seq_nrmjet)) {
                    delete     clust_seq_nrmjet ;
                    set_junked(clust_seq_nrmjet);
                }
            }
            inline void clear_derived () {
                leptons.clear();
                taus.clear();
                jets.clear();
                bjets.clear();
                gamma.clear();
                bquarks.clear();
                bmesons.clear();
                hadrons.clear();
                tojets.clear();
                MET = vector4(0,0,0,0);
                HT = 0 ;
                MHT = 0 ;
                clear_jet_cluster();
            }
            inline void clear_all () {
                clear_derived();
                store.clear();
                event_weight=0;
                event_sigma=0;
            }
            inline double CalcIso (size_t idx, double cone_rad = 0.3) {
                double sum_pt = 0.0;
                if ( (idx>0) && (idx<store.size()) ) {
                    for ( size_t i = 0 ; i < store.size () ; i++ ) {
                        if (
                            (i!=idx) &&
                            (store[i].IsGood ()) &&
                            (store[i](store[idx])<cone_rad)
                        ) {sum_pt=sum_pt+store[i].pt();}
                    }
                    sum_pt = sum_pt / store[idx].pt ();
                }
                return sum_pt;
            }
            inline bool isisolatedlepton (size_t idx) {
                bool ret = ( (store[idx].modeta()<2.5) && (store[idx].IsLepton()) && (store[idx].pt()>20) ) ;
                if (ret) {ret = (CalcIso(idx,0.3)<0.3);}
                return ret ;
            }
            inline bool isisolatedgamma (size_t idx) {
                bool ret = ( (store[idx].modeta()<2.5) && (store[idx].id()==(int)PID::GAMMA) && (store[idx].pt()>20) ) ;
                if (ret) {ret = (CalcIso(idx,0.4)<0.1);}
                return ret ;
            }
            inline void liquidateindex (std::vector <size_t> & indices, vector4s & thevectors) {
                CPPFileIO::deduplicate(indices);
                for (size_t i=0;i<indices.size();i++) { thevectors.push_back (store[indices[i]].getvec()) ; }
            }
            inline void cluster_jets () {
                for(size_t i=0;i<hadrons.size();i++){tojets.push_back(hadrons[i].getpseudojet());}
                clust_seq_nrmjet = new fastjet::ClusterSequence (tojets,jet_def_nrmjet) ;
                jets = sorted_by_pt(clust_seq_nrmjet->inclusive_jets(20.0));
                //printf("DEBUG (cluster_jets): number of jets = %ld\n",jets.size());
            }
            inline bool check_bquarktag (size_t idx) {
                vector4 tmpvec = pseudojet2vector4 (jets[idx]) ;
                if(tmpvec.meta()>2.5) {return false;}
                else for(size_t i=0;i<bquarks.size();i++) if(bquarks[i](tmpvec)<0.3) {return true;}
            }
            inline long getnearestjet2bquark(size_t idx) {
                double nearestdist = 1000 ;
                long ret = -1 ;
                for(size_t i=0;i<jets.size();i++) {
                    vector4 tmpvec = pseudojet2vector4 (jets[i]) ;
                    double tmpdist = bquarks [idx] (tmpvec) ;
                    if (tmpdist<nearestdist) { ret = i ; nearestdist = tmpdist ; }
                }
                if (nearestdist>0.3) {ret=-1;}
                return ret ;
            }
            inline bool check_bmesontag (size_t idx) {
                vector4 tmpvec = pseudojet2vector4 (jets[idx]) ;
                if(tmpvec.meta()>2.5) {return false;}
                else for(size_t i=0;i<bmesons.size();i++) if(bmesons[i](tmpvec)<0.3) {return true;}
            }
            inline void btag_jets () {
                std::vector <size_t> indices ;
                for(size_t i=0;(i<bquarks.size());i++) {
                    long index = getnearestjet2bquark (i) ;
                    if ((index>=0)&&(check_bmesontag((size_t)index))) {
                        indices.push_back(index);
                    }
                }
                pseudojets tmpnjets ;
                for(size_t i=0;i<jets.size();i++) {
                    bool btagged=false ;
                    for(size_t j=0;(j<indices.size())&&(!btagged);j++) if (indices[j]==i) {btagged=true;}
                    if (btagged) {bjets.push_back(jets[i]);}
                    else {tmpnjets.push_back(jets[i]);}
                }
                CPPFileIO::clone_vector(tmpnjets,jets);
                //printf("DEBUG (btag_jets): number of (njets,bjets) = %ld %ld\n",jets.size(),bjets.size());
            }
            inline void Calc_MHT () {
                MHT = 0 ;
                for ( size_t i=0 ; i < jets.size    () ; i++ ) { MHT = MHT + jets[i].pt             () ; }
                for ( size_t i=0 ; i < leptons.size () ; i++ ) { MHT = MHT + leptons[i].getvec().pt () ; }
            }
            inline void QuickConstruct(){set_junked(clust_seq_nrmjet);clear_all();}

        public:

            ssize_t operator >> (CPPFileIO::FileFD & f) {
                ssize_t ret = 0;
                ret = ret + (f.multiwrite2file (event_sigma));
                ret = ret + (f.multiwrite2file (event_weight));
                ret = ret + (f << store);
                return ret;
            }
            ssize_t operator << (CPPFileIO::FileFD & f) {
                clear_all();
                ssize_t ret = 0;
                ret = ret + (f.multiread2file (event_sigma));
                ret = ret + (f.multiread2file (event_weight));
                ret = ret + (f >> store);
                return ret;
            }

            vector4s     bquarks ;
            vector4s     bmesons ;
            vector4s     taus    ;
            vector4s     hadrons ;
            pseudojets   tojets  ;
            pseudojets   jets    ;
            pseudojets   bjets   ;
            pythia_nodes leptons ;
            vector4s     gamma   ;
            vector4      MET     ;
            double       HT      ;
            double       MHT     ;

            inline pythia_node & operator () (size_t i) { return store[i] ; }
            inline size_t operator () () { return store.size() ; }
            inline size_t recurse (size_t idx) {
                if (idx > 0) {
                    int PID = store[idx].id ();
                    size_t d1 = store[idx].daughter1 ();
                    size_t d2 = store[idx].daughter2 ();
                    if ((d1 > idx) && (store[d1].id () == PID)) {
                        return recurse (d1);
                    } else if ((d2 > idx) && (store[d2].id () == PID)) {
                        return recurse (d2);
                    }
                }
                return idx;
            }
            inline pythia_node & operator [] (size_t i) { return store[recurse(i)] ; }
            inline size_t findparticle (long PID, bool takemod=false) {
                size_t ret = 0;
                if(takemod) { PID = CPPFileIO::mymod(PID) ; }
                for (size_t i = 0; i < store.size (); i++) {
                    long ppid = store[i].id ();
                    if(takemod) { ppid = CPPFileIO::mymod(ppid) ; }
                    if (ppid == PID) { return i; }
                }
                return ret;
            }
            inline size_t getdaughter (size_t idx, long pid) {
                size_t ret = 0;
                if (idx > 0) {
                    idx = recurse (idx);
                    size_t d1 = store[idx].daughter1 ();
                    size_t d2 = store[idx].daughter2 ();
                    if ((d1 > 0) || (d2 > 0)) {
                        if (d2 > d1)
                            for (size_t i = d1; (i >= d1) && (i <= d2) && (i > 0); i++) {
                                long PPID = (int) store[i].id ();
                                long MPID = CPPFileIO::mymod (PPID);
                                if (PPID == pid) {
                                    ret = i;
                                    i = d2 + 1;
                                } else if (MPID == pid) {
                                    ret = i;
                                }
                            } else {
                                /* Check the first daughter: */  {
                                    size_t i = d1;
                                    long PPID = (int) store[i].id ();
                                    long MPID = CPPFileIO::mymod (PPID);
                                    if (PPID == pid) {
                                        ret = i;
                                    } else if (MPID == pid) {
                                        ret = i;
                                    }
                                }
                                /* Check the second daughter: */  {
                                    size_t i = d2;
                                    long PPID = (int) store[i].id ();
                                    long MPID = CPPFileIO::mymod (PPID);
                                    if (PPID == pid) {
                                        ret = i;
                                        i = d2 + 1;
                                    } else if (MPID == pid) {
                                        ret = i;
                                    }
                                }
                            }
                    }
                }
                return ret;
            }
            inline void prepare () {
                clear_derived();
                std::vector <size_t> tmp_bmeson ;
                std::vector <size_t> tmp_bquark ;
                std::vector <size_t> tmp_taus   ;
                for(size_t i=0;i<store.size();i++){
                    if (store[i].isFinal()) {
                        if (store[i].isDetectable()) {
                            if(isisolatedlepton(i)) {leptons.push_back(store[i]);}
                            else if (isisolatedgamma(i)) {gamma.push_back(store[i].getvec());}
                            else {hadrons.push_back(store[i].getvec());}
                            MET = MET + store[i].getvec() ;
                            HT = HT + store[i].getvec().pt() ;
                        }
                    }
                    else if (store[i].IsBMeson()) {tmp_bmeson.push_back(recurse(i));}
                    else if (CPPFileIO::mymod(store[i].id())==PID::BOTTOM) {tmp_bquark.push_back(recurse(i));}
                    else if (CPPFileIO::mymod(store[i].id())==PID::TAU) {tmp_taus.push_back(recurse(i));}
                }
                /* De duplicating and transfer: */ {
                    liquidateindex(tmp_bmeson,bmesons);
                    liquidateindex(tmp_bquark,bquarks);
                    liquidateindex(tmp_taus,taus);
                }
                MET = MET * -1  ;
                cluster_jets () ;
                Calc_MHT     () ;
                btag_jets    () ;
            }
            inline void FromPythia (Pythia8::Event & event) {
                clear_all () ;
                for (int i = 0; i < event.size (); i++) { store.push_back (pythia_node(event[i])) ; }
                event_sigma  = 1.0 ;
                event_weight = 1.0 ;
            }
            template < typename T > void ReadFromDelphes (T & delin) {
                clear_all () ;
                event_weight = delin.Event_Weight [0] ;
                event_sigma  = 1.0 ;
                for ( int i = 0 ; i < delin.Particle_ ; i++ ) {
                    float _x = delin.Particle_Px  [i] ;
                    float _y = delin.Particle_Py  [i] ;
                    float _z = delin.Particle_Pz  [i] ;
                    float _t = delin.Particle_E   [i] ;
                    int _d1  = delin.Particle_D1  [i] ;
                    int _d2  = delin.Particle_D2  [i] ;
                    int _pid = delin.Particle_PID [i] ;
                    pythia_node tmp (_x,_y,_z,_t,_d1,_d2,_pid) ;
                    store.push_back (tmp) ;
                }
            }
            EventData(double _jetsize = 0.4) : jet_def_nrmjet(fastjet::antikt_algorithm,_jetsize)
            {QuickConstruct();}
            EventData (Pythia8::Event & event) : jet_def_nrmjet(fastjet::antikt_algorithm,0.4)
            {QuickConstruct();FromPythia(event);}
            ~EventData() {clear_jet_cluster();}
        };
    } typedef MainEventData::EventData EventData ;
    namespace MassDropTagger {
        class HardSubStructureFinder {
        private:
            double max_subjet_mass, mass_drop_threshold, Rfilt, minpt_subjet;
            double mhmax, mhmin, mh; double zcut, rcut_factor;
            size_t nfilt;
            inline void find_structures (const fastjet::PseudoJet & this_jet) {
                fastjet::PseudoJet parent1(0,0,0,0), parent2(0,0,0,0);
                bool haskid=this_jet.validated_cs()->has_parents(this_jet,parent1,parent2);
                if(haskid) {
                    if (parent1.m()<parent2.m()) {std::swap(parent1,parent2);}
                    double kidmass    = parent1.m()  + parent2.m() ;
                    double parentmass = this_jet.m()               ;
                    if(kidmass<parentmass*mass_drop_threshold){
                        t_parts.push_back(parent1);
                        t_parts.push_back(parent2);
                        return;
                    } else if (parent1.m()>mass_drop_threshold*parent2.m()) {find_structures(parent1);return;}
                } else {return;}
                if(this_jet.m()<max_subjet_mass){t_parts.push_back(this_jet);}
                else {
                    fastjet::PseudoJet parent1(0,0,0,0), parent2(0,0,0,0);
                    bool haskid = this_jet.validated_cs()->has_parents(this_jet,parent1,parent2);
                    if (haskid) {
                        if (parent1.m()<parent2.m()) {std::swap(parent1,parent2);}
                        find_structures(parent1);
                        if (parent1.m()<mass_drop_threshold*this_jet.m()) {find_structures(parent2);}
                    }
                }
            }
            inline void run (fastjet::PseudoJet&injet) {
                t_parts.clear(); find_structures(injet);
                if(t_parts.size()>1) {
                    t_parts=sorted_by_pt(t_parts);
                    size_t i=0; size_t j=1;
                    triple = fastjet::join (t_parts[i],t_parts[j]) ;
                    filt_tau_R = std::min ( Rfilt , 0.5*sqrt(t_parts[i].squared_distance(t_parts[j])) ) ;
                    fastjet::JetDefinition filtering_def(fastjet::cambridge_algorithm,filt_tau_R);
                    fastjet::Filter filter(filtering_def,fastjet::SelectorNHardest(nfilt)*fastjet::SelectorPtMin(minpt_subjet));
                    taucandidate=filter(triple); filteredjetmass=taucandidate.m();
                    if((mhmin<filteredjetmass)&&(filteredjetmass<mhmax)&&(taucandidate.pieces().size()>1)){
                        fastjet::JetDefinition   reclustering (fastjet::cambridge_algorithm,10.0)  ;
                        fastjet::ClusterSequence cs_top_sub   (taucandidate.pieces(),reclustering) ;
                        tau_subs=sorted_by_pt(cs_top_sub.exclusive_jets(2));
                        if (tau_subs[1].perp()>minpt_subjet) {
                            HiggsTagged=true;
                            Higgs=tau_subs[0]+tau_subs[1];
                            deltah=CPPFileIO::mymod(taucandidate.m()-mh);
                            tau_hadrons=taucandidate.constituents();
                            double Rprun=injet.validated_cluster_sequence()->jet_def().R();
                            fastjet::JetDefinition jet_def_prune(fastjet::cambridge_algorithm, Rprun);
                            fastjet::Pruner pruner(jet_def_prune,zcut,rcut_factor);
                            prunedjet=pruner(triple);
                            prunedmass=prunedjet.m();
                            unfiltered_mass=triple.m();
                        }
                    }
                }
            }
            inline void initialize () {
                t_parts.clear(); tau_subs.clear(); tau_hadrons.clear();
                max_subjet_mass=30; Rfilt=0.3; minpt_subjet=20;
                mass_drop_threshold=0.7; nfilt=4; filteredjetmass=0.0;
                mh=125.0; mhmax=mh+100.0; mhmin=mh-100.0; filt_tau_R=0;
                zcut=0.1; rcut_factor=0.5; prunedmass=0.0; unfiltered_mass=0.0;
                deltah=10000; HiggsTagged=false;
            }
        public:
            double             filteredjetmass , deltah   , filt_tau_R  , prunedmass   , unfiltered_mass ;
            pseudojets         tau_subs        , t_parts  , tau_hadrons ;
            fastjet::PseudoJet prunedjet       , triple   , Higgs       , taucandidate ;
            bool HiggsTagged ;
            inline void operator () () {initialize();}
            inline void operator () (fastjet::PseudoJet&injet) {run(injet);}
            HardSubStructureFinder(){initialize();}
            ~HardSubStructureFinder(){}
        };
    }
    typedef MassDropTagger::HardSubStructureFinder HardSubStructureFinder;
    namespace DELPHES_DETDATA {
        const size_t BTAG = 1;
        const size_t TAUTAG = 2;
        const size_t NOTHING = 0;
        class DetVector:public vector4 {
        private:
        public:
            int   charge     ;
            float H_fraction ;
            float E_fraction ;
            float M_fraction ;
            DetVector () {
                charge     = 0 ;
                H_fraction = 0 ;
                E_fraction = 0 ;
                M_fraction = 0 ;
            }
            ~DetVector () {}
        }; typedef std::vector < DetVector > DetVectors;
        class JetContainer:public vector4 {
        private:
            inline DetVector & CstGet     (size_t i) {return vectors[0][constituents[i]];}
            inline void        initialize (DetVectors & _vectors, fastjet::PseudoJet & _injet) {
                weight=1.0; TAG=NOTHING;
                /* Get some values: */  {
                    injet      = & _injet       ;
                    vectors    = & _vectors     ;
                    this[0][0] =   _injet.px () ;
                    this[0][1] =   _injet.py () ;
                    this[0][2] =   _injet.pz () ;
                    this[0][3] =   _injet.e  () ;
                }
                /* Get the constituents: */  {
                    constituents.clear ();
                    std::vector <fastjet::PseudoJet> _constituents = _injet.constituents () ;
                    for (size_t j = 0; j < _constituents.size (); j++) {
                        constituents.push_back (_constituents[j].user_index());
                    }
                }
                MainHiggsTauTagger(); MainHiggsTauTagger(_injet);
            }
        public:
            double                   weight       ;
            size_t                   TAG          ;
            std::vector        <int> constituents ;
            DetVectors         *     vectors      ;
            fastjet::PseudoJet *     injet        ;
            MassDropTagger::HardSubStructureFinder MainHiggsTauTagger ;

            int    count_tracks   () {
                int ret = 0;
                for (size_t i = 0; i < constituents.size (); i++) if (CstGet(i).charge != 0) {ret++;}
                return ret;
            }
            int    count_neutrals () {
                int ret = 0;
                for (size_t i = 0; i < constituents.size (); i++) if (CstGet (i).charge == 0) {ret++;}
                return ret;
            }
            int    count_hadrons  () {
                int ret = 0;
                for (size_t i = 0; i < constituents.size (); i++) {
                    if (
                        (CstGet (i).H_fraction > CstGet (i).M_fraction) &&
                        (CstGet (i).H_fraction > CstGet (i).E_fraction)
                    ) {ret++;}
                }
                return ret;
            }
            double Lambda         () {
                std::vector < fastjet::PseudoJet > tmp_const = injet->constituents ();
                fastjet::JetDefinition jet_def (fastjet::antikt_algorithm, 2.0 * 2.0 / pt ());
                fastjet::ClusterSequence clust_seq (tmp_const, jet_def);
                std::vector < fastjet::PseudoJet > inclusive_jets = sorted_by_pt (clust_seq.inclusive_jets (10.0));
                double ret = 0;
                if (inclusive_jets.size () > 0) {
                    ret = log (1.0 - (inclusive_jets[0].pt () / injet->pt ()));
                }
                return ret;
            }
            double project        () {
                NewHEPHeaders::vector3 tmp = xyz.dir ();
                double ret = 0;
                for (size_t i = 0; i < constituents.size (); i++) {ret=ret+(CstGet(i).xyz.dir()*tmp);}
                ret = ret / constituents.size ();
                return ret;
            }
            bool   tau_tag        () {
                if ((TAG&TAUTAG)==TAUTAG) {return true;}
                return false;
            }
            bool   bot_tag        () {
                if ((TAG&BTAG)==BTAG) {return true;}
                return false;
            }
            inline DetVector & operator  () (size_t i) {return CstGet (i);}
            inline size_t      operator  () ()         {return constituents.size ();}
            JetContainer (DetVectors & _vectors, fastjet::PseudoJet & _injet) {initialize (_vectors,_injet);}
            ~JetContainer(){}
        }; typedef std::vector < JetContainer > JetContainers;
        class GenPartContainer {
        private:
        public:
            pythia_nodes store;
            vector4s taus;
            vector4s bots;
            float weight;
            void clear_all () {
                store.clear ();
                taus.clear ();
                bots.clear ();
                weight = 0;
            }
            size_t findparticle (long PID) {
                size_t ret = 0;
                long MPID = CPPFileIO::mymod (PID);

                for (size_t i = 0; i < store.size (); i++) {
                    long ppid = store[i].id ();
                    long mpid = CPPFileIO::mymod (ppid);

                    if (ppid == PID) {
                        return i;
                    } else if (mpid == MPID) {
                        ret = i;
                    }
                }
                return ret;
            }
            size_t recurse (size_t idx) {
                if (idx > 0) {
                    int PID = store[idx].id ();
                    size_t d1 = store[idx].daughter1 ();
                    size_t d2 = store[idx].daughter2 ();
                    if      ((d1>idx)&&(store[d1].id()==PID)) {return recurse (d1);}
                    else if ((d2>idx)&&(store[d2].id()==PID)) {return recurse (d2);}
                }
                return idx;
            }
            size_t getdaughter (size_t idx, long pid) {
                size_t ret = 0;
                if (idx > 0) {
                    idx = recurse (idx);
                    size_t d1 = store[idx].daughter1 ();
                    size_t d2 = store[idx].daughter2 ();
                    if ((d1 > 0) || (d2 > 0)) {
                        if (d2 > d1)
                            for (size_t i = d1; (i >= d1) && (i <= d2) && (i > 0); i++) {
                                long PPID = (int) store[i].id ();
                                long MPID = CPPFileIO::mymod (PPID);

                                if (PPID == pid) {
                                    ret = i;
                                    i = d2 + 1;
                                } else if (MPID == pid) {
                                    ret = i;
                                }
                            } else {
                                /* Check the first daughter: */  {
                                    size_t i = d1;
                                    long PPID = (int) store[i].id ();
                                    long MPID = CPPFileIO::mymod (PPID);

                                    if (PPID == pid) {
                                        ret = i;
                                    } else if (MPID == pid) {
                                        ret = i;
                                    }
                                }
                                /* Check the second daughter: */  {
                                    size_t i = d2;
                                    long PPID = (int) store[i].id ();
                                    long MPID = CPPFileIO::mymod (PPID);

                                    if (PPID == pid) {
                                        ret = i;
                                        i = d2 + 1;
                                    } else if (MPID == pid) {
                                        ret = i;
                                    }
                                }
                            }
                    }
                }
                return ret;
            }
            void GetTaus () {
                taus.clear ();
                std::vector < size_t > tau_indices;
                /* Get genuine taus: */  {
                    size_t Higgs_Index = 0;
                    for (size_t i = 0; (i < store.size ()) && (Higgs_Index == 0); i++) {
                        long ppid = store[i].id ();
                        if ((ppid == NewHEPHeaders::PID::H0)||(ppid == NewHEPHeaders::PID::Z)) {
                            Higgs_Index=i;
                        }
                    }
                    if (Higgs_Index > 0) {
                        Higgs_Index = recurse (Higgs_Index);
                    } else {
                        printf ("SCREWEDDDDDDD!!!! higgs not found\n");
                    }
                    size_t d1 = store[Higgs_Index].daughter1 ();
                    size_t d2 = store[Higgs_Index].daughter2 ();

                    tau_indices.push_back (recurse (d1));
                    tau_indices.push_back (recurse (d2));
                }
                for (size_t i = 0; i < tau_indices.size (); i++) {
                    vector4 tmptauvect = store[tau_indices[i]].getvec();
                    if (CPPFileIO::mymod (tmptauvect.eta ()) < 4.0) {taus.push_back(tmptauvect);}
                }
            }
            void GetBots () {
                bots.clear ();
                std::vector < size_t > tau_indices;
                for (size_t i = 0; i < store.size (); i++)
                    if (CPPFileIO::mymod (store[i].id()) == PID::BOTTOM) {
                        size_t tmpindex = recurse (i);

                        tau_indices.push_back (tmpindex);
                    }
                    CPPFileIO::deduplicate (tau_indices);
                for (size_t i = 0; i < tau_indices.size (); i++) {
                    bots.push_back (store[tau_indices[i]].getvec());
                }
            }
            template < typename T > void ReadFromDelphes (T & delin) {
                clear_all ();
                weight = delin.Event_Weight[0];
                for (int i = 0; i < delin.Particle_; i++) {
                    float _x = delin.Particle_Px[i];
                    float _y = delin.Particle_Py[i];
                    float _z = delin.Particle_Pz[i];
                    float _t = delin.Particle_E[i];
                    int _d1 = delin.Particle_D1[i];
                    int _d2 = delin.Particle_D2[i];
                    int _pid = delin.Particle_PID[i];
                    pythia_node tmp (_x, _y, _z, _t, _d1, _d2, _pid);

                    store.push_back (tmp);
                }
                GetTaus ();
                GetBots ();
            }
            GenPartContainer  () {}
            ~GenPartContainer () {}
        };
        class DelphesContainer {
        private:
        public:
            DetVectors    AllParticles         ;
            DetVectors    ToJets               ;
            JetContainers jets                 ;
            vector4s      Electrons    , Muons ;
            float         weight       , sigma ;
            fastjet::JetDefinition     jet_def   ;
            fastjet::ClusterSequence * clust_seq ;
            std::vector <fastjet::PseudoJet> jetvectors;
            std::vector <fastjet::PseudoJet> inclusive_jets;
            double calciso       (size_t index, double DeltaR) {
                double DeltaR2 = DeltaR * DeltaR;
                double ret = 0;
                if (index < AllParticles.size ()) {
                    for (size_t i = 0; i < index; i++) {
                        if (AllParticles[index].pcone2 (AllParticles[i]) < DeltaR2) {
                            ret = ret + AllParticles[i].pt ();
                        }
                    }
                    for (size_t i = index + 1; i < AllParticles.size (); i++) {
                        if (AllParticles[index].pcone2 (AllParticles[i]) < DeltaR2) {
                            ret = ret + AllParticles[i].pt ();
                        }
                    }
                }
                return ret;
            }
            void   DecideLeptJet () {
                for (size_t i = 0; i < AllParticles.size (); i++) {
                    if ((AllParticles[i].pt () > 15) && (AllParticles[i].charge != 0)) {
                        if (AllParticles[i].M_fraction > 0.9) {
                            if (calciso (i, 0.3) < (0.3 * AllParticles[i].pt ())) {
                                Muons.push_back (NewHEPHeaders::vector4 (AllParticles[i][0], AllParticles[i][1], AllParticles[i][2], AllParticles[i][3]));
                            } else { ToJets.push_back (AllParticles[i]); }
                        } else if (AllParticles[i].E_fraction > 0.9) {
                            if (calciso (i, 0.3) < (0.3 * AllParticles[i].pt ())) {
                                Electrons.push_back (NewHEPHeaders::vector4 (AllParticles[i][0], AllParticles[i][1], AllParticles[i][2], AllParticles[i][3]));
                            } else { ToJets.push_back (AllParticles[i]); }
                        }

                    } else { ToJets.push_back (AllParticles[i]); }
                }
            }
            void   ClusterJets   () {
                jetvectors.clear (); /* Get the vectors: */  {
                    for (size_t i = 0; i < ToJets.size (); i++) {
                        fastjet::PseudoJet tmp (ToJets[i][0], ToJets[i][1], ToJets[i][2], ToJets[i][3]);
                        tmp.set_user_index (i);
                        jetvectors.push_back (tmp);
                    }
                }
                if (jetvectors.size () > 0) {
                    clust_seq = new fastjet::ClusterSequence (jetvectors, jet_def);
                    inclusive_jets = sorted_by_pt (clust_seq->inclusive_jets (100.0));
                    for (size_t i = 0; i < inclusive_jets.size (); i++) {
                        JetContainer tmp (ToJets, inclusive_jets[i]);
                        tmp.weight = weight;
                        jets.push_back (tmp);
                    }
                }
            }
            template < typename T > void ReadFromDelphes (T & indata) {
                weight = indata.Event_Weight[0];
                for (size_t i = 0; i < indata.EFlowTrack_; i++) {
                    DetVector tmpvector; /* Assign the vectors: */  {
                        /* Setting momentum components: */  {
                            TLorentzVector tmp;
                            tmp.SetPtEtaPhiM (indata.EFlowTrack_PT[i], indata.EFlowTrack_Eta[i], indata.EFlowTrack_Phi[i], 0);
                            tmpvector[0] = tmp.Px () ;
                            tmpvector[1] = tmp.Py () ;
                            tmpvector[2] = tmp.Pz () ;
                            tmpvector[3] = tmp.E  () ;
                        }
                        /* Set charge and type: */  {
                            tmpvector.charge = indata.EFlowTrack_Charge[i];
                            if (CPPFileIO::mymod (indata.EFlowTrack_PID[i]) == PID::MUON) {
                                tmpvector.M_fraction = 1.0;
                            } else if (CPPFileIO::mymod (indata.EFlowTrack_PID[i]) == PID::ELECTRON) {
                                tmpvector.E_fraction = 1.0;
                            } else { tmpvector.H_fraction = 1.0; }
                        }
                    }
                    AllParticles.push_back (tmpvector);
                }
                for (size_t i = 0; i < indata.EFlowPhoton_; i++) {
                    DetVector tmpvector; /* Assign the vectors: */  {
                        /* Setting momentum components: */  {
                            TLorentzVector tmp;
                            tmp.SetPtEtaPhiM (indata.EFlowPhoton_ET[i], indata.EFlowPhoton_Eta[i], indata.EFlowPhoton_Phi[i], 0);
                            tmpvector[0] = tmp.Px ();
                            tmpvector[1] = tmp.Py ();
                            tmpvector[2] = tmp.Pz ();
                            tmpvector[3] = tmp.E  ();
                        }
                        /* Set charge and type: */  {
                            tmpvector.charge = 0;
                            tmpvector.E_fraction = indata.EFlowPhoton_Eem[i];
                            tmpvector.H_fraction = indata.EFlowPhoton_Ehad[i];
                        }
                    }
                    AllParticles.push_back (tmpvector);
                }
                for (size_t i = 0; i < indata.EFlowNeutralHadron_; i++) {
                    DetVector tmpvector; /* Assign the vectors: */  {
                        /* Setting momentum components: */  {
                            TLorentzVector tmp;
                            tmp.SetPtEtaPhiM (indata.EFlowNeutralHadron_ET[i], indata.EFlowNeutralHadron_Eta[i], indata.EFlowNeutralHadron_Phi[i], 0);
                            tmpvector[0] = tmp.Px ();
                            tmpvector[1] = tmp.Py ();
                            tmpvector[2] = tmp.Pz ();
                            tmpvector[3] = tmp.E  ();
                        }
                        /* Set charge and type: */  {
                            tmpvector.charge = 0;
                            tmpvector.E_fraction = indata.EFlowNeutralHadron_Eem[i];
                            tmpvector.H_fraction = indata.EFlowNeutralHadron_Ehad[i];
                        }
                    }
                    AllParticles.push_back (tmpvector);
                }
                DecideLeptJet () ;
                ClusterJets   () ;
            }

            DelphesContainer () : jet_def (fastjet::cambridge_aachen_algorithm, 1.0) {
                AllParticles.clear ();
                Electrons.clear ();
                Muons.clear ();
                inclusive_jets.clear ();
                clust_seq = (fastjet::ClusterSequence *) (&junk_address);
            }
            ~DelphesContainer () {
                if (clust_seq!=((fastjet::ClusterSequence*)(&junk_address))) {delete clust_seq;}
            }
        };
        class FullDelphesContainer {
        private:
        public:
            DelphesContainer                 detinfo ;
            GenPartContainer                 geninfo ;
            std::vector      <JetContainer*> taujets ;
            std::vector      <JetContainer*> botjets ;
            template < typename T > inline void ReadFromDelphes (T & indata) {
                taujets.clear ();
                botjets.clear ();
                geninfo.ReadFromDelphes (indata);
                detinfo.ReadFromDelphes (indata);
                for (size_t i = 0; i < geninfo.taus.size (); i++) {
                    double smallest = 10000;
                    long good_index = -10000;
                    for (size_t j = 0; j < detinfo.jets.size (); j++) {
                        double tmpdist = detinfo.jets[j].cone (geninfo.taus[i]);
                        if (tmpdist < smallest) {
                            smallest = tmpdist;
                            good_index = j;
                        }
                    }
                    if (smallest < 0.3) {
                        taujets.push_back (&(detinfo.jets[good_index]));
                        detinfo.jets[good_index].TAG = (detinfo.jets[good_index].TAG | TAUTAG);
                    }
                }
                for (size_t i = 0; i < geninfo.bots.size (); i++) {
                    double smallest = 10000;
                    long good_index = -10000;
                    for (size_t j = 0; j < detinfo.jets.size (); j++) {
                        double tmpdist = detinfo.jets[j].cone (geninfo.bots[i]);
                        if (tmpdist < smallest) {
                            smallest = tmpdist;
                            good_index = j;
                        }
                    }
                    if (smallest < 0.3) {
                        botjets.push_back (&(detinfo.jets[good_index]));
                        detinfo.jets[good_index].TAG = (detinfo.jets[good_index].TAG | BTAG);
                    }
                }
            }
            FullDelphesContainer () {}
            ~FullDelphesContainer () {}
        };
    }
    namespace TopReconstruction {
        class reco_top_had {
        private:
            vector4 *njet[2];
            vector4 *bjet;
        public:
            vector4 top, w;
            bool pass;
            inline void clearthis () {
                top.clearthis ();
                w.clearthis ();
                pass = false;
            }
            inline bool cleared () {
                bool ret = top.cleared () || w.cleared () || (!pass);
                return ret;
            }
            inline double error_w () {
                if (cleared ()) {
                    return 1000000000;
                } else {
                    const double sigma_w = 10;
                    double ret = (w.m () - MASS::W) / sigma_w;
                    ret = ret * ret;
                    return ret;
                }
            }
            inline double error_top () {
                if (cleared ()) {
                    return 1000000000;
                } else {
                    const double sigma_top = 25;
                    double ret = (top.m () - MASS::TOP) / sigma_top;
                    ret = ret * ret;
                    return ret;
                }
            }
            inline double error () {
                if (cleared ()) {
                    return 1000000000;
                } else {
                    return sqrt (error_top () + error_w ());
                }
            }
            inline bool operator > (reco_top_had other) { return error () > other.error () ; }
            inline bool operator < (reco_top_had other) { return error () < other.error () ; }
            inline bool operator > (double other)       { return error () > other          ; }
            inline bool operator < (double other)       { return error () < other          ; }
            inline void clear_parts () {
                njet[0]->clearthis () ;
                njet[1]->clearthis () ;
                bjet->clearthis    () ;
            }
            reco_top_had () {clearthis ();}
            reco_top_had (vector4 & _bjet, vector4 & _njet1, vector4 & _njet0) {
                pass = false;
                top.clearthis ();
                w.clearthis ();
                njet[0] = &_njet0;
                njet[1] = &_njet1;
                bjet = &_bjet;
                if (
                    (!njet[0]->cleared ()) &&
                    (!njet[1]->cleared ()) &&
                    (!bjet->cleared ())
                ) {
                    w = (*njet[0]) + (*njet[1]);
                    top = w + (*bjet);
                    pass = true;
                }
            }
        };
        class reco_top_lept {
        private:
            vector4 * bjet, *lept;
            vector2 * met ;
            void construct ( vector4 & _bjet , vector4 & _lept , vector2 & _met ) {
                clearthis ();
                lept = &_lept;
                met = & _met ;
                bjet = &_bjet;
                const double WMASS2 = MASS::W * MASS::W;
                pnu[0] = NewHEPHeaders::vector4(*met) ;
                pnu[1] = pnu[0];
                double k = (WMASS2 / 2.0) + ((*lept)[0] * (*met)[0]) + ((*lept)[1] * (*met)[1]);
                k = k / (*lept)[3];
                double a = (*lept)[2] / (*lept)[3];
                a = a * a;
                a = a - 1.0;
                double b = 2.0 * k * (*lept)[2] / (*lept)[3];
                double c = (k * k) - ((*met)[0] * (*met)[0]) - ((*met)[1] * (*met)[1]);
                double dis = (b * b) - (4.0 * a * c);
                if (dis > 0) {
                    pnu[0][2] = (-b + sqrt (dis)) / (2.0 * a);
                    pnu[0][3] = pnu[0].p ();
                    pnu[1][2] = (-b - sqrt (dis)) / (2.0 * a);
                    pnu[1][3] = pnu[1].p ();
                    double diff0 = (((*bjet) + pnu[0] + (*lept)).m ()) - MASS::TOP;
                    diff0 = CPPFileIO::mymod (diff0);
                    double diff1 = (((*bjet) + pnu[1] + (*lept)).m ()) - MASS::TOP;
                    diff1 = CPPFileIO::mymod (diff1);
                    if (diff0 > diff1) {
                        CPPFileIO::myswap (pnu[0], pnu[1]);
                    }
                    w = (pnu[0]) + (*lept);
                    top = w + (*bjet);
                    pass = true;
                    if (CPPFileIO::mymod (w.m () - MASS::W) > 1.0) {
                        printf ("ERROR: Something has gone blatently wrong... W Mass error \n");
                    }
                } else {clearthis ();}

            }
        public:
            vector4 top, pnu[2], w;
            bool pass;

            inline void clearthis () {
                top.clearthis ();
                pnu[0].clearthis ();
                pnu[1].clearthis ();
                w.clearthis ();
                pass = false;
            }
            inline bool cleared () {
                bool ret = top.cleared () || w.cleared () || pnu[0].cleared () || (!pass);
                return ret;
            }
            inline double error () {
                if (cleared ()) {
                    return 1000000000;
                } else {
                    const double sigma_top = 15;
                    double ret = (top.m () - MASS::TOP) / sigma_top;
                    return CPPFileIO::mymod (ret);
                }
            }
            inline void clear_parts () {
                bjet->clearthis ();
                met->clearthis  ();
                lept->clearthis ();
            }
            inline bool operator > (reco_top_had other) {
                return error () > other.error ();
            }
            inline bool operator < (reco_top_had other) {
                return error () < other.error ();
            }
            inline bool operator > (reco_top_lept other) {
                return error () > other.error ();
            }
            inline bool operator < (reco_top_lept other) {
                return error () < other.error ();
            }
            inline bool operator > (double other) {
                return error () > other;
            }
            inline bool operator < (double other) {
                return error () < other;
            }

            reco_top_lept () {
                clearthis ();
            }
            reco_top_lept ( vector4 & _bjet , vector4 & _lept , vector2 & _met ) {
                construct (_bjet,_lept,_met) ;
            }
            reco_top_lept ( vector4 & _bjet , vector4 & _lept , vector4 & _met ) {
                construct (_bjet,_lept,_met.xyz.xy) ;
            }
        };
    }
    typedef TopReconstruction::reco_top_had  reco_top_had  ;
    typedef TopReconstruction::reco_top_lept reco_top_lept ;
    class WriteHepmc2Fifo {
      private:
        std::string name;
        HepMC::Pythia8ToHepMC ToHepMC;
        HepMC::IO_GenEvent ascii_io;
      public:
        WriteHepmc2Fifo (std::string _name):ascii_io (&(_name[0]), std::ios::out) {
            name = _name;
        }
         ~WriteHepmc2Fifo () {
        }
        inline void operator  () (Pythia8::Pythia & pythia) {
            HepMC::GenEvent * hepmcevt = new HepMC::GenEvent ();
            ToHepMC.fill_next_event (pythia, hepmcevt);
            ascii_io << hepmcevt;
            delete hepmcevt;
        }
    };
    double LHA2pythia_node (std::string lhafile, std::string outnodefile) {
        CPPFileIO::FileFD outfile (outnodefile) ; outfile.writefile();
        Pythia8::Pythia pythia;
        pythia.readString ("Beams:frameType = 4");
        std::string tmp ("Beams:LHEF = ");
        tmp = tmp + lhafile;
        pythia.readString (&(tmp[0]));
        pythia.init ();
        for (size_t iEvent = 0; (!pythia.info.atEndOfFile ()); iEvent++) if (pythia.next ()) {
            EventData tmpbuf (pythia.event) ;
            tmpbuf >> outfile ;
        }
        pythia.stat ();

        double eventsigma = pythia.info.sigmaGen(); ;
        /* Write cross section: */ {
            std::string sigmafile(outnodefile) ;
            CPPFileIO::FileVector <double> sigmawriter (sigmafile + ".sigma") ;
            sigmawriter.push_back(eventsigma);
        }
        return eventsigma ;
    }
    void LHA2HEPMC (std::string lhafile, std::string hepmcfile) {
        WriteHepmc2Fifo outfile (hepmcfile);

        Pythia8::Pythia pythia;
        pythia.readString ("Beams:frameType = 4");
        std::string tmp ("Beams:LHEF = ");
        tmp = tmp + lhafile;
        pythia.readString (&(tmp[0]));
        pythia.init ();
        for (size_t iEvent = 0; (!pythia.info.atEndOfFile ()); iEvent++)
            if (pythia.next ()) {
                outfile (pythia);
            }
        pythia.stat ();
    }
    class myhist {
      private:
        TH1F * thehist;
        std::vector < float >vals;
          std::vector < float >weights;
          std::string name;
        bool converted;
        bool convert_hist () {
            size_t nums = vals.size ();
            if ((nums > 1) && (!converted)) {
                std::sort (vals.begin (), vals.end ());
                float hist_begin = vals[0];
                float hist_end = vals[nums - 1];
                float diff = hist_end - hist_begin;
                  hist_begin = hist_begin - (diff / 10.0);
                  hist_end = hist_end + (diff / 10.0);
                size_t nbins = CPPFileIO::mymin ((size_t) nums / 4, N_Bins_Max);
                  nbins = CPPFileIO::mymax ((size_t) 8, nbins);
                  thehist = new TH1F (&(name[0]), &(name[0]), nbins, hist_begin, hist_end);
                for (size_t ii = 0; ii < nums; ii++) {
                    thehist->Fill (vals[ii], weights[ii]);
                } vals.clear ();

                weights.clear ();
                converted = true;
            }
            return converted;
        }
        void write_histogram () {
            if (!converted) {
                convert_hist ();
            }
            if (converted) {
                thehist->Write ();
                delete thehist;
            } else {
                printf ("ERROR: Something went horribly wrong... %s\n", &(name[0]));
                vals.clear ();
            }
        }
        void fill (float val, float weight = 1.0) {
            if (converted) {
                thehist->Fill (val, weight);
            } else {
                vals.push_back (val);
                weights.push_back (weight);
                if (vals.size () > 10000) {
                    convert_hist ();
                }
            }
        }
      public:
        size_t N_Bins_Max;
      myhist (std::string _name):name (_name), converted (false) {
            vals.clear ();
            N_Bins_Max = 100;
        }
        ~myhist () {
            write_histogram ();
        }
        inline void Fill (float val, float weight = 1.0) {
            fill (val, weight);
        }
        inline std::string get_name () {
            return name;
        }
        bool force_convert_hist (size_t nbins, float hist_begin, float hist_end) {
            if (!converted) {
                std::sort (vals.begin (), vals.end ());
                thehist = new TH1F (&(name[0]), &(name[0]), nbins, hist_begin, hist_end);
                for (size_t ii = 0; ii < vals.size (); ii++) {
                    thehist->Fill (vals[ii]);
                }
                vals.clear ();
                converted = true;
            }
            return converted;
        }
    };
    class VectorHist {
      private:
        std::string vectorname;
      public:
        myhist pt, eta, phi, m;
        VectorHist (std::string name):pt (name + std::string (":pt")), eta (name + std::string (":eta")), phi (name + std::string (":phi")), m (name + std::string (":mass")), vectorname (name) {
        } ~VectorHist () {
        }
        void Fill (NewHEPHeaders::vector4 invector) {
            if (invector.pass ()) {
                pt.Fill (invector.pt ());
                eta.Fill (invector.eta ());
                phi.Fill (invector.phi ());
                m.Fill (invector.m ());
            }
        }
        void force_convert (double In_pt = 1400, double In_eta = 6.0, double In_phi = 6.3, double In_m = 1200) {
            pt.force_convert_hist (150, -0.1, In_pt);
            eta.force_convert_hist (120, -In_eta, In_eta);
            phi.force_convert_hist (80, -0.1, In_phi);
            m.force_convert_hist (80, -0.1, In_m);
        }
    };
    typedef std::vector < VectorHist > VectorHists;
    class MultiPlot {
      private:
        TFile rootfile;
        std::vector < TH1F * >hists;
        double maxy;
          std::string title;
      public:
          bool normalized;
          MultiPlot (std::string rootfilename, std::string _title):rootfile (&(rootfilename[0])), title (_title) {
            normalized = false;
            maxy = 0;
        } ~MultiPlot () {
            for (size_t i = 0; i < hists.size (); i++) {
                delete hists[i];
            }
        }
        void operator  () (std::string histname) {
            TH1F *tmphist = (TH1F *) rootfile.Get (&(histname[0]));

            if (tmphist != NULL) {
                if (normalized) {
                    double factor = tmphist->Integral ();

                    tmphist->Scale (1.0 / factor);
                }
                int binmax = tmphist->GetMaximumBin ();
                double x = tmphist->GetBinContent (binmax);

                maxy = CPPFileIO::mymax (x, maxy);
                tmphist->SetTitle (&(title[0]));
                hists.push_back (tmphist);
            }
        }
        void Save (std::string dirname = "./") {
            std::vector < int >goodcolors; /* Decide on good colours: */  {
                goodcolors.push_back (kBlack);
                goodcolors.push_back (kRed);
                goodcolors.push_back (kBlue);
                goodcolors.push_back (kGreen + 3);
                goodcolors.push_back (kMagenta + 2);
            }
/* The plotting part: */  {
                TCanvas C;

                for (size_t i = 0; i < hists.size (); i++) {
                    TH1F *thehist = hists[i];

                    thehist->SetLineColor (goodcolors[i]);
                    thehist->SetMaximum (maxy);
                    thehist->Draw ("same:h");
                }
/* Save the canvas: */  {
                    mkdir (&(dirname[0]), 0755);
                    std::string outname = dirname + "/" + title + ".pdf";
                    C.SaveAs (&(outname[0]));
                }
            }
        }
    };
}

#ifdef USE_RNN
#include "MyMVASDTopTagger.cc"
namespace RNN {

    const size_t GFS       =    5 ;
    const size_t GOS       =    7 ;
    const size_t GQS       =   18 ;
    const size_t storesize =  800 ;
    const bool   UseRNN    = true ;

    template <typename T> inline T ProcessAngle (T & th) {
        while ( th >  NewHEPHeaders::CONSTANTS::PI ) { th = th - (NewHEPHeaders::CONSTANTS::PI2) ; }
        while ( th < -NewHEPHeaders::CONSTANTS::PI ) { th = th + (NewHEPHeaders::CONSTANTS::PI2) ; }
    }
    class PreProcessor {
    private:
    public:
        float Pt, Eta, Phi, M ;
        inline void ProcessJet (fastjet::PseudoJet & injet, float * Output,bool doprocess) {
            float pt , eta , phi , m ;
            if (doprocess) {
                pt  = injet.pt  () / Pt  ; Output [0] = pt  ;
                eta = injet.eta () - Eta ; Output [1] = eta ;
                phi = injet.phi () - Phi ; ProcessAngle (phi) ; Output[2] = phi ;
                m   = injet.m   () / M   ; Output [3] = m   ;
            } else if (true)  {
                pt  = injet.pt  () ; Output [0] = pt  ;
                eta = injet.eta () ; Output [1] = eta ;
                phi = injet.phi () ; ProcessAngle (phi) ; Output [2] = phi   ;
                m   = injet.m   () ; Output [3] = m   ;
            } else if (false) {
                pt  = injet.pt  () ; Output [0] = pt  ;
                eta = injet.eta () ; Output [1] = eta ;
                phi = injet.phi () ; Output [2] = phi ;
                m   = injet.m   () ; Output [3] = m   ;
            }
        }
        inline void init (fastjet::PseudoJet & injet) {
            Pt  = injet.pt  () ;
            Eta = injet.eta () ;
            Phi = injet.phi () ;
            M   = injet.m   () + NeuralNetworks::epsG ;
        }
        inline void operator () (fastjet::PseudoJet & injet, float * Output, bool doprocess) {ProcessJet(injet,Output,doprocess);}
        PreProcessor(){}
        ~PreProcessor(){}
    };
    class StoreNode {
    private:
    public:
        bool   lastnode         ;
        float  Input    [ GOS ] ;
        size_t Next     [ 2   ] ;
        inline void debugshow (int index) {
            printf("Showing Node: (%d) {\n",index);
            for (size_t i=0;i<GOS;i++) { printf("\t %ld = %e;\n",i,Input[i]); }
            printf("}\n");
        }
        StoreNode  () {}
        ~StoreNode () {}
    };
    class MainStore {
    private:
    public:
        HEPTopTagger::TopTagger MainTopTagger     ;
        StoreNode               store [storesize] ;
        size_t                  count             ;
        inline StoreNode & operator [] (size_t i) {return store[i];}
        inline void clear(){MainTopTagger();}
        inline void debugshow () {for(size_t i=0;i<storesize;i++){store[i].debugshow(i);}}
        MainStore  () {}
        ~MainStore () {}
    } ;

    class AllParameters {
    private:
        inline void Randomize (CPPFileIO::myrandgen<pcg32_fast>&engine) {
            pars(engine);
            for ( size_t i=0 ; i < weight.size   () ; i++ ) { weight [i] = engine [0] ; }
            for ( size_t i=0 ; i < bias.size     () ; i++ ) { bias   [i] = engine [0] ; }
        }
        inline void WriteToFile (std::string File="./OutParameters") {
            FILE*f=fopen(&(File[0]),"w"); pars(f);
            for ( size_t j=0 ; j<weight.size () ; j++ ) { fprintf ( f , "weight[%ld]=%e;" , j , weight [j] ) ; }
            fprintf (f,"\n") ;
            for ( size_t j=0 ; j<bias.size   () ; j++ ) { fprintf ( f , "bias[%ld]=%e;"   , j , bias   [j] ) ; }
            fprintf (f,"\n") ;
            fclose(f);
        }
        inline void initialize () {
            {OS=GOS;FS=GFS;QS=GQS;}
            {weight.resize(FS,(3*OS)+(2*FS));bias.resize(FS);}
            size_t sizes[4]={QS+FS+OS,10,5,1};
            for(size_t i=0;i<4;i++) {pars.sizes[i]=sizes[i];} pars();
            Randomize();
            #include "OutParameters"
        }
        inline void Randomize () {CPPFileIO::myrandgen<pcg32_fast>engine(1,-1.0,1.0);Randomize(engine);}
    public:
        NeuralNetworks::Parameters <float,3> pars ;
        NeuralNetworks::GoodMatrix weight ; NeuralNetworks::GoodVector bias   ;
        size_t OS, FS, QS;
        inline void operator () (std::string filename) {WriteToFile(filename);}
        inline void operator () () {Randomize();}
        AllParameters  () {initialize();}
        ~AllParameters () {}
    };

    template <typename TP> class JetNN {
    private:
        // The private data elements:
        bool                         lastnode    ;
        NeuralNetworks::GoodVector   Input       ;
        NeuralNetworks::GoodVector   Diffs       ;
        size_t                       Next[2]     ;
        float                      * TrainDeltas ;

        // Important inline accessing methods:
        inline size_t OS () {return parent->OS();}
        inline size_t FS () {return parent->FS();}
        inline float weight (size_t y, size_t x) {return parent->weight(y,x);}
        inline NeuralNetworks::GoodMatrix & weight () {return parent->weight();}
        inline float bias (size_t y) {return parent->bias(y);}
        inline NeuralNetworks::GoodVector & bias() {return parent->bias();}
        inline JetNN <TP> & Next0 () { return parent[0][Next[0]] ; }
        inline JetNN <TP> & Next1 () { return parent[0][Next[1]] ; }
        inline size_t sizeX(){return weight().sizeX();}
        inline size_t sizeY(){return weight().sizeY();}
        size_t list (size_t i) {return parent->Boundary(i);}
        inline float & SelfInput (size_t i) {
            const size_t index = 2 * ( FS() + OS() ) ;
            return Input(index+i);
        }
        inline void SetupInputs () {
            if(!lastnode) {
                Next0().Output=&(Input(0))         ;
                Next1().Output=&(Input(FS()+OS())) ;
            }
        }
        inline void EventInit () {lastnode=true;Input.resize(sizeX());Input.ZeroAll();}
        inline void ReadFromStore (StoreNode & instore) {
            EventInit();
            lastnode=instore.lastnode;
            Next[0]=instore.Next[0]; Next[1]=instore.Next[1];
            SetupInputs();
            for (size_t i=0;i<OS();i++) {SelfInput(i)=instore.Input[i];}
        }

        // The part to evaluate derivatives:
        inline float PI (size_t i) {return Diffs(i)*TrainDeltas[i];}
        inline float DI (size_t k, size_t l, size_t i) {
            if (lastnode) { return 0 ; }
            else {
                if ((list(0)<=i)&&(i<list(1))) { return 0 ; }
                if ((list(1)<=i)&&(i<list(2))) { i=i-list(1); return Next0().BackDer(k,l,i); }
                if ((list(2)<=i)&&(i<list(3))) { return 0 ; }
                if ((list(3)<=i)&&(i<list(4))) { i=i-list(3); return Next1().BackDer(k,l,i); }
                else                           { return 0 ; }
            }
        }
        inline float DI (size_t l, size_t i) {
            if (lastnode) { return 0 ; }
            else {
                if ((list(0)<=i)&&(i<list(1))) { return 0 ; }
                if ((list(1)<=i)&&(i<list(2))) { i=i-list(1); return Next0().BiasDer(l,i); }
                if ((list(2)<=i)&&(i<list(3))) { return 0 ; }
                if ((list(3)<=i)&&(i<list(4))) { i=i-list(3); return Next1().BiasDer(l,i); }
                else                           { return 0 ; }
            }
        }
        inline float DO (size_t k, size_t l, size_t i) {
            float ret=0;
            if (k==i) {ret=Input(l);}
            for (size_t j=0;j<weight().sizeX();j++) { ret = ret + (weight(i,j)*DI(k,l,j)) ; }
            return ret*Diffs(i) ;
        }
        inline float DO (size_t l, size_t i) {
            float ret=0;
            if (l==i) {ret=1.0;}
            for (size_t j=0;j<weight().sizeX();j++) { ret = ret + (weight(i,j)*DI(l,j)) ; }
            return ret*Diffs(i) ;
        }

        // The activation part:
        inline void activate (bool uselru=true) {
            NeuralNetworks::GoodVector OutBuf; OutBuf(weight(),Input,bias());
            if (uselru) { NeuralNetworks::GoodSoftLRU  tmpslave; tmpslave(OutBuf,Diffs); }
            else        { NeuralNetworks::GoodSoftSign tmpslave; tmpslave(OutBuf,Diffs); }
            const size_t mark = 2*(OS()+FS()) ;
            for ( size_t i=0 ; i<OS         () ; i++ ) { Output [i]      = Input  (i+mark) ; }
            for ( size_t i=0 ; i<OutBuf.size() ; i++ ) { Output [i+OS()] = OutBuf (i)      ; }
        }
        inline void SetupStructer (fastjet::PseudoJet & _jet,bool DoPreProcess=true) {
            EventInit () ;
            fastjet::PseudoJet pjet1(0,0,0,0), pjet2(0,0,0,0);
            if (_jet.validated_cs()->has_parents(_jet,pjet1,pjet2)) {
                lastnode=false; if (pjet1.pt()<pjet2.pt()) {std::swap(pjet1,pjet2);}
                Next[0]=parent[0](); Next[1]=parent[0]();
                SetupInputs();
                Next0()(pjet1);Next1()(pjet2);
                SelfInput(4)=pjet1.delta_R(pjet2);
                float tmp_mass = _jet.m() ;
                if (tmp_mass>0.01) {SelfInput(5)=(pjet1.m()+pjet2.m())/tmp_mass;}
                else {SelfInput(5)=0;}
                SelfInput(6)=CPPFileIO::mymax(pjet1.E(),pjet2.E())/_jet.E();
            }
            parent->processer(_jet,&(SelfInput(0)),DoPreProcess);
        }
        inline void ChainActivate (bool uselru=true) {
            if (!lastnode) {Next0().ChainActivate();Next1().ChainActivate();}
            activate(uselru);
        }
    public:

        // The public pointers and derivative data elements:
        float                      * Output      ;
        TP                         * parent      ;
        NeuralNetworks::GoodMatrix   BackDer     ;
        NeuralNetworks::GoodVector   BiasDer     ;

        // Important public functions:
        inline void EvalDerivative (float * _TrainDeltas)                           {
            TrainDeltas = _TrainDeltas ;
            if(!lastnode){
                NeuralNetworks::GoodVector Xi ; Xi.resize(FS()) ;
                Xi.ZeroAll() ;
                for(size_t i=0;i<sizeY();i++) for(size_t jj=list(1);jj<list(2);jj++){
                    size_t j=jj-list(1);
                    Xi(j)=Xi(j)+(PI(i)*weight(i,jj));
                }
                Next0().EvalDerivative(&(Xi(0)));
                Xi.ZeroAll() ;
                for(size_t i=0;i<sizeY();i++) for(size_t jj=list(3);jj<list(4);jj++){
                    size_t j=jj-list(3);
                    Xi(j)=Xi(j)+(PI(i)*weight(i,jj));
                }
                Next1().EvalDerivative(&(Xi(0)));
            }
            BackDer.resize (sizeY(),sizeX()) ; BiasDer.resize (sizeY()) ;
            for(size_t m=0;m<BackDer.sizeY();m++) for(size_t n=0;n<BackDer.sizeX();n++) {
                BackDer(m,n)=PI(m)*Input(n);
                if(!lastnode)
                {BackDer(m,n)=BackDer(m,n)+Next0().BackDer(m,n)+Next1().BackDer(m,n);}
            }
            for(size_t m=0;m<BiasDer.size();m++){BiasDer(m)=PI(m);}
        }
        inline StoreNode GetStore  ()                                               {
            StoreNode ret ;
            ret.lastnode=lastnode;
            for (size_t i=0;i<OS();i++) {ret.Input[i]=SelfInput(i);}
            ret.Next[0]=Next[0] ;
            ret.Next[1]=Next[1] ;
            return ret;
        }
        inline void operator ()    (bool uselru=true)                               {ChainActivate(uselru);}
        inline void operator ()    (StoreNode&instore)                              {ReadFromStore(instore);}
        inline void operator ()    (fastjet::PseudoJet&_jet,bool DoPreProcess=true) {SetupStructer(_jet,DoPreProcess);}

        // The constructors and destructors (Nothing to be done here...) :
        JetNN  () {}
        ~JetNN () {}
    };
    template <typename TDW, typename TDB> class JetStore {
    private:

        // Private Data Elements:
        size_t                                   list[6] ;
        AllParameters                          * pars    ;
        size_t                                   count   ;
        std::vector <JetNN<JetStore<TDW,TDB>>>   Store   ;
        TDW                                      DeltaW  ;
        TDB                                      DeltaB  ;
        size_t                                   EvtCnt  ;

        // Private Access Functions:
        inline void GetStores  (MainStore&stores)     {
            stores.count=count;
            for(size_t i=0;i<count;i++){stores[i]=Store[i].GetStore();}
        }
        inline void initialize (AllParameters&inpars) {
            pars   = & inpars          ;
            count  =   0               ;
            TrainDeltas.resize (FS())  ;
            /* List the boundaries: */ {
                list[0] = 0              ;
                list[1] = OS()           ;
                list[2] = list[1] + FS() ;
                list[3] = list[2] + OS() ;
                list[4] = list[3] + FS() ;
                list[5] = list[4] + OS() ;
            }
            { DeltaW(weight()); DeltaB(bias()); Store.resize(storesize); }
            for (size_t i=0;i<Store.size();i++) {Store[i].parent=this;}
        }
        inline size_t getnew   () {
            size_t ret = count ;
            count++            ;
            return ret         ;
        }
        inline void AnalyzeJet (fastjet::PseudoJet & _jet) {
            Output.resize(OS()+FS());
            processer.init(_jet);
            count=0;
            size_t ret=getnew();
            Store[ret].Output=&(Output[0]);
            Store[ret](_jet,false);
            activate();
        }
        inline void Train (size_t y, size_t x) {DeltaW(y,x,D(y,x));}
        inline void Train (size_t y)           {DeltaB(y,D(y));}
        inline void TrainW () {for(size_t y=0;y<weight().sizeY();y++)for(size_t x=0;x<weight().sizeX();x++){Train(y,x);}}
        inline void TrainB () {for(size_t y=0;y<bias().size();y++){Train(y);}}
        inline void apply (float eta) {
            if(EvtCnt>0){
                DeltaW ((float)eta/EvtCnt) ;
                DeltaB ((float)eta/EvtCnt) ;
                EvtCnt = 0                 ;
            }
        }

    public:

        // The important data elements:
        NeuralNetworks::GoodVector     Output      ;
        NeuralNetworks::GoodVector     TrainDeltas ;
        PreProcessor                   processer   ;

        // Important public functions:
        inline float & weight (size_t x, size_t y)    { return pars->weight (x,y) ; }
        inline float & bias   (size_t x)              { return pars->bias   (x)   ; }
        inline NeuralNetworks::GoodMatrix & weight () { return pars->weight       ; }
        inline NeuralNetworks::GoodVector & bias   () { return pars->bias         ; }
        inline size_t   OS          ()         { return pars->OS               ; }
        inline size_t   FS          ()         { return pars->FS               ; }
        inline size_t   Boundary    (size_t i) { return list   [i]             ; }
        inline size_t   size        ()         { return Output.size ()         ; }
        inline void     activate    ()         { Store[0]           ()         ; }
        inline void     Train       ()         {
            Store[0].EvalDerivative (&(TrainDeltas(0))) ;
            TrainB () ; TrainW () ;
            EvtCnt++ ;
        }
        inline float  & GetD        (size_t i) { return TrainDeltas (i)      ; }
        inline NeuralNetworks::GoodVector & GetD () { return TrainDeltas     ; }
        inline float D (size_t y, size_t x)  { return Store[0].BackDer (y,x) ; }
        inline float D (size_t x)            { return Store[0].BiasDer (x)   ; }
        inline void     ReadMainStore (MainStore&instore) {
            Output.resize (OS()+FS()) ;
            Store[0].Output = & (Output[0]) ;
            count = instore.count ;
            for (size_t i=0;i<count;i++) {Store[i](instore.store[i]);}
            activate();
        }
        inline JetNN<JetStore<TDW,TDB>> & operator [] ( size_t i                  ) {return Store[i];}
        inline void                       operator () ( float eta                 ) {apply(eta);}
        inline void                       operator () ( MainStore & fullstore     ) {GetStores(fullstore);}
        inline void                       operator () ( fastjet::PseudoJet & _jet ) {AnalyzeJet(_jet);}
        inline void                       operator () ( AllParameters & inpars    ) {initialize(inpars);}
        inline size_t                     operator () ()                            {return getnew();}

        JetStore(AllParameters&inpars) {initialize(inpars);}
        JetStore(){}
        ~JetStore(){}
    };

    typedef JetStore <NeuralNetworks::GoodAdaMaxMatrix,NeuralNetworks::GoodAdaMaxVector> JetRNN ;

    class CNNTopTagger {
    private:

        // Important private data elements:
        float                      AvgError ;
        size_t                     count    ;
        NeuralNetworks::GoodFCNN   network  ;
        AllParameters            * pars     ;

        // Important wrapper functions:
        inline size_t QS () {return pars->QS;}
        inline size_t OS () {return pars->OS;}
        inline size_t FS () {return pars->FS;}
        inline void initialize (AllParameters*_pars) { pars = _pars ; network(pars->pars); }
        inline void initialize (AllParameters&_pars) {initialize(&_pars);}
        HEPTopTagger::TopTagger MainTopTagger ;
        inline void Analyze (HEPTopTagger::TopTagger&MainTagger) {
            for (size_t i=0;i<4;i++) {network.Input()(i)=MainTagger.EFCR[i];}
            network.Input()(4)  = MainTagger.tau1tau2       ;
            network.Input()(5)  = MainTagger.tau2tau3       ;
            network.Input()(6)  = MainTagger.tau3tau4       ;
            network.Input()(7)  = MainTagger.DeltaMTop      ;
            network.Input()(8)  = MainTagger.PrunedMass     ;
            network.Input()(9)  = MainTagger.frec           ;
            network.Input()(10) = MainTagger.EnergyFraction ;
            network.Input()(11) = MainTagger.NumBs          ;
            network.Input()(12) = MainTagger.m12            ;
            network.Input()(13) = MainTagger.m23            ;
            network.Input()(14) = MainTagger.m13            ;
            network.Input()(15) = MainTagger.m123           ;
            network.Input()(16) = MainTagger.LSDChi         ;
            network.Input()(17) = MainTagger.ROpt           ;
            MainTopTagger       = MainTagger                ;
            network();
        }
        inline void train (float _ans1) {
            NeuralNetworks::GoodVector ans; ans.resize(1);
            ans(0)=_ans1;
            float delta = CPPFileIO::mymod (network.Output()(0)-_ans1) ;
            network(ans); count++;
            if (CPPFileIO::mymod(delta)>0.9) {AvgError=AvgError+1.0;}
        }

    public:

        // Important public functions:
        inline NeuralNetworks::GoodVector & Output () {return network.Output();}
        inline float   ApplyDelta       (float eta) {
            if(count>0) {
                float ret = (float) AvgError / count ;
                network(eta); count=0; AvgError=0;
                return ret;
            } else {
                printf("No Event Analyzed...\n");
                return 1.0 ;
            }
        }
        inline void    activate         ()          {network();}
        inline void    EvalChain        ()          {network.EvalChain();}
        inline float & getChainDelta    (size_t i)  {
            size_t max = network.Input().size() - 1 - FS() ;
            return network.C(0,max+i);
        }

        // Important accessing functions:
        inline float & operator [] (size_t i)                           {return network.Input()(QS()+i);}
        inline float   operator () ()                                   {return network.Output()(0);}
        inline void    operator () (float _ans1)                        {train(_ans1);}
        inline void    operator () (HEPTopTagger::TopTagger&MainTagger) {Analyze(MainTagger);}
        inline void    operator () (MainStore&store)                    {store.MainTopTagger=MainTopTagger;}
        inline void    operator () (AllParameters*_pars)                {initialize(_pars);}
        inline void    operator () (AllParameters&_pars)                {initialize(_pars);}

        // Constructors and destructorss:
        CNNTopTagger (AllParameters*_pars) { initialize(_pars); }
        CNNTopTagger (AllParameters&_pars) { initialize(_pars); }
        CNNTopTagger(){}
        ~CNNTopTagger(){}
    };

    class FullTopTagger {
    private:
        AllParameters * par2       ;
        CNNTopTagger    toptagger1 ;
        JetRNN          toptagger2 ;
        size_t          EvtCnt     ;
        float           TrueError  ;
        inline void Train (float ans) {
            float nnout=Output()(0);
            float delta = CPPFileIO::mymod(nnout-ans);
            if (delta>0.9) {TrueError=TrueError+1.0;}
            toptagger1(ans);
            if((UseRNN)&&(Rate_RNN>0.0)){
                for(size_t i=0;i<toptagger2.TrainDeltas.size();i++)
                {toptagger2.TrainDeltas(i)=toptagger1.getChainDelta(i);}
                toptagger2.Train();
            }
            EvtCnt++;
        }
        inline void RNN2CNN () {
            for (size_t i=0;i<toptagger2.Output.size();i++) {toptagger1[i]=toptagger2.Output(i);}
            if(false){
                printf("Begin transferring outputs RNN -> CNN {\n");
                for (size_t i=0;i<toptagger2.Output.size();i++)
                {printf("Transferring output %ld = %e\n",i,toptagger2.Output(i));}
                printf("}\n");
            }
        }
        inline void initialize (AllParameters&_par2) {
            toptagger1(&_par2);
            toptagger2(_par2);
            par2=&_par2;
            Rate_RNN=0.01;Rate_CNN=0.1;
            EvtCnt=0;TrueError=0;
        }
    public:
        float Rate_RNN , Rate_CNN ;

        inline NeuralNetworks::GoodVector & Output () {return toptagger1.Output();}
        inline float ApplyDelta (float rate=1.0) {
            float ret;
            if (Rate_CNN>0.0) {toptagger1.ApplyDelta(Rate_CNN*rate);}  else {toptagger1.ApplyDelta(0.0);}
            if ((UseRNN)&&(Rate_RNN>0.0)) {toptagger2(Rate_RNN*rate);} else {toptagger2(0.0);}
            ret=(float)TrueError/EvtCnt;
            EvtCnt=0;TrueError=0;
            return ret;
        }
        inline void operator () (float ans) {Train(ans);}
        inline void operator () (MainStore&store) {
            store.clear(); toptagger1(store);
            if (UseRNN)   {toptagger2(store);}
        }
        inline void putstore    (MainStore&store) {
            if (UseRNN) {
                toptagger2.ReadMainStore(store);
                RNN2CNN();
            }
            toptagger1(store.MainTopTagger);
        }
        inline void operator () (fastjet::PseudoJet & injet, NewHEPHeaders::vector4s & bs) {
            /* The RNN Part: */ if (UseRNN) {toptagger2(injet);RNN2CNN();}
            /* The fully connected part: */ if (true) {
                HEPTopTagger::TopTagger thetagger;
                {thetagger();thetagger(injet);thetagger(bs);}
                toptagger1(thetagger);
            }
        }
        inline void operator () (AllParameters&_par2) {initialize(_par2);}

        FullTopTagger (AllParameters&_par2) {initialize(_par2);}
        FullTopTagger(){}
        ~FullTopTagger(){}
    } ;
}
#endif

#endif
