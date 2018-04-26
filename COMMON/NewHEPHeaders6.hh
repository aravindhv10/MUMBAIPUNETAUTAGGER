#ifndef NewHEPHeaders2_HH
#define NewHEPHeaders2_HH

#include "CPPFileIO2.hh"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include <fastjet/Selector.hh>
#include <fastjet/tools/JHTopTagger.hh>

namespace NewHEPHeaders {
    const bool DEBUG = false;
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
        if (pid == PID::ELECTRON_NU) {return false;}
        if (pid == PID::MUON_NU) {return false;}
        if (pid == PID::TAU_NU) {return false;}
        if (pid == PID::CHI10) {return false;}
        if (pid == PID::CHI20) {return false;}
        return true;
    }
    inline bool islepton (long pid) {
        pid = CPPFileIO::mymod (pid);
        if (pid == PID::ELECTRON) {return true;}
        if (pid == PID::MUON) {return true;}
        if (pid == PID::TAU) {return true;}
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
        template < typename TR=double > class plane2vector   {
        private:
            void construct (TR _x = 0, TR _y = 0) {x[0] = _x;x[1] = _y;}
        public:
            TR x[2];
            inline void SetPtPhi (const TR _pt=0,const TR _phi=0) {x[0]=_pt*cos(_phi);x[1]=_pt*sin(_phi);}
            inline TR pt2 () const { return (x[0] * x[0]) + (x[1] * x[1]); }
            inline TR pt () const { return (sqrt (pt2())); }
            inline TR safenorm2 () const {
                TR mag = pt2 ();
                if (CPPFileIO::mymod (mag) < 0.0000000001) {mag = CPPFileIO::mysign (mag) * 0.0000000001;}
                return mag;
            }
            inline TR phi () const {
                TR ret = acos (x[0] / sqrt (safenorm2 ()));
                if (x[1] < 0) {ret = CONSTANTS::PI2 - ret;}
                return ret;
            }
            inline TR dphi (const plane2vector < TR > b) const {
                TR ret = CPPFileIO::mymod (b.phi () - phi ());
                if (ret > CONSTANTS::PI) { ret = CONSTANTS::PI2 - ret; }
                return ret;
            }
            inline plane2vector < TR > operator + (const plane2vector < TR > b) const
            { return plane2vector < TR > (x[0] + b.x[0], x[1] + b.x[1]); }

            inline plane2vector < TR > operator - (const plane2vector < TR > b) const
            { return plane2vector < TR > (x[0] - b.x[0], x[1] - b.x[1]); }

            inline TR operator * (const plane2vector < TR > b) const
            { return (x[0] * b.x[0]) + (x[1] * b.x[1]); }

            inline plane2vector <TR> operator * (const TR b) const
            { return plane2vector <TR> (x[0] * b, x[1] * b); }

            inline plane2vector <TR> operator / (const TR b) const
            { return plane2vector <TR> (x[0] / b, x[1] / b); }

            inline TR operator  () (const plane2vector <TR> b) const { return dphi (b); }

            inline plane2vector <TR> flip () const
            { return plane2vector <TR> (-x[0], -x[1]); }

            inline plane2vector <TR> dir () const {
                plane2vector <TR> ret (x[0], x[1]);
                TR mag = sqrt (ret.safenorm2 ());
                ret = ret / mag;
                return ret;
            }

            inline bool operator > (const plane2vector < TR > b) const { return pt2 () > b.pt2 (); }
            inline bool operator < (const plane2vector < TR > b) const { return pt2 () < b.pt2 (); }
            inline ssize_t operator >> (CPPFileIO::FileFD & f) const { return f.multiwrite2file (*this); }
            inline ssize_t operator << (CPPFileIO::FileFD & f) const { return f.multiread2file (*this); }

            inline TR & operator [] (size_t i)       { return x[i]; }
            inline TR   operator [] (size_t i) const { return x[i]; }

            inline void clearthis () { x[0] = 0; x[1] = 0; }
            inline bool operator == (const plane2vector < TR > b) const {
                plane2vector < TR > tmp (x[0], x[1]);
                tmp = tmp - b;
                TR diff = tmp.pt2 ();
                diff = CPPFileIO::mymod (diff);
                return diff < VECTOR_EQUALITY_LIMIT;
            }

            plane2vector (const TR _x = 0, const TR _y = 0) {construct (_x, _y);}
            plane2vector (const plane2vector < TR > &c) {construct(c.x[0], c.x[1]);}
            ~plane2vector () {}
        };
        template < typename TR=double > class euclid3vector  {
        private:
        public:
            plane2vector < TR > xy; TR z;
            inline void SetPtEtaPhi(const TR _pt=0, const TR _eta=0, const TR _phi=0){
                xy.SetPtPhi (_pt,_phi) ;
                TR K = 2.0 * _eta ;
                K = exp (K) ;
                K = (K-1.0) / (K+1.0) ;
                K = K * K ;
                K = K / ( 1.0 - K ) ;
                z = sqrt(_pt*_pt*K) * CPPFileIO::mysign(_eta) ;
            }
            inline TR phi () const {return xy.phi ();}
            inline TR pt2 () const {return xy.pt2 ();}
            inline TR pt () const {return xy.pt ();}
            inline euclid3vector < TR > operator + (const euclid3vector < TR > a) const
            { return euclid3vector < TR > (xy + a.xy, z + a.z); }

            inline euclid3vector < TR > operator - (const euclid3vector < TR > a) const
            { return euclid3vector < TR > (xy - a.xy, z - a.z); }

            inline TR operator * (const euclid3vector <TR> a) const
            { return (xy * a.xy) + (z * a.z); }

            inline euclid3vector <TR> operator * (const TR a) const
            { return euclid3vector < TR > (xy * a, z * a); }

            inline euclid3vector <TR> operator / (const TR a) const
            { return euclid3vector < TR > (xy / a, z / a); }

            inline TR p2 () const { return xy.pt2 () + (z * z); }
            inline TR p () const { return sqrt (p2 ()); }

            inline TR & operator [] (size_t ret) {
                if (ret > 1) { return z; }
                else { return xy[ret]; }
            }
            inline TR operator [] (size_t ret) const {
                if (ret > 1) { return z; }
                else { return xy[ret]; }
            }

            inline TR eta () const
            { TR tmp_p = p (); return 0.5 * log ((tmp_p + z) / (tmp_p - z)); }

            inline TR meta () const { return CPPFileIO::mymod (eta ()); }

            inline TR cone2 (const euclid3vector < TR > b) const {
                TR tphi = xy.dphi (b.xy);
                tphi = tphi * tphi;
                TR teta = eta () - b.eta ();
                teta = teta * teta;
                TR ret = teta + tphi;
                return ret;
            }
            inline TR cone (const euclid3vector < TR > b) const
            { return sqrt (cone2 (b)); }

            inline TR dphi (const euclid3vector < TR > b) const
            { TR tphi = xy.dphi (b.xy); return tphi; }

            inline TR operator  () (const euclid3vector <TR> b) const { return cone (b); }
            inline TR safenorm2 () const {
                TR mag = xy.pt2 () + (z * z);
                if (CPPFileIO::mymod (mag) < 0.0000000001)
                { mag = CPPFileIO::mysign (mag) * 0.0000000001; }
                return mag;
            }
            inline euclid3vector < TR > flip () const { return euclid3vector < TR > (xy.flip (), -z); }
            inline euclid3vector < TR > trans () const { return euclid3vector < TR > (xy, 0); }
            inline euclid3vector < TR > dir () const {
                euclid3vector < TR > ret (*this);
                TR mag = sqrt (ret.safenorm2 ());
                ret = ret / mag;
                return ret;
            }
            inline bool operator > (const euclid3vector < TR > b) const { return pt2 () > b.pt2 (); }
            inline bool operator < (const euclid3vector < TR > b) const { return pt2 () < b.pt2 (); }
            inline ssize_t operator >> (CPPFileIO::FileFD & f) const { return f.multiwrite2file (*this); }
            inline ssize_t operator << (CPPFileIO::FileFD & f) const { return f.multiread2file (*this); }
            inline void clearthis () { xy = plane2vector < TR > (0, 0); z = 0; }
            inline bool operator == (const euclid3vector < TR > b) const {
                euclid3vector < TR > tmp = (*this) - b;
                TR diff = tmp.p2 ();
                diff = CPPFileIO::mymod (diff);
                return diff < VECTOR_EQUALITY_LIMIT;
            }
            euclid3vector (const TR _x = 0, const TR _y = 0, const TR _z = 0):xy (_x, _y) {z = _z;}
            euclid3vector (const plane2vector < TR > a, const TR _z = 0):xy (a) {z = _z;}
            euclid3vector (const euclid3vector < TR > &a):xy (a.xy) {z = a.z;}
        };
        template < typename TR=double > class lorentz4vector {
        private:
        public:
            euclid3vector < TR > xyz; TR t;
            inline void SetPtEtaPhiM(const TR _pt, const TR _eta, const TR _phi, const TR _m){
                xyz.xy.SetPtPhi(_pt,_phi);
                TR K = exp(2.0*_eta) ;
                K=(K+1.0)/(K-1.0);
                TR _pt2 = _pt * _pt ;
                TR _m2 = _m*_m ;
                TR _z2 = (_pt2+_m2)/((K*K)-1.0) ;
                TR _E2 = _pt2 + _z2 + _m2 ;
                t = sqrt(_E2);
                xyz[2]=sqrt(_z2)*CPPFileIO::mysign(_eta);
            }

            inline TR & operator [] (const size_t ref)
            { if (ref > 2) { return t; } else { return xyz[ref]; } }

            inline TR operator [] (const size_t ref) const
            { if (ref > 2) { return t; } else { return xyz[ref]; } }

            inline TR pt2 () const { return xyz.pt2 (); }
            inline TR pt () const { return xyz.pt (); }
            inline TR p2 () const { return xyz.p2 (); }
            inline TR p () const { return xyz.p (); }
            inline TR phi () const { return xyz.phi (); }
            inline TR m2 () const { return CPPFileIO::mymod ((t * t) - p2 ()); }
            inline TR n2 () const { return CPPFileIO::mymod ((t * t) + p2 ()); }
            inline TR eta () const { return (0.5 * log ((t + xyz.z) / (t - xyz.z))); }
            inline TR meta () const { return CPPFileIO::mymod (eta ()); }
            inline TR peta () const { return xyz.eta (); }
            inline TR pmeta () const { return xyz.meta (); }
            inline TR m () const { return sqrt (m2 ()); }
            inline TR n () const { return sqrt (n2 ()); }
            inline TR dphi (const lorentz4vector < TR > b) const { return xyz.dphi (b.xyz); }
            inline lorentz4vector <TR> operator + (const lorentz4vector < TR > b) const
            { return lorentz4vector <TR> (xyz + b.xyz, t + b.t); }

            inline lorentz4vector <TR> operator - (const lorentz4vector < TR > b) const
            { return lorentz4vector <TR> (xyz - b.xyz, t - b.t); }

            inline TR operator * (const lorentz4vector <TR> b) const
            { return (t * b.t) - (xyz * b.xyz); }

            inline lorentz4vector < TR > operator * (const TR b) const
            { return lorentz4vector < TR > (xyz * b, t * b); }

            inline lorentz4vector < TR > operator / (const TR b) const
            { return lorentz4vector < TR > (xyz / b, t / b); }

            inline void operator = ( const euclid3vector < TR > other ) {
                xyz = other     ;
                t   = other.p() ;
            }
            inline void operator = ( const lorentz4vector < TR > other ) {
                xyz = other.xyz ;
                t   = other.t ;
            }

            inline TR pcone2 (const lorentz4vector < TR > b) const {
                TR tphi = xyz.dphi (b.xyz);
                tphi = tphi * tphi;
                TR teta = peta () - b.peta ();
                teta = teta * teta;
                TR ret = teta + tphi;
                return ret;
            }
            inline TR pcone (const lorentz4vector < TR > b) const { return sqrt (pcone2 (b)); }
            inline TR cone2 (const lorentz4vector < TR > b) const {
                TR tphi = xyz.dphi (b.xyz);
                tphi = tphi * tphi;
                TR teta = eta () - b.eta ();
                teta = teta * teta;
                TR ret = teta + tphi;
                return ret;
            }
            inline TR cone (const lorentz4vector < TR > b) const { return sqrt (cone2 (b)); }
            inline TR operator  () (const lorentz4vector < TR > b) const { return cone (b); }
            inline lorentz4vector < TR > flip () const
            { return lorentz4vector < TR > (xyz.flip (), t); }

            inline lorentz4vector < TR > trans () const
            { return lorentz4vector < TR > (xyz.trans (), t); }

            inline lorentz4vector < TR > dir () const
            { return lorentz4vector < TR > (xyz.dir (), t); }

            inline bool operator > (const lorentz4vector < TR > b) const { return pt2 () > b.pt2 (); }
            inline bool operator < (const lorentz4vector < TR > b) const { return pt2 () < b.pt2 (); }
            inline ssize_t operator >> (CPPFileIO::FileFD & f) const { return f.multiwrite2file (*this); }
            inline ssize_t operator << (CPPFileIO::FileFD & f) const { return f.multiread2file (*this); }
            inline bool cleared () const { return (t < 0); }
            inline bool pass () const { return (t > 0); }
            inline void clearthis () { t = -1; xyz = euclid3vector < TR > (0, 0, 0); }
            inline bool operator == (const lorentz4vector < TR > b) const {
                lorentz4vector < TR > tmp = (*this) - b;
                TR diff = tmp.n2 ();
                diff = CPPFileIO::mymod (diff);
                return diff < VECTOR_EQUALITY_LIMIT;
            }
            inline TR gamma () const { return (TR) t / m (); }
            inline TR gamma2 () const { return t * t / m2 (); }
            inline TR beta () const {
                TR g = (TR) gamma2 ();
                g = 1.0 / g;
                g = 1.0 - g;
                return sqrt (g);
            }
            inline euclid3vector <TR> Velocity () const { return xyz.dir () * beta (); }
            inline lorentz4vector <TR> LorentzBoost (const euclid3vector < TR > booster) const {
                euclid3vector < TR > parallel = booster * ((xyz * booster) / booster.p2 ());
                euclid3vector < TR > perpendicular = xyz - parallel;
                TR gm = booster.p2 ();
                gm = 1.0 - gm;
                gm = (1.0 / gm);
                gm = sqrt (gm);
                lorentz4vector < TR > ret;
                ret.t = gm * (t - (parallel * booster));
                parallel = (parallel - (booster * t)) * (TR) gm;
                ret.xyz = parallel + perpendicular;
                return ret;
            }
            inline lorentz4vector < TR > LorentzBoostGamma (const euclid3vector <TR> booster) const {
                TR gm2   = (TR) booster.p2 ()        ;
                TR gm    = (TR) sqrt       (gm2)     ;
                TR beta2 = (TR) 1.0    -   (1.0/gm2) ;
                TR beta  = (TR) sqrt       (beta2)   ;
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
            inline fastjet::PseudoJet getpseudojet () const { return fastjet::PseudoJet (xyz[0], xyz[1], xyz[2], t); }
            inline void operator = (const fastjet::PseudoJet & injet) {this[0]=lorentz4vector<TR>(injet);}
            lorentz4vector (const TR _x = 0, const TR _y = 0, const TR _z = 0, const TR _t = 0):xyz (_x, _y, _z) { t = _t; }
            lorentz4vector (const euclid3vector < TR > a, const TR _t = -1) :xyz (a) { if(_t<0) {t=a.p();} else {t=_t;} }
            lorentz4vector (const plane2vector  < TR > a, const TR _z = 0 , const TR _t = -1) :xyz (a)
            { if (_t<0) {t=a.pt();} else {t=_t;} xyz.z = _z ; }
            lorentz4vector (const lorentz4vector < TR > &a):xyz (a.xyz) { t = a.t; }
            lorentz4vector (const fastjet::PseudoJet & injet) {
                t      = injet.e  () ;
                xyz[2] = injet.pz () ;
                xyz[1] = injet.py () ;
                xyz[0] = injet.px () ;
            }
            ~lorentz4vector(){}
        };
        template < typename TR=double, typename TI=int > class DelphesVectors {
        private:
        public:
            lorentz4vector <TR> momentum      ;
            TI                  Charge        ;
            TR                  Eem, Ehad, Emu;
            inline void Set_Ehad_Fraction (const TR InFrac) { Ehad = momentum[3] * InFrac ; }
            inline void Set_Eem_Fraction  (const TR InFrac) { Eem  = momentum[3] * InFrac ; }
            inline void Set_Emu_Fraction  (const TR InFrac) { Emu  = momentum[3] * InFrac ; }
            inline void SetPtEtaPhiM      (const TR _pt, const TR _eta, const TR _phi, const TR _m)
            { momentum.SetPtEtaPhiM (_pt,_eta,_phi,_m) ; }

            inline TR & operator [] (const size_t ref) {
                if      ( ref <= 3 ) { return momentum[ref] ; }
                else if ( ref == 4 ) { return Eem           ; }
                else if ( ref == 5 ) { return Ehad          ; }
                else if ( ref == 6 ) { return Emu           ; }
            }
            inline TR operator [] (const size_t ref) const {
                if      ( ref <= 3 ) { return momentum[ref] ; }
                else if ( ref == 4 ) { return Eem           ; }
                else if ( ref == 5 ) { return Ehad          ; }
                else if ( ref == 6 ) { return Emu           ; }
            }

            inline TR pt2   () const { return momentum.pt2  () ; }
            inline TR pt    () const { return momentum.pt   () ; }
            inline TR p2    () const { return momentum.p2   () ; }
            inline TR p     () const { return momentum.p    () ; }
            inline TR phi   () const { return momentum.phi  () ; }
            inline TR m2    () const { return momentum.m2   () ; }
            inline TR n2    () const { return momentum.n2   () ; }
            inline TR eta   () const { return momentum.eta  () ; }
            inline TR meta  () const { return momentum.meta () ; }
            inline TR peta  () const { return momentum.peta () ; }
            inline TR pmeta () const { return momentum.meta () ; }
            inline TR m     () const { return momentum.m    () ; }
            inline TR n     () const { return momentum.n    () ; }
            inline TR dphi       (const DelphesVectors <TR,TI> b) const { return momentum.dphi (b.momentum) ; }
            inline TR operator * (const DelphesVectors <TR,TI> b) const { return momentum * b.momentum      ; }

            inline DelphesVectors <TR,TI> operator * (const TR b) const
            { return DelphesVectors <TR,TI> (momentum*b,Eem*b,Ehad*b,Emu*b,Charge); }

            inline lorentz4vector < TR > operator / (const TR b) const
            { return DelphesVectors <TR,TI> (momentum/b,Eem/b,Ehad/b,Emu/b,Charge); }

            inline void operator = ( const DelphesVectors <TR,TI> other )
            { momentum=other.momentum; Eem=other.Eem; Ehad=other.Ehad; Emu=other.Emu; Charge=other.Charge;}

            inline DelphesVectors <TR,TI> operator + (const DelphesVectors <TR,TI> b) const {
                return DelphesVectors <TR,TI> (
                    momentum+b.momentum,
                    Eem+b.Eem,
                    Ehad+b.Ehad,
                    Emu+b.Emu,
                    Charge+b.Charge
                );
            }
            inline DelphesVectors <TR,TI> operator - (const DelphesVectors <TR,TI> b) const {
                return DelphesVectors <TR,TI> (
                    momentum-b.momentum,
                    Eem-b.Eem,
                    Ehad-b.Ehad,
                    Emu-b.Emu,
                    Charge-b.Charge
                );
            }

            inline TR cone (const DelphesVectors <TR,TI> b) const { return momentum.cone(b.momentum); }
            inline TR operator () (const DelphesVectors <TR,TI> b) const { return momentum(b.momentum); }
            inline bool operator > (const DelphesVectors <TR,TI> b) const { return momentum.pt2 () > b.momentum.pt2 (); }
            inline bool operator < (const DelphesVectors <TR,TI> b) const { return momentum.pt2 () < b.momentum.pt2 (); }
            inline ssize_t operator >> (CPPFileIO::FileFD & f) const { return f.multiwrite2file (*this); }
            inline ssize_t operator << (CPPFileIO::FileFD & f) const { return f.multiread2file (*this); }
            inline bool cleared () const { return (momentum[3] < 0); }
            inline bool pass () const { return (momentum[3] > 0); }
            inline void clearthis () { momentum.clearthis(); Eem=-10000; Ehad=-10000; Emu=-10000; Charge=0; }
            inline TR gamma () const { return (TR) momentum[3] / momentum.m(); }
            inline TR gamma2 () const { return (TR) momentum[3] * momentum[3] / momentum.m2(); }
            inline TR beta () const {return momentum.beta();}
            inline euclid3vector < TR > Velocity () const { return momentum.Velocity(); }
            inline DelphesVectors <TR,TI> LorentzBoost (const euclid3vector < TR > booster) const {
                lorentz4vector<TR>ret=momentum.LorentzBoost(booster);
                TR ratio = ret[3] / momentum[3] ;
                TR _Eem=Eem*ratio, _Ehad=Ehad*ratio, _Emu=Emu*ratio;
                return DelphesVectors <TR,TI> (ret,_Eem,_Ehad,_Emu,Charge) ;
            }
            inline DelphesVectors <TR,TI> LorentzBoostGamma (const euclid3vector < TR > booster) const {
                lorentz4vector<TR>ret=momentum.LorentzBoostGamma(booster);
                TR ratio=ret[3]/momentum[3];
                TR _Eem=Eem*ratio, _Ehad=Ehad*ratio, _Emu=Emu*ratio;
                return DelphesVectors <TR,TI> (ret,_Eem,_Ehad,_Emu,Charge) ;
            }

            inline fastjet::PseudoJet getpseudojet () const
            { return fastjet::PseudoJet (momentum[0],momentum[1],momentum[2],momentum[3]) ; }

            DelphesVectors
            (const lorentz4vector<TR>_momentum, const TR _Eem=0, const TR _Ehad=0, const TR _Emu=0, const TI _Charge=0)
            { momentum=_momentum; Eem=_Eem; Ehad=_Ehad; Emu=_Emu; Charge=_Charge; }

            DelphesVectors
            (const TR _x=0, const TR _y=0, const TR _z=0, const TR _t=0, const TR _Eem=0, const TR _Ehad=0, const TR _Emu=0, const TI _Charge=0)
            { momentum=lorentz4vector<TR>(_x,_y,_z,_t); Eem=_Eem; Ehad=_Ehad; Emu=_Emu; Charge=_Charge; }

            ~DelphesVectors () {}
        } ;
        template < typename TRF=double, typename TRI=long > class ParticleNode   {
        private:
            lorentz4vector <TRF> momentum;
            TRI d1, d2, pid;
        public:

            inline TRI id () { return pid; }
            inline bool isFinal () { return (d1 == -1) && (d2 == -1); }
            inline TRI daughter1 () { return d1; }
            inline TRI daughter2 () { return d2; }
            inline TRF px () { return momentum[0]; }
            inline TRF py () { return momentum[1]; }
            inline TRF pz () { return momentum[2]; }
            inline TRF e () { return momentum[3]; }
            inline TRF pt () { return momentum.pt (); }
            inline TRF eta () { return momentum.eta (); }
            inline TRF modeta () { return CPPFileIO::mymod (eta()); }
            inline bool isDetectable () {
                bool ret = (pt () > 0.5) && (modeta () < 6.0);
                ret = ret && detectable (pid);
                return ret;
            }
            inline bool IsGood () { return isFinal () && isDetectable (); }
            inline bool IsLepton () { return islepton (pid); }
            inline bool IsBLike () { return isblike(pid); }
            inline bool IsBMeson () {
                TRI tmppid = CPPFileIO::mymod(pid) ;
                return ((tmppid>100)&&isblike(tmppid));
            }
            inline bool IsBQuakr () {
                TRI tmppid = CPPFileIO::mymod(pid) ;
                return (tmppid==PID::BOTTOM);
            }
            inline TRF operator [] (size_t i) { return momentum[i] ; }
            inline lorentz4vector <TRF> & getvec () { return momentum; }
            inline lorentz4vector <TRF> & operator () () { return momentum; }
            inline TRF operator  () (ParticleNode b) { return momentum (b.momentum); }
            inline TRF operator  () (lorentz4vector <TRF> b) { return momentum (b); }
            inline TRF pcone (ParticleNode b) { return momentum.pcone (b.momentum); }
            inline TRF pcone (lorentz4vector <TRF> b) { return momentum.pcone (b); }
            inline fastjet::PseudoJet getpseudojet () { return momentum.getpseudojet (); }

            ParticleNode () {d1=-1;d2= -1;pid = 0;momentum.clearthis ();}

            ParticleNode (TRF _x, TRF _y, TRF _z, TRF _t, TRI _d1, TRI _d2, TRI _pid) :
            momentum (_x, _y, _z, _t) { d1 = _d1;d2 = _d2;pid = _pid; }

            ParticleNode (const Pythia8::Particle & part) {
                momentum = lorentz4vector <TRF> (part.px (), part.py (), part.pz (), part.e ());
                pid = part.id ();
                if (part.isFinal ()) { d1 = -1; d2 = -1; }
                else { d1 = part.daughter1 (); d2 = part.daughter2 (); }
            }
            ~ParticleNode () {}
        };
        template < typename TR=double > class JetContainer : public lorentz4vector <TR> {
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
                    for (size_t j = 0; j < _constituents.size (); j++)
                    { constituents.push_back (_constituents[j].user_index ()); }
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
                if (cleared ()) { return 1000000000; }
                else { return sqrt (error_top () + error_w ()); }
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
            inline bool operator > (reco_top_had other) { return error () > other.error (); }
            inline bool operator < (reco_top_had other) { return error () < other.error (); }
            inline bool operator > (reco_top_lept other) { return error () > other.error (); }
            inline bool operator < (reco_top_lept other) { return error () < other.error (); }
            inline bool operator > (double other) { return error () > other; }
            inline bool operator < (double other) { return error () < other; }
            reco_top_lept () { clearthis (); }
            reco_top_lept ( vector4 & _bjet , vector4 & _lept , vector2 & _met )
            { construct (_bjet,_lept,_met) ; }

            reco_top_lept ( vector4 & _bjet , vector4 & _lept , vector4 & _met )
            { construct (_bjet,_lept,_met.xyz.xy) ; }
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
        WriteHepmc2Fifo (std::string _name):ascii_io (&(_name[0]), std::ios::out) { name = _name; }
        ~WriteHepmc2Fifo () {}
        inline void operator  () (Pythia8::Pythia & pythia) {
            HepMC::GenEvent * hepmcevt = new HepMC::GenEvent ();
            ToHepMC.fill_next_event (pythia, hepmcevt);
            ascii_io << hepmcevt;
            delete hepmcevt;
        }
    };
    void LHA2HEPMC (std::string lhafile, std::string hepmcfile) {
        WriteHepmc2Fifo outfile (hepmcfile);
        Pythia8::Pythia pythia;
        pythia.readString ("Beams:frameType = 4");
        std::string tmp ("Beams:LHEF = ");
        tmp = tmp + lhafile;
        pythia.readString (&(tmp[0]));
        pythia.init ();
        for (size_t iEvent = 0; (!pythia.info.atEndOfFile ()); iEvent++) if (pythia.next ())
        {outfile (pythia);}
        pythia.stat ();
    }
}
#endif
