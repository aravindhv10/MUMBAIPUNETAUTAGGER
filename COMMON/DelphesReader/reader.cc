#include "./NewHEPHeaders6.hh"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TString.h"
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "./Nsubjettiness/main.hh"
#include "./EnergyCorrelator/EnergyCorrelator.cc"
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
#include <fastjet/Selector.hh>
#include <fastjet/tools/JHTopTagger.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "TColor.h"

namespace Step1 {

    constexpr double epsilon = 0.0000001 ;
    typedef std::vector <size_t> indices;

    class DelphesReader {
    private:
        std::string        ListOfFiles        ;
        TChain             chain              ;
        ExRootTreeReader * treeReader         ;
        size_t             numberOfEntries    ;
        TClonesArray     * EFlowTrack         ;
        TClonesArray     * EFlowPhoton        ;
        TClonesArray     * EFlowNeutralHadron ;
        TClonesArray     * GenParticles       ;
    public:
        NewHEPHeaders::pseudojets                                 jetvectors ;
        std::vector < NewHEPHeaders::VECTORS::DelphesVectors <> > Vectors    ;
        std::vector < NewHEPHeaders::VECTORS::ParticleNode   <> > GenVectors ;
        NewHEPHeaders::VECTORS::lorentz4vector               <>   ZVector    ;
    private:
        inline void clear () {Vectors.clear();GenVectors.clear();jetvectors.clear();ZVector.clearthis();}
        inline NewHEPHeaders::pseudojets & Analyze   (size_t entry) {
            clear(); treeReader->ReadEntry(entry);
            size_t limit_EFlowTrack         = EFlowTrack->GetEntries         () ;
            size_t limit_EFlowPhoton        = EFlowPhoton->GetEntries        () ;
            size_t limit_EFlowNeutralHadron = EFlowNeutralHadron->GetEntries () ;
            size_t limit_GenParticles       = GenParticles->GetEntries       () ;
            /* Read gen level info: */ {
                GenVectors.resize(limit_GenParticles);
                for(size_t j=0;j<limit_GenParticles;j++){
                    GenParticle*tmp=(GenParticle*)GenParticles->At(j);
                    GenVectors[j]=NewHEPHeaders::VECTORS::ParticleNode<>(
                        tmp->Px,tmp->Py,tmp->Pz,tmp->E,tmp->D1,tmp->D2,tmp->PID
                    );
                }
                for(size_t j=0;(j<GenVectors.size())&&(ZVector.cleared());j++)
                {if(GenVectors[j].id()==NewHEPHeaders::PID::Z){ZVector=GenVectors[j].getvec();}}
            }
            /* Read tracks: */ {
                for(size_t i=0;i<limit_EFlowTrack;i++){
                    Track*tmp=(Track*)EFlowTrack->At(i);
                    NewHEPHeaders::VECTORS::DelphesVectors<>tmp2;
                    tmp2.SetPtEtaPhiM(tmp->PT,tmp->Eta,tmp->Phi,0);
                    if(CPPFileIO::mymod(tmp->PID)==NewHEPHeaders::PID::MUON){
                        tmp2.Eem  = tmp2[3] * 0.0 ;
                        tmp2.Ehad = tmp2[3] * 0.0 ;
                        tmp2.Emu  = tmp2[3] * 1.0 ;
                    } else if(CPPFileIO::mymod(tmp->PID)==NewHEPHeaders::PID::ELECTRON) {
                        tmp2.Eem  = tmp2[3] * 1.0 ;
                        tmp2.Ehad = tmp2[3] * 0.0 ;
                        tmp2.Emu  = tmp2[3] * 0.0 ;
                    } else {
                        tmp2.Eem  = tmp2[3] * 0.1 ;
                        tmp2.Ehad = tmp2[3] * 0.9 ;
                        tmp2.Emu  = tmp2[3] * 0.0 ;
                    }
                    tmp2.Charge=(CPPFileIO::mymod(tmp->Charge));
                    fastjet::PseudoJet tmpjet = tmp2.getpseudojet();
                    tmpjet.set_user_index (Vectors.size()) ;
                    Vectors.push_back (tmp2) ; jetvectors.push_back (tmpjet) ;
                }
            }
            /* Read Photons: */ {
                for(size_t i=0;i<limit_EFlowPhoton;i++){
                    Tower*tmp=(Tower*)EFlowPhoton->At(i);
                    NewHEPHeaders::VECTORS::DelphesVectors<>tmp2;
                    tmp2.SetPtEtaPhiM(tmp->ET,tmp->Eta,tmp->Phi,0);
                    tmp2.Eem=tmp->Eem; tmp2.Ehad=tmp->Ehad;
                    tmp2.Emu=0; tmp2.Charge=0;
                    fastjet::PseudoJet tmpjet = tmp2.getpseudojet();
                    tmpjet.set_user_index (Vectors.size()) ;
                    Vectors.push_back (tmp2) ; jetvectors.push_back (tmpjet) ;
                }
            }
            /* Read Neutral Hadrons: */ {
                for(size_t i=0;i<limit_EFlowNeutralHadron;i++){
                    Tower*tmp=(Tower*)EFlowNeutralHadron->At(i);
                    NewHEPHeaders::VECTORS::DelphesVectors<>tmp2;
                    tmp2.SetPtEtaPhiM(tmp->ET,tmp->Eta,tmp->Phi,0);
                    tmp2.Eem=tmp->Eem; tmp2.Ehad=tmp->Ehad;
                    tmp2.Emu=0; tmp2.Charge=0;
                    fastjet::PseudoJet tmpjet=tmp2.getpseudojet();
                    tmpjet.set_user_index(Vectors.size());
                    Vectors.push_back(tmp2); jetvectors.push_back(tmpjet);
                }
            }
            return jetvectors;
        }
        inline void                        construct ()             {
            chain.Add            (&(ListOfFiles[0]))                             ;
            treeReader         = new ExRootTreeReader   ( &chain               ) ;
            numberOfEntries    = treeReader->GetEntries (                      ) ;
            EFlowTrack         = treeReader->UseBranch  ( "EFlowTrack"         ) ;
            EFlowPhoton        = treeReader->UseBranch  ( "EFlowPhoton"        ) ;
            EFlowNeutralHadron = treeReader->UseBranch  ( "EFlowNeutralHadron" ) ;
            GenParticles       = treeReader->UseBranch  ( "Particle"           ) ;
        }
        inline void                        destroy   ()             {delete treeReader;}
    public:
        inline NewHEPHeaders::pseudojets & operator () (size_t entry) {return Analyze(entry);}
        inline size_t count_tracks           (std::vector <size_t> & in_indices) const {
            size_t ret=0;
            for(size_t i=0;i<in_indices.size();i++){ret=ret+CPPFileIO::mymod(Vectors[in_indices[i]].Charge);}
            return ret;
        }
        inline double hadronic_energy        (std::vector <size_t> & in_indices) const {
            double ret=0;
            for(size_t i=0;i<in_indices.size();i++){ret=ret+Vectors[in_indices[i]].Ehad;}
            return ret;
        }
        inline double electromagnetic_energy (std::vector <size_t> & in_indices) const {
            double ret=0;
            for(size_t i=0;i<in_indices.size();i++){ret=ret+Vectors[in_indices[i]].Eem;}
            return ret;
        }
        inline size_t operator ()            ()                                        {return numberOfEntries;}
        DelphesReader (std::string&_ListOfFiles): ListOfFiles(_ListOfFiles), chain("Delphes") {construct();}
        ~DelphesReader () {destroy();}
    };

    class HardSubStructureFinder {
    public:
        const double max_subjet_mass, mass_drop_threshold, Rfilt, minpt_subjet, mh, mhmin, mhmax, zcut, rcut_factor;
        const size_t nfilt;
        double filteredjetmass, deltah, filt_tau_R, prunedmass, unfiltered_mass,
        EFC[6], EFCDR[4], frac_em, frac_had, nsub[5], nsub_ratio[4];
        bool    HiggsTagged        ;
        size_t  n_tracks           ;
        indices index_constituents ;
        NewHEPHeaders::pseudojets tau_subs  , t_parts  , tau_hadrons ;
        fastjet::PseudoJet        prunedjet , triple   , Higgs       , taucandidate ;
    private:
        inline void clear                   (                                     ) {
            t_parts.clear ()      ; tau_subs.clear ()  ; tau_hadrons.clear () ; index_constituents.clear() ;
            filteredjetmass = 0.0 ; filt_tau_R =     0 ; prunedmass  = 0.0    ; n_tracks = 0 ;
            unfiltered_mass = 0.0 ; deltah     = 10000 ; HiggsTagged = false  ;
            for (size_t i=0;i<6;i++) { EFC        [i] = -10000.0 ; }
            for (size_t i=0;i<4;i++) { EFCDR      [i] = -10000.0 ; }
            for (size_t i=0;i<5;i++) { nsub       [i] = -10000.0 ; }
            for (size_t i=0;i<4;i++) { nsub_ratio [i] = -10000.0 ; }
        }
        inline void read_extra_variables    ( const DelphesReader      & in       ) {
            n_tracks = in.count_tracks           (index_constituents) ;
            frac_em  = in.electromagnetic_energy (index_constituents) / taucandidate.E () ;
            frac_had = in.hadronic_energy        (index_constituents) / taucandidate.E () ;
        }
        inline void get_constituent_indices ( const fastjet::PseudoJet & this_jet ) {
            NewHEPHeaders::pseudojets vectors=this_jet.constituents();
            const size_t limit=vectors.size(); index_constituents.resize(limit);
            for (size_t i=0;i<limit;i++) {index_constituents[i]=vectors[i].user_index();}
        }
        inline void EvalEnergyCorrelation   ( fastjet::PseudoJet & this_jet       ) {
            using namespace fastjet          ;
            using namespace fastjet::contrib ;
            const double beta    = 2.0                    ;
            if(this_jet.constituents().size()>0){
                /* The energy correlation part: */ {
                    const auto   measure = EnergyCorrelator::pt_R ;
                    EnergyCorrelator ECF0 ( 0, beta, measure ) ; EFC[0] = ECF0 (this_jet) ;
                    EnergyCorrelator ECF1 ( 1, beta, measure ) ; EFC[1] = ECF1 (this_jet) ;
                    EnergyCorrelator ECF2 ( 2, beta, measure ) ; EFC[2] = ECF2 (this_jet) ;
                    EnergyCorrelator ECF3 ( 3, beta, measure ) ; EFC[3] = ECF3 (this_jet) ;
                    EnergyCorrelator ECF4 ( 4, beta, measure ) ; EFC[4] = ECF4 (this_jet) ;
                    EnergyCorrelator ECF5 ( 5, beta, measure ) ; EFC[5] = ECF5 (this_jet) ;
                    for (size_t i=0;i<4;i++) if (EFC[i+1]>epsilon) {EFCDR[i]=(EFC[i]+EFC[i+2])/(EFC[i+1]*EFC[i+1]);}
                }
                /* The NSubjettiness part: */ {
                    Nsubjettiness nSub1 ( 1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    Nsubjettiness nSub2 ( 2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    Nsubjettiness nSub3 ( 3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    Nsubjettiness nSub4 ( 4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    Nsubjettiness nSub5 ( 5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    nsub[0]=nSub1(this_jet); nsub[1]=nSub2(this_jet); nsub[2]=nSub3(this_jet);
                    nsub[3]=nSub4(this_jet); nsub[4]=nSub5(this_jet);
                    for(size_t i=0;i<4;i++) if(nsub[i]>epsilon) {nsub_ratio[i]=nsub[i+1]/nsub[i];}
                }
            }
        }
        inline void find_structures         ( const fastjet::PseudoJet & this_jet ) {
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
                } else {find_structures(parent1);return;}
            } else {return;}
            // The part below is never reached...
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
        inline void run_filter              (                                     ) {
            t_parts    = sorted_by_pt  (t_parts)               ;
            triple     = fastjet::join (t_parts[0],t_parts[1]) ;
            filt_tau_R = std::min ( Rfilt , 0.5 * sqrt (t_parts[0].squared_distance(t_parts[1])) ) ;
            fastjet::JetDefinition filtering_def (fastjet::cambridge_algorithm,filt_tau_R) ;
            fastjet::Filter filter (filtering_def,fastjet::SelectorNHardest(nfilt)*fastjet::SelectorPtMin(minpt_subjet)) ;
            taucandidate    = filter         (triple) ;
            filteredjetmass = taucandidate.m ()       ;
            EvalEnergyCorrelation   ( taucandidate )  ;
            get_constituent_indices ( taucandidate )  ;
        }
        inline void run_recluster           (                                     ) {
            fastjet::JetDefinition   reclustering (fastjet::cambridge_algorithm,10.0)  ;
            fastjet::ClusterSequence cs_top_sub   (taucandidate.pieces(),reclustering) ;
            tau_subs=sorted_by_pt(cs_top_sub.exclusive_jets(2));
        }
        inline void run_variable_evaluater  ( fastjet::PseudoJet       & injet    ) {
            HiggsTagged  = true;
            Higgs        = tau_subs[0]+tau_subs[1];
            deltah       = CPPFileIO::mymod(taucandidate.m()-mh);
            tau_hadrons  = taucandidate.constituents();
            double Rprun = injet.validated_cluster_sequence()->jet_def().R();
            fastjet::JetDefinition jet_def_prune(fastjet::cambridge_algorithm, Rprun);
            fastjet::Pruner pruner(jet_def_prune,zcut,rcut_factor);
            prunedjet       = pruner      (triple) ;
            prunedmass      = prunedjet.m ()       ;
            unfiltered_mass = triple.m    ()       ;
        }
        inline void run                     ( fastjet::PseudoJet       & injet    ) {
            clear(); find_structures(injet);
            if (t_parts.size()>1) {
                run_filter();
                if((mhmin<filteredjetmass)&&(filteredjetmass<mhmax)&&(taucandidate.pieces().size()>1)){
                    run_recluster();
                    if (tau_subs[1].perp()>minpt_subjet)
                    {run_variable_evaluater(injet);}
                }
            }
        }
    public:
        inline void   operator () (const DelphesReader&in)   {read_extra_variables(in);}
        inline void   operator () ()                         {clear();}
        inline double operator () (fastjet::PseudoJet&injet) {
            run(injet);
            if (HiggsTagged) {return filteredjetmass;}
            else             {return -1000;}
        }
        HardSubStructureFinder() :
        max_subjet_mass ( 30.0 ) , mass_drop_threshold (  0.5 ) , minpt_subjet (      20.0 ) ,
        Rfilt           (  0.3 ) , mh                  ( 90.0 ) , mhmin        ( mh - 50.0 ) ,
        zcut            (  0.1 ) , rcut_factor         (  0.5 ) , mhmax        ( mh + 50.0 ) ,
        nfilt           (  4   ) { clear               (      ) ; }
        ~HardSubStructureFinder(){}
    };

    class OutPutVariables {
    public:
        double filteredjetmass, deltah, filt_tau_R, prunedmass, unfiltered_mass,
        EFC[6], EFCDR[4], frac_em, frac_had, nsub[5], nsub_ratio[4];
        bool    HiggsTagged;
        size_t  n_tracks;
        NewHEPHeaders::VECTORS::lorentz4vector <> tau_subs[2] , t_parts[2]  , tau_hadrons[2] ;
        NewHEPHeaders::VECTORS::lorentz4vector <> prunedjet   , triple      , taucandidate   , Higgs ;
    private:
        inline void clear () {
            t_parts[0].clearthis () ; tau_subs[0].clearthis () ; tau_hadrons[0].clearthis () ;
            t_parts[1].clearthis () ; tau_subs[1].clearthis () ; tau_hadrons[1].clearthis () ;
            filteredjetmass = 0.0 ; filt_tau_R =     0 ; prunedmass  = 0.0    ; n_tracks = 0 ;
            unfiltered_mass = 0.0 ; deltah     = 10000 ; HiggsTagged = false  ;
            for (size_t i=0;i<6;i++) { EFC        [i] = -10000.0 ; }
            for (size_t i=0;i<4;i++) { EFCDR      [i] = -10000.0 ; }
            for (size_t i=0;i<5;i++) { nsub       [i] = -10000.0 ; }
            for (size_t i=0;i<4;i++) { nsub_ratio [i] = -10000.0 ; }
        }
        inline void ReadFrom (HardSubStructureFinder&other) {
            frac_em         = other.frac_em         ;
            frac_had        = other.frac_had        ;
            filteredjetmass = other.filteredjetmass ;
            filt_tau_R      = other.filt_tau_R      ;
            prunedmass      = other.prunedmass      ;
            n_tracks        = other.n_tracks        ;
            unfiltered_mass = 0.0                   ;
            deltah          = 10000                 ;
            HiggsTagged     = false                 ;
            for (size_t i=0;i<6;i++) { EFC        [i] = other.EFC        [i] ; }
            for (size_t i=0;i<4;i++) { EFCDR      [i] = other.EFCDR      [i] ; }
            for (size_t i=0;i<5;i++) { nsub       [i] = other.nsub       [i] ; }
            for (size_t i=0;i<4;i++) { nsub_ratio [i] = other.nsub_ratio [i] ; }
            for ( size_t i=0 ; i < CPPFileIO::mymin ( other.t_parts.size     () , (size_t)2 ) ; i++ )
            { t_parts     [i] = other.t_parts     [i] ; }
            for ( size_t i=0 ; i < CPPFileIO::mymin ( other.tau_subs.size    () , (size_t)2 ) ; i++ )
            { tau_subs    [i] = other.tau_subs    [i] ; }
            for ( size_t i=0 ; i < CPPFileIO::mymin ( other.tau_hadrons.size () , (size_t)2 ) ; i++ )
            { tau_hadrons [i] = other.tau_hadrons [i] ; }
        }
    public:
        inline void operator () () {clear();}
        inline void operator = (HardSubStructureFinder&other) {ReadFrom(other);}
        OutPutVariables  () {}
        ~OutPutVariables () {}
    };

    class MainAnalyzer {
    private:
        DelphesReader * MainReader;
        std::string OutFileName;
        CPPFileIO::FileVector<OutPutVariables>Writer;
        inline void Analyze(size_t i){
            NewHEPHeaders::pseudojets&jetvectors=MainReader[0](i);
            fastjet::JetAlgorithm algorithm=fastjet::antikt_algorithm;
            const double jet_rad=1.00;
            fastjet::JetDefinition jetDef(algorithm,jet_rad);
            fastjet::ClusterSequence clust_seq(jetvectors,jetDef);
            NewHEPHeaders::pseudojets antikt_jets=sorted_by_pt(clust_seq.inclusive_jets());
            for (size_t j=0;(j<antikt_jets.size())&&(j<2);j++) if(antikt_jets[j].perp()>200.0) {
                fastjet::PseudoJet&this_jet=antikt_jets[j];
                if(this_jet.m()>40){
                    HardSubStructureFinder tmpslave ;
                    double mass = tmpslave(this_jet);
                    if(tmpslave.HiggsTagged){
                        tmpslave(MainReader[0]);
                        OutPutVariables tmp; tmp = tmpslave ;
                        Writer.push_back(tmp);
                    }
                }
            }
        }
        inline void Analyze(){
            size_t limit=MainReader[0]();
            for(size_t i=0;i<limit;i++){
                if((i%500)==0){printf("Analyzing %ld event:\n",i);}
                Analyze(i);
            }
        }
        inline void AnalyzeNewFile(std::string _delphesfilename){
            DelphesReader tmpreader(_delphesfilename);
            MainReader=&tmpreader; Analyze();
        }

    public:
        inline void operator()(std::string _delphesfilename){AnalyzeNewFile(_delphesfilename);}
        MainAnalyzer (std::string _OutFileName): OutFileName (_OutFileName), Writer(OutFileName) {}
        ~MainAnalyzer(){}
    };

    inline void ProcessType (std::string _name) {
        std::string name(_name)      ;
        std::string NoISRList   [16] ;
        std::string WithISRList [16] ;
        /* Get the names of root files: */ {
            for(size_t i=0;i<16;i++){
                char tmp[512] ;
                sprintf (tmp,"./%s/%ld/NoISRout.root",&(name[0]),i+1) ; NoISRList[i]   = std::string(tmp);
                sprintf (tmp,"./%s/%ld/out.root",&(name[0]),i+1)      ; WithISRList[i] = std::string(tmp);
            }
        }
        /* Prepare the output: */ {
            mkdir("./SKIM_DATA/",(mode_t)0755);
            char tmp[512] ;
            sprintf(tmp,"./SKIM_DATA/%s",&(name[0])); mkdir(tmp,(mode_t)0755);
            CPPFileIO::ForkMe forker;
            if(forker.InKid()){
                sprintf(tmp,"./SKIM_DATA/%s/NoMPI",&(name[0])); MainAnalyzer NoMPI(tmp);
                for(size_t i=0;i<16;i++){NoMPI(NoISRList[i]);}
            }
            if(forker.InKid()){
                sprintf(tmp,"./SKIM_DATA/%s/WithMPI",&(name[0])); MainAnalyzer WithMPI(tmp);
                for(size_t i=0;i<16;i++){WithMPI(WithISRList[i]);}
            }
        }
    }

    inline void ProcessAll () {
        CPPFileIO::ForkMe forker;
        if(forker.InKid()){ProcessType("BoostedZToNuNuBar");}
        if(forker.InKid()){ProcessType("BoostedZ");}
        if(forker.InKid()){ProcessType("UnBoostedZ");}
    }
}

namespace Step2 {

    class MyHist1D {
    private:
        std::string histname      ;
        std::string lefthistname  ;
        TH1F        LH            ;
    public:
        inline void Fill  (double a) {if(a>-90.0){LH.Fill(a);}}
        inline void Write () {
            TCanvas C;
            LH.Scale(1.0/LH.Integral());
            LH.SetLineWidth(3);
            int binmaxL = LH.GetMaximumBin ()        ;
            double xL   = LH.GetBinContent (binmaxL) ;
            double maxy = xL;
            LH.SetMaximum(maxy);
            LH.SetLineColor(TColor::GetColor("#990000"));
            LH.Draw("hist same");
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            std::string outname = "./GRAPHS/" + histname + ".pdf" ;
            C.SaveAs(&(outname[0]));
        }
        MyHist1D  (std::string _histname, size_t nbins, double min, double max):
        histname(_histname),
        lefthistname(_histname),
        LH(&(lefthistname[0]),&(histname[0]),nbins,min,max) {}
        ~MyHist1D (){Write();}
    };

    class MyHist2D {
    private:
        std::string histname  ;
        std::string histname1 ;
        std::string histname2 ;
        TH1F H1, H2;
    public:
        inline void Fill1 (double a) {if(a>-90.0){H1.Fill(a);}}
        inline void Fill2 (double a) {if(a>-90.0){H2.Fill(a);}}
        inline void Write () {
            TCanvas C;
            H1.Scale(1.0/H1.Integral());
            H2.Scale(1.0/H2.Integral());
            H1.SetLineWidth(3); H2.SetLineWidth(3);
            int binmax1 = H1.GetMaximumBin ()        ;
            int binmax2 = H2.GetMaximumBin ()        ;
            double x1   = H1.GetBinContent (binmax1) ;
            double x2   = H2.GetBinContent (binmax2) ;
            double maxy = CPPFileIO::mymax (x1,x2)   ;
            H1.SetMaximum(maxy); H2.SetMaximum(maxy);
            H1.SetLineColor(TColor::GetColor("#990000"));
            H2.SetLineColor(TColor::GetColor("#000099"));
            H1.Draw("hist same"); H2.Draw("hist same");
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            std::string outname = "./GRAPHS/" + histname + ".pdf" ;
            C.SaveAs(&(outname[0]));
        }
        MyHist2D  (std::string _histname, size_t nbins, double min, double max):
        histname(_histname),
        histname1("H1"+_histname), histname2("H2"+_histname),
        H1(&(histname1[0]),&(histname[0]),nbins,min,max),
        H2(&(histname2[0]),&(histname[0]),nbins,min,max) {}
        ~MyHist2D () {Write();}
    };

    class MyHist3D {
    private:
        std::string histname ;
        std::string H1name, H2name, H3name  ;
        TH1F H1, H2, H3 ;
    public:
        bool WriteStacked;
        inline void H1Fill (double a) {if(a>-90.0){H1.Fill(a);}}
        inline void H2Fill (double a) {if(a>-90.0){H2.Fill(a);}}
        inline void H3Fill (double a) {if(a>-90.0){H3.Fill(a);}}

        inline void Write () {
            H1.Scale(1.0/H1.Integral());
            H2.Scale(1.0/H2.Integral());
            H3.Scale(1.0/H3.Integral());

            int binmax1 = H1.GetMaximumBin () ;
            int binmax2 = H2.GetMaximumBin () ;
            int binmax3 = H3.GetMaximumBin () ;

            double x=CPPFileIO::mymax
            (H3.GetBinContent(binmax3),CPPFileIO::mymax(H1.GetBinContent(binmax1),H2.GetBinContent(binmax2)));

            H1.SetLineColor(TColor::GetColor("#990000"));
            H2.SetLineColor(TColor::GetColor("#009900"));
            H3.SetLineColor(TColor::GetColor("#000099"));
            H1.SetMaximum     (x) ; H2.SetMaximum     (x) ; H3.SetMaximum     (x) ;
            H1.SetLineWidth   (3) ; H2.SetLineWidth   (3) ; H3.SetLineWidth   (3) ;
            TCanvas C;
            H1.Draw ("hist same") ; H2.Draw ("hist same") ; H3.Draw ("hist same") ;
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            std::string outname = "./GRAPHS/" + histname + ".pdf" ;
            C.SaveAs(&(outname[0]));
        }

        inline void WriteStack () {
            H1.SetLineColor(TColor::GetColor("#990000")); H1.SetFillColor(TColor::GetColor("#990000"));
            H2.SetLineColor(TColor::GetColor("#009900")); H2.SetFillColor(TColor::GetColor("#009900"));
            H3.SetLineColor(TColor::GetColor("#000099")); H3.SetFillColor(TColor::GetColor("#000099"));
            H1.SetLineWidth(3); H2.SetLineWidth(3); H3.SetLineWidth(3);
            THStack hs ("hs","Stacked 1D histograms") ;
            hs.Add(&H1);
            hs.Add(&H2);
            hs.Add(&H3);
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            TCanvas C; hs.Draw(); std::string outname = "./GRAPHS/" + histname + ".pdf" ;
            C.SaveAs(&(outname[0]));
        }

        MyHist3D  (std::string _histname, size_t nbins, double min, double max):
        histname(_histname),
        H1name("H1"+histname), H2name("H2"+histname), H3name("H3"+histname),
        H1(&(H1name[0]),&(histname[0]),nbins,min,max),
        H2(&(H2name[0]),&(histname[0]),nbins,min,max),
        H3(&(H3name[0]),&(histname[0]),nbins,min,max) {WriteStacked=false;}
        ~MyHist3D (){
            if(WriteStacked){WriteStack();}
            else{Write();}
        }
    } ;

    inline void PlotHist (std::string name, std::vector<float>&vals) {
        size_t limit=vals.size();
        if(limit>0){
            std::sort(vals.begin(),vals.end());
            TH1F hist(&(name[0]),&(name[0]),100,vals[0],vals[limit-1]);
            for(size_t i=0;i<limit;i++){hist.Fill(vals[i]);}
            TCanvas C;
            hist.Scale(1.0/hist.Integral());
            hist.SetLineWidth(3);
            int binmaxL = hist.GetMaximumBin ()        ;
            double xL   = hist.GetBinContent (binmaxL) ;
            double maxy = xL;
            hist.SetMaximum(maxy);
            hist.SetLineColor(TColor::GetColor("#990000"));
            hist.Draw("hist same");
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            std::string outname = "./GRAPHS/" + name + ".pdf" ;
            C.SaveAs(&(outname[0]));
        }
    }

    inline void PlotHist (std::string name, std::vector<float>&vals, std::vector<float>&vals2) {
        TCanvas C                      ;
        size_t limit  = vals.size  ()  ;
        size_t limit2 = vals2.size ()  ;
        float xmin , xmax              ;
        std::string name2 = name + "2" ;
        /* Check for limits of histograms: */ {
            if(limit>0){
                std::sort ( vals.begin() , vals.end() ) ;
                xmin = vals [0]       ;
                xmax = vals [limit-1] ;
            }
            if(limit2>0){
                std::sort ( vals2.begin() , vals2.begin() ) ;
                xmin = CPPFileIO::mymin ( vals2 [0]       , xmin ) ;
                xmax = CPPFileIO::mymax ( vals2 [limit-1] , xmax ) ;
            }
        }
        TH1F hist  ( & ( name  [0] ) , & ( name  [0] ) , 100 , xmin , xmax ) ;
        TH1F hist2 ( & ( name2 [0] ) , & ( name2 [0] ) , 100 , xmin , xmax ) ;
        /* prepare the histograms: */ {
            /* Fill the histograms: */ {
                for ( size_t i = 0 ; i < vals.size  () ; i++ ) { hist.Fill  ( vals  [i] ) ; }
                for ( size_t i = 0 ; i < vals2.size () ; i++ ) { hist2.Fill ( vals2 [i] ) ; }
            }
            /* Rescale and color the histograms: */ {
                hist.Scale  ( 1.0 / hist.Integral  () ) ; hist.SetLineWidth  (3) ; hist.SetLineColor  (TColor::GetColor("#990000"));
                hist2.Scale ( 1.0 / hist2.Integral () ) ; hist2.SetLineWidth (3) ; hist2.SetLineColor (TColor::GetColor("#000099"));
            }
            /* Set the maximum */ {
                float x1 = hist.GetBinContent  ( hist.GetMaximumBin  () ) ;
                float x2 = hist2.GetBinContent ( hist2.GetMaximumBin () ) ;
                float x = CPPFileIO::mymax(x1,x2);
                hist.SetMaximum(x); hist2.SetMaximum(x);
            }
        }
        /* Draw and save the histograms: */ {
            TCanvas C;
            hist.Draw  ("hist same") ;
            hist2.Draw ("hist same") ;
            name = name + ".pdf";
            C.SaveAs(&(name[0]));
        }
    }

    inline void PlotHist (std::string name, std::vector<float>&vals, std::vector<float>&vals2, std::vector<float>&vals3) {

        float xmin , xmax              ;
        std::string name2 = name + "2" ;
        std::string name3 = name + "3" ;

        size_t limit  = vals.size  ()  ;
        size_t limit2 = vals2.size ()  ;
        size_t limit3 = vals3.size ()  ;

        /* Check for limits of histograms: */ {
            if(limit>0){
                std::sort ( vals.begin() , vals.end() ) ;
                xmin = vals [0]       ;
                xmax = vals [limit-1] ;
            }
            if(limit2>0){
                std::sort ( vals2.begin() , vals2.begin() ) ;
                xmin = CPPFileIO::mymin ( vals2 [0]       , xmin ) ;
                xmax = CPPFileIO::mymax ( vals2 [limit-1] , xmax ) ;
            }
            if(limit3>0){
                std::sort ( vals3.begin() , vals3.begin() ) ;
                xmin = CPPFileIO::mymin ( vals3 [0]       , xmin ) ;
                xmax = CPPFileIO::mymax ( vals3 [limit-1] , xmax ) ;
            }
        }

        TH1F hist  ( & ( name  [0] ) , & ( name  [0] ) , 100 , xmin , xmax ) ;
        TH1F hist2 ( & ( name2 [0] ) , & ( name2 [0] ) , 100 , xmin , xmax ) ;
        TH1F hist3 ( & ( name3 [0] ) , & ( name3 [0] ) , 100 , xmin , xmax ) ;

        /* prepare the histograms: */ {
            /* Fill the histograms: */ {
                for ( size_t i = 0 ; i < vals.size  () ; i++ ) { hist.Fill  ( vals  [i] ) ; }
                for ( size_t i = 0 ; i < vals2.size () ; i++ ) { hist2.Fill ( vals2 [i] ) ; }
                for ( size_t i = 0 ; i < vals3.size () ; i++ ) { hist3.Fill ( vals3 [i] ) ; }
            }
            /* Rescale and color the histograms: */ {
                hist.Scale  ( 1.0 / hist.Integral  () ) ; hist.SetLineWidth  (3) ; hist.SetLineColor  (TColor::GetColor("#990000"));
                hist2.Scale ( 1.0 / hist2.Integral () ) ; hist2.SetLineWidth (3) ; hist2.SetLineColor (TColor::GetColor("#009900"));
                hist3.Scale ( 1.0 / hist3.Integral () ) ; hist3.SetLineWidth (3) ; hist3.SetLineColor (TColor::GetColor("#000099"));
            }
            /* Set the maximum */ {
                float x1 = hist.GetBinContent  ( hist.GetMaximumBin  () ) ;
                float x2 = hist2.GetBinContent ( hist2.GetMaximumBin () ) ;
                float x3 = hist3.GetBinContent ( hist3.GetMaximumBin () ) ;
                float x  = CPPFileIO::mymax ( CPPFileIO::mymax (x1,x2) , x3 ) ;
                hist.SetMaximum(x); hist2.SetMaximum(x); hist3.SetMaximum(x);
            }
        }

        /* Draw and save the histograms: */ {
            TCanvas C;
            hist.Draw  ("hist same") ;
            hist2.Draw ("hist same") ;
            hist3.Draw ("hist same") ;
            name = name + ".pdf";
            C.SaveAs(&(name[0]));
        }

    }

    class PlotAll {
    private:
        CPPFileIO::FileArray <Step1::OutPutVariables> reader1, reader2;
        size_t Limit1, Limit2;
        Step1::OutPutVariables *element1, *element2 ;
    public:
        inline void Plot_Masses () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].frac_had;}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].frac_had;}
            PlotHist("Masses",Masses1,Masses2);
        }

        PlotAll():
        reader1 ( "./SKIM_DATA/BoostedZ/NoMPI"   ) , Limit1 (reader1.size()) , element1(&(reader1(0,Limit1))) ,
        reader2 ( "./SKIM_DATA/BoostedZ/WithMPI" ) , Limit2 (reader2.size()) , element2(&(reader2(0,Limit2)))
        {Plot_Masses();}

        ~PlotAll(){}
    } ;

}
