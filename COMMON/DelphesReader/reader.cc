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

constexpr double epsilon = 0.0000001 ;
typedef std::vector <size_t> indices;

class DelphesReader {
private:
    std::string               & ListOfFiles        ;
    TChain                      chain              ;
    ExRootTreeReader          * treeReader         ;
    size_t                      numberOfEntries    ;
    TClonesArray              * EFlowTrack         ;
    TClonesArray              * EFlowPhoton        ;
    TClonesArray              * EFlowNeutralHadron ;

    inline NewHEPHeaders::pseudojets & Analyze (size_t entry) {
        treeReader->ReadEntry (entry) ;
        jetvectors.clear () ; Vectors.clear () ;
        size_t limit_EFlowTrack         = EFlowTrack->GetEntries         () ;
        size_t limit_EFlowPhoton        = EFlowPhoton->GetEntries        () ;
        size_t limit_EFlowNeutralHadron = EFlowNeutralHadron->GetEntries () ;
        for(size_t i=0;i<limit_EFlowTrack;i++){
            Track*tmp=(Track*)EFlowTrack->At(i);
            NewHEPHeaders::VECTORS::DelphesVectors<>tmp2;
            tmp2.SetPtEtaPhiM(tmp->PT,tmp->Eta,tmp->Phi,0);
            if(CPPFileIO::mymod(tmp->PID)==NewHEPHeaders::PID::MUON){
                tmp2.Eem=tmp2[3]*0.0;
                tmp2.Ehad=tmp2[3]*0.0;
                tmp2.Emu=tmp2[3]*1.0;
            } else if(CPPFileIO::mymod(tmp->PID)==NewHEPHeaders::PID::ELECTRON) {
                tmp2.Eem=tmp2[3]*1.0;
                tmp2.Ehad=tmp2[3]*0.0;
                tmp2.Emu=tmp2[3]*0.0;
            } else {
                tmp2.Eem=tmp2[3]*0.2;
                tmp2.Ehad=tmp2[3]*0.8;
                tmp2.Emu=tmp2[3]*0.0;
            }
            tmp2.Charge=(CPPFileIO::mymod(tmp->Charge));
            fastjet::PseudoJet tmpjet = tmp2.getpseudojet();
            tmpjet.set_user_index (Vectors.size()) ;
            Vectors.push_back (tmp2) ; jetvectors.push_back (tmpjet) ;
        }
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
        for(size_t i=0;i<limit_EFlowNeutralHadron;i++){
            Tower*tmp=(Tower*)EFlowNeutralHadron->At(i);
            NewHEPHeaders::VECTORS::DelphesVectors<>tmp2;
            tmp2.SetPtEtaPhiM(tmp->ET,tmp->Eta,tmp->Phi,0);
            tmp2.Eem=tmp->Eem; tmp2.Ehad=tmp->Ehad;
            tmp2.Emu=0; tmp2.Charge=0;
            fastjet::PseudoJet tmpjet = tmp2.getpseudojet();
            tmpjet.set_user_index (Vectors.size()) ;
            Vectors.push_back (tmp2) ; jetvectors.push_back (tmpjet) ;
        }
        return jetvectors;
    }
    inline void construct(){
        chain.Add            (&(ListOfFiles[0]))                             ;
        treeReader         = new ExRootTreeReader   ( &chain               ) ;
        numberOfEntries    = treeReader->GetEntries (                      ) ;
        EFlowTrack         = treeReader->UseBranch  ( "EFlowTrack"         ) ;
        EFlowPhoton        = treeReader->UseBranch  ( "EFlowPhoton"        ) ;
        EFlowNeutralHadron = treeReader->UseBranch  ( "EFlowNeutralHadron" ) ;
    }
    inline void destroy () {delete treeReader;}
public:
    NewHEPHeaders::pseudojets jetvectors ;
    std::vector <NewHEPHeaders::VECTORS::DelphesVectors<>> Vectors ;
    inline NewHEPHeaders::pseudojets & operator () (size_t entry) {return Analyze(entry);}

    inline size_t count_tracks (std::vector <size_t> & in_indices) const {
        size_t ret=0;
        for(size_t i=0;i<in_indices.size();i++){ret=ret+Vectors[in_indices[i]].Charge;}
        return ret;
    }
    inline double hadronic_energy (std::vector <size_t> & in_indices) const {
        double ret=0;
        for(size_t i=0;i<in_indices.size();i++){ret=ret+Vectors[in_indices[i]].Ehad;}
        return ret;
    }
    inline double electromagnetic_energy (std::vector <size_t> & in_indices) const {
        double ret=0;
        for(size_t i=0;i<in_indices.size();i++){ret=ret+Vectors[in_indices[i]].Eem;}
        return ret;
    }

    inline size_t operator () () {return numberOfEntries;}
    DelphesReader (std::string&_ListOfFiles) : ListOfFiles(_ListOfFiles), chain("Delphes") {construct();}
    ~DelphesReader () {destroy();}
};

class HardSubStructureFinder {
public:
    const double max_subjet_mass, mass_drop_threshold, Rfilt, minpt_subjet, mh, mhmin, mhmax, zcut, rcut_factor;
    const size_t nfilt;
    double filteredjetmass, deltah, filt_tau_R, prunedmass, unfiltered_mass, EFC[6], EFCDR[4], frac_em, frac_had;
    bool HiggsTagged ;
    size_t n_tracks ;
    indices index_constituents ;
    NewHEPHeaders::pseudojets tau_subs  , t_parts  , tau_hadrons ;
    fastjet::PseudoJet        prunedjet , triple   , Higgs       , taucandidate ;
private:
    void read_extra_variables (const DelphesReader & in) {
        n_tracks = in.count_tracks           (index_constituents) ;
        frac_em  = in.electromagnetic_energy (index_constituents) ;
        frac_had = in.hadronic_energy        (index_constituents) ;
    }
    inline void get_constituent_indices (const fastjet::PseudoJet & this_jet) {
        NewHEPHeaders::pseudojets vectors = this_jet.constituents();
        const size_t limit = vectors.size(); index_constituents.resize(limit);
        for (size_t i=0;i<limit;i++) {index_constituents[i]=vectors[i].user_index();}
    }
    inline void EvalEnergyCorrelation (const fastjet::PseudoJet & this_jet) {
        using namespace fastjet          ;
        using namespace fastjet::contrib ;
        const double beta    = 2.0                    ;
        const auto   measure = EnergyCorrelator::pt_R ;
        EnergyCorrelator ECF0(0,beta,measure); EFC[0]=ECF0(this_jet);
        EnergyCorrelator ECF1(1,beta,measure); EFC[1]=ECF1(this_jet);
        EnergyCorrelator ECF2(2,beta,measure); EFC[2]=ECF2(this_jet);
        EnergyCorrelator ECF3(3,beta,measure); EFC[3]=ECF3(this_jet);
        EnergyCorrelator ECF4(4,beta,measure); EFC[4]=ECF4(this_jet);
        EnergyCorrelator ECF5(5,beta,measure); EFC[5]=ECF5(this_jet);
        for (size_t i=0;i<4;i++) if (EFC[i+1]>epsilon) {EFCDR[i]=(EFC[i]+EFC[i+2])/EFC[i+1];}
    }
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
    inline void run_filter () {
        t_parts    = sorted_by_pt  (t_parts)               ;
        triple     = fastjet::join (t_parts[0],t_parts[1]) ;
        filt_tau_R = std::min ( Rfilt , 0.5 * sqrt (t_parts[0].squared_distance(t_parts[1])) ) ;
        fastjet::JetDefinition filtering_def (fastjet::cambridge_algorithm,filt_tau_R) ;
        fastjet::Filter filter (filtering_def,fastjet::SelectorNHardest(nfilt)*fastjet::SelectorPtMin(minpt_subjet)) ;
        taucandidate    = filter         (triple) ;
        filteredjetmass = taucandidate.m ()       ;
        EvalEnergyCorrelation   ( taucandidate ) ;
        get_constituent_indices ( taucandidate ) ;
    }
    inline void run_recluster () {
        fastjet::JetDefinition   reclustering (fastjet::cambridge_algorithm,10.0)  ;
        fastjet::ClusterSequence cs_top_sub   (taucandidate.pieces(),reclustering) ;
        tau_subs=sorted_by_pt(cs_top_sub.exclusive_jets(2));
    }
    inline void run_variable_evaluater (fastjet::PseudoJet&injet) {
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
    inline void run (fastjet::PseudoJet&injet) {
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
    inline void clear () {
        t_parts.clear ()      ; tau_subs.clear ()  ; tau_hadrons.clear () ; index_constituents.clear() ;
        filteredjetmass = 0.0 ; filt_tau_R =     0 ; prunedmass  = 0.0    ; n_tracks = 0 ;
        unfiltered_mass = 0.0 ; deltah     = 10000 ; HiggsTagged = false  ;
        for (size_t i=0;i<6;i++) { EFC   [i] = -10000.0 ; }
        for (size_t i=0;i<4;i++) { EFCDR [i] = -10000.0 ; }
    }
public:
    inline void   operator () () {clear();}
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



class MainAnalyzer {
private:
    DelphesReader * MainReader;
    std::string OutFileName;
    std::string delphesfilename ;
    TH1F OutMassHist, nsubjettinesshist1, nsubjettinesshist2, nsubjettinesshist3, NumTracks ;
    inline void Analyze(size_t i){
        NewHEPHeaders::pseudojets&jetvectors=MainReader[0](i);
        fastjet::JetAlgorithm algorithm=fastjet::antikt_algorithm;
        const double jet_rad=1.00;
        fastjet::JetDefinition jetDef(algorithm,jet_rad);
        fastjet::ClusterSequence clust_seq(jetvectors,jetDef);
        NewHEPHeaders::pseudojets antikt_jets=sorted_by_pt(clust_seq.inclusive_jets());
        for (size_t j=0;(j<antikt_jets.size())&&(j<2);j++) if(antikt_jets[j].perp()>200.0) {
            fastjet::PseudoJet&this_jet = antikt_jets[j];
            using namespace fastjet::contrib;
            double beta = 1.0;
            Nsubjettiness         nSub1_beta1(1,   OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
            Nsubjettiness         nSub2_beta1(2,   OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
            Nsubjettiness         nSub3_beta1(3,   OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
            NsubjettinessRatio   nSub21_beta1(2,1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
            NsubjettinessRatio   nSub32_beta1(3,2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
            double  tau1_beta1 =  nSub1_beta1(this_jet);
            double  tau2_beta1 =  nSub2_beta1(this_jet);
            double  tau3_beta1 =  nSub3_beta1(this_jet);
            double tau21_beta1 = nSub21_beta1(this_jet);
            double tau32_beta1 = nSub32_beta1(this_jet);
            beta = 2.0;
            Nsubjettiness         nSub1_beta2(1,   OnePass_KT_Axes(), UnnormalizedMeasure(beta));
            Nsubjettiness         nSub2_beta2(2,   OnePass_KT_Axes(), UnnormalizedMeasure(beta));
            Nsubjettiness         nSub3_beta2(3,   OnePass_KT_Axes(), UnnormalizedMeasure(beta));
            NsubjettinessRatio   nSub21_beta2(2,1, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
            NsubjettinessRatio   nSub32_beta2(3,2, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
            double  tau1_beta2 =  nSub1_beta2(this_jet);
            double  tau2_beta2 =  nSub2_beta2(this_jet);
            double  tau3_beta2 =  nSub3_beta2(this_jet);
            double tau21_beta2 = nSub21_beta2(this_jet);
            double tau32_beta2 = nSub32_beta2(this_jet);
            if(this_jet.m()>40){
                HardSubStructureFinder tmpslave ;
                OutMassHist.Fill        ( tmpslave(this_jet) ) ;
                if(tmpslave.HiggsTagged){
                    NewHEPHeaders::pseudojets constituents = this_jet.constituents();
                    std::vector <size_t> indices ;
                    for(size_t ii=0;ii<constituents.size();ii++)
                    {indices.push_back((size_t)constituents[ii].user_index());}
                    NumTracks.Fill((float)MainReader->count_tracks(indices));
                }
                nsubjettinesshist1.Fill ( tau1_beta1         ) ;
                nsubjettinesshist2.Fill ( tau2_beta1         ) ;
                nsubjettinesshist3.Fill ( tau3_beta1         ) ;
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
    inline void WriteHistograms(){
        CPPFileIO::ToDir circhanger(true); circhanger(OutFileName);
        /*OutMassHist*/ {
            TCanvas C;
            OutMassHist.Draw("hist");
            C.SaveAs("OutMassHist.pdf");
        }
        /*nsubjettinesshist1*/ {
            TCanvas C;
            nsubjettinesshist1.Draw("hist");
            C.SaveAs("nsubjettinesshist1.pdf");
        }
        /*nsubjettinesshist2*/ {
            TCanvas C;
            nsubjettinesshist2.Draw("hist");
            C.SaveAs("nsubjettinesshist2.pdf");
        }
        /*nsubjettinesshist3*/ {
            TCanvas C;
            nsubjettinesshist3.Draw("hist");
            C.SaveAs("nsubjettinesshist3.pdf");
        }
        /*nNumber of tracks*/ {
            TCanvas C;
            NumTracks.Draw("hist");
            C.SaveAs("NumTracks.pdf");
        }
    }
public:
    inline void operator()(std::string _delphesfilename){AnalyzeNewFile(_delphesfilename);}

    MainAnalyzer (std::string _OutFileName):
    OutFileName        (_OutFileName) ,
    OutMassHist        ( "OutMassHist"        , "OutMassHist"        , 100 , 20.0 , 150.0 ) ,
    nsubjettinesshist1 ( "nsubjettinesshist1" , "nsubjettinesshist1" , 100 , -0.1 ,  30.1 ) ,
    nsubjettinesshist2 ( "nsubjettinesshist2" , "nsubjettinesshist2" , 100 , -0.1 ,  30.1 ) ,
    nsubjettinesshist3 ( "nsubjettinesshist3" , "nsubjettinesshist3" , 100 , -0.1 ,  30.1 ) ,
    NumTracks          ( "NumTracks"          , "NumTracks"          , 120 , -0.1 ,  60.1 )
    {}

    ~MainAnalyzer(){WriteHistograms();}
};
