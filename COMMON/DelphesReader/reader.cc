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
    double                    filteredjetmass , deltah   , filt_tau_R  , prunedmass   , unfiltered_mass ;
    NewHEPHeaders::pseudojets tau_subs        , t_parts  , tau_hadrons ;
    fastjet::PseudoJet        prunedjet       , triple   , Higgs       , taucandidate ;
    bool HiggsTagged ;
    inline void operator () () {initialize();}
    inline double operator () (fastjet::PseudoJet&injet) {
        run(injet);
        if(HiggsTagged){return filteredjetmass;}
        else {return -1000;}
    }
    HardSubStructureFinder(){initialize();}
    ~HardSubStructureFinder(){}
};

class DelphesReader {
private:
    std::string               & ListOfFiles        ;
    TChain                      chain              ;
    ExRootTreeReader          * treeReader         ;
    size_t                      numberOfEntries    ;
    TClonesArray              * EFlowTrack         ;
    TClonesArray              * EFlowPhoton        ;
    TClonesArray              * EFlowNeutralHadron ;
    NewHEPHeaders::pseudojets   jetvectors         ;
    inline NewHEPHeaders::pseudojets & Analyze (size_t entry) {
        treeReader->ReadEntry (entry) ; jetvectors.clear () ;
        size_t limit_EFlowTrack         = EFlowTrack->GetEntries         () ;
        size_t limit_EFlowPhoton        = EFlowPhoton->GetEntries        () ;
        size_t limit_EFlowNeutralHadron = EFlowNeutralHadron->GetEntries () ;
        for(size_t i=0;i<limit_EFlowTrack;i++){
            Track*tmp=(Track*)EFlowTrack->At(i);
            NewHEPHeaders::VECTORS::lorentz4vector<>tmp2;
            tmp2.SetPtEtaPhiM(tmp->PT,tmp->Eta,tmp->Phi,0);
            jetvectors.push_back(tmp2.getpseudojet());
        }
        for(size_t i=0;i<limit_EFlowPhoton;i++){
            Tower*tmp=(Tower*)EFlowPhoton->At(i);
            NewHEPHeaders::VECTORS::lorentz4vector<>tmp2;
            tmp2.SetPtEtaPhiM(tmp->ET,tmp->Eta,tmp->Phi,0);
            jetvectors.push_back(tmp2.getpseudojet());
        }
        for(size_t i=0;i<limit_EFlowNeutralHadron;i++){
            Tower*tmp=(Tower*)EFlowNeutralHadron->At(i);
            NewHEPHeaders::VECTORS::lorentz4vector<>tmp2;
            tmp2.SetPtEtaPhiM(tmp->ET,tmp->Eta,tmp->Phi,0);
            jetvectors.push_back(tmp2.getpseudojet());
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
    inline NewHEPHeaders::pseudojets & operator () (size_t entry) {return Analyze(entry);}
    inline size_t operator () () {return numberOfEntries;}
    DelphesReader (std::string&_ListOfFiles) : ListOfFiles(_ListOfFiles), chain("Delphes") {construct();}
    ~DelphesReader () {destroy();}
};
class MainAnalyzer {
private:
    DelphesReader * MainReader;
    std::string OutFileName;
    std::string delphesfilename ;
    TH1F OutMassHist, nsubjettinesshist1, nsubjettinesshist2, nsubjettinesshist3 ;
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
    }
public:
    inline void operator()(std::string _delphesfilename){AnalyzeNewFile(_delphesfilename);}

    MainAnalyzer (std::string _OutFileName):
    OutFileName        (_OutFileName)                                            ,
    OutMassHist        ("OutMassHist","OutMassHist",100,20,150)                  ,
    nsubjettinesshist1 ("nsubjettinesshist1","nsubjettinesshist1",100,-0.1,30.1) ,
    nsubjettinesshist2 ("nsubjettinesshist2","nsubjettinesshist2",100,-0.1,30.1) ,
    nsubjettinesshist3 ("nsubjettinesshist3","nsubjettinesshist3",100,-0.1,30.1) {}

    ~MainAnalyzer(){WriteHistograms();}
};
