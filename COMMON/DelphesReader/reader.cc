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
    inline NewHEPHeaders::pseudojets&Analyze(size_t entry){
        treeReader->ReadEntry (entry) ;
        size_t limit_EFlowTrack         = EFlowTrack->GetEntries         () ;
        size_t limit_EFlowPhoton        = EFlowPhoton->GetEntries        () ;
        size_t limit_EFlowNeutralHadron = EFlowNeutralHadron->GetEntries () ;
        jetvectors.clear();
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
    inline void destroy(){delete treeReader;}
public:
    inline NewHEPHeaders::pseudojets&operator()(size_t entry){return Analyze(entry);}
    inline size_t operator()(){return numberOfEntries;}
    DelphesReader(std::string&_ListOfFiles):ListOfFiles(_ListOfFiles),chain("Delphes"){construct();}
    ~DelphesReader(){destroy();}
};
class MainAnalyzer {
private:
    DelphesReader*MainReader;
    std::string delphesfilename;
    TH1F OutMassHist ;
    inline void construct(std::string _delphesfilename)
    {delphesfilename=_delphesfilename;MainReader=new DelphesReader(delphesfilename);}
    inline void destroy(){delete MainReader;}
    inline void Analyze(size_t i){
        NewHEPHeaders::pseudojets&jetvectors=MainReader[0](i);
        fastjet::JetAlgorithm algorithm=fastjet::antikt_algorithm;
        const double jet_rad=1.00;
        fastjet::JetDefinition jetDef(algorithm,jet_rad);
        fastjet::ClusterSequence clust_seq(jetvectors,jetDef);
        NewHEPHeaders::pseudojets antikt_jets=sorted_by_pt(clust_seq.inclusive_jets());
        for (size_t j=0;(j<antikt_jets.size())&&(j<2);j++) if(antikt_jets[j].perp()>200.0) {
            fastjet::PseudoJet&this_jet = antikt_jets[j];
            if(this_jet.m()>40){OutMassHist.Fill(this_jet.m());}
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
        }
    }
    inline void Analyze(){
        size_t limit=MainReader[0]();
        for(size_t i=0;i<limit;i++){
            if((i%100)==0){printf("Analyzing %ld event:\n",i);}
            Analyze(i);
        }
    }
    inline void AnalyzeNewFile(std::string _delphesfilename){destroy();construct(_delphesfilename);Analyze();}
public:
    inline void operator()(std::string _delphesfilename){AnalyzeNewFile(_delphesfilename);}
    MainAnalyzer (std::string _delphesfilename):
    OutMassHist("OutMassHist","OutMassHist",100,40,200)
    {construct(_delphesfilename);Analyze();}
    ~MainAnalyzer(){
        destroy();
        TCanvas C;
        OutMassHist.Draw("hist");
        C.SaveAs("MassHist.pdf");
    }
};
