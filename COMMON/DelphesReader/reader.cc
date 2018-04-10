#include <iostream>
#include <utility>
#include <vector>
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
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
#include "./NewHEPHeaders4.hh"

using namespace NewHEPHeaders;

class DelphesReader {
private:
    std::string      & ListOfFiles        ;
    TChain             chain              ;
    ExRootTreeReader * treeReader         ;
    size_t             numberOfEntries    ;
    TClonesArray     * EFlowTrack         ;
    TClonesArray     * EFlowPhoton        ;
    TClonesArray     * EFlowNeutralHadron ;
    inline void Analyze (size_t entry) {
        treeReader->ReadEntry(entry);
        size_t limit_EFlowTrack = EFlowTrack->GetEntries() ;
        size_t limit_EFlowPhoton = EFlowPhoton->GetEntries() ;
        size_t limit_EFlowNeutralHadron = EFlowNeutralHadron->GetEntries();
        for(size_t i=0;i<limit_EFlowTrack;i++){
            Track*tmp=(Track*)EFlowTrack->At(i);
        }
        for(size_t i=0;i<limit_EFlowPhoton;i++){
            Tower*tmp=(Tower*)EFlowPhoton->At(i);
        }
        for(size_t i=0;i<limit_EFlowNeutralHadron;i++){
            Tower*tmp=(Tower*)EFlowNeutralHadron->At(i);
        }
    }
    inline void construct(){
        chain.Add(&(ListOfFiles[0]));
        numberOfEntries=treeReader->GetEntries();
        EFlowTrack         = treeReader->UseBranch ( "EFlowTrack"         ) ;
        EFlowPhoton        = treeReader->UseBranch ( "EFlowPhoton"        ) ;
        EFlowNeutralHadron = treeReader->UseBranch ( "EFlowNeutralHadron" ) ;
    }
    inline void destroy(){}
public:

    DelphesReader (std::string&_ListOfFiles): ListOfFiles(_ListOfFiles), chain("Delphes") {construct();}
    ~DelphesReader() {destroy();}
};
