#include <sstream>
#include <iomanip>
#include <fstream>
#include <exception>
#include "Pythia8/Pythia.h"
//#include "Pythia8Plugins/HepMC2.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

template <typename T> inline T GetMod (T a) {
    if (a<0) {return -a;}
    else {return a;}
}

int main () {
    Pythia8::Pythia pythia ;
    pythia.readString("HiggsSM:ffbar2HZ = on");
    pythia.readString("PhaseSpace:pTHatMin = 500");
    pythia.readString("15:onMode = off");
    pythia.readString("15:onIfAny = 111 211");
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 12");
    pythia.readString("25:onMode = off");
    pythia.readString("25:onIfAny = 15");
    pythia.init();
    fastjet::JetDefinition jetDef (fastjet::antikt_algorithm,0.8) ;
    FILE *f = fopen("./pts","w");
    FILE *mass = fopen("./mass","w");
    for(size_t i=0;i<10000;i++) if (pythia.next()) {
        std::vector <fastjet::PseudoJet> pseudojets ;
        for(size_t j=0;j<pythia.event.size();j++) if(pythia.event[j].isFinal()) {
            size_t mpid = GetMod(pythia.event[j].id()) ;
            if((mpid!=12)&&(mpid!=14)&&(mpid!=16)) {
                pseudojets.push_back(
                    fastjet::PseudoJet(pythia.event[j].px(),pythia.event[j].py(),pythia.event[j].pz(),pythia.event[j].e())
                );
            }
        }
        fastjet::ClusterSequence clustSeq (pseudojets,jetDef) ;
        std::vector <fastjet::PseudoJet> sortedJets;
        sortedJets = sorted_by_pt(clustSeq.inclusive_jets(200.0));
        if(sortedJets.size()>0) {
            fastjet::PseudoJet hig=sortedJets[0];
            fprintf(mass,"%e\n",hig.m());
        }
        for(size_t k=0;k<sortedJets.size();k++) {
            fprintf(f,"%e ",sortedJets[k].pt());
        }
        fprintf(f,"\n");
    }
    fclose(f);
    fclose(mass);
    return 0;
}
