#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
using namespace Pythia8;
int main(int argc, char ** argv) {
    if(argc!=3){
        printf("Wrong usage...\n");
        printf("Usage: %s <LHE File> <HepMC File>\n",argv[0]);
    } else {
        Pythia pythia; pythia.readString("Beams:frameType = 4"); /**/ {
            char buf[1024]; sprintf(buf,"Beams:LHEF = %s",argv[1]);
            pythia.readString(buf); pythia.init();
        }
        HepMC::Pythia8ToHepMC ToHepMC; mkfifo((const char*)argv[2],(mode_t)0755);
        HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);
        for (size_t iEvent=0;!pythia.info.atEndOfFile();iEvent++) if (pythia.next()) {
            HepMC::GenEvent*hepmcevt=new HepMC::GenEvent();
            ToHepMC.fill_next_event(pythia,hepmcevt);
            ascii_io<<hepmcevt; delete hepmcevt;
        }
        pythia.stat();
    }
    return 0;
}
