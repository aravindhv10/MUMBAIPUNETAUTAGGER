#include "./NewHEPHeaders4.hh"
int generate ( int ThID ) {
    Pythia8::Pythia pythia ; /**/ {
        /* The random number part: */ {
            pythia.readString("Random:setSeed = on");
            char tmp[256];
            sprintf(tmp,"Random:seed = %d",ThID);
            pythia.readString(tmp);
        }
        /* The event part: */ {
            pythia.readString("HiggsSM:ffbar2HZ = on");
            pythia.readString("PhaseSpace:pTHatMin = 500");
        }
        /* The decay part: */ {
            pythia.readString("15:onMode = off");
            pythia.readString("15:onIfAny = 111 211");
            pythia.readString("23:onMode = off");
            pythia.readString("23:onIfAny = 12");
            pythia.readString("25:onMode = off");
            pythia.readString("25:onIfAny = 15");
        }
        pythia.init();
    }
    char tmp[1024] ;
    sprintf(tmp,"./hepmc/%d/hepmc.fifo",ThID);
    NewHEPHeaders::WriteHepmc2Fifo writer(tmp);
    for(size_t i=0;i<10000;i++) if (pythia.next()) {writer(pythia);}
    return 0;
}
int RunDelphes ( int ThID ) {
    char tmp[256] ;
    sprintf(tmp,"%d",ThID);
    execl("./RunDelphes","./RunDelphes",tmp,NULL);
    return 0;
}
int GenSim ( int ThID ) {
    /* Create the directory structure: */ {
        CPPFileIO::ToDir dirchanger (true) ;
        dirchanger("hepmc"); dirchanger(ThID);
        mkfifo("./hepmc.fifo",(mode_t)0755);
    }
    CPPFileIO::ForkMe forker ;
    if(forker.InKid()) {generate(ThID);}
    if(forker.InKid()) {RunDelphes(ThID);}
    return 0;
}
