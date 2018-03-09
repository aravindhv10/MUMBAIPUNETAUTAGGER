#include "./NewHEPHeaders4.hh"
class GenerateHiggs {
private:
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
public:
    int GenSim ( int ThID ) {
        /* Create the directory structure: */ {
            char tmp[256]; sprintf(tmp,"hepmc/%d",ThID);
            mkdir("hepmc",0755); mkdir(tmp,0755);
            sprintf(tmp,"hepmc/%d/hepmc.fifo",ThID);
            mkfifo(tmp,(mode_t)0755);
        }
        CPPFileIO::ForkMe forker ;
        //printf("Now Starting generate...\n");
        if(forker.InKid()) {generate(ThID);}
        if(forker.InKid()) {RunDelphes(ThID);}
        return 0;
    }
    GenerateHiggs(){}
    ~GenerateHiggs(){}
};
class GenerateQCD {
private:
    int generate ( int ThID ) {
        Pythia8::Pythia pythia ; /**/ {
            /* The random number part: */ {
                pythia.readString("Random:setSeed = on");
                char tmp[256];
                sprintf(tmp,"Random:seed = %d",ThID);
                pythia.readString(tmp);
            }
            /* The event part: */ {
                pythia.readString("HardQCD:all = on");
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
public:
    int GenSim ( int ThID ) {
        /* Create the directory structure: */ {
            char tmp[256]; sprintf(tmp,"hepmc/%d",ThID);
            mkdir("hepmc",0755); mkdir(tmp,0755);
            sprintf(tmp,"hepmc/%d/hepmc.fifo",ThID);
            mkfifo(tmp,(mode_t)0755);
        }
        CPPFileIO::ForkMe forker ;
        if(forker.InKid()) {generate(ThID);}
        if(forker.InKid()) {RunDelphes(ThID);}
        return 0;
    }
    GenerateQCD(){}
    ~GenerateQCD(){}
};

class Analyzer {
private:
    TH1F MassHist ;
    inline void ProcessData (NewHEPHeaders::DELPHES_DETDATA::FullDelphesContainer&indata) {
        if(indata.detinfo.jets.size()>0) {
            if(indata.detinfo.jets[0].MainHiggsTauTagger.HiggsTagged){
                MassHist.Fill(indata.detinfo.jets[0].MainHiggsTauTagger.filteredjetmass);
                //MassHist.Fill(indata.detinfo.jets[0].MainHiggsTauTagger.prunedmass);
                NewHEPHeaders::vector4 tmpbuf(indata.detinfo.jets[0]);
                //if(indata.geninfo.Higgs(tmpbuf)<0.3){printf("Real Jet...\n");}
                //else{printf("Wrong jet...\n");}
            }
        }
    }
public:
    template <typename T> inline void operator () (T&intree) {
        if (intree.fChain == 0) {return;}
        long nentries = intree.fChain->GetEntriesFast();
        long nbytes = 0, nb = 0;
        for (long jentry=0;jentry<nentries;jentry++) {
            long ientry = intree.LoadTree(jentry);
            if (ientry < 0) {break;}
            nb = intree.fChain->GetEntry (jentry) ;
            nbytes = nbytes + nb ;
            NewHEPHeaders::DELPHES_DETDATA::FullDelphesContainer datareader ;
            datareader.ReadFromDelphes(intree);
            ProcessData(datareader);
        }
    }
    Analyzer() : MassHist("MassHist","MassHist",100,20,160) {}
    ~Analyzer() {
        TCanvas C;
        MassHist.Draw();
        C.SaveAs("MassHist.pdf");
    }
};
