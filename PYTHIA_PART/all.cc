#include "./NewHEPHeaders4.hh"
const size_t MemTotal = 16321488 ;
const size_t processor = 7 ;
class Generator {
private:
    inline void MakeStruct ( int ThID ) {
        char tmp[256]; sprintf(tmp,"hepmc/%d",ThID);
        mkdir("hepmc",0755); mkdir(tmp,0755);
        sprintf(tmp,"hepmc/%d/hepmc.fifo",ThID);
        mkfifo(tmp,(mode_t)0755);
    }
    inline void RunDelphes ( int ThID ) {
        char tmp[256] ;
        sprintf(tmp,"%d",ThID);
        execl("./RunDelphes","./RunDelphes",tmp,NULL);
    }
    inline void generate ( int ThID , Pythia8::Pythia&pythia ) {
        char tmp[1024] ;
        /* Common Phase Space part: */ {
            if(pthatmin>0){
                char tmp[1024]; sprintf(tmp,"PhaseSpace:pTHatMin = %ld",pthatmin);
                pythia.readString(tmp);
            }
        }
        /* The random number part: */ {
            pythia.readString("Random:setSeed = on");
            sprintf(tmp,"Random:seed = %d",ThID);
            pythia.readString(tmp);
        }
        /* The Common decay part: */ {
            pythia.readString("15:onMode = off");
            pythia.readString("15:onIfAny = 111 211");
        }
        pythia.init();
        sprintf(tmp,"./hepmc/%d/hepmc.fifo",ThID);
        NewHEPHeaders::WriteHepmc2Fifo writer(tmp);
        for(size_t i=0;i<10000;i++) if (pythia.next()) {writer(pythia);}
    }
    inline int GenSim ( int ThID , Pythia8::Pythia&pythia ) {
        MakeStruct(ThID); CPPFileIO::ForkMe forker ;
        if(forker.InKid()) {generate(ThID,pythia);}
        if(forker.InKid()) {RunDelphes(ThID);}
        return 0;
    }
    inline void GenerateHiggs ( int ThID ) {
        Pythia8::Pythia pythia ; /* Configure pythia: */ {
            pythia.readString("HiggsSM:ffbar2HZ = on");
            /* The Common decay part: */ {
                pythia.readString("23:onMode = off");
                pythia.readString("23:onIfAny = 12");
                pythia.readString("25:onMode = off");
                pythia.readString("25:onIfAny = 15");
            }
        }
        GenSim(ThID,pythia);
    }
    inline void GenerateQCD ( int ThID ) {
        Pythia8::Pythia pythia ; /* Configure pythia: */ {
            pythia.readString("HardQCD:all = on");
        }
        GenSim(ThID,pythia);
    }
    inline void GenerateZ ( int ThID ) {
        Pythia8::Pythia pythia ; /* Configure pythia: */ {
            pythia.readString("WeakZ0:gmZmode = 2");
            pythia.readString("WeakBosonAndParton:ffbar2gmZgm = on");
            /* The Common decay part: */ {
                pythia.readString("23:onMode = off");
                pythia.readString("23:onIfAny = 15");
            }
        }
        GenSim(ThID,pythia);
    }
public:
    long pthatmin ;
    inline void GenerateHiggs(){
        #pragma omp parallel for
        for(size_t i=1;i<=(processor+1);i++){GenerateHiggs(i);}
    }
    inline void GenerateQCD(){
        #pragma omp parallel for
        for(size_t i=1;i<=(processor+1);i++){GenerateQCD(i);}
    }
    inline void GenerateZ(){
        #pragma omp parallel for
        for(size_t i=1;i<=(processor+1);i++){GenerateZ(i);}
    }
    Generator  () {pthatmin=-1;}
    ~Generator () {}
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
class AnalyzerZ{
private:
    TH1F MassHist ;
    template <typename T> inline void ProcessData (T&indata) {
        NewHEPHeaders::EventData reader ; reader.ReadFromDelphes(indata); reader.prepare();
        fastjet::JetDefinition jet_def_fat_jet(fastjet::cambridge_aachen_algorithm,1.0) ;
        fastjet::ClusterSequence clust_seq_nrmjet(reader.tojets,jet_def_fat_jet);
        NewHEPHeaders::pseudojets jets = sorted_by_pt(clust_seq_nrmjet.inclusive_jets(100.0));
        if(jets.size()>0){
            NewHEPHeaders::HardSubStructureFinder tagger; tagger(jets[0]);
            MassHist.Fill(tagger.filteredjetmass);
        }
        return;
        NewHEPHeaders::vector4s taus;
        for(size_t i=0;i<indata.Jet_;i++){
            if(indata.Jet_TauTag[i]==1){
                TLorentzVector tmp; /* Read the jet vectors: */ {
                    tmp.SetPtEtaPhiM (
                        indata.Jet_PT[i], indata.Jet_Eta[i], indata.Jet_Phi[i], 0
                    );
                }
                NewHEPHeaders::vector4 tmpvector; /* Read it into something better: */ {
                    tmpvector[0] = tmp.Px () ;
                    tmpvector[1] = tmp.Py () ;
                    tmpvector[2] = tmp.Pz () ;
                    tmpvector[3] = tmp.E  () ;
                }
                taus.push_back(tmpvector);
            }
        }
        if(taus.size()>1){MassHist.Fill((taus[0]+taus[1]).m());}
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
            ProcessData(intree);
        }
    }
    AnalyzerZ  () : MassHist("MassHist","MassHist",100,20,160) {}
    ~AnalyzerZ () { TCanvas C; MassHist.Draw(); C.SaveAs("MassHist.pdf"); }
};
