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
                    for (size_t i=0;i<4;i++) if (EFC[i+1]>epsilon) {EFCDR[i]=(EFC[i]+EFC[i+2])/EFC[i+1];}
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
            for (size_t i=0;i<CPPFileIO::mymin(other.t_parts.size(),(size_t)2);i++)     {t_parts[i]=other.t_parts[i];}
            for (size_t i=0;i<CPPFileIO::mymin(other.tau_subs.size(),(size_t)2);i++)    {tau_subs[i]=other.tau_subs[i];}
            for (size_t i=0;i<CPPFileIO::mymin(other.tau_hadrons.size(),(size_t)2);i++) {tau_hadrons[i]=other.tau_hadrons[i];}
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
                    tmpslave(this_jet);
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
        inline void WriteHistograms(){
        }
    public:
        inline void operator()(std::string _delphesfilename){AnalyzeNewFile(_delphesfilename);}
        MainAnalyzer (std::string _OutFileName): OutFileName (_OutFileName), Writer(OutFileName) {}
        ~MainAnalyzer(){WriteHistograms();}
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
            sprintf(tmp,"./SKIM_DATA/%s/NoMPI",&(name[0])); MainAnalyzer NoMPI(tmp);
            sprintf(tmp,"./SKIM_DATA/%s/WithMPI",&(name[0])); MainAnalyzer WithMPI(tmp);
            CPPFileIO::ForkMe forker;
            if(forker.InKid()){for(size_t i=0;i<16;i++){NoMPI(NoISRList[i]);}}
            if(forker.InKid()){for(size_t i=0;i<16;i++){WithMPI(WithISRList[i]);}}
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
    inline void PlotAll (std::string filename) {
        CPPFileIO::FileArray <Step1::OutPutVariables> reader(filename);
        Step1::OutPutVariables*element;
        size_t Limit = reader.size();
        element=&(reader(0,Limit));
        TH1F Masses("Masses","Masses",150,-0.1,150.1);
        for(size_t i=0;i<Limit;i++){Masses.Fill(element[i].filteredjetmass);}
        TCanvas C;
        Masses.Draw();
        filename=filename+".pdf";
        C.SaveAs(&(filename[0]));
    }
}
