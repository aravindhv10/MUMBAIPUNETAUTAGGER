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
#include "TColor.h"
#include "TH2F.h"

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
    private:
        std::vector <NewHEPHeaders::VECTORS::DelphesVectors<>> * Vectors ;
    public:
        const double max_subjet_mass  , mass_drop_threshold , Rfilt      , minpt_subjet ,
        mh                            , mhmin               , mhmax      , zcut         , rcut_factor     ;
        const size_t nfilt            ;
        double       filteredjetmass  , deltah              , filt_tau_R , prunedmass   , unfiltered_mass ,
        Planar_Flow[3]                , EFC[3][6]           , frac_em    , frac_had     , nsub[3][5]      ,
        Spread[3]                     ;
        bool    HiggsTagged           ;
        size_t  n_tracks              ;
        indices index_constituents[3] ;
        NewHEPHeaders::pseudojets                 tau_subs  , t_parts  , tau_hadrons , tau_pieces   ;
        fastjet::PseudoJet                        prunedjet , triple   , Higgs       , taucandidate ;
        NewHEPHeaders::VECTORS::lorentz4vector <> ZVector   ;
    private:
        inline void   clear                      (                                     )       {
            t_parts.clear ()      ; tau_subs.clear ()  ; tau_hadrons.clear () ; tau_pieces.clear () ;
            filteredjetmass = 0.0 ; filt_tau_R =     0 ; prunedmass  = 0.0    ; n_tracks = 0.0      ;
            unfiltered_mass = 0.0 ; deltah     = 10000 ; HiggsTagged = false  ;
            Planar_Flow[0] = -10000.0 ; Planar_Flow[1] = -10000.0 ; Planar_Flow[2] = -10000.0 ;
            Spread[0]      = -10000.0 ; Spread[1]      = -10000.0 ; Spread[2]      = -10000.0 ;
            index_constituents[0].clear() ; index_constituents[1].clear() ; index_constituents[2].clear() ;
            ZVector.clearthis();
            for (size_t i=0;i<6;i++) {
                EFC[0][i]  = -10000.0 ;
                EFC[1][i]  = -10000.0 ;
                EFC[2][i]  = -10000.0 ;
            }
            for (size_t i=0;i<5;i++) {
                nsub[0][i] = -10000.0 ;
                nsub[1][i] = -10000.0 ;
                nsub[2][i] = -10000.0 ;
            }
            CPPFileIO::set_junked(Vectors);
        }
        inline size_t count_tracks               ( std::vector <size_t> & in_indices   ) const {
            size_t ret=0;
            for(size_t i=0;i<in_indices.size();i++){ret=ret+CPPFileIO::mymod(Vectors[0][in_indices[i]].Charge);}
            return ret;
        }
        inline double hadronic_energy            ( std::vector <size_t> & in_indices   ) const {
            double ret=0;
            for(size_t i=0;i<in_indices.size();i++){ret=ret+Vectors[0][in_indices[i]].Ehad;}
            return ret;
        }
        inline double electromagnetic_energy     ( std::vector <size_t> & in_indices   ) const {
            double ret=0;
            for(size_t i=0;i<in_indices.size();i++){ret=ret+Vectors[0][in_indices[i]].Eem;}
            return ret;
        }
        inline void   read_extra_variables       ( const DelphesReader & in            )       {
            n_tracks       = in.count_tracks           (index_constituents[0]) ;
            frac_em        = in.electromagnetic_energy (index_constituents[0]) / taucandidate.E () ;
            frac_had       = in.hadronic_energy        (index_constituents[0]) / taucandidate.E () ;
        }
        inline void   read_extra_variables       ( std::vector <NewHEPHeaders::VECTORS::DelphesVectors<>> & _Vectors ) {
            Vectors        = & _Vectors                                     ;
            n_tracks       = count_tracks           (index_constituents[0]) ;
            frac_em        = electromagnetic_energy (index_constituents[0]) / taucandidate.E () ;
            frac_had       = hadronic_energy        (index_constituents[0]) / taucandidate.E () ;
        }
        inline double Get_Spread                 ( const fastjet::PseudoJet & injet    )       {
            double ret = 0.0 ;
            NewHEPHeaders::pseudojets constituents = injet.constituents();
            size_t limit = constituents.size();
            for(size_t i=0;i<limit;i++){ret+=constituents[i].delta_R(injet)*constituents[i].perp();}
            return ret ;
        }
        inline double Get_Planar_Flow            ( const fastjet::PseudoJet & injet    )       {
            double ret = -10000.0 ;
            NewHEPHeaders::pseudojets constituents = injet.constituents();
            size_t limit = constituents.size();
            if(limit>2){
                std::vector <NewHEPHeaders::VECTORS::lorentz4vector<>> momentas ; /* Read in the 4 vectors: */ {
                    momentas.resize(limit);
                    for(size_t i=0;i<limit;i++)
                    {momentas[i]=NewHEPHeaders::VECTORS::lorentz4vector<>(constituents[i]);}
                }
                NewHEPHeaders::VECTORS::lorentz4vector <> jetvector(injet);
                double mJ = jetvector.m();
                NewHEPHeaders::VECTORS::euclid3vector<>JetAxis[3]; /* Get the 3 independent axis: */ {
                    /* Read the independent vectors: */ {
                        JetAxis[0] = jetvector.xyz.dir   () ;
                        JetAxis[1] = momentas[0].xyz.dir () ;
                        JetAxis[2] = momentas[1].xyz.dir () ;
                    }
                    /* Use gram-smidt procedure: */ {
                        JetAxis[1] = ( JetAxis[1] - (JetAxis[0]*(JetAxis[1]*JetAxis[0])) ).dir () ;
                        JetAxis[2] = ( JetAxis[2] - (JetAxis[0]*(JetAxis[2]*JetAxis[0])) ).dir () ;
                        JetAxis[2] = ( JetAxis[2] - (JetAxis[1]*(JetAxis[2]*JetAxis[1])) ).dir () ;
                    }
                }
                double matrix[2][2]; {matrix[0][0]=0;matrix[0][1]=0;matrix[1][0]=0;matrix[1][1]=0;}
                for(size_t i=0;i<momentas.size();i++){
                    NewHEPHeaders::VECTORS::euclid3vector<>tmp=momentas[i].xyz;
                    matrix[0][1] += (tmp*JetAxis[1]) * (tmp*JetAxis[2]) / momentas[i][3] ;
                    matrix[0][0] += (tmp*JetAxis[1]) * (tmp*JetAxis[1]) / momentas[i][3] ;
                    matrix[1][1] += (tmp*JetAxis[2]) * (tmp*JetAxis[2]) / momentas[i][3] ;
                }
                matrix[0][0]/=mJ; matrix[1][1]/=mJ; matrix[0][1]/=mJ;
                matrix[1][0]=matrix[0][1];
                ret = 4.0 * ( (matrix[0][0]*matrix[1][1]) - (matrix[1][0]*matrix[0][1]) ) ;
                double trace = (matrix[0][0]+matrix[1][1]) ;
                ret/=(trace*trace);
            }
            return ret;
        }
        inline void   get_constituent_indices    ( const fastjet::PseudoJet & this_jet )       {
            /* The Full Jet Part: */ {
                NewHEPHeaders::pseudojets vectors=this_jet.constituents();
                const size_t limit=vectors.size(); index_constituents[0].resize(limit);
                for (size_t i=0;i<limit;i++) {index_constituents[0][i]=vectors[i].user_index();}
            }
            /* The First SubJet part: */ if(tau_subs.size()>0) {
                NewHEPHeaders::pseudojets vectors=tau_subs[0].constituents();
                const size_t limit=vectors.size(); index_constituents[1].resize(limit);
                for (size_t i=0;i<limit;i++) {index_constituents[1][i]=vectors[i].user_index();}
            }
            /* The First SubJet part: */ if(tau_subs.size()>1) {
                NewHEPHeaders::pseudojets vectors=tau_subs[1].constituents();
                const size_t limit=vectors.size(); index_constituents[2].resize(limit);
                for (size_t i=0;i<limit;i++) {index_constituents[2][i]=vectors[i].user_index();}
            }
        }
        inline void   EvalEnergyCorrelation      ( fastjet::PseudoJet       & this_jet )       {
            using namespace fastjet          ;
            using namespace fastjet::contrib ;
            const double beta = 2.0          ;
            if(this_jet.constituents().size()>0){
                /* The energy correlation part: */ {
                    const auto   measure = EnergyCorrelator::pt_R ;
                    EnergyCorrelator ECF0 ( 0, beta, measure ) ; EFC[0][0] = ECF0 (this_jet) ;
                    EnergyCorrelator ECF1 ( 1, beta, measure ) ; EFC[0][1] = ECF1 (this_jet) ;
                    EnergyCorrelator ECF2 ( 2, beta, measure ) ; EFC[0][2] = ECF2 (this_jet) ;
                    EnergyCorrelator ECF3 ( 3, beta, measure ) ; EFC[0][3] = ECF3 (this_jet) ;
                    EnergyCorrelator ECF4 ( 4, beta, measure ) ; EFC[0][4] = ECF4 (this_jet) ;
                    EnergyCorrelator ECF5 ( 5, beta, measure ) ; EFC[0][5] = ECF5 (this_jet) ;
                }
                /* The NSubjettiness part: */ {
                    Nsubjettiness nSub1 ( 1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    Nsubjettiness nSub2 ( 2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    Nsubjettiness nSub3 ( 3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    Nsubjettiness nSub4 ( 4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    Nsubjettiness nSub5 ( 5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                    nsub[0][0]=nSub1(this_jet); nsub[0][1]=nSub2(this_jet); nsub[0][2]=nSub3(this_jet);
                    nsub[0][3]=nSub4(this_jet); nsub[0][4]=nSub5(this_jet);
                }
                Planar_Flow[0]=Get_Planar_Flow(this_jet);
                Spread[0]=Get_Spread(this_jet);
            }
        }
        inline void   read_reclustered_variables () {
            if(tau_subs.size()>0){
                Planar_Flow[1]=Get_Planar_Flow(tau_subs[0]);
                Spread[1]=Get_Spread(tau_subs[0]);
                using namespace fastjet ; using namespace fastjet::contrib ; const double beta = 2.0 ;
                Nsubjettiness nSub1 ( 1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                Nsubjettiness nSub2 ( 2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                Nsubjettiness nSub3 ( 3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                Nsubjettiness nSub4 ( 4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                Nsubjettiness nSub5 ( 5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                nsub[1][0]=nSub1(tau_subs[0]); nsub[1][1]=nSub2(tau_subs[0]); nsub[1][2]=nSub3(tau_subs[0]);
                nsub[1][3]=nSub4(tau_subs[0]); nsub[1][4]=nSub5(tau_subs[0]);
                const auto   measure = EnergyCorrelator::pt_R ;
                EnergyCorrelator ECF0 ( 0, beta, measure ) ; EFC[1][0] = ECF0 (tau_subs[0]) ;
                EnergyCorrelator ECF1 ( 1, beta, measure ) ; EFC[1][1] = ECF1 (tau_subs[0]) ;
                EnergyCorrelator ECF2 ( 2, beta, measure ) ; EFC[1][2] = ECF2 (tau_subs[0]) ;
                EnergyCorrelator ECF3 ( 3, beta, measure ) ; EFC[1][3] = ECF3 (tau_subs[0]) ;
                EnergyCorrelator ECF4 ( 4, beta, measure ) ; EFC[1][4] = ECF4 (tau_subs[0]) ;
                EnergyCorrelator ECF5 ( 5, beta, measure ) ; EFC[1][5] = ECF5 (tau_subs[0]) ;
            }
            if(tau_subs.size()>1){
                Planar_Flow[2]=Get_Planar_Flow(tau_subs[1]);
                Spread[2]=Get_Spread(tau_subs[1]);
                using namespace fastjet ; using namespace fastjet::contrib ; const double beta = 2.0 ;
                Nsubjettiness nSub1 ( 1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                Nsubjettiness nSub2 ( 2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                Nsubjettiness nSub3 ( 3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                Nsubjettiness nSub4 ( 4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                Nsubjettiness nSub5 ( 5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta) ) ;
                nsub[2][0]=nSub1(tau_subs[1]); nsub[2][1]=nSub2(tau_subs[1]); nsub[2][2]=nSub3(tau_subs[1]);
                nsub[2][3]=nSub4(tau_subs[1]); nsub[2][4]=nSub5(tau_subs[1]);
                const auto   measure = EnergyCorrelator::pt_R ;
                EnergyCorrelator ECF0 ( 0, beta, measure ) ; EFC[2][0] = ECF0 (tau_subs[1]) ;
                EnergyCorrelator ECF1 ( 1, beta, measure ) ; EFC[2][1] = ECF1 (tau_subs[1]) ;
                EnergyCorrelator ECF2 ( 2, beta, measure ) ; EFC[2][2] = ECF2 (tau_subs[1]) ;
                EnergyCorrelator ECF3 ( 3, beta, measure ) ; EFC[2][3] = ECF3 (tau_subs[1]) ;
                EnergyCorrelator ECF4 ( 4, beta, measure ) ; EFC[2][4] = ECF4 (tau_subs[1]) ;
                EnergyCorrelator ECF5 ( 5, beta, measure ) ; EFC[2][5] = ECF5 (tau_subs[1]) ;
            }
        }
        inline void   find_structures            ( const fastjet::PseudoJet & this_jet )       {
            fastjet::PseudoJet parent1(0,0,0,0), parent2(0,0,0,0);
            bool haskid=this_jet.validated_cs()->has_parents(this_jet,parent1,parent2);
            if(haskid) {
                if (parent1.m()<parent2.m()) {std::swap(parent1,parent2);}
                double kidmass    = parent1.m()  + parent2.m() ;
                double parentmass = this_jet.m()               ;
                if(kidmass<parentmass*mass_drop_threshold)
                { t_parts.push_back(parent1); t_parts.push_back(parent2); return; }
                else { find_structures(parent1); return; }
            } else {return;}
        }
        inline void   run_filter                 (                                     )       {
            t_parts    = sorted_by_pt  (t_parts)               ;
            triple     = fastjet::join (t_parts[0],t_parts[1]) ;
            filt_tau_R = std::min ( Rfilt , 0.5 * sqrt (t_parts[0].squared_distance(t_parts[1])) ) ;
            fastjet::JetDefinition filtering_def (fastjet::cambridge_algorithm,filt_tau_R) ;
            fastjet::Filter filter (filtering_def,fastjet::SelectorNHardest(nfilt)*fastjet::SelectorPtMin(minpt_subjet)) ;
            taucandidate    = filter         (triple) ;
            filteredjetmass = taucandidate.m ()       ;
            EvalEnergyCorrelation ( taucandidate )    ;
        }
        inline void   run_recluster              (                                     )       {
            get_constituent_indices ( taucandidate )  ;
            /* The pieces part: */ {
                fastjet::JetDefinition   reclustering (fastjet::cambridge_algorithm,10.0)  ;
                fastjet::ClusterSequence cs_top_sub   (taucandidate.pieces(),reclustering) ;
                tau_pieces=sorted_by_pt(cs_top_sub.exclusive_jets(2));
            }
            /* The Full SubJets part: */ {
                fastjet::JetDefinition   reclustering (fastjet::cambridge_algorithm,10.0)  ;
                fastjet::ClusterSequence cs_top_sub   (taucandidate.constituents(),reclustering) ;
                tau_subs=sorted_by_pt(cs_top_sub.exclusive_jets(2));
                read_reclustered_variables();
            }
        }
        inline void   run_variable_evaluater     ( fastjet::PseudoJet       & injet    )       {
            HiggsTagged  = true;
            Higgs        = tau_pieces[0]+tau_pieces[1];
            deltah       = CPPFileIO::mymod(taucandidate.m()-mh);
            tau_hadrons  = taucandidate.constituents();
            double Rprun = injet.validated_cluster_sequence()->jet_def().R();
            fastjet::JetDefinition jet_def_prune(fastjet::cambridge_algorithm, Rprun);
            fastjet::Pruner pruner(jet_def_prune,zcut,rcut_factor);
            prunedjet       = pruner      (triple) ;
            prunedmass      = prunedjet.m ()       ;
            unfiltered_mass = triple.m    ()       ;
        }
        inline void   run                        ( fastjet::PseudoJet       & injet    )       {
            clear(); find_structures(injet);
            if (t_parts.size()>1) {
                run_filter();
                if((mhmin<filteredjetmass)&&(filteredjetmass<mhmax)&&(taucandidate.pieces().size()>1)){
                    run_recluster();
                    if((tau_pieces[1].perp()>minpt_subjet)&&(tau_subs[1].perp()>minpt_subjet))
                    {run_variable_evaluater(injet);}
                }
            }
        }
    public:
        inline void   operator () (std::vector <NewHEPHeaders::VECTORS::DelphesVectors<>> & _Vectors)
        {read_extra_variables(_Vectors);}
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
        double filteredjetmass , deltah    , filt_tau_R , prunedmass , unfiltered_mass ,
        Planar_Flow[3]         , EFC[3][6] , frac_em    , frac_had   , nsub[3][5]      ,
        Spread[3]              ;
        bool    HiggsTagged ;
        size_t  n_tracks    ;
        NewHEPHeaders::VECTORS::lorentz4vector <> tau_subs[2] , t_parts[2]  , tau_hadrons[2] , tau_pieces[2] ;
        NewHEPHeaders::VECTORS::lorentz4vector <> prunedjet   , triple      , taucandidate   , Higgs         , ZVector ;
    private:
        inline void clear () {
            t_parts[0].clearthis () ; tau_subs[0].clearthis () ; tau_hadrons[0].clearthis () ; tau_pieces[0].clearthis () ;
            t_parts[1].clearthis () ; tau_subs[1].clearthis () ; tau_hadrons[1].clearthis () ; tau_pieces[1].clearthis () ;
            prunedjet.clearthis  () ; triple.clearthis      () ; taucandidate.clearthis   () ;
            Higgs.clearthis      () ; ZVector.clearthis     () ;
            Planar_Flow[0]    = 0.0 ; Planar_Flow[1] =     0.0 ; Planar_Flow[2]      =   0.0 ;
            Spread[0]         = 0.0 ; Spread[1]      =     0.0 ; Spread[2]           =   0.0 ;
            filteredjetmass   = 0.0 ; filt_tau_R     =     0.0 ; prunedmass          =   0.0 ; n_tracks = 0 ;
            unfiltered_mass   = 0.0 ; deltah         = 10000.0 ; HiggsTagged         = false ;
            for (size_t i=0;i<6;i++) {
                EFC[0][i]  = -10000.0 ;
                EFC[1][i]  = -10000.0 ;
                EFC[2][i]  = -10000.0 ;
            }
            for (size_t i=0;i<5;i++) {
                nsub[0][i] = -10000.0 ;
                nsub[1][i] = -10000.0 ;
                nsub[2][i] = -10000.0 ;
            }
        }
        inline void ReadFrom (HardSubStructureFinder&other) {
            Higgs           = other.Higgs           ;
            ZVector         = other.ZVector         ;
            prunedjet       = other.prunedjet       ;
            triple          = other.triple          ;
            taucandidate    = other.taucandidate    ;
            frac_em         = other.frac_em         ;
            frac_had        = other.frac_had        ;
            filteredjetmass = other.filteredjetmass ;
            filt_tau_R      = other.filt_tau_R      ;
            prunedmass      = other.prunedmass      ;
            n_tracks        = other.n_tracks        ;
            unfiltered_mass = other.unfiltered_mass ;
            deltah          = other.deltah          ;
            HiggsTagged     = other.HiggsTagged     ;
            Planar_Flow[0]  = other.Planar_Flow[0]  ;
            Planar_Flow[1]  = other.Planar_Flow[1]  ;
            Planar_Flow[2]  = other.Planar_Flow[2]  ;
            Spread[0]       = other.Spread[0]       ;
            Spread[1]       = other.Spread[1]       ;
            Spread[2]       = other.Spread[2]       ;
            for (size_t i=0;i<6;i++) {
                EFC[0][i]  = other.EFC[0][i]  ;
                EFC[1][i]  = other.EFC[1][i]  ;
                EFC[2][i]  = other.EFC[2][i]  ;
            }
            for (size_t i=0;i<5;i++) {
                nsub[0][i] = other.nsub[0][i] ;
                nsub[1][i] = other.nsub[1][i] ;
                nsub[2][i] = other.nsub[2][i] ;
            }
            for ( size_t i=0 ; i < CPPFileIO::mymin ( other.t_parts.size     () , (size_t)2 ) ; i++ )
            { t_parts     [i] = other.t_parts     [i] ; }
            for ( size_t i=0 ; i < CPPFileIO::mymin ( other.tau_subs.size    () , (size_t)2 ) ; i++ )
            { tau_subs    [i] = other.tau_subs    [i] ; }
            for ( size_t i=0 ; i < CPPFileIO::mymin ( other.tau_hadrons.size () , (size_t)2 ) ; i++ )
            { tau_hadrons [i] = other.tau_hadrons [i] ; }
            for ( size_t i=0 ; i < CPPFileIO::mymin ( other.tau_pieces.size  () , (size_t)2 ) ; i++ )
            { tau_pieces  [i] = other.tau_pieces  [i] ; }
        }
    public:
        inline double nsub_ratio (size_t i, size_t j=0) {
            double ret=0;
            if(nsub[j][i]>epsilon){ret=nsub[j][i+1]/nsub[j][i];}
            return ret;
        }
        inline double EFCDR (size_t i, size_t j=0) {
            double ret=0;
            if(EFC[j][i+1]>epsilon){ret=EFC[j][i+2]*EFC[j][i]/(EFC[j][i+1]*EFC[j][i+1]);}
            return ret;
        }
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
    private:
        inline void Analyze(size_t i){
            NewHEPHeaders::pseudojets&jetvectors=MainReader[0](i);
            fastjet::JetAlgorithm algorithm=fastjet::antikt_algorithm;
            const double jet_rad=1.00;
            fastjet::JetDefinition jetDef(algorithm,jet_rad);
            fastjet::ClusterSequence clust_seq(jetvectors,jetDef);
            NewHEPHeaders::pseudojets antikt_jets=sorted_by_pt(clust_seq.inclusive_jets());
            for (size_t j=0;(j<antikt_jets.size())&&(j<2);j++) if(antikt_jets[j].perp()>200.0) {
                // Recluster the jets:
                NewHEPHeaders::pseudojets constituents = antikt_jets[j].constituents();
                const double tmp_jet_rad=100.00;
                fastjet::JetDefinition tmp_jetDef(fastjet::cambridge_aachen_algorithm,tmp_jet_rad);
                fastjet::ClusterSequence tmp_clust_seq(constituents,tmp_jetDef);
                NewHEPHeaders::pseudojets CAJets_jets=sorted_by_pt(tmp_clust_seq.inclusive_jets());
                fastjet::PseudoJet&this_jet=CAJets_jets[0];
                if(this_jet.m()>40){
                    HardSubStructureFinder tmpslave ;
                    double mass = tmpslave(this_jet);
                    if(tmpslave.HiggsTagged){
                        tmpslave(MainReader[0].Vectors);
                        tmpslave.ZVector=MainReader[0].ZVector;
                        OutPutVariables tmp; tmp = tmpslave ;
                        Writer.push_back(tmp);
                    }
                }
            }
        }
        inline void Analyze(){
            size_t limit=MainReader[0]();
            for(size_t i=0;i<limit;i++){if((i%500)==0){printf("Analyzing %ld event:\n",i);}Analyze(i);}
        }
        inline void AnalyzeNewFile(std::string _delphesfilename){
            DelphesReader tmpreader(_delphesfilename);
            MainReader=&tmpreader; Analyze();
        }
    public:
        inline void operator()(std::string _delphesfilename){AnalyzeNewFile(_delphesfilename);}
        MainAnalyzer (std::string _OutFileName): OutFileName (_OutFileName), Writer(OutFileName) {}
        ~MainAnalyzer(){}
    };

    inline void ProcessType (std::string _name) {
        std::string name(_name)      ;
        std::string NoISRList   [16] ;
        std::string WithISRList [16] ;
        /* Get the names of root files: */ {
            for(size_t i=0;i<16;i++){
                char tmp[512] ;
                sprintf (tmp,"./DATA/%s/%ld/NoISRout.root",&(name[0]),i+1) ; NoISRList[i]   = std::string(tmp);
                sprintf (tmp,"./DATA/%s/%ld/out.root",&(name[0]),i+1)      ; WithISRList[i] = std::string(tmp);
            }
        }
        /* Prepare the output: */ {
            mkdir("./SKIM_DATA/",(mode_t)0755);
            char tmp[512] ;
            sprintf(tmp,"./SKIM_DATA/%s",&(name[0])); mkdir(tmp,(mode_t)0755);
            CPPFileIO::ForkMe forker;
            if(forker.InKid()){
                sprintf(tmp,"./SKIM_DATA/%s/NoMPI",&(name[0])); MainAnalyzer NoMPI(tmp);
                for(size_t i=0;i<16;i++){NoMPI(NoISRList[i]);}
            }
            if(forker.InKid()){
                sprintf(tmp,"./SKIM_DATA/%s/WithMPI",&(name[0])); MainAnalyzer WithMPI(tmp);
                for(size_t i=0;i<16;i++){WithMPI(WithISRList[i]);}
            }
        }
    }
    inline void ProcessAll () {
        CPPFileIO::ForkMe forker;
        if (forker.InKid()) { ProcessType ( "BoostedZToNuNuBar" ) ; }
        if (forker.InKid()) { ProcessType ( "BoostedZ"          ) ; }
        if (forker.InKid()) { ProcessType ( "UnBoostedZ"        ) ; }
        if (forker.InKid()) { ProcessType ( "BoostedZToBBbar"   ) ; }
    }
}

namespace Step2 {

    inline void PlotHist (std::string name, std::vector<float>&vals) {
        size_t limit=vals.size();
        if(limit>0){
            std::sort(vals.begin(),vals.end());
            TH1F hist(&(name[0]),&(name[0]),100,vals[0],vals[limit-1]);
            for(size_t i=0;i<limit;i++){hist.Fill(vals[i]);}
            TCanvas C;
            hist.Scale(1.0/hist.Integral());
            hist.SetLineWidth(3);
            int binmaxL = hist.GetMaximumBin ()        ;
            double xL   = hist.GetBinContent (binmaxL) ;
            double maxy = xL;
            hist.SetMaximum(maxy);
            hist.SetLineColor(TColor::GetColor("#990000"));
            hist.Draw("hist same");
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            std::string outname = "./GRAPHS/" + name + ".pdf" ;
            C.SaveAs(&(outname[0]));
        }
    }

    inline void PlotHist (
        std::string name,
        std::vector<float>&vals, std::vector<float>&vals2
    ) {
        TCanvas C                      ;
        size_t limit  = vals.size  ()  ;
        size_t limit2 = vals2.size ()  ;
        float xmin , xmax              ;
        std::string name2 = name + "2" ;
        /* Check for limits of histograms: */ {
            if(limit>0){
                std::sort ( vals.begin() , vals.end() ) ;
                xmin = vals [0]       ;
                xmax = vals [limit-1] ;
            }
            if(limit2>0){
                std::sort ( vals2.begin() , vals2.begin() ) ;
                xmin = CPPFileIO::mymin ( vals2 [0]       , xmin ) ;
                xmax = CPPFileIO::mymax ( vals2 [limit-1] , xmax ) ;
            }
        }
        TH1F hist  ( & ( name  [0] ) , & ( name  [0] ) , 100 , xmin , xmax ) ;
        TH1F hist2 ( & ( name2 [0] ) , & ( name2 [0] ) , 100 , xmin , xmax ) ;
        /* prepare the histograms: */ {
            /* Fill the histograms: */ {
                for ( size_t i = 0 ; i < vals.size  () ; i++ ) { hist.Fill  ( vals  [i] ) ; }
                for ( size_t i = 0 ; i < vals2.size () ; i++ ) { hist2.Fill ( vals2 [i] ) ; }
            }
            /* Rescale and color the histograms: */ {
                hist.Scale  ( 1.0 / hist.Integral  () ) ; hist.SetLineWidth  (3) ; hist.SetLineColor  (TColor::GetColor("#990000"));
                hist2.Scale ( 1.0 / hist2.Integral () ) ; hist2.SetLineWidth (3) ; hist2.SetLineColor (TColor::GetColor("#000099"));
            }
            /* Set the maximum */ {
                float x1 = hist.GetBinContent  ( hist.GetMaximumBin  () ) ;
                float x2 = hist2.GetBinContent ( hist2.GetMaximumBin () ) ;
                float x = CPPFileIO::mymax(x1,x2);
                hist.SetMaximum(x); hist2.SetMaximum(x);
            }
        }
        /* Draw and save the histograms: */ {
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            TCanvas C;
            hist.Draw  ("hist same") ;
            hist2.Draw ("hist same") ;
            name = "./GRAPHS/" + name + ".pdf";
            C.SaveAs(&(name[0]));
        }
    }

    inline void PlotHist (
        std::string name,
        std::vector<float>&vals, std::vector<float>&vals2, std::vector<float>&vals3
    ) {
        float xmin , xmax              ;
        std::string name2 = name + "2" ;
        std::string name3 = name + "3" ;
        size_t limit  = vals.size  ()  ;
        size_t limit2 = vals2.size ()  ;
        size_t limit3 = vals3.size ()  ;
        /* Check for limits of histograms: */ {
            if(limit>0){
                std::sort ( vals.begin() , vals.end() ) ;
                xmin = vals [0]       ;
                xmax = vals [limit-1] ;
            }
            if(limit2>0){
                std::sort ( vals2.begin() , vals2.begin() ) ;
                xmin = CPPFileIO::mymin ( vals2 [0]        , xmin ) ;
                xmax = CPPFileIO::mymax ( vals2 [limit2-1] , xmax ) ;
            }
            if(limit3>0){
                std::sort ( vals3.begin() , vals3.begin() ) ;
                xmin = CPPFileIO::mymin ( vals3 [0]        , xmin ) ;
                xmax = CPPFileIO::mymax ( vals3 [limit3-1] , xmax ) ;
            }
        }
        TH1F hist  ( & ( name  [0] ) , & ( name  [0] ) , 100 , xmin , xmax ) ;
        TH1F hist2 ( & ( name2 [0] ) , & ( name2 [0] ) , 100 , xmin , xmax ) ;
        TH1F hist3 ( & ( name3 [0] ) , & ( name3 [0] ) , 100 , xmin , xmax ) ;
        /* prepare the histograms: */ {
            /* Fill the histograms: */ {
                for ( size_t i = 0 ; i < vals.size  () ; i++ ) { hist.Fill  ( vals  [i] ) ; }
                for ( size_t i = 0 ; i < vals2.size () ; i++ ) { hist2.Fill ( vals2 [i] ) ; }
                for ( size_t i = 0 ; i < vals3.size () ; i++ ) { hist3.Fill ( vals3 [i] ) ; }
            }
            /* Rescale and color the histograms: */ {
                hist.Scale  ( 1.0 / hist.Integral  () ) ; hist.SetLineWidth  (3) ; hist.SetLineColor  (TColor::GetColor("#990000"));
                hist2.Scale ( 1.0 / hist2.Integral () ) ; hist2.SetLineWidth (3) ; hist2.SetLineColor (TColor::GetColor("#009900"));
                hist3.Scale ( 1.0 / hist3.Integral () ) ; hist3.SetLineWidth (3) ; hist3.SetLineColor (TColor::GetColor("#000099"));
            }
            /* Set the maximum */ {
                float x1 = hist.GetBinContent  ( hist.GetMaximumBin  () ) ;
                float x2 = hist2.GetBinContent ( hist2.GetMaximumBin () ) ;
                float x3 = hist3.GetBinContent ( hist3.GetMaximumBin () ) ;
                float x  = CPPFileIO::mymax ( CPPFileIO::mymax (x1,x2) , x3 ) ;
                hist.SetMaximum(x); hist2.SetMaximum(x); hist3.SetMaximum(x);
            }
        }
        /* Draw and save the histograms: */ {
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            TCanvas C;
            hist.Draw  ("hist same") ;
            hist2.Draw ("hist same") ;
            hist3.Draw ("hist same") ;
            name = "./GRAPHS/" + name + ".pdf";
            C.SaveAs(&(name[0]));
        }
    }

    inline void PlotHist (
        std::string name,
        std::vector<float>&vals, std::vector<float>&vals2, std::vector<float>&vals3, std::vector<float>&vals4,
        float xmin=0 , float xmax=-1
    ) {
        std::string name2 = name + "2" ;
        std::string name3 = name + "3" ;
        std::string name4 = name + "4" ;
        size_t limit  = vals.size  ()  ;
        size_t limit2 = vals2.size ()  ;
        size_t limit3 = vals3.size ()  ;
        size_t limit4 = vals4.size ()  ;
        /* Check for limits of histograms: */ {
            if(xmin>xmax) {
                if(limit>0){
                    std::sort ( vals.begin() , vals.end() ) ;
                    xmin = vals [0]       ;
                    xmax = vals [limit-1] ;
                }
                if(limit2>0){
                    std::sort ( vals2.begin() , vals2.begin() ) ;
                    xmin = CPPFileIO::mymin ( vals2 [0]        , xmin ) ;
                    xmax = CPPFileIO::mymax ( vals2 [limit2-1] , xmax ) ;
                }
                if(limit3>0){
                    std::sort ( vals3.begin() , vals3.begin() ) ;
                    xmin = CPPFileIO::mymin ( vals3 [0]        , xmin ) ;
                    xmax = CPPFileIO::mymax ( vals3 [limit3-1] , xmax ) ;
                }
                if(limit4>0){
                    std::sort ( vals4.begin() , vals4.begin() ) ;
                    xmin = CPPFileIO::mymin ( vals4 [0]        , xmin ) ;
                    xmax = CPPFileIO::mymax ( vals4 [limit4-1] , xmax ) ;
                }
            }
        }
        TH1F hist  ( & ( name  [0] ) , & ( name  [0] ) , 100 , xmin , xmax ) ;
        TH1F hist2 ( & ( name2 [0] ) , & ( name2 [0] ) , 100 , xmin , xmax ) ;
        TH1F hist3 ( & ( name3 [0] ) , & ( name3 [0] ) , 100 , xmin , xmax ) ;
        TH1F hist4 ( & ( name4 [0] ) , & ( name4 [0] ) , 100 , xmin , xmax ) ;
        /* prepare the histograms: */ {
            /* Fill the histograms: */ {
                for ( size_t i = 0 ; i < vals.size  () ; i++ ) { hist.Fill  ( vals  [i] ) ; }
                for ( size_t i = 0 ; i < vals2.size () ; i++ ) { hist2.Fill ( vals2 [i] ) ; }
                for ( size_t i = 0 ; i < vals3.size () ; i++ ) { hist3.Fill ( vals3 [i] ) ; }
                for ( size_t i = 0 ; i < vals4.size () ; i++ ) { hist4.Fill ( vals4 [i] ) ; }
            }
            /* Rescale and color the histograms: */ {
                hist.Scale  ( 1.0 / hist.Integral  () ) ; hist.SetLineWidth  (3) ; hist.SetLineColor  (TColor::GetColor("#990000"));
                hist2.Scale ( 1.0 / hist2.Integral () ) ; hist2.SetLineWidth (3) ; hist2.SetLineColor (TColor::GetColor("#009900"));
                hist3.Scale ( 1.0 / hist3.Integral () ) ; hist3.SetLineWidth (3) ; hist3.SetLineColor (TColor::GetColor("#000099"));
                hist4.Scale ( 1.0 / hist4.Integral () ) ; hist4.SetLineWidth (3) ; hist4.SetLineColor (TColor::GetColor("#000000"));
            }
            /* Set the maximum */ {
                float x1 = hist.GetBinContent  ( hist.GetMaximumBin  () ) ;
                float x2 = hist2.GetBinContent ( hist2.GetMaximumBin () ) ;
                float x3 = hist3.GetBinContent ( hist3.GetMaximumBin () ) ;
                float x4 = hist4.GetBinContent ( hist4.GetMaximumBin () ) ;
                float x  = CPPFileIO::mymax ( CPPFileIO::mymax ( CPPFileIO::mymax ( x1 , x2 ) , x3 ) , x4 ) ;
                hist.SetMaximum(x); hist2.SetMaximum(x); hist3.SetMaximum(x); hist4.SetMaximum(x);
            }
        }
        /* Draw and save the histograms: */ {
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            TCanvas C;
            hist.Draw  ("hist same") ;
            hist2.Draw ("hist same") ;
            hist3.Draw ("hist same") ;
            hist4.Draw ("hist same") ;
            name = "./GRAPHS/" + name + ".pdf";
            C.SaveAs(&(name[0]));
        }
    }

    inline void PlotHistLog (
        std::string name,
        std::vector<float>&vals, std::vector<float>&vals2, std::vector<float>&vals3, std::vector<float>&vals4,
        float xmin=0 , float xmax=-1
    ) {
        std::string name2 = name + "2" ;
        std::string name3 = name + "3" ;
        std::string name4 = name + "4" ;
        size_t limit  = vals.size  ()  ;
        size_t limit2 = vals2.size ()  ;
        size_t limit3 = vals3.size ()  ;
        size_t limit4 = vals4.size ()  ;
        /* Check for limits of histograms: */ {
            if(xmin>xmax) {
                if(limit>0){
                    std::sort ( vals.begin() , vals.end() ) ;
                    xmin = vals [0]       ;
                    xmax = vals [limit-1] ;
                }
                if(limit2>0){
                    std::sort ( vals2.begin() , vals2.begin() ) ;
                    xmin = CPPFileIO::mymin ( vals2 [0]        , xmin ) ;
                    xmax = CPPFileIO::mymax ( vals2 [limit2-1] , xmax ) ;
                }
                if(limit3>0){
                    std::sort ( vals3.begin() , vals3.begin() ) ;
                    xmin = CPPFileIO::mymin ( vals3 [0]        , xmin ) ;
                    xmax = CPPFileIO::mymax ( vals3 [limit3-1] , xmax ) ;
                }
                if(limit4>0){
                    std::sort ( vals4.begin() , vals4.begin() ) ;
                    xmin = CPPFileIO::mymin ( vals4 [0]        , xmin ) ;
                    xmax = CPPFileIO::mymax ( vals4 [limit4-1] , xmax ) ;
                }
            }
        }
        TH1F hist  ( & ( name  [0] ) , & ( name  [0] ) , 100 , xmin , xmax ) ;
        TH1F hist2 ( & ( name2 [0] ) , & ( name2 [0] ) , 100 , xmin , xmax ) ;
        TH1F hist3 ( & ( name3 [0] ) , & ( name3 [0] ) , 100 , xmin , xmax ) ;
        TH1F hist4 ( & ( name4 [0] ) , & ( name4 [0] ) , 100 , xmin , xmax ) ;
        /* prepare the histograms: */ {
            /* Fill the histograms: */ {
                for ( size_t i = 0 ; i < vals.size  () ; i++ ) if ( vals  [i] > -100.0 ) { hist.Fill  ( vals  [i] ) ; }
                for ( size_t i = 0 ; i < vals2.size () ; i++ ) if ( vals2 [i] > -100.0 ) { hist2.Fill ( vals2 [i] ) ; }
                for ( size_t i = 0 ; i < vals3.size () ; i++ ) if ( vals3 [i] > -100.0 ) { hist3.Fill ( vals3 [i] ) ; }
                for ( size_t i = 0 ; i < vals4.size () ; i++ ) if ( vals4 [i] > -100.0 ) { hist4.Fill ( vals4 [i] ) ; }
            }
            /* Rescale and color the histograms: */ {
                hist.Scale  ( 1.0 / hist.Integral  () ) ; hist.SetLineWidth  (3) ; hist.SetLineColor  (TColor::GetColor("#990000"));
                hist2.Scale ( 1.0 / hist2.Integral () ) ; hist2.SetLineWidth (3) ; hist2.SetLineColor (TColor::GetColor("#009900"));
                hist3.Scale ( 1.0 / hist3.Integral () ) ; hist3.SetLineWidth (3) ; hist3.SetLineColor (TColor::GetColor("#000099"));
                hist4.Scale ( 1.0 / hist4.Integral () ) ; hist4.SetLineWidth (3) ; hist4.SetLineColor (TColor::GetColor("#000000"));
            }
            /* Set the maximum */ {
                float x1 = hist.GetBinContent  ( hist.GetMaximumBin  () ) ;
                float x2 = hist2.GetBinContent ( hist2.GetMaximumBin () ) ;
                float x3 = hist3.GetBinContent ( hist3.GetMaximumBin () ) ;
                float x4 = hist4.GetBinContent ( hist4.GetMaximumBin () ) ;
                float x  = CPPFileIO::mymax ( CPPFileIO::mymax ( CPPFileIO::mymax ( x1 , x2 ) , x3 ) , x4 ) ;
                hist.SetMaximum(x); hist2.SetMaximum(x); hist3.SetMaximum(x); hist4.SetMaximum(x);
            }
        }
        /* Draw and save the histograms: */ {
            mkdir((const char*)"./GRAPHS",(mode_t)0755);
            TCanvas C; C.SetLogy(1);
            hist.Draw  ("hist same") ; hist2.Draw ("hist same") ;
            hist3.Draw ("hist same") ; hist4.Draw ("hist same") ;
            name = "./GRAPHS/" + name + ".pdf";
            C.SaveAs(&(name[0]));
        }
    }

    template <typename T> class Normalizer  {
    private:
        std::vector <T> values ;
        T L2Norm2, L2Norm ;
    public:
    private:
        inline T DotProduct (const Normalizer<T>&other) const {
            T ret = 0 ;
            if(values.size()!=other()){printf("U'r Retarded... sizes are different...\n");}
            else {for(size_t i=0;i<values.size();i++){ret+=values[i]*other(i);}}
            return ret;
        }
    public:
        inline T        operator *  ( const  Normalizer <T> & other ) const { return DotProduct(other) ; }
        inline T      & operator () ( size_t i                      )       { return values[i]         ; }
        inline T        operator () ( size_t i                      ) const { return values[i]         ; }
        inline size_t   operator () (                               ) const { return values.size()     ; }
        Normalizer(const std::vector<T>&_values){
            T mean=0;
            for(size_t i=0;i<_values.size();i++){mean+=_values[i];}
            mean/=(T)_values.size();
            values.resize(_values.size());
            for(size_t i=0;i<_values.size();i++){values[i]=_values[i]-mean;}
            L2Norm2=0;
            for(size_t i=0;i<values.size();i++){L2Norm2+=values[i]*values[i];}
            L2Norm=sqrt(L2Norm2);
            for(size_t i=0;i<values.size();i++){values[i]/=L2Norm;}
        }
        ~Normalizer(){}
    } ;

    class PlotAll2 {
    private:
        CPPFileIO::FileArray   <Step1::OutPutVariables>   reader1  ,   reader2  ,   reader3  ,   reader4  ;
        size_t                                            Limit1   ,   Limit2   ,   Limit3   ,   Limit4   ;
        Step1::OutPutVariables                          * element1 , * element2 , * element3 , * element4 ;

        inline void CalculateCorrelationSpread (size_t index) {
            std::vector <float> spreads, pfs;
            spreads.resize(Limit1); pfs.resize(Limit1);
            for(size_t i=0;i<Limit1;i++){
                spreads [i] = element1[i].Spread      [index] ;
                pfs     [i] = element1[i].Planar_Flow [index] ;
            }
            Normalizer <float> Nspreads ( spreads ) ;
            Normalizer <float> Npfs     ( pfs     ) ;
            printf("Spread Values(%ld): %e %e %e\n",index,Nspreads*Nspreads,Nspreads*Npfs,Npfs*Npfs);
        }
        inline void CalculateCorrelationNSub   (size_t index, size_t j=0) {
            std::vector <float> nsubs, pfs;
            nsubs.resize(Limit1); pfs.resize(Limit1);
            for(size_t i=0;i<Limit1;i++){nsubs[i]=element1[i].nsub[j][index];pfs[i]=element1[i].Planar_Flow[j];}
            Normalizer <float> Nnsubs ( nsubs ) ;
            Normalizer <float> Npfs   ( pfs   ) ;
            printf("Values(%ld,%ld): %e %e %e\n",index,j,Nnsubs*Nnsubs,Nnsubs*Npfs,Npfs*Npfs);
        }
        inline void CalculateCorrelation () {
            for(size_t j=0;j<3;j++){
                for(size_t i=1;i<=4;i++)
                {CalculateCorrelationNSub(i,j);}
                CalculateCorrelationSpread(j);
            }
        }

        inline void Plot2D (size_t index) {
            double max=20;
            if(index==2){max=10;}
            if(index>2){max=5;}
            /* Element-1 */ {
                char name[512]; sprintf(name,"T1NSub%ldvsPf",index);
                TH2F NSub1vsPf(name,name,110,-0.1,max,110,-0.1,1.1);
                for(size_t i=0;i<Limit1;i++)
                {NSub1vsPf.Fill((float)element1[i].nsub[0][index],(float)element1[i].Planar_Flow[0]);}
                std::string outname(name);
                mkdir((const char*)"./GRAPHS",(mode_t)0755);
                TCanvas C; NSub1vsPf.Draw("colz");
                outname = "./GRAPHS/" + outname + ".pdf";
                C.SaveAs(&(outname[0]));
            }
            /* Element-2 */ {
                char name[512]; sprintf(name,"T2NSub%ldvsPf",index);
                TH2F NSub1vsPf(name,name,110,-0.1,max,110,-0.1,1.1);
                for(size_t i=0;i<Limit2;i++)
                {NSub1vsPf.Fill((float)element2[i].nsub[0][index],(float)element2[i].Planar_Flow[0]);}
                std::string outname(name);
                mkdir((const char*)"./GRAPHS",(mode_t)0755);
                TCanvas C; NSub1vsPf.Draw("colz");
                outname = "./GRAPHS/" + outname + ".pdf";
                C.SaveAs(&(outname[0]));
            }
            /* Element-3 */ {
                char name[512]; sprintf(name,"T3NSub%ldvsPf",index);
                TH2F NSub1vsPf(name,name,110,-0.1,max,110,-0.1,1.1);
                for(size_t i=0;i<Limit3;i++)
                {NSub1vsPf.Fill((float)element3[i].nsub[0][index],(float)element3[i].Planar_Flow[0]);}
                std::string outname(name);
                mkdir((const char*)"./GRAPHS",(mode_t)0755);
                TCanvas C; NSub1vsPf.Draw("colz");
                outname = "./GRAPHS/" + outname + ".pdf";
                C.SaveAs(&(outname[0]));
            }
            /* Element-3 */ {
                char name[512]; sprintf(name,"T4NSub%ldvsPf",index);
                TH2F NSub1vsPf(name,name,110,-0.1,max,110,-0.1,1.1);
                for(size_t i=0;i<Limit4;i++)
                {NSub1vsPf.Fill((float)element4[i].nsub[0][index],(float)element4[i].Planar_Flow[0]);}
                std::string outname(name);
                mkdir((const char*)"./GRAPHS",(mode_t)0755);
                TCanvas C; NSub1vsPf.Draw("colz");
                outname = "./GRAPHS/" + outname + ".pdf";
                C.SaveAs(&(outname[0]));
            }
        }
        inline void Plot2D () {Plot2D(1);Plot2D(2);Plot2D(3);Plot2D(4);Plot2D(5);}

        inline void Plot_Masses () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for (size_t i=0;i<Limit1;i++) {Masses1[i]=element1[i].filteredjetmass;}
            for (size_t i=0;i<Limit2;i++) {Masses2[i]=element2[i].filteredjetmass;}
            for (size_t i=0;i<Limit3;i++) {Masses3[i]=element3[i].filteredjetmass;}
            for (size_t i=0;i<Limit4;i++) {Masses4[i]=element4[i].filteredjetmass;}
            PlotHist("Masses",Masses1,Masses2,Masses3,Masses4);
        }

        inline void Plot_Spread (size_t index) {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].Spread[index];}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].Spread[index];}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].Spread[index];}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].Spread[index];}
            char tmp[512]; sprintf(tmp,"Spread%ld",index);
            PlotHist(tmp,Masses1,Masses2,Masses3,Masses4);
        }
        inline void Plot_Spread () {Plot_Spread(0);Plot_Spread(1);Plot_Spread(2);}

        inline void Plot_EFrac () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].frac_em;}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].frac_em;}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].frac_em;}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].frac_em;}
            PlotHist("EFrac",Masses1,Masses2,Masses3,Masses4);
        }

        inline void Plot_HFrac () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].frac_had;}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].frac_had;}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].frac_had;}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].frac_had;}
            PlotHist("HFrac",Masses1,Masses2,Masses3,Masses4);
        }

        inline void Plot_NTracks () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].n_tracks;}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].n_tracks;}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].n_tracks;}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].n_tracks;}
            PlotHist("NTracks",Masses1,Masses2,Masses3,Masses4);
        }

        inline void Plot_nsub (size_t j, size_t k=0) {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].nsub[k][j];}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].nsub[k][j];}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].nsub[k][j];}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].nsub[k][j];}
            char tmp[512];
            if (k==0) { sprintf ( tmp , "NSub%ld"            , j+1     ) ; }
            else      { sprintf ( tmp , "NSub%ldOnSubjet%ld" , j+1 , k ) ; }
            PlotHist(tmp,Masses1,Masses2,Masses3,Masses4);
        }
        inline void Plot_nsub () {for(size_t i=0;i<5;i++){Plot_nsub(i,0);Plot_nsub(i,1);Plot_nsub(i,2);}}

        inline void Plot_nsub_ratio (size_t j) {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].nsub_ratio(j);}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].nsub_ratio(j);}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].nsub_ratio(j);}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].nsub_ratio(j);}
            char tmp[512];
            sprintf(tmp,"NSubRatio%ld",j+1);
            PlotHist(tmp,Masses1,Masses2,Masses3,Masses4,0.0,2.0);
        }
        inline void Plot_nsub_ratio ()
        {Plot_nsub_ratio(0);Plot_nsub_ratio(1);Plot_nsub_ratio(2);Plot_nsub_ratio(3);}

        inline void Plot_ECorrDR (size_t j) {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].EFCDR(j);}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].EFCDR(j);}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].EFCDR(j);}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].EFCDR(j);}
            char tmp[512]; sprintf(tmp,"ECorrDR%ld",j+1);
            PlotHist(tmp,Masses1,Masses2,Masses3,Masses4,0.0,2.0);
        }
        inline void Plot_ECorrDR ()
        {Plot_ECorrDR(0);Plot_ECorrDR(1);Plot_ECorrDR(2);Plot_ECorrDR(3);}

        inline void Plot_ECorr (size_t j, size_t k=0) {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].EFC[k][j];}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].EFC[k][j];}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].EFC[k][j];}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].EFC[k][j];}
            char tmp[512];
            if   (k==0) { sprintf ( tmp , "ECorr%ld"            , j+1     ) ; }
            else        { sprintf ( tmp , "ECorr%ldOnSubjet%ld" , j+1 , k ) ; }
            PlotHist(tmp,Masses1,Masses2,Masses3,Masses4);
        }
        inline void Plot_ECorr () {for(size_t i=0;i<6;i++){Plot_ECorr(i,0);Plot_ECorr(i,1);Plot_ECorr(i,2);}}

        inline void Plot_PlanarFlow0 () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].Planar_Flow[0];}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].Planar_Flow[0];}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].Planar_Flow[0];}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].Planar_Flow[0];}
            PlotHistLog("PlanarFlow",Masses1,Masses2,Masses3,Masses4,-0.0001,0.5);
        }
        inline void Plot_PlanarFlow1 () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].Planar_Flow[1];}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].Planar_Flow[1];}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].Planar_Flow[1];}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].Planar_Flow[1];}
            PlotHistLog("PlanarFlow1",Masses1,Masses2,Masses3,Masses4,-0.0001,1.5);
        }
        inline void Plot_PlanarFlow2 () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].Planar_Flow[1];}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].Planar_Flow[1];}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].Planar_Flow[1];}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].Planar_Flow[1];}
            PlotHistLog("PlanarFlow2",Masses1,Masses2,Masses3,Masses4,-0.0001,1.5);
        }
        inline void Plot_PlanarFlow () {Plot_PlanarFlow0();Plot_PlanarFlow1();Plot_PlanarFlow2();}

        inline void Plot_Z_Pt () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            std::vector <float> Masses3; Masses3.resize(Limit3);
            std::vector <float> Masses4; Masses4.resize(Limit4);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].ZVector.pt();}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].ZVector.pt();}
            for(size_t i=0;i<Limit3;i++){Masses3[i]=element3[i].ZVector.pt();}
            for(size_t i=0;i<Limit4;i++){Masses4[i]=element4[i].ZVector.pt();}
            PlotHistLog("GenLevelZPt",Masses1,Masses2,Masses3,Masses4);
        }

    public:

        PlotAll2():
        reader1 ( "./SKIM_DATA/BoostedZ/WithMPI"          ) ,
        reader2 ( "./SKIM_DATA/BoostedZToNuNuBar/WithMPI" ) ,
        reader3 ( "./SKIM_DATA/BoostedZToBBbar/WithMPI"   ) ,
        reader4 ( "./SKIM_DATA/UnBoostedZ/WithMPI"        ) ,
        Limit1(reader1.size()) , Limit2(reader2.size()) ,
        Limit3(reader3.size()) , Limit4(reader4.size()) ,
        element1(&(reader1(0,Limit1))) , element2(&(reader2(0,Limit2))) ,
        element3(&(reader3(0,Limit3))) , element4(&(reader4(0,Limit4))) {
            Plot_Masses     () ; Plot_EFrac       () ; Plot_HFrac       () ; Plot_NTracks () ;
            Plot_ECorrDR    () ; Plot_ECorr       () ; Plot_nsub_ratio  () ; Plot_nsub    () ;
            Plot_PlanarFlow () ; Plot_PlanarFlow1 () ; Plot_PlanarFlow2 () ; Plot_Z_Pt    () ;
            Plot_Spread     () ;
            Plot2D      (1) ; Plot2D      (2) ; Plot2D      (3) ; Plot2D (4) ; Plot2D (5) ;
            CalculateCorrelation();
        }
        ~PlotAll2(){}
    } ;

    class PlotAll {
    private:
        CPPFileIO::FileArray <Step1::OutPutVariables> reader1, reader2;
        size_t Limit1, Limit2;
        Step1::OutPutVariables *element1, *element2 ;
    public:
        inline void Plot_Masses () {
            std::vector <float> Masses1; Masses1.resize(Limit1);
            std::vector <float> Masses2; Masses2.resize(Limit2);
            for(size_t i=0;i<Limit1;i++){Masses1[i]=element1[i].frac_had;}
            for(size_t i=0;i<Limit2;i++){Masses2[i]=element2[i].frac_had;}
            PlotHist("Masses",Masses1,Masses2);
        }

        PlotAll():
        reader1 ( "./SKIM_DATA/BoostedZ/NoMPI"   ) , Limit1 (reader1.size()) , element1(&(reader1(0,Limit1))) ,
        reader2 ( "./SKIM_DATA/BoostedZ/WithMPI" ) , Limit2 (reader2.size()) , element2(&(reader2(0,Limit2)))
        {Plot_Masses();}

        ~PlotAll(){}
    } ;

}
