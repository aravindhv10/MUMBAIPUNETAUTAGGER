#include "./NewHEPHeaders4.hh"

// Copied form hep top tagger source file:
class HardSubStructureFinder {
private:
    inline void find_structures (const fastjet::PseudoJet&this_jet) {
        fastjet::PseudoJet parent1(0,0,0,0), parent2(0,0,0,0);
        if (this_jet.m() < max_subjet_mass || !this_jet.validated_cs()->has_parents(this_jet, parent1, parent2)) {t_parts.push_back(this_jet);}
        else {
            if (parent1.m()<parent2.m()) {std::swap(parent1, parent2);}
            find_structures(parent1);
            if (parent1.m() < mass_drop_threshold * this_jet.m())
                find_structures(parent2);
        }
    }
    inline void initialize () {
        max_subjet_mass=30;
        mass_drop_threshold=0.7;
    }
public:
    NewHEPHeaders::pseudojets t_parts;
    double max_subjet_mass, mass_drop_threshold;
    HardSubStructureFinder(){initialize();}
    ~HardSubStructureFinder(){}
};
