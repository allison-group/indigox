#include <iostream>

#include <indigox/indigox.hpp>

int main() {
    using namespace indigox;
  Options::DATA_DIRECTORY = "/Users/iwelsh/GitHub/indigox/data";
    // Set the options to use
    Options::AssignElectrons::ALGORITHM = Options::AssignElectrons::Algorithm::FPT;
    Options::AssignElectrons::FPT::ADD_EDGES_TO_TD = false;
    Options::AssignElectrons::FPT::MINIMUM_PROPAGATION_DEPTH = 1;
    Options::AssignElectrons::USE_ELECTRON_PAIRS = false;

    // Prepare elements
    PeriodicTable PT = IXPeriodicTable::GetInstance();
    Element H = PT->GetElement("H");
    Element C = PT->GetElement("C");
    Element O = PT->GetElement("O");
    Element N = PT->GetElement("N");

    // Build the molecule
    Molecule m = Molecule(new IXMolecule());
    m->SetTotalCharge(-1);

  Atom a1 = m->NewAtom();
  std::cout << a1->ToString() << std::endl;
  
    // Add the atoms
    Atom c1  = m->NewAtom(C);  c1->SetName("C1") ;
    Atom c2  = m->NewAtom(C);  c2->SetName("C2") ;
    Atom c3  = m->NewAtom(C);  c3->SetName("C3") ;
    Atom c4  = m->NewAtom(C);  c4->SetName("C4") ;
    Atom c5  = m->NewAtom(C);  c5->SetName("C5") ;
    Atom c6  = m->NewAtom(C);  c6->SetName("C6") ;
    Atom h7  = m->NewAtom(H);  h7->SetName("H7") ;
    Atom h8  = m->NewAtom(H);  h8->SetName("H8") ;
    Atom n9  = m->NewAtom(N);  n9->SetName("N9") ;
    Atom h10 = m->NewAtom(H); h10->SetName("H10");
    Atom h11 = m->NewAtom(H); h11->SetName("H11");
    Atom c12 = m->NewAtom(C); c12->SetName("C12");
    Atom o13 = m->NewAtom(O); o13->SetName("O13");
    Atom o14 = m->NewAtom(O); o14->SetName("O14");
    Atom o15 = m->NewAtom(O); o15->SetName("O15");
    Atom o16 = m->NewAtom(O); o16->SetName("O16");

    // Add the bonds
    m->NewBond(c1,c2);
    m->NewBond(c1,c6);
    m->NewBond(c1,h7);
    m->NewBond(c2,c3);
    m->NewBond(c2,h8);
    m->NewBond(c3,c4);
    m->NewBond(c3,n9);
    m->NewBond(c4,c5);
    m->NewBond(c4,h10);
    m->NewBond(c5,c6);
    m->NewBond(c5,h11);
    m->NewBond(c6,c12);
    m->NewBond(n9,o13);
    m->NewBond(n9,o14);
    m->NewBond(c12,o15);
    m->NewBond(c12,o16);
    
    // Print out some information
    std::cout << "Number of atoms: " << m->NumAtoms() << ", number of bonds: " << m->NumBonds() << "\n";

    // Calculate bond orders and formal charges
    size_t count = m->AssignElectrons();
    std::cout << count << " resonance structure(s) calculated with a score of " << m->GetMinimumElectronAssignmentScore() << ".\n";

    // Print out each of the structures
    for (size_t i = 0; i < count; ++i) {
        m->ApplyElectronAssignment(i);
        for (MolAtomIterator it = m->BeginAtom(); it != m->EndAtom(); ++it) {
            Atom at = *it;
            if (at->GetFormalCharge() != 0)
                std::cout << "Atom " << at->GetName() << " has a formal charge of " << at->GetFormalCharge() << ".\n";
        }
        for (MolBondIterator it = m->BeginBond(); it != m->EndBond(); ++it) {
            Bond bt = *it;
            if (bt->GetOrder() != 1)
                std::cout << "Bond between " << bt->GetSourceAtom()->GetName() << " and " << bt->GetTargetAtom()->GetName() << " has an order of " << bt->GetOrder() << ".\n";
        }
        std::cout << std::endl;
    }

    return 0;
}

