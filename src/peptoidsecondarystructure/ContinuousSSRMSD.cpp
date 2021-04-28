/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ContinuousSSRMSD.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"

namespace PLMD {
namespace peptoidsecondarystructure {

void ContinuousSSRMSD::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the alpha-beta variables you wish to calculate. "
                              "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one alpha-beta values will be "
                              "calculated for each ATOM keyword you specify (all ATOM keywords should "
                              "specify the indices of continuous backbone atoms).  The eventual number of quantities calculated by this "
                              "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
  keys.add("compulsory","TYPE","DRMSD","the manner in which RMSD alignment is performed. Should be OPTIMAL, SIMPLE or DRMSD. "
                                       "For more details on the OPTIMAL and SIMPLE methods see \\ref RMSD. For more details on the "
                                       "DRMSD method see \\ref DRMSD.");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions");
  keys.add("compulsory","R_0","0.08","The r_0 parameter of the switching function.");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","NN","8","The n parameter of the switching function");
  keys.add("compulsory","MM","12","The m parameter of the switching function");
  keys.reserve("optional","STRANDS_CUTOFF","If in a segment of protein the two strands are further apart then the calculation "
               "of the actual RMSD is skipped as the structure is very far from being beta-sheet like. "
               "This keyword speeds up the calculation enormously when you are using the LESS_THAN option. "
               "However, if you are using some other option, then this cannot be used");
  keys.addFlag("VERBOSE",false,"write a more detailed output");
  keys.add("hidden","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities "
           "that contributed less than TOL at the previous neighbor list update step are ignored.");
  ActionWithVessel::registerKeywords( keys );
  keys.use("LESS_THAN"); keys.use("MIN"); keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
  keys.setComponentsIntroduction("By default this Action calculates the number of structural units that are within a certain "
                                 "distance of a idealized secondary structure element. This quantity can then be referenced "
                                 "elsewhere in the input by using the label of the action. However, this Action can also be used to "
                                 "calculate the following quantities by using the keywords as described below.  The quantities then "
                                 "calculated can be referenced using the label of the action followed by a dot and then the name "
                                 "from the table below.  Please note that you can use the LESS_THAN keyword more than once.  The resulting "
                                 "components will be labelled <em>label</em>.lessthan-1, <em>label</em>.lessthan-2 and so on unless you "
                                 "exploit the fact that these labels can be given custom labels by using the LABEL keyword in the "
                                 "description of you LESS_THAN function that you are computing");
}


ContinuousSSRMSD::ContinuousSSRMSD(const ActionOptions&ao):
  Action(ao),
  SecondaryStructureRMSD(ao) { }

void ContinuousSSRMSD::init() {
  // read in the backbone atoms
  std::vector<unsigned> chains;
  getBackboneChains("peptoid", chains);

  // get reference structure from whatever class is deriving us.
  std::vector<Vector> refstructure = getRefStructure();
  plumed_massert(refstructure.size() != 0, "reference structure must not be empty");
  plumed_massert(refstructure.size() % ATOMS_IN_BB_RES == 0,
                 "reference structure be a multiple of " + std::to_string(ATOMS_IN_BB_RES) + " residues");


  // This constructs all conceivable continuous sections with the same length as the
  // reference in the backbone of the chains
  unsigned nprevious=0; std::vector<unsigned> nlist(refstructure.size());
  for(unsigned i=0; i<chains.size(); ++i) {
    if( chains[i]< refstructure.size() ) error("segment of backbone defined is not long enough to form this reference. "
                                               "Each backbone fragment must contain a minimum of " + std::to_string(refstructure.size()) + " residues");
    unsigned nres = chains[i] / ATOMS_IN_BB_RES;
    if(chains[i] % ATOMS_IN_BB_RES != 0 ) error("backbone segment received does not contain a multiple of " +
                                                std::to_string(ATOMS_IN_BB_RES) + " residues");

    // Loop over all segments of reference length in this chain
    for(unsigned ires=0; ires <= nres - refstructure.size()/ATOMS_IN_BB_RES; ires++) {
      unsigned accum = nprevious + ATOMS_IN_BB_RES * ires;
      for(unsigned k=0; k < refstructure.size(); ++k) nlist[k] = accum + k; //indices of atoms in all_atoms
      addColvar( nlist );
    }
    nprevious+=chains[i];
  }

  // Store the secondary structure ( last number makes sure we convert to internal units nm )
  setSecondaryStructure( refstructure, 0.17/atoms.getUnits().getLength(), 0.1/atoms.getUnits().getLength() );
}

void ContinuousSSRMSD::readBackboneAtoms(const std::string& moltype, std::vector< std::vector<AtomNumber> > &backatoms ) {
  std::vector<AtomNumber> t;
  for(int i=1;; ++i) {
    parseAtomList("ATOMS", i, t);
    if (t.empty()) break;
    backatoms.push_back(t);
    log << "  Backbone " << i << " is calculated from atoms : ";
    for (unsigned j = 0; j < t.size(); ++j) log << t[j].serial() <<" ";
    log << "\n";
    t.resize(0);
  }
}

}
}
