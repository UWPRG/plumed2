/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "SecondaryStructureRMSD.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace peptoidsecondarystructure {

//+PLUMEDOC COLVAR ALPHARMSD
/*
Probe the alpha helical content of a protein structure.

Any chain of six contiguous residues in a protein chain can form an alpha helix. This
colvar thus generates the set of all possible six residue sections and calculates
the RMSD distance between the configuration in which the residues find themselves
and an idealized alpha helical structure. These distances can be calculated by either
aligning the instantaneous structure with the reference structure and measuring each
atomic displacement or by calculating differences between the set of inter-atomic
distances in the reference and instantaneous structures.

This colvar is based on the following reference \cite pietrucci09jctc.  The authors of
this paper use the set of distances from the alpha helix configurations to measure
the number of segments that have an alpha helical configuration. This is done by calculating
the following sum of functions of the rmsd distances:

\f[
s = \sum_i \frac{ 1 - \left(\frac{r_i-d_0}{r_0}\right)^n } { 1 - \left(\frac{r_i-d_0}{r_0}\right)^m }
\f]

where the sum runs over all possible segments of alpha helix.  By default the
NN, MM and D_0 parameters are set equal to those used in \cite pietrucci09jctc.  The R_0
parameter must be set by the user - the value used in \cite pietrucci09jctc was 0.08 nm.

If you change the function in the above sum you can calculate quantities such as the average
distance from a purely the alpha helical configuration or the distance between the set of
residues that is closest to an alpha helix and the reference configuration. To do these sorts of
calculations you can use the AVERAGE and MIN keywords. In addition you can use the LESS_THAN
keyword if you would like to change the form of the switching function. If you use any of these
options you no longer need to specify NN, R_0, MM and D_0.

Please be aware that for codes like gromacs you must ensure that plumed
reconstructs the chains involved in your CV when you calculate this CV using
anything other than TYPE=DRMSD.  For more details as to how to do this see \ref WHOLEMOLECULES.

\par Examples

The following input calculates the number of six residue segments of
protein that are in an alpha helical configuration.

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
alpha: ALPHARMSD RESIDUES=all
\endplumedfile

Here the same is done use RMSD instead of DRMSD

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
WHOLEMOLECULES ENTITY0=1-100
alpha: ALPHARMSD RESIDUES=all TYPE=OPTIMAL R_0=0.1
\endplumedfile

*/
//+ENDPLUMEDOC

class AlphaDPlusCisRMSD : public SecondaryStructureRMSD {
public:
  static void registerKeywords( Keywords& keys );
  explicit AlphaDPlusCisRMSD(const ActionOptions&);
  static constexpr int RES_IN_REF = 3; // # of residues we're comparing against
};

PLUMED_REGISTER_ACTION(AlphaDPlusCisRMSD, "ALPHADPLUSCISRMSD")

void AlphaDPlusCisRMSD::registerKeywords(Keywords& keys ) {
  SecondaryStructureRMSD::registerKeywords( keys );
}

AlphaDPlusCisRMSD::AlphaDPlusCisRMSD(const ActionOptions&ao):
  Action(ao),
  SecondaryStructureRMSD(ao)
{
  // read in the backbone atoms
  std::vector<unsigned> chains; readBackboneAtoms( "protein", chains);

  // This constructs all conceivable sections of alpha helix in the backbone of the chains
  int natoms = RES_IN_REF * ATOMS_IN_BB_RES;
  unsigned nprevious=0; std::vector<unsigned> nlist(natoms);
  for(unsigned i=0; i<chains.size(); ++i) {
    if( chains[i]< natoms ) error("segment of backbone defined is not long enough to form an alpha helix. "
                                  "Each backbone fragment must contain a minimum of " + std::to_string(RES_IN_REF) + " residues");
    unsigned nres= chains[i] / ATOMS_IN_BB_RES;
    if(chains[i] % ATOMS_IN_BB_RES != 0 ) error("backbone segment received does not contain a multiple of " +
                                                std::to_string(ATOMS_IN_BB_RES) + " residues");
    for(unsigned ires=0; ires< nres - ATOMS_IN_BB_RES; ires++) {
      unsigned accum= nprevious + ATOMS_IN_BB_RES * ires;
      for(unsigned k=0; k< natoms; ++k) nlist[k] = accum + k; //indices of atoms in all_atoms
      addColvar( nlist );
    }
    nprevious+=chains[i];
  }

  // Build the reference structure ( in angstroms )
  std::vector<Vector> reference(natoms);
  reference[0] = Vector(-0.046,  1.837, -0.547); // CLP    i
  reference[1] = Vector(-0.408,  1.267, -1.572); // OL
  reference[2] = Vector(-0.806,  2.807,  0.058); // NL
  reference[3] = Vector(-0.297,  3.697,  1.034); // CA
  reference[4] = Vector(-2.084,  3.170, -0.401); // CB1
  reference[5] = Vector( 1.245, -0.942,  0.140); // CLP    i+1
  reference[6] = Vector( 1.830, -1.992,  0.061); // OL
  reference[7] = Vector( 1.860,  0.187, -0.332); // NL
  reference[8] = Vector( 1.310,  1.481,  0.025); // CA
  reference[9] = Vector( 3.241,  0.171, -0.772); // CB1
  reference[10] = Vector(-1.554, -2.433, -0.485 ); // CLP    i+2
  reference[11] = Vector(-2.229, -3.435, -0.491 ); // OL
  reference[12] = Vector(-0.879, -2.087,  0.638 ); // NL
  reference[13] = Vector(-0.164, -0.855,  0.765 ); // CA
  reference[14] = Vector(-1.017, -2.874,  1.886 ); // CB1
  // Store the secondary structure ( last number makes sure we convert to internal units nm )
  setSecondaryStructure( reference, 0.17/atoms.getUnits().getLength(), 0.1/atoms.getUnits().getLength() );
}

}
}
