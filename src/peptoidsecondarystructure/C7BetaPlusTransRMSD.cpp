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
#include "ContinuousSSRMSD.h"
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

class C7BetaPlusTransRMSD : public ContinuousSSRMSD {
public:
  static void registerKeywords( Keywords& keys );
  explicit C7BetaPlusTransRMSD(const ActionOptions&);
  std::vector<Vector> getRefStructure() const override;
};

PLUMED_REGISTER_ACTION(C7BetaPlusTransRMSD, "C7BETAPLUSTRANSRMSD")

void C7BetaPlusTransRMSD::registerKeywords(Keywords& keys ) {
  ContinuousSSRMSD::registerKeywords( keys );
}

std::vector<Vector> C7BetaPlusTransRMSD::getRefStructure() const {
  return  { Vector(-2.424, -0.383, 1.022), // CLP    i
            Vector(-2.834, -1.493, 1.432), // OL
            Vector(-3.294, 0.447, 0.242), // NL
            Vector(-4.704, -0.033, 0.142), // CA
            Vector(-2.814, 1.597, -0.538), // CB1
            Vector(0.856, 0.167, -0.268), // CLP    i+1
            Vector(0.646, 1.387, -0.138), // OL
            Vector(0.056, -0.653, 0.532), // NL
            Vector(-0.984, -0.103, 1.332), // CA
            Vector(0.026, -2.113, 0.422), // CB1
            Vector(3.906, 0.117, -0.328), // CLP    i+2
            Vector(3.516, -0.533, 0.642), // OL
            Vector(3.026, 0.397, -1.248), // NL
            Vector(1.856, -0.413, -1.208), // CA
            Vector(3.166, 1.617, -2.038) };
}

C7BetaPlusTransRMSD::C7BetaPlusTransRMSD(const ActionOptions&ao):
  Action(ao),
  ContinuousSSRMSD(ao)
{
  init();
}

}
}
