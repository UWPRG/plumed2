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

#ifndef __PLUMED_peptoidsecondarystructure_ContinuousSSRMSD_h
#define __PLUMED_peptoidsecondarystructure_ContinuousSSRMSD_h

#include "secondarystructure/SecondaryStructureRMSD.h"

namespace PLMD {

namespace peptoidsecondarystructure {

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc

class ContinuousSSRMSD :
  public secondarystructure::SecondaryStructureRMSD
{
public:
  explicit ContinuousSSRMSD(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  /// Reference str
  virtual std::vector<Vector> getRefStructure() const = 0;
protected:
  void readBackboneAtoms(const std::string& moltype, std::vector< std::vector<AtomNumber> > &backatoms ) override;
  /// Call this function in the constructor of the derived class.
  void init();
private:
  static constexpr int ATOMS_IN_BB_RES = 5; // # of backbone atoms per residue
};

}
}

#endif
