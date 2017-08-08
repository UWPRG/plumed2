/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION COMBINE
/*
Calculate a polynomial combination of a set of other variables.

The functional form of this function is
\f[
C=\sum_{i=1}^{N_{arg}} c_i (x_i-a_i)^{p_i}
\f]

The coefficients c, the parameters a and the powers p are provided as vectors.

Notice that COMBINE is not able to predict which will be periodic domain
of the computed value automatically. The user is thus forced to specify it
explicitly. Use PERIODIC=NO if the resulting variable is not periodic,
and PERIODIC=A,B where A and B are the two boundaries if the resulting variable
is periodic.



\par Examples

The following input tells plumed to print the distance between atoms 3 and 5
its square (as computed from the x,y,z components) and the distance
again as computed from the square root of the square.
\plumedfile
DISTANCE LABEL=dist      ATOMS=3,5 COMPONENTS
COMBINE  LABEL=distance2 ARG=dist.x,dist.y,dist.z POWERS=2,2,2 PERIODIC=NO
COMBINE  LABEL=distance  ARG=distance2 POWERS=0.5 PERIODIC=NO
PRINT ARG=distance,distance2
\endplumedfile
(See also \ref PRINT and \ref DISTANCE).

The following input tells plumed to add a restraint on the
cube of a dihedral angle. Notice that since the angle has a
periodic domain
-pi,pi its cube has a domain -pi**3,pi**3.
\plumedfile
t: TORSION ATOMS=1,3,5,7
c: COMBINE ARG=t POWERS=3 PERIODIC=-31.0062766802998,31.0062766802998
RESTRAINT ARG=c KAPPA=10 AT=0
\endplumedfile



*/
//+ENDPLUMEDOC


class Combine :
  public Function
{
  bool normalize;
  std::vector<double> coefficients;
  std::vector<double> parameters;
  std::vector<double> powers;
public:
  explicit Combine(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Combine,"COMBINE")

void Combine::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory","COEFFICIENTS","1.0","the coefficients of the arguments in your function");
  keys.add("compulsory","PARAMETERS","0.0","the parameters of the arguments in your function");
  keys.add("compulsory","POWERS","1.0","the powers to which you are raising each of the arguments in your function");
  keys.addFlag("NORMALIZE",false,"normalize all the coefficents so that in total they are equal to one");
}

Combine::Combine(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  normalize(false),
  coefficients(getNumberOfArguments(),1.0),
  parameters(getNumberOfArguments(),0.0),
  powers(getNumberOfArguments(),1.0)
{
  if( !numberedkeys ) {
     if( getPntrToArgument(0)->getRank()>0 && arg_ends.size()!=2 ) error("should only specify one non-scalar argument in input to ARG keyword");
 
     coefficients.resize( getNumberOfScalarArguments() ); parseVector("COEFFICIENTS",coefficients);
     if(coefficients.size()!=static_cast<unsigned>(getNumberOfScalarArguments()))
       error("Size of COEFFICIENTS array should be the same as number for arguments");
  
     parameters.resize(getNumberOfScalarArguments()); parseVector("PARAMETERS",parameters);
     if(parameters.size()!=static_cast<unsigned>(getNumberOfScalarArguments()))
       error("Size of PARAMETERS array should be the same as number for arguments");
  
     powers.resize(getNumberOfScalarArguments()); parseVector("POWERS",powers);
     if(powers.size()!=static_cast<unsigned>(getNumberOfScalarArguments()))
       error("Size of POWERS array should be the same as number for arguments");
  } else {
     parseVector("COEFFICIENTS",coefficients);
     if(coefficients.size()!=static_cast<unsigned>(getNumberOfArguments()))
       error("Size of COEFFICIENTS array should be the same as number for arguments");

     parseVector("PARAMETERS",parameters);
     if(parameters.size()!=static_cast<unsigned>(getNumberOfArguments()))
       error("Size of PARAMETERS array should be the same as number for arguments");

     parseVector("POWERS",powers);
     if(powers.size()!=static_cast<unsigned>(getNumberOfArguments()))
       error("Size of POWERS array should be the same as number for arguments");
  }
  parseFlag("NORMALIZE",normalize);

  if(normalize) {
    double n=0.0;
    for(unsigned i=0; i<coefficients.size(); i++) n+=coefficients[i];
    for(unsigned i=0; i<coefficients.size(); i++) coefficients[i]*=(1.0/n);
  }

  addValueWithDerivatives();
  checkRead();

  bool allsame=true; double coeff=coefficients[0];
  for(unsigned i=1;i<coefficients.size(); i++){
      if( coefficients[i]!=coeff ){ allsame=false; break; }
  }
  if( allsame ){
      log.printf("  with all coefficients equal to %f\n",coeff);
  } else {
      log.printf("  with coefficients:");
      for(unsigned i=0; i<coefficients.size(); i++) log.printf(" %f",coefficients[i]);
      log.printf("\n");
  }
  allsame=true; double param=parameters[0];
  for(unsigned i=1;i<parameters.size(); i++){
      if( parameters[i]!=param ){ allsame=false; break; }
  }
  if( allsame ){
      log.printf("  with all parameters equal to %f\n",param);
  } else {
      log.printf("  with parameters:");
      for(unsigned i=0; i<parameters.size(); i++) log.printf(" %f",parameters[i]);
      log.printf("\n");
  }
  allsame=true; double power=powers[0];
  for(unsigned i=1;i<powers.size(); i++){
      if( powers[i]!=power ){ allsame=false; break; }
  }
  if( allsame ){
      log.printf("  with all powers equal to %f\n",power);
  } else {
      log.printf("  and powers:");
      for(unsigned i=0; i<powers.size(); i++) log.printf(" %f",powers[i]);
      log.printf("\n");
  }
}

void Combine::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  double combine=0.0;
  if( args.size()==1 && !numberedkeys ) {
      unsigned ind = myvals.getTaskIndex(); plumed_dbg_assert( ind<parameters.size() );
      double cv = getPntrToArgument(0)->difference( parameters[ind], args[0] );
      combine = coefficients[ind]*pow(cv,powers[ind]); 
      addDerivative( 0, 0, coefficients[ind]*powers[ind]*pow(cv,powers[ind]-1.0), myvals ); 
  } else {
      for(unsigned i=0; i<coefficients.size(); ++i) {
        double cv = difference(i,parameters[i],args[i]);   
        combine+=coefficients[i]*pow(cv,powers[i]);
        addDerivative(0, i, coefficients[i]*powers[i]*pow(cv,powers[i]-1.0), myvals );
      };
  }
  addValue( 0, combine, myvals );
}

}
}


