// this file was auto-generated by wrapit v0.1.0-54-g4322429
#include <type_traits>
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

#include "cpp/jlGeant4.h"


#ifdef VERBOSE_IMPORT
#  define DEBUG_MSG(a) std::cerr << a << "\n"
#else
#  define DEBUG_MSG(a)
#endif
#define __HERE__  __FILE__ ":" QUOTE2(__LINE__)
#define QUOTE(arg) #arg
#define QUOTE2(arg) QUOTE(arg)
void add_methods_for_CLHEP_HepLorentzRotation_HepLorentzRotation_row(jlcxx::Module& types, jlcxx::TypeWrapper<CLHEP::HepLorentzRotation::HepLorentzRotation_row>& t176) {


  /**********************************************************************/
  /* Wrappers for the methods of class CLHEP::HepLorentzRotation::HepLorentzRotation_row
   */


  DEBUG_MSG("Adding wrapper for void CLHEP::HepLorentzRotation::HepLorentzRotation_row::HepLorentzRotation_row(const CLHEP::HepLorentzRotation &, int) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/CLHEP/Vector/LorentzRotation.h:178:12
  t176.constructor<const CLHEP::HepLorentzRotation &, int>(/*finalize=*/true);
  types.set_override_module(jl_base_module);


  DEBUG_MSG("Adding getindex method to wrap double CLHEP::HepLorentzRotation::HepLorentzRotation_row::operator[](int) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/CLHEP/Vector/LorentzRotation.h:179:19
  t176.method("getindex",
    [](CLHEP::HepLorentzRotation::HepLorentzRotation_row& a, int i){
    return a[i];
  });

  /* End of CLHEP::HepLorentzRotation::HepLorentzRotation_row class method wrappers
   **********************************************************************/

}
