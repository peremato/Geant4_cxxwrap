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
void add_methods_for_HepGeom_Reflect3D(jlcxx::Module& types, jlcxx::TypeWrapper<HepGeom::Reflect3D>& t77) {


  /**********************************************************************/
  /* Wrappers for the methods of class HepGeom::Reflect3D
   */


  DEBUG_MSG("Adding wrapper for void HepGeom::Reflect3D::Reflect3D(double, double, double, double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Geometry/Transform3D.h:634:5
  t77.constructor<double, double, double, double>(/*finalize=*/true);

  /* End of HepGeom::Reflect3D class method wrappers
   **********************************************************************/

}
