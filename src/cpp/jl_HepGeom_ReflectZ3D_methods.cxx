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
void add_methods_for_HepGeom_ReflectZ3D(jlcxx::Module& types, jlcxx::TypeWrapper<HepGeom::ReflectZ3D>& t78) {


  /**********************************************************************/
  /* Wrappers for the methods of class HepGeom::ReflectZ3D
   */


  DEBUG_MSG("Adding wrapper for void HepGeom::ReflectZ3D::ReflectZ3D(double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/CLHEP/Geometry/Transform3D.h:702:5
  t78.constructor<double>(/*finalize=*/true);

  /* End of HepGeom::ReflectZ3D class method wrappers
   **********************************************************************/

}
