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
void add_methods_for_HepGeom_Rotate3D(jlcxx::Module& types, jlcxx::TypeWrapper<HepGeom::Rotate3D>& t70) {


  /**********************************************************************/
  /* Wrappers for the methods of class HepGeom::Rotate3D
   */


  DEBUG_MSG("Adding wrapper for void HepGeom::Rotate3D::Rotate3D(const CLHEP::HepRotation &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Geometry/Transform3D.h:378:12
  t70.constructor<const CLHEP::HepRotation &>(/*finalize=*/true);

  /* End of HepGeom::Rotate3D class method wrappers
   **********************************************************************/

}
