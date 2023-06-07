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
void add_methods_for_G4Isotope(jlcxx::Module& types, jlcxx::TypeWrapper<G4Isotope>& t36) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Isotope
   */


  DEBUG_MSG("Adding wrapper for void G4Isotope::G4Isotope(const G4String &, G4int, G4int, G4double, G4int) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:77:5
  t36.constructor<const G4String &, G4int, G4int>(/*finalize=*/false);
  t36.constructor<const G4String &, G4int, G4int, G4double>(/*finalize=*/false);
  t36.constructor<const G4String &, G4int, G4int, G4double, G4int>(/*finalize=*/false);

  DEBUG_MSG("Adding wrapper for const G4String & G4Isotope::GetName() (" __HERE__ ")");
  // signature to use in the veto list: const G4String & G4Isotope::GetName()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:87:21
  t36.method("GetName", static_cast<const G4String & (G4Isotope::*)()  const>(&G4Isotope::GetName));

  DEBUG_MSG("Adding wrapper for G4int G4Isotope::GetZ() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4Isotope::GetZ()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:90:14
  t36.method("GetZ", static_cast<G4int (G4Isotope::*)()  const>(&G4Isotope::GetZ));

  DEBUG_MSG("Adding wrapper for G4int G4Isotope::GetN() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4Isotope::GetN()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:93:14
  t36.method("GetN", static_cast<G4int (G4Isotope::*)()  const>(&G4Isotope::GetN));

  DEBUG_MSG("Adding wrapper for G4double G4Isotope::GetA() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Isotope::GetA()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:96:14
  t36.method("GetA", static_cast<G4double (G4Isotope::*)()  const>(&G4Isotope::GetA));

  DEBUG_MSG("Adding wrapper for G4int G4Isotope::Getm() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4Isotope::Getm()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:99:14
  t36.method("Getm", static_cast<G4int (G4Isotope::*)()  const>(&G4Isotope::Getm));

  DEBUG_MSG("Adding wrapper for G4Isotope * G4Isotope::GetIsotope(const G4String &, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4Isotope * G4Isotope::GetIsotope(const G4String &, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:102:16
  t36.method("G4Isotope!GetIsotope", static_cast<G4Isotope * (*)(const G4String &, G4bool) >(&G4Isotope::GetIsotope));
  t36.method("G4Isotope!GetIsotope", [](const G4String & arg0)->G4Isotope *{ return G4Isotope::GetIsotope(arg0); });

  DEBUG_MSG("Adding wrapper for size_t G4Isotope::GetNumberOfIsotopes() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4Isotope::GetNumberOfIsotopes()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:108:12
  t36.method("G4Isotope!GetNumberOfIsotopes", static_cast<size_t (*)() >(&G4Isotope::GetNumberOfIsotopes));

  DEBUG_MSG("Adding wrapper for size_t G4Isotope::GetIndex() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4Isotope::GetIndex()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:110:12
  t36.method("GetIndex", static_cast<size_t (G4Isotope::*)()  const>(&G4Isotope::GetIndex));
  types.set_override_module(jl_base_module);

  DEBUG_MSG("Adding wrapper for G4bool G4Isotope::operator==(const G4Isotope &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4Isotope::operator==(const G4Isotope &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:123:12
  t36.method("==", static_cast<G4bool (G4Isotope::*)(const G4Isotope &)  const>(&G4Isotope::operator==));

  DEBUG_MSG("Adding wrapper for G4bool G4Isotope::operator!=(const G4Isotope &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4Isotope::operator!=(const G4Isotope &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:124:12
  t36.method("!=", static_cast<G4bool (G4Isotope::*)(const G4Isotope &)  const>(&G4Isotope::operator!=));

  types.unset_override_module();

  DEBUG_MSG("Adding wrapper for void G4Isotope::SetName(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Isotope::SetName(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Isotope.hh:131:10
  t36.method("SetName", static_cast<void (G4Isotope::*)(const G4String &) >(&G4Isotope::SetName));

  /* End of G4Isotope class method wrappers
   **********************************************************************/

}
