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
void add_methods_for_G4ExtrudedSolid_ZSection(jlcxx::Module& types, jlcxx::TypeWrapper<G4ExtrudedSolid::ZSection>& t187) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4ExtrudedSolid::ZSection
   */


  DEBUG_MSG("Adding wrapper for void G4ExtrudedSolid::ZSection::ZSection(G4double, const G4TwoVector &, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ExtrudedSolid.hh:79:7
  t187.constructor<G4double, const G4TwoVector &, G4double>(/*finalize=*/true);

  DEBUG_MSG("Adding fZ methods  to provide read access to the field fZ (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ExtrudedSolid.hh:82:19
  // signature to use in the veto list: G4ExtrudedSolid::ZSection::fZ
  t187.method("fZ", [](const G4ExtrudedSolid::ZSection& a) -> G4double { return a.fZ; });
  t187.method("fZ", [](G4ExtrudedSolid::ZSection& a) -> G4double { return a.fZ; });
  t187.method("fZ", [](const G4ExtrudedSolid::ZSection* a) -> G4double { return a->fZ; });
  t187.method("fZ", [](G4ExtrudedSolid::ZSection* a) -> G4double { return a->fZ; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ExtrudedSolid.hh:82:19
  // signature to use in the veto list: G4ExtrudedSolid::ZSection::fZ
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fZ! methods to provide write access to the field fZ (" __HERE__ ")");
  t187.method("fZ!", [](G4ExtrudedSolid::ZSection& a, G4double val) -> G4double { return a.fZ = val; });

  DEBUG_MSG("Adding fZ! methods to provide write access to the field fZ (" __HERE__ ")");
  t187.method("fZ!", [](G4ExtrudedSolid::ZSection* a, G4double val) -> G4double { return a->fZ = val; });

  DEBUG_MSG("Adding fOffset methods  to provide read access to the field fOffset (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ExtrudedSolid.hh:83:19
  // signature to use in the veto list: G4ExtrudedSolid::ZSection::fOffset
  t187.method("fOffset", [](const G4ExtrudedSolid::ZSection& a) -> const G4TwoVector& { return a.fOffset; });
  t187.method("fOffset", [](G4ExtrudedSolid::ZSection& a) -> G4TwoVector& { return a.fOffset; });
  t187.method("fOffset", [](const G4ExtrudedSolid::ZSection* a) -> const G4TwoVector& { return a->fOffset; });
  t187.method("fOffset", [](G4ExtrudedSolid::ZSection* a) -> G4TwoVector& { return a->fOffset; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ExtrudedSolid.hh:83:19
  // signature to use in the veto list: G4ExtrudedSolid::ZSection::fOffset
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fOffset! methods to provide write access to the field fOffset (" __HERE__ ")");
  t187.method("fOffset!", [](G4ExtrudedSolid::ZSection& a, const G4TwoVector& val) -> G4TwoVector& { return a.fOffset = val; });

  DEBUG_MSG("Adding fOffset! methods to provide write access to the field fOffset (" __HERE__ ")");
  t187.method("fOffset!", [](G4ExtrudedSolid::ZSection* a, const G4TwoVector& val) -> G4TwoVector& { return a->fOffset = val; });

  DEBUG_MSG("Adding fScale methods  to provide read access to the field fScale (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ExtrudedSolid.hh:84:19
  // signature to use in the veto list: G4ExtrudedSolid::ZSection::fScale
  t187.method("fScale", [](const G4ExtrudedSolid::ZSection& a) -> G4double { return a.fScale; });
  t187.method("fScale", [](G4ExtrudedSolid::ZSection& a) -> G4double { return a.fScale; });
  t187.method("fScale", [](const G4ExtrudedSolid::ZSection* a) -> G4double { return a->fScale; });
  t187.method("fScale", [](G4ExtrudedSolid::ZSection* a) -> G4double { return a->fScale; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4ExtrudedSolid.hh:84:19
  // signature to use in the veto list: G4ExtrudedSolid::ZSection::fScale
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fScale! methods to provide write access to the field fScale (" __HERE__ ")");
  t187.method("fScale!", [](G4ExtrudedSolid::ZSection& a, G4double val) -> G4double { return a.fScale = val; });

  DEBUG_MSG("Adding fScale! methods to provide write access to the field fScale (" __HERE__ ")");
  t187.method("fScale!", [](G4ExtrudedSolid::ZSection* a, G4double val) -> G4double { return a->fScale = val; });

  /* End of G4ExtrudedSolid::ZSection class method wrappers
   **********************************************************************/

}
