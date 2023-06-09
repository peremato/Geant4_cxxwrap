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
void add_methods_for_G4UniformMagField(jlcxx::Module& types, jlcxx::TypeWrapper<G4UniformMagField>& t212) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4UniformMagField
   */


  DEBUG_MSG("Adding wrapper for void G4UniformMagField::G4UniformMagField(const G4ThreeVector &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UniformMagField.hh:45:5
  t212.constructor<const G4ThreeVector &>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4UniformMagField::G4UniformMagField(G4double, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UniformMagField.hh:48:5
  t212.constructor<G4double, G4double, G4double>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4UniformMagField::G4UniformMagField(const G4UniformMagField &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UniformMagField.hh:54:5
  t212.constructor<const G4UniformMagField &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4UniformMagField & G4UniformMagField::operator=(const G4UniformMagField &) (" __HERE__ ")");
  // signature to use in the veto list: G4UniformMagField & G4UniformMagField::operator=(const G4UniformMagField &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UniformMagField.hh:55:24
  t212.method("assign", static_cast<G4UniformMagField & (G4UniformMagField::*)(const G4UniformMagField &) >(&G4UniformMagField::operator=));

  DEBUG_MSG("Adding wrapper for void G4UniformMagField::SetFieldValue(const G4ThreeVector &) (" __HERE__ ")");
  // signature to use in the veto list: void G4UniformMagField::SetFieldValue(const G4ThreeVector &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UniformMagField.hh:61:10
  t212.method("SetFieldValue", static_cast<void (G4UniformMagField::*)(const G4ThreeVector &) >(&G4UniformMagField::SetFieldValue));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4UniformMagField::GetConstantFieldValue() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4UniformMagField::GetConstantFieldValue()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UniformMagField.hh:63:19
  t212.method("GetConstantFieldValue", static_cast<G4ThreeVector (G4UniformMagField::*)()  const>(&G4UniformMagField::GetConstantFieldValue));

  DEBUG_MSG("Adding wrapper for G4Field * G4UniformMagField::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4Field * G4UniformMagField::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UniformMagField.hh:66:22
  t212.method("Clone", static_cast<G4Field * (G4UniformMagField::*)()  const>(&G4UniformMagField::Clone));

  /* End of G4UniformMagField class method wrappers
   **********************************************************************/

}
