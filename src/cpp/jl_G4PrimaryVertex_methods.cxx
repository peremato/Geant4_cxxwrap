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
void add_methods_for_G4PrimaryVertex(jlcxx::Module& types, jlcxx::TypeWrapper<G4PrimaryVertex>& t96) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4PrimaryVertex
   */


  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::G4PrimaryVertex(G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:54:5
  t96.constructor<G4double, G4double, G4double, G4double>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::G4PrimaryVertex(G4ThreeVector, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:55:5
  t96.constructor<G4ThreeVector, G4double>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::G4PrimaryVertex(const G4PrimaryVertex &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:61:5
  t96.constructor<const G4PrimaryVertex &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4PrimaryVertex & G4PrimaryVertex::operator=(const G4PrimaryVertex &) (" __HERE__ ")");
  // signature to use in the veto list: G4PrimaryVertex & G4PrimaryVertex::operator=(const G4PrimaryVertex &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:62:22
  t96.method("assign", static_cast<G4PrimaryVertex & (G4PrimaryVertex::*)(const G4PrimaryVertex &) >(&G4PrimaryVertex::operator=));
  types.set_override_module(jl_base_module);

  DEBUG_MSG("Adding wrapper for G4bool G4PrimaryVertex::operator==(const G4PrimaryVertex &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PrimaryVertex::operator==(const G4PrimaryVertex &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:65:12
  t96.method("==", static_cast<G4bool (G4PrimaryVertex::*)(const G4PrimaryVertex &)  const>(&G4PrimaryVertex::operator==));

  DEBUG_MSG("Adding wrapper for G4bool G4PrimaryVertex::operator!=(const G4PrimaryVertex &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PrimaryVertex::operator!=(const G4PrimaryVertex &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:66:12
  t96.method("!=", static_cast<G4bool (G4PrimaryVertex::*)(const G4PrimaryVertex &)  const>(&G4PrimaryVertex::operator!=));

  types.unset_override_module();

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4PrimaryVertex::GetPosition() (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4PrimaryVertex::GetPosition()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:75:26
  t96.method("GetPosition", static_cast<G4ThreeVector (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::GetPosition));

  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetPosition(G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4PrimaryVertex::SetPosition(G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:76:17
  t96.method("SetPosition", static_cast<void (G4PrimaryVertex::*)(G4double, G4double, G4double) >(&G4PrimaryVertex::SetPosition));

  DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetX0() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4PrimaryVertex::GetX0()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:77:21
  t96.method("GetX0", static_cast<G4double (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::GetX0));

  DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetY0() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4PrimaryVertex::GetY0()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:78:21
  t96.method("GetY0", static_cast<G4double (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::GetY0));

  DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetZ0() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4PrimaryVertex::GetZ0()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:79:21
  t96.method("GetZ0", static_cast<G4double (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::GetZ0));

  DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetT0() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4PrimaryVertex::GetT0()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:80:21
  t96.method("GetT0", static_cast<G4double (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::GetT0));

  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetT0(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4PrimaryVertex::SetT0(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:81:17
  t96.method("SetT0", static_cast<void (G4PrimaryVertex::*)(G4double) >(&G4PrimaryVertex::SetT0));

  DEBUG_MSG("Adding wrapper for G4int G4PrimaryVertex::GetNumberOfParticle() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4PrimaryVertex::GetNumberOfParticle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:82:18
  t96.method("GetNumberOfParticle", static_cast<G4int (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::GetNumberOfParticle));

  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetPrimary(G4PrimaryParticle *) (" __HERE__ ")");
  // signature to use in the veto list: void G4PrimaryVertex::SetPrimary(G4PrimaryParticle *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:83:17
  t96.method("SetPrimary", static_cast<void (G4PrimaryVertex::*)(G4PrimaryParticle *) >(&G4PrimaryVertex::SetPrimary));

  DEBUG_MSG("Adding wrapper for G4PrimaryParticle * G4PrimaryVertex::GetPrimary(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4PrimaryParticle * G4PrimaryVertex::GetPrimary(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:84:24
  t96.method("GetPrimary", static_cast<G4PrimaryParticle * (G4PrimaryVertex::*)(G4int)  const>(&G4PrimaryVertex::GetPrimary));
  t96.method("GetPrimary", [](G4PrimaryVertex const& a)->G4PrimaryParticle *{ return a.GetPrimary(); });
  t96.method("GetPrimary", [](G4PrimaryVertex const* a)->G4PrimaryParticle *{ return a->GetPrimary(); });

  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetNext(G4PrimaryVertex *) (" __HERE__ ")");
  // signature to use in the veto list: void G4PrimaryVertex::SetNext(G4PrimaryVertex *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:85:17
  t96.method("SetNext", static_cast<void (G4PrimaryVertex::*)(G4PrimaryVertex *) >(&G4PrimaryVertex::SetNext));

  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::ClearNext() (" __HERE__ ")");
  // signature to use in the veto list: void G4PrimaryVertex::ClearNext()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:86:17
  t96.method("ClearNext", static_cast<void (G4PrimaryVertex::*)() >(&G4PrimaryVertex::ClearNext));

  DEBUG_MSG("Adding wrapper for G4PrimaryVertex * G4PrimaryVertex::GetNext() (" __HERE__ ")");
  // signature to use in the veto list: G4PrimaryVertex * G4PrimaryVertex::GetNext()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:87:29
  t96.method("GetNext", static_cast<G4PrimaryVertex * (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::GetNext));

  DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetWeight() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4PrimaryVertex::GetWeight()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:88:21
  t96.method("GetWeight", static_cast<G4double (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::GetWeight));

  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetWeight(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4PrimaryVertex::SetWeight(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:89:17
  t96.method("SetWeight", static_cast<void (G4PrimaryVertex::*)(G4double) >(&G4PrimaryVertex::SetWeight));

  DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::Print() (" __HERE__ ")");
  // signature to use in the veto list: void G4PrimaryVertex::Print()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PrimaryVertex.hh:93:10
  t96.method("Print", static_cast<void (G4PrimaryVertex::*)()  const>(&G4PrimaryVertex::Print));

  /* End of G4PrimaryVertex class method wrappers
   **********************************************************************/

}
