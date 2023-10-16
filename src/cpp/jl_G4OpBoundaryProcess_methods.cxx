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
void add_methods_for_G4OpBoundaryProcess(jlcxx::Module& types, jlcxx::TypeWrapper<G4OpBoundaryProcess>& t232) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4OpBoundaryProcess
   */


  DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::G4OpBoundaryProcess(const G4String &, G4ProcessType) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:122:12
  t232.constructor<const G4String &>(/*finalize=*/true);
  t232.constructor<const G4String &, G4ProcessType>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4bool G4OpBoundaryProcess::IsApplicable(const G4ParticleDefinition &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4OpBoundaryProcess::IsApplicable(const G4ParticleDefinition &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:126:18
  t232.method("IsApplicable", static_cast<G4bool (G4OpBoundaryProcess::*)(const G4ParticleDefinition &) >(&G4OpBoundaryProcess::IsApplicable));

  DEBUG_MSG("Adding wrapper for G4double G4OpBoundaryProcess::GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4OpBoundaryProcess::GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:130:20
  t232.method("GetMeanFreePath", static_cast<G4double (G4OpBoundaryProcess::*)(const G4Track &, G4double, G4ForceCondition *) >(&G4OpBoundaryProcess::GetMeanFreePath));

  DEBUG_MSG("Adding wrapper for G4VParticleChange * G4OpBoundaryProcess::PostStepDoIt(const G4Track &, const G4Step &) (" __HERE__ ")");
  // signature to use in the veto list: G4VParticleChange * G4OpBoundaryProcess::PostStepDoIt(const G4Track &, const G4Step &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:136:22
  t232.method("PostStepDoIt", static_cast<G4VParticleChange * (G4OpBoundaryProcess::*)(const G4Track &, const G4Step &) >(&G4OpBoundaryProcess::PostStepDoIt));

  DEBUG_MSG("Adding wrapper for G4OpBoundaryProcessStatus G4OpBoundaryProcess::GetStatus() (" __HERE__ ")");
  // signature to use in the veto list: G4OpBoundaryProcessStatus G4OpBoundaryProcess::GetStatus()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:140:37
  t232.method("GetStatus", static_cast<G4OpBoundaryProcessStatus (G4OpBoundaryProcess::*)()  const>(&G4OpBoundaryProcess::GetStatus));

  DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::SetInvokeSD(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4OpBoundaryProcess::SetInvokeSD(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:143:16
  t232.method("SetInvokeSD", static_cast<void (G4OpBoundaryProcess::*)(G4bool) >(&G4OpBoundaryProcess::SetInvokeSD));

  DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::PreparePhysicsTable(const G4ParticleDefinition &) (" __HERE__ ")");
  // signature to use in the veto list: void G4OpBoundaryProcess::PreparePhysicsTable(const G4ParticleDefinition &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:146:16
  t232.method("PreparePhysicsTable", static_cast<void (G4OpBoundaryProcess::*)(const G4ParticleDefinition &) >(&G4OpBoundaryProcess::PreparePhysicsTable));

  DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::Initialise() (" __HERE__ ")");
  // signature to use in the veto list: void G4OpBoundaryProcess::Initialise()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:148:16
  t232.method("Initialise", static_cast<void (G4OpBoundaryProcess::*)() >(&G4OpBoundaryProcess::Initialise));

  DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::SetVerboseLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4OpBoundaryProcess::SetVerboseLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4OpBoundaryProcess.hh:150:8
  t232.method("SetVerboseLevel", static_cast<void (G4OpBoundaryProcess::*)(G4int) >(&G4OpBoundaryProcess::SetVerboseLevel));

  /* End of G4OpBoundaryProcess class method wrappers
   **********************************************************************/

}
