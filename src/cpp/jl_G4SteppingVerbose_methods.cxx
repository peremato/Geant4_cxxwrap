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
void add_methods_for_G4SteppingVerbose(jlcxx::Module& types, jlcxx::TypeWrapper<G4SteppingVerbose>& t148) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4SteppingVerbose
   */

  DEBUG_MSG("Adding wrapper for G4VSteppingVerbose * G4SteppingVerbose::Clone() (" __HERE__ ")");
  // signature to use in the veto list: G4VSteppingVerbose * G4SteppingVerbose::Clone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:51:33
  t148.method("Clone", static_cast<G4VSteppingVerbose * (G4SteppingVerbose::*)() >(&G4SteppingVerbose::Clone));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::NewStep() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::NewStep()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:56:18
  t148.method("NewStep", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::NewStep));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::AtRestDoItInvoked() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::AtRestDoItInvoked()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:57:18
  t148.method("AtRestDoItInvoked", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::AtRestDoItInvoked));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::AlongStepDoItAllDone() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::AlongStepDoItAllDone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:58:18
  t148.method("AlongStepDoItAllDone", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::AlongStepDoItAllDone));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::PostStepDoItAllDone() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::PostStepDoItAllDone()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:59:18
  t148.method("PostStepDoItAllDone", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::PostStepDoItAllDone));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::AlongStepDoItOneByOne() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::AlongStepDoItOneByOne()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:60:18
  t148.method("AlongStepDoItOneByOne", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::AlongStepDoItOneByOne));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::PostStepDoItOneByOne() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::PostStepDoItOneByOne()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:61:18
  t148.method("PostStepDoItOneByOne", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::PostStepDoItOneByOne));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::StepInfo() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::StepInfo()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:62:18
  t148.method("StepInfo", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::StepInfo));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::TrackingStarted() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::TrackingStarted()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:63:18
  t148.method("TrackingStarted", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::TrackingStarted));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::DPSLStarted() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::DPSLStarted()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:64:18
  t148.method("DPSLStarted", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::DPSLStarted));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::DPSLUserLimit() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::DPSLUserLimit()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:65:18
  t148.method("DPSLUserLimit", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::DPSLUserLimit));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::DPSLPostStep() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::DPSLPostStep()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:66:18
  t148.method("DPSLPostStep", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::DPSLPostStep));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::DPSLAlongStep() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::DPSLAlongStep()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:67:18
  t148.method("DPSLAlongStep", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::DPSLAlongStep));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::VerboseTrack() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::VerboseTrack()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:68:18
  t148.method("VerboseTrack", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::VerboseTrack));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::VerboseParticleChange() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::VerboseParticleChange()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:69:18
  t148.method("VerboseParticleChange", static_cast<void (G4SteppingVerbose::*)() >(&G4SteppingVerbose::VerboseParticleChange));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::ShowStep() (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::ShowStep()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:70:18
  t148.method("ShowStep", static_cast<void (G4SteppingVerbose::*)()  const>(&G4SteppingVerbose::ShowStep));

  DEBUG_MSG("Adding wrapper for void G4SteppingVerbose::UseBestUnit(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4SteppingVerbose::UseBestUnit(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:76:17
  t148.method("G4SteppingVerbose!UseBestUnit", static_cast<void (*)(G4int) >(&G4SteppingVerbose::UseBestUnit));
  t148.method("G4SteppingVerbose!UseBestUnit", []()->void{ G4SteppingVerbose::UseBestUnit(); });

  DEBUG_MSG("Adding wrapper for G4int G4SteppingVerbose::BestUnitPrecision() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4SteppingVerbose::BestUnitPrecision()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4SteppingVerbose.hh:77:18
  t148.method("G4SteppingVerbose!BestUnitPrecision", static_cast<G4int (*)() >(&G4SteppingVerbose::BestUnitPrecision));

  /* End of G4SteppingVerbose class method wrappers
   **********************************************************************/

}
