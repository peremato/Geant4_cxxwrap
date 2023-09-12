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
void add_methods_for_G4PVReplica(jlcxx::Module& types, jlcxx::TypeWrapper<G4PVReplica>& t166) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4PVReplica
   */


  DEBUG_MSG("Adding wrapper for void G4PVReplica::G4PVReplica(const G4String &, G4LogicalVolume *, G4LogicalVolume *, const EAxis, const G4int, const G4double, const G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:94:5
  t166.constructor<const G4String &, G4LogicalVolume *, G4LogicalVolume *, const EAxis, const G4int, const G4double>(/*finalize=*/false);
  t166.constructor<const G4String &, G4LogicalVolume *, G4LogicalVolume *, const EAxis, const G4int, const G4double, const G4double>(/*finalize=*/false);


  DEBUG_MSG("Adding wrapper for void G4PVReplica::G4PVReplica(const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, const EAxis, const G4int, const G4double, const G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:102:5
  t166.constructor<const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, const EAxis, const G4int, const G4double>(/*finalize=*/false);
  t166.constructor<const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, const EAxis, const G4int, const G4double, const G4double>(/*finalize=*/false);

  DEBUG_MSG("Adding wrapper for EVolume G4PVReplica::VolumeType() (" __HERE__ ")");
  // signature to use in the veto list: EVolume G4PVReplica::VolumeType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:121:21
  t166.method("VolumeType", static_cast<EVolume (G4PVReplica::*)()  const>(&G4PVReplica::VolumeType));

  DEBUG_MSG("Adding wrapper for G4bool G4PVReplica::IsMany() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVReplica::IsMany()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:123:12
  t166.method("IsMany", static_cast<G4bool (G4PVReplica::*)()  const>(&G4PVReplica::IsMany));

  DEBUG_MSG("Adding wrapper for G4bool G4PVReplica::IsReplicated() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVReplica::IsReplicated()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:124:12
  t166.method("IsReplicated", static_cast<G4bool (G4PVReplica::*)()  const>(&G4PVReplica::IsReplicated));

  DEBUG_MSG("Adding wrapper for G4int G4PVReplica::GetCopyNo() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4PVReplica::GetCopyNo()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:126:19
  t166.method("GetCopyNo", static_cast<G4int (G4PVReplica::*)()  const>(&G4PVReplica::GetCopyNo));

  DEBUG_MSG("Adding wrapper for void G4PVReplica::SetCopyNo(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4PVReplica::SetCopyNo(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:127:19
  t166.method("SetCopyNo", static_cast<void (G4PVReplica::*)(G4int) >(&G4PVReplica::SetCopyNo));

  DEBUG_MSG("Adding wrapper for G4bool G4PVReplica::IsParameterised() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVReplica::IsParameterised()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:128:20
  t166.method("IsParameterised", static_cast<G4bool (G4PVReplica::*)()  const>(&G4PVReplica::IsParameterised));

  DEBUG_MSG("Adding wrapper for G4VPVParameterisation * G4PVReplica::GetParameterisation() (" __HERE__ ")");
  // signature to use in the veto list: G4VPVParameterisation * G4PVReplica::GetParameterisation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:129:36
  t166.method("GetParameterisation", static_cast<G4VPVParameterisation * (G4PVReplica::*)()  const>(&G4PVReplica::GetParameterisation));

  DEBUG_MSG("Adding wrapper for G4int G4PVReplica::GetMultiplicity() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4PVReplica::GetMultiplicity()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:130:19
  t166.method("GetMultiplicity", static_cast<G4int (G4PVReplica::*)()  const>(&G4PVReplica::GetMultiplicity));

  DEBUG_MSG("Adding wrapper for void G4PVReplica::GetReplicationData(EAxis &, G4int &, G4double &, G4double &, G4bool &) (" __HERE__ ")");
  // signature to use in the veto list: void G4PVReplica::GetReplicationData(EAxis &, G4int &, G4double &, G4double &, G4bool &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:131:18
  t166.method("GetReplicationData", static_cast<void (G4PVReplica::*)(EAxis &, G4int &, G4double &, G4double &, G4bool &)  const>(&G4PVReplica::GetReplicationData));

  DEBUG_MSG("Adding wrapper for void G4PVReplica::SetRegularStructureId(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4PVReplica::SetRegularStructureId(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:137:18
  t166.method("SetRegularStructureId", static_cast<void (G4PVReplica::*)(G4int) >(&G4PVReplica::SetRegularStructureId));

  DEBUG_MSG("Adding wrapper for G4bool G4PVReplica::IsRegularStructure() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4PVReplica::IsRegularStructure()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:143:12
  t166.method("IsRegularStructure", static_cast<G4bool (G4PVReplica::*)()  const>(&G4PVReplica::IsRegularStructure));

  DEBUG_MSG("Adding wrapper for G4int G4PVReplica::GetRegularStructureId() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4PVReplica::GetRegularStructureId()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:144:11
  t166.method("GetRegularStructureId", static_cast<G4int (G4PVReplica::*)()  const>(&G4PVReplica::GetRegularStructureId));

  DEBUG_MSG("Adding wrapper for G4int G4PVReplica::GetInstanceID() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4PVReplica::GetInstanceID()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:149:18
  t166.method("GetInstanceID", static_cast<G4int (G4PVReplica::*)()  const>(&G4PVReplica::GetInstanceID));

  DEBUG_MSG("Adding wrapper for void G4PVReplica::InitialiseWorker(G4PVReplica *) (" __HERE__ ")");
  // signature to use in the veto list: void G4PVReplica::InitialiseWorker(G4PVReplica *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:155:10
  t166.method("InitialiseWorker", static_cast<void (G4PVReplica::*)(G4PVReplica *) >(&G4PVReplica::InitialiseWorker));

  DEBUG_MSG("Adding wrapper for void G4PVReplica::TerminateWorker(G4PVReplica *) (" __HERE__ ")");
  // signature to use in the veto list: void G4PVReplica::TerminateWorker(G4PVReplica *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:159:10
  t166.method("TerminateWorker", static_cast<void (G4PVReplica::*)(G4PVReplica *) >(&G4PVReplica::TerminateWorker));

  /* End of G4PVReplica class method wrappers
   **********************************************************************/

}
