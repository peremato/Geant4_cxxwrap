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
void add_methods_for_G4NavigationHistory(jlcxx::Module& types, jlcxx::TypeWrapper<G4NavigationHistory>& t35) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4NavigationHistory
   */


  DEBUG_MSG("Adding wrapper for void G4NavigationHistory::G4NavigationHistory(const G4NavigationHistory &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:66:3
  t35.constructor<const G4NavigationHistory &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4NavigationHistory & G4NavigationHistory::operator=(const G4NavigationHistory &) (" __HERE__ ")");
  // signature to use in the veto list: G4NavigationHistory & G4NavigationHistory::operator=(const G4NavigationHistory &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:69:31
  t35.method("assign", static_cast<G4NavigationHistory & (G4NavigationHistory::*)(const G4NavigationHistory &) >(&G4NavigationHistory::operator=));

  DEBUG_MSG("Adding wrapper for void G4NavigationHistory::Reset() (" __HERE__ ")");
  // signature to use in the veto list: void G4NavigationHistory::Reset()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:72:15
  t35.method("Reset", static_cast<void (G4NavigationHistory::*)() >(&G4NavigationHistory::Reset));

  DEBUG_MSG("Adding wrapper for void G4NavigationHistory::Clear() (" __HERE__ ")");
  // signature to use in the veto list: void G4NavigationHistory::Clear()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:76:15
  t35.method("Clear", static_cast<void (G4NavigationHistory::*)() >(&G4NavigationHistory::Clear));

  DEBUG_MSG("Adding wrapper for void G4NavigationHistory::SetFirstEntry(G4VPhysicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: void G4NavigationHistory::SetFirstEntry(G4VPhysicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:80:15
  t35.method("SetFirstEntry", static_cast<void (G4NavigationHistory::*)(G4VPhysicalVolume *) >(&G4NavigationHistory::SetFirstEntry));

  DEBUG_MSG("Adding wrapper for const G4AffineTransform & G4NavigationHistory::GetTopTransform() (" __HERE__ ")");
  // signature to use in the veto list: const G4AffineTransform & G4NavigationHistory::GetTopTransform()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:84:35
  t35.method("GetTopTransform", static_cast<const G4AffineTransform & (G4NavigationHistory::*)()  const>(&G4NavigationHistory::GetTopTransform));

  DEBUG_MSG("Adding wrapper for const G4AffineTransform * G4NavigationHistory::GetPtrTopTransform() (" __HERE__ ")");
  // signature to use in the veto list: const G4AffineTransform * G4NavigationHistory::GetPtrTopTransform()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:87:35
  t35.method("GetPtrTopTransform", static_cast<const G4AffineTransform * (G4NavigationHistory::*)()  const>(&G4NavigationHistory::GetPtrTopTransform));

  DEBUG_MSG("Adding wrapper for G4int G4NavigationHistory::GetTopReplicaNo() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4NavigationHistory::GetTopReplicaNo()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:90:16
  t35.method("GetTopReplicaNo", static_cast<G4int (G4NavigationHistory::*)()  const>(&G4NavigationHistory::GetTopReplicaNo));

  DEBUG_MSG("Adding wrapper for EVolume G4NavigationHistory::GetTopVolumeType() (" __HERE__ ")");
  // signature to use in the veto list: EVolume G4NavigationHistory::GetTopVolumeType()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:93:18
  t35.method("GetTopVolumeType", static_cast<EVolume (G4NavigationHistory::*)()  const>(&G4NavigationHistory::GetTopVolumeType));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4NavigationHistory::GetTopVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4NavigationHistory::GetTopVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:96:29
  t35.method("GetTopVolume", static_cast<G4VPhysicalVolume * (G4NavigationHistory::*)()  const>(&G4NavigationHistory::GetTopVolume));

  DEBUG_MSG("Adding wrapper for size_t G4NavigationHistory::GetDepth() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4NavigationHistory::GetDepth()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:99:22
  t35.method("GetDepth", static_cast<size_t (G4NavigationHistory::*)()  const>(&G4NavigationHistory::GetDepth));

  DEBUG_MSG("Adding wrapper for size_t G4NavigationHistory::GetMaxDepth() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4NavigationHistory::GetMaxDepth()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:102:22
  t35.method("GetMaxDepth", static_cast<size_t (G4NavigationHistory::*)()  const>(&G4NavigationHistory::GetMaxDepth));

  DEBUG_MSG("Adding wrapper for const G4AffineTransform & G4NavigationHistory::GetTransform(G4int) (" __HERE__ ")");
  // signature to use in the veto list: const G4AffineTransform & G4NavigationHistory::GetTransform(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:106:35
  t35.method("GetTransform", static_cast<const G4AffineTransform & (G4NavigationHistory::*)(G4int)  const>(&G4NavigationHistory::GetTransform));

  DEBUG_MSG("Adding wrapper for G4int G4NavigationHistory::GetReplicaNo(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4NavigationHistory::GetReplicaNo(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:109:16
  t35.method("GetReplicaNo", static_cast<G4int (G4NavigationHistory::*)(G4int)  const>(&G4NavigationHistory::GetReplicaNo));

  DEBUG_MSG("Adding wrapper for EVolume G4NavigationHistory::GetVolumeType(G4int) (" __HERE__ ")");
  // signature to use in the veto list: EVolume G4NavigationHistory::GetVolumeType(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:112:18
  t35.method("GetVolumeType", static_cast<EVolume (G4NavigationHistory::*)(G4int)  const>(&G4NavigationHistory::GetVolumeType));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4NavigationHistory::GetVolume(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4NavigationHistory::GetVolume(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:115:29
  t35.method("GetVolume", static_cast<G4VPhysicalVolume * (G4NavigationHistory::*)(G4int)  const>(&G4NavigationHistory::GetVolume));

  DEBUG_MSG("Adding wrapper for void G4NavigationHistory::NewLevel(G4VPhysicalVolume *, EVolume, G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4NavigationHistory::NewLevel(G4VPhysicalVolume *, EVolume, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:118:15
  t35.method("NewLevel", static_cast<void (G4NavigationHistory::*)(G4VPhysicalVolume *, EVolume, G4int) >(&G4NavigationHistory::NewLevel));
  t35.method("NewLevel", [](G4NavigationHistory& a, G4VPhysicalVolume * arg0)->void{ a.NewLevel(arg0); });
  t35.method("NewLevel", [](G4NavigationHistory& a, G4VPhysicalVolume * arg0, EVolume arg1)->void{ a.NewLevel(arg0, arg1); });
  t35.method("NewLevel", [](G4NavigationHistory* a, G4VPhysicalVolume * arg0)->void{ a->NewLevel(arg0); });
  t35.method("NewLevel", [](G4NavigationHistory* a, G4VPhysicalVolume * arg0, EVolume arg1)->void{ a->NewLevel(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4NavigationHistory::BackLevel() (" __HERE__ ")");
  // signature to use in the veto list: void G4NavigationHistory::BackLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:123:15
  t35.method("BackLevel", static_cast<void (G4NavigationHistory::*)() >(&G4NavigationHistory::BackLevel));

  DEBUG_MSG("Adding wrapper for void G4NavigationHistory::BackLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4NavigationHistory::BackLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NavigationHistory.hh:127:15
  t35.method("BackLevel", static_cast<void (G4NavigationHistory::*)(G4int) >(&G4NavigationHistory::BackLevel));

  /* End of G4NavigationHistory class method wrappers
   **********************************************************************/

}
