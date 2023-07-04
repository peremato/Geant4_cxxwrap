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
void add_methods_for_G4VisAttributes(jlcxx::Module& types, jlcxx::TypeWrapper<G4VisAttributes>& t51) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4VisAttributes
   */


  DEBUG_MSG("Adding wrapper for void G4VisAttributes::G4VisAttributes(G4bool) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:79:3
  t51.constructor<G4bool>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4VisAttributes::G4VisAttributes(const G4Colour &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:80:3
  t51.constructor<const G4Colour &>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4VisAttributes::G4VisAttributes(G4bool, const G4Colour &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:81:3
  t51.constructor<G4bool, const G4Colour &>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4VisAttributes::G4VisAttributes(const G4VisAttributes &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:82:3
  t51.constructor<const G4VisAttributes &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for G4VisAttributes & G4VisAttributes::operator=(const G4VisAttributes &) (" __HERE__ ")");
  // signature to use in the veto list: G4VisAttributes & G4VisAttributes::operator=(const G4VisAttributes &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:84:20
  t51.method("assign", static_cast<G4VisAttributes & (G4VisAttributes::*)(const G4VisAttributes &) >(&G4VisAttributes::operator=));

  DEBUG_MSG("Adding wrapper for const G4VisAttributes & G4VisAttributes::GetInvisible() (" __HERE__ ")");
  // signature to use in the veto list: const G4VisAttributes & G4VisAttributes::GetInvisible()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:86:33
  t51.method("G4VisAttributes!GetInvisible", static_cast<const G4VisAttributes & (*)() >(&G4VisAttributes::GetInvisible));
  types.set_override_module(jl_base_module);

  DEBUG_MSG("Adding wrapper for G4bool G4VisAttributes::operator!=(const G4VisAttributes &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VisAttributes::operator!=(const G4VisAttributes &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:88:10
  t51.method("!=", static_cast<G4bool (G4VisAttributes::*)(const G4VisAttributes &)  const>(&G4VisAttributes::operator!=));

  DEBUG_MSG("Adding wrapper for G4bool G4VisAttributes::operator==(const G4VisAttributes &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VisAttributes::operator==(const G4VisAttributes &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:89:10
  t51.method("==", static_cast<G4bool (G4VisAttributes::*)(const G4VisAttributes &)  const>(&G4VisAttributes::operator==));

  types.unset_override_module();

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetVisibility(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetVisibility(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:91:8
  t51.method("SetVisibility", static_cast<void (G4VisAttributes::*)(G4bool) >(&G4VisAttributes::SetVisibility));
  t51.method("SetVisibility", [](G4VisAttributes& a)->void{ a.SetVisibility(); });
  t51.method("SetVisibility", [](G4VisAttributes* a)->void{ a->SetVisibility(); });

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetDaughtersInvisible(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetDaughtersInvisible(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:92:8
  t51.method("SetDaughtersInvisible", static_cast<void (G4VisAttributes::*)(G4bool) >(&G4VisAttributes::SetDaughtersInvisible));
  t51.method("SetDaughtersInvisible", [](G4VisAttributes& a)->void{ a.SetDaughtersInvisible(); });
  t51.method("SetDaughtersInvisible", [](G4VisAttributes* a)->void{ a->SetDaughtersInvisible(); });

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetColour(const G4Colour &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetColour(const G4Colour &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:93:8
  t51.method("SetColour", static_cast<void (G4VisAttributes::*)(const G4Colour &) >(&G4VisAttributes::SetColour));

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetColor(const G4Color &) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetColor(const G4Color &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:94:8
  t51.method("SetColor", static_cast<void (G4VisAttributes::*)(const G4Color &) >(&G4VisAttributes::SetColor));

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetColour(G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetColour(G4double, G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:95:8
  t51.method("SetColour", static_cast<void (G4VisAttributes::*)(G4double, G4double, G4double, G4double) >(&G4VisAttributes::SetColour));
  t51.method("SetColour", [](G4VisAttributes& a, G4double arg0, G4double arg1, G4double arg2)->void{ a.SetColour(arg0, arg1, arg2); });
  t51.method("SetColour", [](G4VisAttributes* a, G4double arg0, G4double arg1, G4double arg2)->void{ a->SetColour(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetColor(G4double, G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetColor(G4double, G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:97:8
  t51.method("SetColor", static_cast<void (G4VisAttributes::*)(G4double, G4double, G4double, G4double) >(&G4VisAttributes::SetColor));
  t51.method("SetColor", [](G4VisAttributes& a, G4double arg0, G4double arg1, G4double arg2)->void{ a.SetColor(arg0, arg1, arg2); });
  t51.method("SetColor", [](G4VisAttributes* a, G4double arg0, G4double arg1, G4double arg2)->void{ a->SetColor(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetLineStyle(G4VisAttributes::LineStyle) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetLineStyle(G4VisAttributes::LineStyle)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:99:8
  t51.method("SetLineStyle", static_cast<void (G4VisAttributes::*)(G4VisAttributes::LineStyle) >(&G4VisAttributes::SetLineStyle));

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetLineWidth(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetLineWidth(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:100:8
  t51.method("SetLineWidth", static_cast<void (G4VisAttributes::*)(G4double) >(&G4VisAttributes::SetLineWidth));

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetForceWireframe(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetForceWireframe(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:101:8
  t51.method("SetForceWireframe", static_cast<void (G4VisAttributes::*)(G4bool) >(&G4VisAttributes::SetForceWireframe));
  t51.method("SetForceWireframe", [](G4VisAttributes& a)->void{ a.SetForceWireframe(); });
  t51.method("SetForceWireframe", [](G4VisAttributes* a)->void{ a->SetForceWireframe(); });

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetForceSolid(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetForceSolid(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:102:8
  t51.method("SetForceSolid", static_cast<void (G4VisAttributes::*)(G4bool) >(&G4VisAttributes::SetForceSolid));
  t51.method("SetForceSolid", [](G4VisAttributes& a)->void{ a.SetForceSolid(); });
  t51.method("SetForceSolid", [](G4VisAttributes* a)->void{ a->SetForceSolid(); });

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetForceCloud(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetForceCloud(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:103:8
  t51.method("SetForceCloud", static_cast<void (G4VisAttributes::*)(G4bool) >(&G4VisAttributes::SetForceCloud));
  t51.method("SetForceCloud", [](G4VisAttributes& a)->void{ a.SetForceCloud(); });
  t51.method("SetForceCloud", [](G4VisAttributes* a)->void{ a->SetForceCloud(); });

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetForceNumberOfCloudPoints(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetForceNumberOfCloudPoints(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:104:8
  t51.method("SetForceNumberOfCloudPoints", static_cast<void (G4VisAttributes::*)(G4int) >(&G4VisAttributes::SetForceNumberOfCloudPoints));

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetForceAuxEdgeVisible(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetForceAuxEdgeVisible(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:106:8
  t51.method("SetForceAuxEdgeVisible", static_cast<void (G4VisAttributes::*)(G4bool) >(&G4VisAttributes::SetForceAuxEdgeVisible));
  t51.method("SetForceAuxEdgeVisible", [](G4VisAttributes& a)->void{ a.SetForceAuxEdgeVisible(); });
  t51.method("SetForceAuxEdgeVisible", [](G4VisAttributes* a)->void{ a->SetForceAuxEdgeVisible(); });

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetForceLineSegmentsPerCircle(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetForceLineSegmentsPerCircle(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:107:8
  t51.method("SetForceLineSegmentsPerCircle", static_cast<void (G4VisAttributes::*)(G4int) >(&G4VisAttributes::SetForceLineSegmentsPerCircle));

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetStartTime(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetStartTime(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:112:8
  t51.method("SetStartTime", static_cast<void (G4VisAttributes::*)(G4double) >(&G4VisAttributes::SetStartTime));

  DEBUG_MSG("Adding wrapper for void G4VisAttributes::SetEndTime(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4VisAttributes::SetEndTime(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:113:8
  t51.method("SetEndTime", static_cast<void (G4VisAttributes::*)(G4double) >(&G4VisAttributes::SetEndTime));

  DEBUG_MSG("Adding wrapper for G4bool G4VisAttributes::IsVisible() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VisAttributes::IsVisible()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:117:19
  t51.method("IsVisible", static_cast<G4bool (G4VisAttributes::*)()  const>(&G4VisAttributes::IsVisible));

  DEBUG_MSG("Adding wrapper for G4bool G4VisAttributes::IsDaughtersInvisible() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VisAttributes::IsDaughtersInvisible()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:118:19
  t51.method("IsDaughtersInvisible", static_cast<G4bool (G4VisAttributes::*)()  const>(&G4VisAttributes::IsDaughtersInvisible));

  DEBUG_MSG("Adding wrapper for const G4Colour & G4VisAttributes::GetColour() (" __HERE__ ")");
  // signature to use in the veto list: const G4Colour & G4VisAttributes::GetColour()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:119:19
  t51.method("GetColour", static_cast<const G4Colour & (G4VisAttributes::*)()  const>(&G4VisAttributes::GetColour));

  DEBUG_MSG("Adding wrapper for const G4Color & G4VisAttributes::GetColor() (" __HERE__ ")");
  // signature to use in the veto list: const G4Color & G4VisAttributes::GetColor()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:120:19
  t51.method("GetColor", static_cast<const G4Color & (G4VisAttributes::*)()  const>(&G4VisAttributes::GetColor));

  DEBUG_MSG("Adding wrapper for G4VisAttributes::LineStyle G4VisAttributes::GetLineStyle() (" __HERE__ ")");
  // signature to use in the veto list: G4VisAttributes::LineStyle G4VisAttributes::GetLineStyle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:121:19
  t51.method("GetLineStyle", static_cast<G4VisAttributes::LineStyle (G4VisAttributes::*)()  const>(&G4VisAttributes::GetLineStyle));

  DEBUG_MSG("Adding wrapper for G4double G4VisAttributes::GetLineWidth() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VisAttributes::GetLineWidth()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:122:19
  t51.method("GetLineWidth", static_cast<G4double (G4VisAttributes::*)()  const>(&G4VisAttributes::GetLineWidth));

  DEBUG_MSG("Adding wrapper for G4bool G4VisAttributes::IsForceDrawingStyle() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VisAttributes::IsForceDrawingStyle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:123:19
  t51.method("IsForceDrawingStyle", static_cast<G4bool (G4VisAttributes::*)()  const>(&G4VisAttributes::IsForceDrawingStyle));

  DEBUG_MSG("Adding wrapper for G4VisAttributes::ForcedDrawingStyle G4VisAttributes::GetForcedDrawingStyle() (" __HERE__ ")");
  // signature to use in the veto list: G4VisAttributes::ForcedDrawingStyle G4VisAttributes::GetForcedDrawingStyle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:124:22
  t51.method("GetForcedDrawingStyle", static_cast<G4VisAttributes::ForcedDrawingStyle (G4VisAttributes::*)()  const>(&G4VisAttributes::GetForcedDrawingStyle));

  DEBUG_MSG("Adding wrapper for G4int G4VisAttributes::GetForcedNumberOfCloudPoints() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VisAttributes::GetForcedNumberOfCloudPoints()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:125:19
  t51.method("GetForcedNumberOfCloudPoints", static_cast<G4int (G4VisAttributes::*)()  const>(&G4VisAttributes::GetForcedNumberOfCloudPoints));

  DEBUG_MSG("Adding wrapper for G4bool G4VisAttributes::IsForceAuxEdgeVisible() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VisAttributes::IsForceAuxEdgeVisible()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:126:19
  t51.method("IsForceAuxEdgeVisible", static_cast<G4bool (G4VisAttributes::*)()  const>(&G4VisAttributes::IsForceAuxEdgeVisible));

  DEBUG_MSG("Adding wrapper for G4bool G4VisAttributes::IsForcedAuxEdgeVisible() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VisAttributes::IsForcedAuxEdgeVisible()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:127:19
  t51.method("IsForcedAuxEdgeVisible", static_cast<G4bool (G4VisAttributes::*)()  const>(&G4VisAttributes::IsForcedAuxEdgeVisible));

  DEBUG_MSG("Adding wrapper for G4bool G4VisAttributes::IsForceLineSegmentsPerCircle() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4VisAttributes::IsForceLineSegmentsPerCircle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:128:19
  t51.method("IsForceLineSegmentsPerCircle", static_cast<G4bool (G4VisAttributes::*)()  const>(&G4VisAttributes::IsForceLineSegmentsPerCircle));

  DEBUG_MSG("Adding wrapper for G4int G4VisAttributes::GetForcedLineSegmentsPerCircle() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VisAttributes::GetForcedLineSegmentsPerCircle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:129:19
  t51.method("GetForcedLineSegmentsPerCircle", static_cast<G4int (G4VisAttributes::*)()  const>(&G4VisAttributes::GetForcedLineSegmentsPerCircle));

  DEBUG_MSG("Adding wrapper for G4double G4VisAttributes::GetStartTime() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VisAttributes::GetStartTime()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:130:19
  t51.method("GetStartTime", static_cast<G4double (G4VisAttributes::*)()  const>(&G4VisAttributes::GetStartTime));

  DEBUG_MSG("Adding wrapper for G4double G4VisAttributes::GetEndTime() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4VisAttributes::GetEndTime()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:131:19
  t51.method("GetEndTime", static_cast<G4double (G4VisAttributes::*)()  const>(&G4VisAttributes::GetEndTime));

  DEBUG_MSG("Adding wrapper for G4int G4VisAttributes::GetMinLineSegmentsPerCircle() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4VisAttributes::GetMinLineSegmentsPerCircle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VisAttributes.hh:132:19
  t51.method("G4VisAttributes!GetMinLineSegmentsPerCircle", static_cast<G4int (*)() >(&G4VisAttributes::GetMinLineSegmentsPerCircle));

  /* End of G4VisAttributes class method wrappers
   **********************************************************************/

}