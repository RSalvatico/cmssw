// Hook for setting shower scale in top and W resonances
// for Powheg ttb_NLO_dec and b_bbar_4l processes
// C++ port of algorithm by Jezo et. al. (arXiv:1607.04538, Appendix B.2)
#include "GeneratorInterface/Pythia8Interface/plugins/PowhegHooksBB4L.h"

class PowhegResHook : public Pythia8::UserHooks {
public:
  // Constructor and destructor.
  PowhegResHook() {}
  ~PowhegResHook() override {}

  bool canSetResonanceScale() override { return true; }

  double scaleResonance(const int iRes, const Pythia8::Event& event) override;

  //--------------------------------------------------------------------------

private:
  bool calcScales_;
};
