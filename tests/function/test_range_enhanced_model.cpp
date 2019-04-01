#include "function/range_enhanced_model.hpp"


using namespace Ariadne;

int main() {
    typedef TaylorModel<ValidatedTag,FloatDP> TaylorModelType;
    typedef RangeEnhancedModel<TaylorModel<ValidatedTag,FloatDP>> RangeEnhancedTaylorModelType;

    typedef ScaledFunctionPatch<TaylorModel<ValidatedTag,FloatDP>> TaylorFunctionModel;
    typedef ScaledFunctionPatch<RangeEnhancedModel<TaylorModel<ValidatedTag,FloatDP>>> RangeEnhancedTaylorFunctionModel;

//    TaylorModelType tm;
//    RangeEnhancedTaylorModelType rtm = tm;

    TaylorFunctionModel tf;
    RangeEnhancedTaylorFunctionModel rtf;

    rtf=rtf+rtf;
    foo(rtf);
    rtf=tf;
}
