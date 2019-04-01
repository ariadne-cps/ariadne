#include "numeric/numeric.hpp"
#include "algebra/algebra.hpp"
#include "algebra/series.hpp"
#include "algebra/operations.hpp"

#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"
#include "function/scaled_function_patch.hpp"

namespace Ariadne {

template<class M> class RangeEnhancedModel
    : DispatchAlgebraOperations<RangeEnhancedModel<M>,NumericType<M>>
    , public M
{
    using X = typename M::NumericType;
    using Y = typename M::GenericNumericType;
  public:
    typedef M ModelType;
    typedef typename M::RangeType RangeType;
    //typedef typename M::NumericType NumericType;
/*
    typedef typename M::Paradigm Paradigm;
    typedef typename M::PrecisionType PrecisionType;
    typedef typename M::ErrorPrecisionType ErrorPrecisionType;
    typedef typename M::RawFloatType RawFloatType;
    typedef typename M::CoefficientType CoefficientType;

*/
    ModelType _model; mutable RangeType _range; mutable bool _up_to_date;

  public:
    //RangeEnhancedModel(ModelType model, RangeType range) : _model(std::move(model)), _range(std::move(range)), _up_to_date(false) { }
    //RangeType range() const { if (not _up_to_date) { _range=intersection(_model.range(),_range); _up_to_date=true; } return _range; }

    RangeEnhancedModel(ModelType model, RangeType range) : ModelType(std::move(model)), _range(intersection(range,_model.range())) { }
    ModelType const& model() const { return *this; }
    RangeType const& range() const { return _range; }

//    RangeEnhancedModel<M>& operator=(M const&) = delete;
//    RangeEnhancedModel<M>& operator=(RangeEnhancedModel<M> const&);
    RangeEnhancedModel<M>& operator=(Y const& c) { this->M::operator=(c); _range=c; return *this; }
};

template<class M> struct AlgebraOperations<RangeEnhancedModel<M>> {
  public:
    template<class OP> static RangeEnhancedModel<M> apply(OP op, RangeEnhancedModel<M> const& m1, RangeEnhancedModel<M> const& m2) {
        return RangeEnhancedModel<M>(op(m1.model(),m2.model()),op(m1.range(),m2.range())); }
    template<class OP> static RangeEnhancedModel<M> apply(OP op, RangeEnhancedModel<M> const& m1, NumericType<M> const& c2) {
        return RangeEnhancedModel<M>(op(m1.model(),c2),op(m1.range(),c2)); }
    template<class OP> static RangeEnhancedModel<M> apply(OP op, RangeEnhancedModel<M> const& m) {
        return RangeEnhancedModel<M>(op(m.model()),op(m.range())); }
};


using BaseModelType = TaylorModel<ValidatedTag,FloatDP>;
template class ScaledFunctionPatch<RangeEnhancedModel<BaseModelType>>;

template<class... TS> void foo(FunctionModelInterface<TS...> const& f) { }

}
