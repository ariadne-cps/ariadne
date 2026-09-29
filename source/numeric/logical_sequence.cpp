/***************************************************************************
 *            numeric/logical_sequence.cpp
 *
 *  Copyright  2026  Ariadne contributors
 *
 ****************************************************************************/
#include "numeric/integer.hpp"
#include "numeric/sequence.hpp"
#include "numeric/logical_sequence.hpp"
namespace Ariadne {
namespace {
class SequenceDisjunction : public LogicalInterface {
    Sequence<LowerKleenean> _sequence;
  public: explicit SequenceDisjunction(Sequence<LowerKleenean> s):_sequence(s){}
  private:
    LogicalInterface* _copy() const override { return new SequenceDisjunction(*this); }
    LogicalValue _check(Effort e) const override {
        for(Natural k=0u;k!=e.work();++k) if(definitely(_sequence[k].check(e))) return LogicalValue::TRUE;
        return LogicalValue::INDETERMINATE;
    }
    OutputStream& _write(OutputStream& os) const override { return os<<"disjunction("<<_sequence[0u]<<","<<_sequence[1u]<<","<<_sequence[2u]<<",...)"; }
};
class SequenceConjunction : public LogicalInterface {
    Sequence<UpperKleenean> _sequence;
  public: explicit SequenceConjunction(Sequence<UpperKleenean> s):_sequence(s){}
  private:
    LogicalInterface* _copy() const override { return new SequenceConjunction(*this); }
    LogicalValue _check(Effort e) const override {
        for(Natural k=0u;k!=e.work();++k) if(definitely(not _sequence[k].check(e))) return LogicalValue::FALSE;
        return LogicalValue::INDETERMINATE;
    }
    OutputStream& _write(OutputStream& os) const override { return os<<"conjunction("<<_sequence[0u]<<","<<_sequence[1u]<<","<<_sequence[2u]<<",...)"; }
};
}
LowerKleenean disjunction(Sequence<LowerKleenean> const& s) {
    return LowerKleenean(LogicalHandle(make_handle<const SequenceDisjunction>(s)));
}
UpperKleenean conjunction(Sequence<UpperKleenean> const& s) {
    return UpperKleenean(LogicalHandle(make_handle<const SequenceConjunction>(s)));
}
}
