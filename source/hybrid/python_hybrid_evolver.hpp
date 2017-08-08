/***************************************************************************
 *            python_hybrid_evolver.hpp
 *
 *  Copyright  2009  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file python_hybrid_evolver.hpp
 *  \brief Wrapper for an evolver for hybrid systems written in Python
 */

#ifndef ARIADNE_PYTHON_HYBRID_EVOLVER_HPP
#define ARIADNE_PYTHON_HYBRID_EVOLVER_HPP

#ifdef HAVE_BOOST_PYTHON_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/python.hpp>

#include "utility/tuple.hpp"

#include "function/taylor_function.hpp"
#include "taylor_set.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_time.hpp"

#include "dynamics/orbit.hpp"

#include "hybrid/hybrid_automata.hpp"
#include "dynamics/evolver_interface.hpp"
#include "dynamics/evolver_base.hpp"
#include "evolution_parameters.hpp"

#include "utility/logging.hpp"

namespace Ariadne {

template<class Sys, class BS> class Evolver;

class ValidatedScalarTaylorFunctionModelDP;
class ValidatedVectorTaylorFunctionModelDP;
class TaylorConstrainedImageSet;
typedef Pair<DiscreteLocation,TaylorConstrainedImageSet> HybridTaylorConstrainedImageSet;
template<class ES> class Orbit;

class EvolutionParameters;
template<class X> class TaylorModel;
template<class MDL> class CalculusInterface;

class EvolutionProfiler;

class HybridTime;



/*! \brief A class for computing the evolution of a hybrid system.
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class PythonHybridEvolver
    : public EvolverBase<HybridAutomaton,HybridTaylorConstrainedImageSet>
    , public Loggable
{
    typedef VectorFunction FunctionType;
    typedef Vector<ExactIntervalType> BoxType;
    typedef ValidatedVectorTaylorFunctionModelDP FunctionModelType;
    typedef ValidatedVectorTaylorFunctionModelDP MapModelType;
    typedef ValidatedVectorTaylorFunctionModelDP FlowModelType;
    typedef ValidatedScalarTaylorFunctionModelDP ConstraintModelType;
    typedef TaylorModel TimeModelType;
    typedef TaylorConstrainedImageSet SetModelType;
    typedef TaylorConstrainedImageSet TimedSetModelType;
  public:
    typedef ContinuousEvolutionParameters EvolutionParametersType;
    typedef HybridAutomaton::TimeType TimeType;
    typedef Int IntegerType;
    typedef FloatDP RealType;
    typedef List<DiscreteEvent> EventListType;
    typedef HybridAutomaton SystemType;
    typedef TaylorConstrainedImageSet ContinuousEnclosureType;
    typedef Pair<DiscreteLocation,TaylorConstrainedImageSet> HybridEnclosureType;
    typedef HybridEnclosureType EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
    typedef FloatDP ContinuousTimeType;
  public:

    //! \brief Default constructor.
    PythonHybridEvolver() : _parameters() { this->initialise_python(); }

    //! \brief Construct from parameters using a default integrator.
    PythonHybridEvolver(const EvolutionParametersType& parameters) : _parameters(parameters) { this->initialise_python(); }

    /*! \brief Make a dynamically-allocated copy. */
    PythonHybridEvolver* clone() const { return new PythonHybridEvolver(*this); }

    //@{
    //! \name Parameters controlling the evolution.
    //! \brief A reference to the parameters controlling the evolution.
    EvolutionParametersType& parameters() { return this->_parameters; }
    const EvolutionParametersType& parameters() const { return this->_parameters; }

    //@}


    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using upper semantics.
    Orbit<EnclosureType> orbit(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;


    //! \brief Compute an approximation to the evolution set using upper semantics.
    EnclosureListType evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,UPPER_SEMANTICS,false);
        return final; }

    //! \brief Compute an approximation to the evolution set under upper semantics.
    EnclosureListType reach(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,UPPER_SEMANTICS,true);
        return intermediate; }

  protected:
    typedef Tuple<DiscreteLocation, EventListType, SetModelType, TimeModelType> HybridTimedSetType;

    // This is the only method which is called in Python
    virtual Void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const SystemType& system, const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, Bool reach) const;

    Void initialise_python();

  private:
    EvolutionParametersType _parameters;
};


Orbit<PythonHybridEvolver::EnclosureType>
PythonHybridEvolver::
orbit(const SystemType& system,
      const EnclosureType& initial_set,
      const TimeType& time,
      Semantics semantics) const
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,
                     system,initial_set,time,semantics,false);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}

enum Bool { False, True };

template<class SET>
ListSet< Pair<DiscreteLocation,SET> >*
make_hybrid_list_set(const boost::python::list& pylst)
//make_hybrid_list_set(const boost::python::object& pyobj)
{
    ListSet< Pair<DiscreteLocation,SET> >* result=new ListSet< Pair<DiscreteLocation,SET> >();
    //boost::python::list pylst=boost::python::extract<boost::python::list>(pyobj);
    for(Int i=0; i!=len(pylst); ++i) {
        boost::python::Tuple pytup=boost::python::extract<boost::python::Tuple>(pylst[i]);
        Ariadne::DiscreteLocation q(boost::python::extract<Int>(pytup[0]));
        SET s(boost::python::extract<SET>(pytup[1]));
        //Pair<Ariadne::DiscreteLocation,SET> pr=std::make_pair(q,s);
        //result->adjoin(std::make_pair(q,s));
    }
    return result;
}

Void
PythonHybridEvolver::initialise_python()
{
        Py_Initialize();

        boost::python::object main_module(boost::python::handle<>(PyImport_ImportModule("__main__")));
        boost::python::object main_namespace = main_module.attr("__dict__");

        boost::python::object evolution_module(boost::python::handle<>(PyImport_ImportModule("hybrid_evolver")));
        boost::python::object evolution_namespace = evolution_module.attr("__dict__");
        main_namespace["hybrid_evolver"]=evolution_module;

        boost::python::class_< ListSet< Pair<DiscreteLocation,TaylorConstrainedImageSet> > >
            enclosure_list_class("HybridTaylorConstrainedImageSetList",boost::python::no_init);
        enclosure_list_class.def("__init__",boost::python::make_constructor(&make_hybrid_list_set<TaylorConstrainedImageSet>));

        boost::python::class_<HybridAutomaton>("HybridAutomaton",boost::python::no_init);
        boost::python::class_<EnclosureType>("HybridTaylorConstrainedImageSet",boost::python::no_init);
        boost::python::class_<HybridTime>("HybridTime",boost::python::no_init);
        //boost::python::enum_<Semantics>("Semantics")
        //    .value("UPPER_SEMANTICS", UPPER_SEMANTICS).value("LOWER_SEMANTICS", LOWER_SEMANTICS);

        evolution_namespace["HybridTaylorConstrainedImageSetList"]=enclosure_list_class;
        //main_namespace["TaylorConstrainedImageSetList"]=enclosure_list_class;
}


Void
PythonHybridEvolver::_evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                const SystemType& system, const EnclosureType& initial, const TimeType& time,
                                Semantics semantics, Bool reach) const
{
    try {
        boost::python::object main_module(boost::python::handle<>(PyImport_ImportModule("__main__")));
        boost::python::object main_namespace = main_module.attr("__dict__");

        main_namespace["system"]=system;
        main_namespace["initial"]=initial;
        main_namespace["time"]=time;
        main_namespace["semantics"]=Int(semantics);

        //PyRun_SimpleString("print\nprint \"Starting Python evaluate\"\nprint dir()\nprint\nprint dir(hybrid_evolver)\nprint\n");
        PyRun_SimpleString("evolver=hybrid_evolver.HybridEvolverPrototype()");
        PyRun_SimpleString("(final,reachable,intermediate)=evolver.orbit(system,initial,time,semantics)");

        final=boost::python::extract<EnclosureListType>(main_namespace["final"]);
        reachable=boost::python::extract<EnclosureListType>(main_namespace["reachable"]);
        intermediate=boost::python::extract<EnclosureListType>(main_namespace["intermediate"]);
        ARIADNE_LOG(2,"final.size()="<<final.size()<<", "<<
                      "reachable.size()="<<reachable.size()<<", "<<
                      "intermediate.size()="<<intermediate.size()<<"\n");
    }
    catch( boost::python::error_already_set ) {
        PyErr_Print();
    }
    return;
}


} // namespace Ariadne

#endif // HAVE_BOOST_PYTHON_HPP

#endif // ARIADNE_PYTHON_HYBRID_EVOLVER_HPP
