#ifndef ARIADNE_STANDARD_BOUNDER_H
#define ARIADNE_STANDARD_BOUNDER_H

#include <boost/shared_ptr.hpp>

#include "bounder_interface.h"

namespace Ariadne {
  

    /*! \ingroup Integrators
     *  \brief A class for bounding the flow of a vector field.
     */
    class StandardBounder
      : public BounderInterface
    {
     public:
      //@{ 
      //! \name Destructors, constructors and cloning operations.
 
      /*! \brief Default constructor. */
      StandardBounder();

      /*! \brief Construct with a given maximum step size. */
      StandardBounder(const Rational& max_step_size);

      /*! \brief Make a dynamically-allocated copy. */
      StandardBounder* clone() const;
     
      //@}

      //@{ 
      //! \name Data access.
      Rational maximum_step_size() const;
      //@}

      //@{ 
      //! \name Methods for computing bounding set for a flow. */

      /*! \brief Computes an integration time and a bounding box for the flow of \a vector_field starting in \a initial_set. */
      virtual std::pair< Rational, Box >
      flow_bounds(const FunctionInterface& vector_field,
                  const Box& initial_set) const;

      /*! \brief Computes an integration time and a bounding box for the flow of \a vector_field starting in \a initial_set. */
      virtual std::pair< Rational, Box >
      flow_bounds(const FunctionInterface& vector_field,
                  const Box& initial_set,
                  const Rational& maximum_step_size) const;


     /*! \brief Verifies that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. 
       *
       *  This method may return \a false even if the flow remains in \a bound. 
       */
      virtual bool check_flow_bounds(const FunctionInterface& vector_field,
                                     const Box& initial_set,
                                     const Rational& integration_time,
                                     const Box& bound) const;
      
      /*! \brief Compute a set \a bound such that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to \a integration_time, given a bound \a estimated_bound. */
      virtual Box refine_flow_bounds(const FunctionInterface& vector_field,
                                     const Box& initial_set,
                                     const Box& estimated_bound,
                                     const Rational& integration_time) const;

     private:
      Rational _maximum_step_size;
     public:
      uint verbosity;
    };





  
} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_H */
