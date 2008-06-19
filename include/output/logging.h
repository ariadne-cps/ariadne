/***************************************************************************
 *            output/logging.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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



#ifndef _ARIADNE_LOGGING_H
#define _ARIADNE_LOGGING_H

#include <iostream>
#include <fstream>

/*! Send a message to the global logging stream. */
#define ARIADNE_LOG(level,msg)                  \
  if(verbosity >= level) { std::clog << msg << std::flush; }

namespace Ariadne {
  
  
    // Global log output file
    extern std::ofstream log_file_stream;
  
    /*! \brief Redirect logging output to file \a filename. */
    void redirect_log(const char* filename);

   
    extern int 
      ariadne_verbosity,
      numeric_verbosity,
      linear_algebra_verbosity,
      function_verbosity,
      geometry_verbosity,
      solver_verbosity,
      applicator_verbosity, 
      integrator_verbosity, 
      detector_verbosity, 
      hybrid_evolver_verbosity
    ; 

    extern int verbosity;
   
    /*! \brief Set the verbosity level for the %Numeric module. */
    inline void set_numeric_verbosity(int v) { verbosity=v; }
    /*! \brief Set the verbosity level for the Linear Algebra module. */
    inline void set_linear_algebra_verbosity(int v) { verbosity=v; }
    /*! \brief Set the verbosity level for the Linear Programming module. */
    inline void set_linear_programming_verbosity(int v) { verbosity=v; }
    /*! \brief Set the verbosity level for the %Combinatoric module. */
    inline void set_combinatoric_verbosity(int v) { verbosity=v; }
    /*! \brief Set the verbosity level for the %Function module. */
    inline void set_function_verbosity(int v) { verbosity=v; }
    /*! \brief Set the verbosity level for the %Geometry module. */
    inline void set_geometry_verbosity(int v) { verbosity=v; }
    /*! \brief Set the verbosity level for the %System module. */
    inline void set_system_verbosity(int v) { verbosity=v; }
    /*! \brief Set the verbosity level for the %Evaluation module. */
    inline void set_evaluation_verbosity(int v) { 
      solver_verbosity=v;
      applicator_verbosity=v;
      integrator_verbosity=v; 
      detector_verbosity=v; 
      hybrid_evolver_verbosity=v;
    }
    /*! \brief Set the verbosity level for the %Solver class. */
    inline void set_solver_verbosity(int v) { solver_verbosity=v; }
    /*! \brief Set the verbosity level for the %Applicator class. */
    inline void set_applicator_verbosity(int v) { applicator_verbosity=v; }
    /*! \brief Set the verbosity level for the %Integrator class. */
    inline void set_integrator_verbosity(int v) { integrator_verbosity=v; }
    /*! \brief Set the verbosity level for the %Detector class. */
    inline void set_detector_verbosity(int v) { detector_verbosity=v; }
    /*! \brief Set the verbosity level for the %HybridEvolver class. */
    inline void set_hybrid_evolver_verbosity(int v) { hybrid_evolver_verbosity=v; }

    /*! \brief Set the verbosity level for the %Input class. */
    inline void set_input_verbosity(int v) { verbosity=v; }


} // namespace Ariadne

#endif // _ARIADNE_LOGGING_H
