/***************************************************************************
 *            test_gnuplot.cpp
 *
 *  Copyright  2008-21  Mirko Albanese
 *
 ****************************************************************************
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */



#include "../test.hpp"
#include "ariadne.hpp"
#include "numeric/numeric.hpp"
#include "numeric/float_bounds.hpp"
#include "algebra/tensor.hpp"
#include "algebra/spec_def.t.hpp"

#include "dynamics/1D_pde.hpp"
#include "dynamics/2D_pde.hpp"

using namespace Ariadne;
using namespace std;

class TestPdes
{
    public:

        TestPdes(){ };
        ~TestPdes(){ };

        template<class PR>
        void test(PR pr)
        {
            ARIADNE_TEST_CALL(sinusoidal_string_animation(pr));
            ARIADNE_TEST_CALL(triangular_string_animation(pr));
            ARIADNE_TEST_CALL(gauss_3d_animation(pr));
        }

        template<class PR>
        void sinusoidal_string_animation(PR pr)
        {
            SizeType Nx = 200;             // # point in space

            Parameter1D<PR> stringModel(pr);

            stringModel.length = 200;  // Length
            stringModel.tension = 10000;
            stringModel.mass = cast_exact(ApproximateDouble(0.9));

            FloatValue<PR> c = sqrt((stringModel.tension/(stringModel.mass/stringModel.length))).value();

            FloatValue<PR> wavelength = (stringModel.length*2).value();

            FloatValue<PR> frequency = (c/wavelength).value(); // Frequency

            stringModel.amp = cast_exact(ApproximateDouble(0.8));
            stringModel.damping = 100;

            stringModel.CourantNumber = cast_exact(ApproximateDouble(0.8));

            stringModel.k = ((2*pi)/wavelength).value();
            stringModel.omega = (2*pi*frequency).value();

            stringModel.is_triangular = false;

            Tensor<2, FloatValue<PR>> data = pde_1d_solver(stringModel, Nx, pr);

            FloatValue<PR> tollerance(cast_exact(ApproximateDouble(1.2)),pr);

            for (SizeType time = 0; time < data.size(1); time++)
            {
                for (SizeType x1 = 0; x1 < Nx; x1++)
                {
                    if ( (data[{x1, time}]) > (stringModel.amp+tollerance).value() )
                    {
                        ARIADNE_TEST_ASSERT(abs(data[{x1, time}]) <= (stringModel.amp+tollerance));
                        return;
                    }
                    
                }
                
            }
            

        }//String Evolution over time

        template<class PR>
        void triangular_string_animation(PR pr)
        {
            SizeType Nx = 200;             // # point in space

            Parameter1D<PR> stringModel(pr);

            stringModel.length = 200;  // Length
            stringModel.tension = 10000;
            stringModel.mass = cast_exact(ApproximateDouble(0.9));

            FloatValue<PR> c = sqrt((stringModel.tension/(stringModel.mass/stringModel.length))).value();

            stringModel.x0 = (0.85_q*stringModel.length).value();  // Point of max amplitube - Triangular IC
            FloatValue<PR> wavelength = (stringModel.length*2).value();

            FloatValue<PR> frequency = (c/wavelength).value(); // Frequency

            stringModel.amp = cast_exact(ApproximateDouble(0.8));
            stringModel.damping = 100;

            stringModel.CourantNumber = cast_exact(ApproximateDouble(0.8));

            stringModel.k = ((2*pi)/wavelength).value();
            stringModel.omega = (2*pi*frequency).value();

            stringModel.is_triangular = true;

            FloatValue<PR> tollerance(cast_exact(ApproximateDouble(1.2)),pr);

            Tensor<2, FloatValue<PR>> data = pde_1d_solver(stringModel, Nx, pr);

            for (SizeType time = 0; time < data.size(1); time++)
            {
                for (SizeType x1 = 0; x1 < Nx; x1++)
                {
                    if ( (data[{x1, time}]) > (stringModel.amp+tollerance).value() )
                    {
                        ARIADNE_TEST_ASSERT(abs(data[{x1, time}]) <= (stringModel.amp+tollerance));
                        return;
                    }
                    
                }
                
            }

        }//String Evolution over time

        template<class PR>
        void gauss_3d_animation(PR pr)
        {
            int Nx = 25;   //Mesh 1° dim
            int Ny = 25;   //Mesh 2° dim

            Parameter2D<PR> firstDim(pr), secondDim(pr);
            firstDim.length = 25;
            secondDim.length = 25;

            firstDim.x0 = 15;
            secondDim.x0 = 15;

            FloatValue<PR> tollerance(cast_exact(ApproximateDouble(1.2)),pr);

            Tensor<3, FloatValue<PR>> data = pde_2d_solver(firstDim, secondDim, SizeType(Nx), SizeType(Ny), pr);

            for (SizeType time = 0; time < data.size(1); time++)
            {
                for (SizeType x2 = 0; x2 < SizeType(Ny); x2++)
                {
                    for(SizeType x1 = 0; x1 < SizeType(Nx); x1++)
                    {
                        if ( (data[{x1, x2, time}]) > tollerance )
                        {
                            ARIADNE_TEST_ASSERT(abs(data[{x1, x2, time}]) <= tollerance);
                            return;
                        }
                    }                        
                }                
            }

        }
};

int main(int argc, const char** argv) {

    TestPdes testPdes;

    testPdes.test(double_precision);

    return ARIADNE_TEST_FAILURES;
}



