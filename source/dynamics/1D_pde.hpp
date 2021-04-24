#include "numeric/numeric.hpp"
#include "utility/array.hpp"
#include "algebra/tensor.hpp"

namespace Ariadne
{
    template<class PR>
    class Parameter1D
        {
            public:
                Parameter1D(PR pr) :    length(cast_exact(ApproximateDouble(0.0)), pr),
                                        tension(cast_exact(ApproximateDouble(0.0)), pr),
                                        mass(cast_exact(ApproximateDouble(0.0)), pr),
                                        damping(cast_exact(ApproximateDouble(0.0)), pr),
                                        CourantNumber(cast_exact(ApproximateDouble(0.0)), pr),
                                        x0(cast_exact(ApproximateDouble(0.0)), pr),
                                        amp(cast_exact(ApproximateDouble(0.0)), pr),
                                        k(cast_exact(ApproximateDouble(0.0)), pr),
                                        omega(cast_exact(ApproximateDouble(0.0)), pr),
                                        is_triangular(false)
                {}
                FloatValue<PR> length;
                FloatValue<PR> tension;
                FloatValue<PR> mass;
                FloatValue<PR> damping;
                FloatValue<PR> CourantNumber;
                FloatValue<PR> x0;
                FloatValue<PR> amp;
                FloatValue<PR> k;
                FloatValue<PR> omega;
                Bool is_triangular;

        };

    template<class X>
    Array<FloatValue<X>> linspace1d(FloatValue<X> L, SizeType n, X pr)
    {
        Array<FloatValue<X>> linspaced(n, FloatValue<X>(cast_exact(ApproximateDouble(0.0)), pr));
        if (n == 0)
            return linspaced;
        if (n == 1)
        {
            linspaced[0] = L;
            return linspaced;
        }
        FloatValue<X> delta = (L/n).value();
        for (SizeType i = 0; i < (n - 1); i++)
        {
            linspaced[i] = (0 + delta*i).value();
        }
        linspaced[n - 1] = L;
        
        return linspaced;
    }

    // Set initial condition
    template<class X>
    Tensor<2, FloatValue<X>> set_ic(Tensor<2, FloatValue<X>>& uts, SizeType Nx, Array<FloatValue<X>> spacePoint, Parameter1D<X>& stringModel, bool isTriangular)
    {
        if (!isTriangular)
        {   //Sinusoind init condition
            for (SizeType i = 0; i < Nx; i++){
                uts[{i, 0}] = (stringModel.amp*sin(stringModel.k*spacePoint[i])).value(); 
            }
        }
        else{
            //Triangular init condition
            for (SizeType i = 0; i < Nx; i++){
                if (definitely(spacePoint[i] < stringModel.x0))
                    uts[{i, 0}] = (stringModel.amp*(spacePoint[i]/stringModel.x0)).value();
                else
                    uts[{i, 0}] = (stringModel.amp*(stringModel.length - spacePoint[i])/(stringModel.length - stringModel.x0)).value();
            }
        }

        return uts;
    }

    //Solving one dimensional pde
    template<class PR>
    Tensor<2, FloatValue<PR>> pde_1d_solver(Parameter1D<PR>& stringParameter, SizeType Nx, PR pr)
    {
        FloatValue<PR> c = sqrt((stringParameter.tension/(stringParameter.mass/stringParameter.length))).value();

        FloatValue<PR> T(cast_exact(ApproximateDouble(0.03)), pr);
        FloatValue<PR> zb(cast_exact(ApproximateDouble(0.0)), pr);

        FloatValue<PR> C2 = pow(stringParameter.CourantNumber, 2).value();

        Array<FloatValue<PR>> space = linspace1d(stringParameter.length, Nx, pr);

        FloatValue<PR> dx = (space[1] - space[0]).value();
        FloatValue<PR> dt = (stringParameter.CourantNumber*dx/c).value();

        FloatValue<PR> Nt = round(T/dt).value();
        SizeType Ntime = Nt.get_d();

        Array<FloatValue<PR>> time = linspace1d(T, Ntime, pr);

        Tensor<2, FloatValue<PR>> uts({Nx, Ntime}, zb);

        uts = set_ic(uts, Nx, space, stringParameter, true);

        // Set first Time step
        SizeType n = 0;
        uts[{0, n}] = 0;
        uts[{Nx - 1, n}] = 0;
        for (SizeType i = 1; i < Nx-1; i++)
        {
            uts[{i, n+1}] = (uts[{i, n}] - C2*(uts[{i-1, n}] - 2u*uts[{i, n}] + uts[{i+1, n}])/2u
                + pow(dt,2)*sin(stringParameter.k*space[i])*cos(stringParameter.omega*time[n])).value();
        }// first time step

        // For each time step - main loop
        for (n = 1; n < Nt.get_d(); n++)
        {
            uts[{0, n}] = 0;
            uts[{Nx - 1, n}] = 0;
            // Update inner points
            for (SizeType i = 1; i < Nx-1; i++)
            {
                uts[{i, n+1}] = ((1/(stringParameter.damping*dt/2u+1))*(2*uts[{i, n}] - uts[{i, n-1}] + (stringParameter.damping*dt/2)*uts[{i, n-1}] +
                    C2*(uts[{i-1, n}] - 2*uts[{i, n}] + uts[{i+1, n}])) + pow(dt, 2)*sin(stringParameter.k*space[i])*cos(stringParameter.omega*time[n])).value();
            }
        }//main loop
        return uts;
    }




}   // namespace Ariadne

