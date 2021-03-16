#include "1D_pde.hpp"

//using namespace std;
namespace Ariadne{

    template<class X>
    Array<FloatValue<X>> linspace1d(FloatValue<X> L, SizeType n)
    {
        Array<FloatValue<X>> linspaced(n, FloatValue<X>(0.0));
        if (n == 0)
            return linspaced;
        if (n == 1)
        {
            linspaced[0] = L;
            return linspaced;
        }
        FloatValue<X> delta = (L/(n - 1));
        for (SizeType i = 0; i < (n - 1); i++)
        {
            linspaced[i] = (0 + delta*i);
        }
        linspaced[n - 1] = L;
        
        return linspaced;
    }

    // Set initial condition
    template<class X>
    Tensor<2, FloatValue<X>> set_ic(Tensor<2, FloatValue<X>>& uts,/* std::function<FloatValue<X>(FloatValue<X>)> &phi0, */SizeType Nx, Array<FloatValue<X>> spacePoint, Parameter1D<X>& stringModel, bool isTriangular)
    {
        if (!isTriangular)
        {
            for (SizeType i = 0; i < Nx; i++){
            //uts[{i, 0}] = phi0(spacePoint[i]);
                uts[{i, 0}] = sin(stringModel.k*spacePoint[i]); //Sinusoind init condition
            }
        }
        else{
            for (SizeType i = 0; i < Nx; i++){
                if (definitely(spacePoint[i] < stringModel.x0))
                    uts[{i, 0}] = (stringModel.amp*(spacePoint[i]/stringModel.x0));//.value();
                else
                    uts[{i, 0}] (stringModel.amp*(stringModel.length - spacePoint[i])/(stringModel.length - stringModel.x0));//.value();
            }
        }

        return uts;
    }

    //Solving one dimensional pde
    template<class PR>
    Tensor<2, FloatValue<PR>> pde_1d_solver(/*std::function<FloatValue<PR>(FloatValue<PR>)>& phi0, std::function<FloatValue<PR>(FloatValue<PR>, FloatValue<PR>)>& source, */Parameter1D<PR>& stringParameter, SizeType Nx)
    {
        FloatValue<PR> c = sqrt((stringParameter.tension/(stringParameter.mass/stringParameter.length)));
        //FloatDPValue c = sqrt((stringParameter.tension/(stringParameter.mass/stringParameter.length))).value();

        FloatValue<PR> T = 0.060;
        FloatValue<PR> zb = 0;

        FloatValue<PR> C2 = pow(stringParameter.CourantNumber, 2);
        //FloatDPValue C2 = pow(stringParameter.CourantNumber, 2).value();

        Array<FloatValue<PR>> space = linspace1d(stringParameter.length, Nx);
        //Array<FloatDPValue> space = linspace1D(stringParameter.length, Nx);

        FloatValue<PR> dx = (space[1] - space[0]);
        //FloatDPValue dx = (space[1] - space[0]).value();
        FloatValue<PR> dt = (stringParameter.CourantNumber*dx/c);
        //FloatDPValue dt = (stringParameter.CourantNumber*dx/c).value();
        //FFloatDPValue Nt = round(T/dt).value();
        FloatValue<PR> Nt = round(T/dt);
        //FloatDPValue Nt = round(T/dt).value();
        SizeType Ntime = Nt.get_d();

        Array<FloatValue<PR>> time = linspace1d(T, Ntime);

        Tensor<2, FloatValue<PR>> uts({Nx, Ntime}, zb);

        uts = set_ic(uts, /*phi0, */Nx, space);

        // Set first Time step
        SizeType n = 0;
        uts[{0, n}] = 0;
        uts[{Nx - 1, n}] = 0;
        for (SizeType i = 1; i < Nx - 1; i++)
        {
            uts[{i, n+1}] = (uts[{i, n}] - C2*(uts[{i-1, n}] - 2u*uts[{i, n}] + uts[{i+1, n}])/2u
                + pow(dt,2)*sin(stringParameter.k*space[i])*cos(stringParameter.omega*time[n]));//source(space[i], time[n]));//.value();
        }// first time step

        // For each time step - main loop
        for (n = 1; n < Nt.get_d(); n++)
        {
            uts[{0, n}] = 0;
            uts[{Nx - 1, n}] = 0;
            // Update inner points
            for (SizeType i = 1; i < Nx - 1; i++)
            {
                uts[{i, n+1}] = ((1/(stringParameter.damping*dt/2u+1))*(2*uts[{i, n}] - uts[{i, n-1}] + (stringParameter.damping*dt/2)*uts[{i, n-1}] +
                    C2*(uts[{i-1, n}] - 2*uts[{i, n}] + uts[{i+1, n}])) + 1000u*pow(dt, 2)*sin(stringParameter.k*space[i])*cos(stringParameter.omega*time[n]));//source(space[i], time[n]));//.value();
            }
        }//main loop
        return uts;
    }
        
}// namespace Ariadne