#include "2D_pde.hpp"

namespace Ariadne
{
    template<class PR>
    Array<FloatValue<PR>> linspace2d(FloatValue<PR> L, SizeType n)
    {
        Array<FloatValue<PR>> linspaced(n,FloatValue<PR>(0.0_exact,DoublePrecision()));
        if (n == 0)
            return linspaced;
        if (n == 1)
        {
            linspaced[0] = L;
            return linspaced;
        }
        FloatValue<PR> delta = (L/(n - 1)).value();
        for (SizeType i = 0; i < (n - 1); i++)
        {
            linspaced[i] = ((0 + delta*i)).value();
        }
        linspaced[n - 1] = L;
        
        return linspaced;
    }
template<class PR>
    Tensor<2, FloatValue<PR>> set_ic_2d(std::function<FloatValue<PR>(FloatValue<PR>, FloatValue<PR>)>& phi0, SizeType Nx, SizeType Ny, Array<FloatValue<PR>> spacePointX, Array<FloatValue<PR>> spacePointY)
    {
        FloatValue<PR> zb(0, DoublePrecision());
        Tensor<2, FloatValue<PR>> u({Nx, Ny}, zb);
        for (SizeType i = 0; i < Nx; i++)
        {
            for (SizeType j = 0; j < Ny; j++)
            {
                u[{i, j}] = phi0(spacePointX[i], spacePointY[j]);
            }
        }
        return u;
    }

template<class PR>
Tensor<3, FloatValue<PR>> pde_2d_solver(std::function<FloatValue<PR>(FloatValue<PR>, FloatValue<PR>)>& phi0, std::function<FloatValue<PR>(FloatValue<PR>, FloatValue<PR>, FloatValue<PR>)>& source, Parameter2D<PR>& firstDim, Parameter2D<PR>& secondDim, SizeType Nx, SizeType Ny)
    {
        auto spaceX = linspace2d(FloatValue<PR>(firstDim.length, DoublePrecision()), Nx);    //Mesh point x
        auto spaceY = linspace2d(FloatValue<PR>(secondDim.length, DoublePrecision()), Ny);   //Mesh point y

        auto dx = (spaceX[1] - spaceX[0]).value();
        auto dy = (spaceY[1] - spaceY[0]).value();

        FloatValue<PR> c = 10.0_exact;
        FloatValue<PR> T = 1.0_exact;

        FloatValue<PR> zb(0, DoublePrecision());

        FloatValue<PR> dt = 0.015625_exact;
        
        FloatValue<PR> Nt = round(T/dt).value();
        SizeType Ntime = Nt.get_d();
        auto time = linspace2d(T, Ntime);                   //Mesh point t

        FloatValue<PR> Cx = (c*dt/dx).value();                       //Courant Numbers
        FloatValue<PR> Cy = (c*dt/dy).value();
        FloatValue<PR> Cx2 = pow(Cx, 2).value();
        FloatValue<PR> Cy2 = pow(Cy, 2).value();

        Tensor<2, FloatValue<PR>> u({Nx, Ny}, zb);
        Tensor<3, FloatValue<PR>> uts({Nx, Ny, Ntime}, zb);

        u = set_ic_2d(phi0, Nx, Ny, spaceX, spaceY);
        for (SizeType x = 0; x < Nx; x++)
        {
            for (SizeType y = 0; y < Ny; y++)
            {
                uts[{x, y, 0}] = u[{x, y}];
            }   
        }

        //Compute first time step
        SizeType n = 0;
        for (SizeType i = 1; i < Nx - 1; i++)
        {
            for (SizeType j = 1; j < Ny - 1; j++)
            {
                auto u_xx = (uts[{i-1, j, n}] - 2*uts[{i, j, n}] + uts[{i+1, j, n}]).value();
                auto u_yy = (uts[{i, j-1, n}] - 2*uts[{i, j, n}] + uts[{i, j+1, n}]).value();
                uts[{i, j, n+1}] = (uts[{i, j, n}] + Cx2*u_xx/2u + Cy2*u_yy/2u + pow(dt, 2)*source(spaceX[i], spaceY[j], time[n])).value();
                //uts[{i, j, n+1}] += dt*VelIC(spaceX[i], spaceY[j]);
            }
        }
        //Compute each time step
        for (n = 1; n < Ntime; n++)
        {
            for (SizeType i = 1; i < Nx - 1; i++)
            {
                for (SizeType j = 1; j < Ny - 1; j++)
                {
                    auto u_xx = (uts[{i-1, j, n}] - 2*uts[{i, j, n}] + uts[{i+1, j, n}]).value();
                    auto u_yy = (uts[{i, j-1, n}] - 2*uts[{i, j, n}] + uts[{i, j+1, n}]).value();
                    uts[{i, j, n+1}] = (2*uts[{i, j, n}] - uts[{i, j, n-1}] + Cx2*u_xx + Cy2*u_yy + pow(dt, 2)*source(spaceX[i], spaceY[j], time[n])).value();
                }
                
            }
        }
    return uts;
    }


}//namespace Ariadne