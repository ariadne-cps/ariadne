/***************************************************************************
 *            numeric/double.code.h
 *
 *  Copyright  2008  Pieter Collins
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
 
namespace Ariadne {
namespace Numeric {

double exp(const double& x) {
  if(x>1.0) { double r=x/2; r=exp(r); return r*r; }
  else if(x<0) { return exp(1./x); }
  else if(x==0) { return 1.0; }
  else {
    static const int N=24;
    double r=0.0;
    double p=1.0;
    double t[N+1];
    t[0]=p;
    for(int i=1; i<=N; ++i) {
      p*=x;
      p/=i;
      t[i]=p;
    }
    for(int i=N; i>=0; --i) {
      r+=t[i];
    }
    return r;
  }
}

}
}
