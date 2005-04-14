/***************************************************************************
 *            model2.cc
 *
 *  Thu Feb  3 14:06:15 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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

#include <vector>

#include <linear_algebra.h>
#include <polyhedron.h>
#include <poly_map.h>
#include <ptree.h>
#include <vector_field.h>
#include <poly_vf.h>
#include <automaton.h>
#include <hybrid_set.h>
#include <solver.h>

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::Geometry::IO_Operators;
using namespace Ariadne::Map;
using namespace Ariadne::Map::Affine;
using namespace Ariadne::Evaluation;
using namespace Ariadne::VectorField;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::VectorField::Affine;
using namespace Ariadne::HybridDefinitions;
using namespace Ariadne::HybridDefinitions::IO_Operators;

namespace UBLAS = boost::numeric::ublas;

using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;

typedef AriadneRationalType Real;
typedef AriadneState<> State;
typedef AriadneRectangle<State> Rectangle;
typedef AriadnePolyhedron<State> Poly;
typedef AriadneDenotableSet<Poly> DenotableSet;
typedef AriadnePTree<Rectangle> PTree;
typedef UBLAS::matrix<Real> MatrixT;
typedef UBLAS::vector<Real> VectorT;
typedef AriadnePolyAffineMap< State > PMap;

typedef AriadneAffineMap< PMap > AMap;

typedef PolyAffineIntegrator< State > PInt;
typedef AriadneAffineIntegrator< PInt > AInt;
typedef AriadneIntegrator< AInt > Integrator;

typedef AriadneAffineVectorField< AMap > AVectF;

typedef AriadneDiscreteLocation< AVectF > DiscreteLocation;

typedef AriadneLeavingDiscreteTransition< DiscreteLocation, AMap > LeavingTrans;
typedef AriadneHybridAutomaton< LeavingTrans > Automaton; 

typedef AriadneLocationDenotableSet<DiscreteLocation> LocationDenotableSet;
typedef AriadneHybridDenotableSet<LocationDenotableSet> HybridDenotableSet;

typedef AriadneLocationTrace<DiscreteLocation> LocationTrace;
typedef AriadneTrace<LocationTrace> Trace;

typedef PolyhedronMatlabExporter<State> BSExporter;
typedef DenotableSetExporter<BSExporter> Exporter;

typedef AriadneSolver<Automaton, HybridDenotableSet, Trace , 
	PTree, Integrator, Exporter > Solver;

int main() {

	Real max_flow(5,100);
	Real delta_step(1,1000);
	
	uint dim=4;
	
	/* Define Constants */
	Real a(-3,2), b(955, 10);
	Real a_p=-21, b_p=21840;
	Real a_phi(4,350), b_phi(27,35);
	Real T_m=15;
	Real G=390000;
	Real a_L(-168,100), b_L=105;
	Real k(45,100000); k=k/100000;
	Real k_alpha(259,10000), k_n=6;
	
	Real min_n=740,max_n=860, n_error=10; 
	Real min_theta=0,max_theta=180, theta_error(1,10); 
	Real min_p=5000, max_p=12000, p_error=10;
	Real min_m=k*min_p, max_m=k*max_p, m_error=abs((max_m+min_m)*(delta_step*10));
	Real min_T=0.6*G*min_m, max_T=G*max_m, T_error=abs((max_T+min_T)*(delta_step*10));
	Real min_phi=-15, max_phi=20, phi_error(2,10); 
	Real min_alpha=0, max_alpha=20, alpha_error(1,10);
	
	Real n_0=800;
	Real theta_0=0;
	Real p_0=7300;
	Real m_0=k*p_0;
	Real T_0=G*m_0;
	Real phi_0=20;
	Real alpha0=-(a_p*p_0)/b_p;;

	/* Define Reachability Set Bounding Box */
	State uc(dim),lc(dim), d_rec(dim), init_val(dim);
	
	lc[0]=min_n; 		uc[0]=max_n; 		init_val[0]=n_0;		d_rec[0]=n_error; 
	lc[1]=min_theta; 	uc[1]=max_theta;	init_val[1]=theta_0;		d_rec[1]=theta_error; 
	lc[2]=min_T; 		uc[2]=max_T;		init_val[2]=T_0;		d_rec[2]=T_error; 
	lc[3]=min_phi; 		uc[3]=max_phi;		init_val[3]=phi_0;		d_rec[3]=phi_error; 
	
	/* Define Automaton */
	Automaton H("Model2");
	
	/* Define Variables (Remember previous order) */
	Variable n(0);
	Variable theta(1);
	Variable T(2);
	Variable phi(3);
	Variable last_space_dim(3);
	
	/* Define Vector Field */
	MatrixT A=zero_matrix<Real>(dim);
	
	A(0,0)=a;A(0,2)=b; 		/* dot{n} = a * n + b * T */
	A(1,0)=k_n;			/* dot{theta} = k_n * n */
	
	AVectF vf(A);
	
	/* Define Location S */
	Constraint_System cs_inv1;	
	cs_inv1.insert( denumerator(theta_error) * theta <= denumerator(theta_error) * 180 + numerator(theta_error) );
	cs_inv1.insert( 0*last_space_dim >= -1 );

	Poly inv1(cs_inv1);

	DiscreteLocation l0("S", vf , inv1);
	
	/* Define Location S+ */
	Constraint_System cs_inv2;	
	cs_inv1.insert( denumerator(theta_error) * theta <= denumerator(theta_error) * 180 + numerator(theta_error) );
	cs_inv2.insert( theta >= 160 );
	cs_inv2.insert( 0*last_space_dim >= -1 );

	Poly inv2(cs_inv2);

	DiscreteLocation l1("S+", vf , inv2);

	/* Define activation from S to S+ */
	Constraint_System cs_activ1;	
	cs_activ1.insert( theta >= 160 );
	cs_activ1.insert( 0*last_space_dim >= -1 );
	Poly activ1(cs_activ1);
	
	/* Define Reset from S to S+ */
	MatrixT M1=identity_matrix<Real>(dim);
	VectorT v1=null_vector<Real>(dim);
	
	M1(3,1)=-1; M1(3,3)=0; v1(3)=180;	/* phi= 180 - theta */
	
	AMap reset1(M1,v1);
	
	/* Add to H the arc from S to S+ */
	H.add_arc(l0, l1, activ1, reset1);

	/* Define activation from S+ to S */
	Constraint_System cs_activ2;	
	cs_activ2.insert( theta >= 180 );
	cs_activ2.insert( 0*last_space_dim >= -1 );
	Poly activ2(cs_activ2);
	
	/* Define Reset from S to S+ */
	MatrixT M2=identity_matrix<Real>(dim);
	VectorT v2=null_vector<Real>(dim);
	
	v2(1)=-180;				/* theta= theta - 180 */
	M2(2,3)=T_m * a_phi; M2(2,2)=0; 
	v2(2)=T_m * b_phi;	/* T = T_m * a_phi * phi + T_m * b_phi */
	
	AMap reset2(M2,v2);
	
	/* Add to H the arc from S to S+ */
	H.add_arc(l1, l0, activ2, reset2);
	
	/* Define activation from S to S */
	Constraint_System cs_activ3;	
	cs_activ3.insert( theta >= 180 );
	cs_activ3.insert( 0*last_space_dim >= -1 );
	Poly activ3(cs_activ3);
	
	/* Define Reset from S to S */
	MatrixT M3=identity_matrix<Real>(dim);
	VectorT v3=null_vector<Real>(dim);
	
	v3(1)=-180;			/* theta= theta - 180 */
	M3(2,2)=0; v3(2)=T_m * b_phi;	/* T = T_m * b_phi */
	
	AMap reset3(M3,v3);
	
	/* Add to H the arc from S to S+ */
	H.add_arc(l0, l0, activ3, reset3);
	
	/* Define initial region */
	Constraint_System cs_init;	
	Real var_coef, known_term,delta_term;
	
	for (uint i=0; i< dim; i++) {
		
		var_coef=denumerator(init_val[i])*denumerator(d_rec[i]);
		known_term=denumerator(d_rec[i]) * numerator(init_val[i]);
		delta_term=denumerator(init_val[i]) * numerator(d_rec[i]);
		
		cs_init.insert( var_coef * Variable(i) >= known_term - delta_term );
		cs_init.insert( var_coef * Variable(i) <= known_term + delta_term );
	}

	Poly bs_init(cs_init);
	
	LocationDenotableSet init(l0, bs_init);
	
	Solver solver;
	
	Rectangle bounding_box(lc,uc);

	std::vector<uint> dims;
	
	dims.resize(3);
	
	dims[0]=0;
	dims[1]=1;
	dims[2]=2;
	
	Exporter exporter(H.name(),bounding_box,dims);
	//PTree ptree(bounding_box,9);	
	
	try {
		solver.reach(H, init, max_flow,
				15, delta_step, exporter);
		
	} catch (std::exception &e) {
		
		std::cout << __LINE__ << ":" << __FILE__ <<  " "<< e.what()<< std::endl;
	}
	
	/* Export the automaton */
	
	dot_print(H);

	return 0;
}
