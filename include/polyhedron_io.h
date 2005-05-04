/***************************************************************************
 *            polyhedron_io.h
 *
 *  Wed Feb 16 19:37:07 2005
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
 
#ifndef _POLYHEDRON_IO_H
#define _POLYHEDRON_IO_H

#include <matlabexp.h>
#include <polyhedron.h>
#include <denotable_set.h>

inline void _print_polyhedron_(std::ostream &os,
				const Parma_Polyhedra_Library::NNC_Polyhedron &p) {
	
	using namespace Parma_Polyhedra_Library::IO_Operators;
	
	os << p;
}

namespace Ariadne {	
namespace Geometry {

namespace IO_Operators{	

template <typename S>
std::ostream& operator<<(std::ostream &os, const Polyhedron<S> &r) {
	
	_print_polyhedron_(os,r._poly);
	
	return os;
}

template <typename BSE>
class PolyhedronSetExporter {
	
	public:
		typedef BSE Exporter;
		typedef typename Exporter::BasicSet BasicSet;
		typedef typename BasicSet::State State;
		typedef typename BasicSet::Real Real;
	    	typedef DenotableSet<BasicSet> DenotableSet;
		typedef Rectangle<State> Rectangle;
	
	private:
		
		Exporter _e;
	
	public:
	
		PolyhedronSetExporter() {}
			
		PolyhedronSetExporter(const Exporter &e): _e(e) {}
			
		PolyhedronSetExporter(const std::string &file_name, 
			const Rectangle &bbox): _e(file_name,bbox) {}
			
		PolyhedronSetExporter(const std::string &file_name, 
				const Rectangle &bbox, std::vector<uint> dims): 
				_e(file_name, bbox, dims) {}
			 
		PolyhedronSetExporter(std::ofstream os, Rectangle bbox): _e(os,bbox) {}
			
		PolyhedronSetExporter(const std::string &file_name): _e(bbox) {}
			
		PolyhedronSetExporter(const std::string &file_name,
					std::vector<uint> dims): 
				_e(file_name, dims) {}
					
		PolyhedronSetExporter(std::ofstream os): _e(os) {}
			
		void export_denotableset(const DenotableSet &A) {
			
			BasicSet set=convex_hull(A[0], A[A.size()-1]);
		
			(this->_e).export_basicset(set);
		}
					
		void export_denotableset(const DenotableSet &A, std::vector<uint> dims) {
		
			BasicSet set=convex_hull(A[0], A[A.size()-1]);
		
			(this->_e).export_basicset(set,dims);
		}
};	

template <typename S>
class PolyhedronMatlabExporter {
		
	public:
		typedef S State;
		typedef typename S::Real Real;
		typedef Polyhedron<State> BasicSet;
		typedef Rectangle<State> Rectangle;
	
	private:
		
		bool _os_given;
	
		std::ofstream _os;
	
		std::string _f_name;
	
		bool _with_dims;
	
		std::vector<uint> _dims;
	
		uint _polyhedron_n;
	
		void _open_stream_if_not_opened() {
		
			if (_os_given) return;
			
			_os_given=true;
			
			open_matlab_stream(this->_f_name,this->_os);
		}
	
	public:
		PolyhedronMatlabExporter(): _os_given(false), _with_dims(false), 
					_polyhedron_n(0) {}
		
		PolyhedronMatlabExporter(const std::string &file_name, 
						const Rectangle &bbox): _os_given(false), 
				 _f_name(file_name), _with_dims(false), _polyhedron_n(0) {}
					
		PolyhedronMatlabExporter(const std::string &file_name, 
					 const Rectangle &bbox, const std::vector<uint> &dims): 
				_os_given(false), 
				_f_name(file_name), _with_dims(true),
				_dims(dims), _polyhedron_n(0) {}
					
		PolyhedronMatlabExporter(const std::string &file_name, 
					const std::vector<uint> &dims): 
				_os_given(false), 
				_f_name(file_name), _with_dims(true),
				_dims(dims), _polyhedron_n(0) {}
					
		PolyhedronMatlabExporter(std::ofstream os, 
					const Rectangle &bbox): _os_given(true), 
				_os(os), _with_dims(false), _polyhedron_n(0){}
		
		PolyhedronMatlabExporter(std::ofstream os, const Rectangle &bbox, 
					const std::vector<uint> &dims): _os_given(true), 
				_os(os), _with_dims(true), _dims(dims){}
			
		~PolyhedronMatlabExporter() {
			if (_os_given) {
				
				if (this->_polyhedron_n>0) {
					
					this->_os << "hold on" <<std::endl
							<< "options.edgecolor='b'" <<std::endl;	
				
					for (size_t i=0; i< this->_polyhedron_n; i++) {
					
						this->_os << "plot(P"<< i<< ", options)" <<std::endl;
					}
				}
				
				close_mathlab_stream(this->_os);	
			}
		}
		
		void export_basicset(const BasicSet &p) {
			
			if (!this->_with_dims) {
				
				if (p.dim()!=3) 
					throw std::invalid_argument("No dimensions to project");
			
				std::vector<uint> dims;

				dims.resize(3);

				for (size_t i=0; i< 3; i++) {
					dims[i]=i;
				}
							
			
				this->export_basicset(p,dims);
				return;	
			}
			
			this->export_basicset(p,_dims);
		}
		
		void export_basicset(const BasicSet &p,
			const std::vector<uint> &dims) {
				
			this->_open_stream_if_not_opened();
	
			Real known_term, error(1,100);	

			BasicSet proj_p=p;	
			proj_p=proj_p.project_on_dimentions(dims);
			proj_p.set_precision_to_upperapproximating_for_output(error);
			
			Parma_Polyhedra_Library::Constraint_System::const_iterator j_cs, begin, end;

			std::vector<Real> coeff;
				
			coeff.resize(dims.size());
				
			begin=((proj_p._poly).constraints()).begin();
			end=((proj_p._poly).constraints()).end();
	
			this->_os << "C"<< this->_polyhedron_n << " = [ " ;
				
			bool first=true;
				
			for (j_cs=begin; j_cs!=end; j_cs++) {
				
				if (!first) {
					this->_os << "; ";
				} 
				
				bool not_zero=false;
				
				for (size_t i=0; i< coeff.size(); i++) {
					coeff[i]=j_cs->coefficient(Parma_Polyhedra_Library::Variable(dims[i]));
					
					not_zero=not_zero || (coeff[i]!=0);
				}
				
				if (not_zero) {
				
					first=false;
					
					for (size_t i=0; i< coeff.size(); i++) {
						this->_os << -coeff[i] << " ";
					}
				}
					
			}
	
			this->_os << " ]" <<std::endl;
			
			this->_os << "d"<< this->_polyhedron_n << " = [ ";
			
			first=true;
			
			for (j_cs=begin; j_cs!=end; j_cs++) {
		
				if (!first) {
					this->_os << " ; ";
				} 
				
				known_term=j_cs->inhomogeneous_term();
		
				bool not_zero=false;
				
				for (size_t i=0; i< coeff.size(); i++) {
					coeff[i]=j_cs->coefficient(Parma_Polyhedra_Library::Variable(dims[i]));
					
					not_zero=not_zero || (coeff[i]!=0);
				}
				
				if (not_zero) {
				
					first=false;
				
					this->_os << known_term;
				}
					
			}
	
			this->_os << " ]" <<std::endl;
			
			this->_os << "P" << this->_polyhedron_n << 
				" = polytope( C"<< this->_polyhedron_n <<  
				" , d"<< this->_polyhedron_n << " )"<<std::endl;
			
			this->_polyhedron_n++;
		}

		inline std::string &file_name() { return this->_f_name; }
		
		inline std::ofstream &stream() { return this->_os; }
		
		inline bool &os_given() { return this->_os_given; }
};

template <typename S>
class PolyhedronVTKOctavePlotExporter {
		
	public:
		typedef S State;
		typedef typename S::Real Real;
		typedef Polyhedron<State> BasicSet;
		typedef Rectangle<State> Rectangle;
	
	private:
		
		bool _os_given;
	
		std::ofstream _os;
	
		std::string _f_name;
	
		bool _with_dims;
	
		std::vector<uint> _dims;
	
		uint _polyhedron_n;
	
		void _open_stream_if_not_opened() {
		
			if (_os_given) return;
			
			_os_given=true;
			
			open_matlab_stream(this->_f_name,this->_os);
		}
	
	public:
		PolyhedronVTKOctavePlotExporter(): _os_given(false), _with_dims(false), 
					_polyhedron_n(0) {}
		
		PolyhedronVTKOctavePlotExporter(const std::string &file_name, 
						const Rectangle &bbox): _os_given(false), 
				 _f_name(file_name), _with_dims(false), _polyhedron_n(0) {}
					
		PolyhedronVTKOctavePlotExporter(const std::string &file_name, 
					 const Rectangle &bbox, const std::vector<uint> &dims): 
				_os_given(false), 
				_f_name(file_name), _with_dims(true),
				_dims(dims), _polyhedron_n(0) {}
					
		PolyhedronVTKOctavePlotExporter(const std::string &file_name, 
					const std::vector<uint> &dims): 
				_os_given(false), 
				_f_name(file_name), _with_dims(true),
				_dims(dims), _polyhedron_n(0) {}
					
		PolyhedronVTKOctavePlotExporter(std::ofstream os, 
					const Rectangle &bbox): _os_given(true), 
				_os(os), _with_dims(false), _polyhedron_n(0){}
		
		PolyhedronVTKOctavePlotExporter(std::ofstream os, const Rectangle &bbox, 
					const std::vector<uint> &dims): _os_given(true), 
				_os(os), _with_dims(true), _dims(dims){}
			
		~PolyhedronVTKOctavePlotExporter() {
			if (_os_given) {
				close_mathlab_stream(this->_os);	
			}
		}
		
		void export_basicset(const BasicSet &p) {
			
			if (!this->_with_dims) {
				
				if (p.dim()!=3) 
					throw std::invalid_argument("No dimensions to project");
			
				std::vector<uint> dims;

				dims.resize(3);

				for (size_t i=0; i< 3; i++) {
					dims[i]=i;
				}
							
			
				this->export_basicset(p,dims);
				return;	
			}
			
			this->export_basicset(p,_dims);
		}
		
		void export_basicset(const BasicSet &p,
			const std::vector<uint> &dims) {
				
			this->_open_stream_if_not_opened();
	
			Real known_term, error(1,100);	

			BasicSet proj_p=p;	
			proj_p=proj_p.project_on_dimentions(dims);
			proj_p.set_precision_to_upperapproximating_for_output(error);
			
			Parma_Polyhedra_Library::Constraint_System::const_iterator j_cs, begin, end;

			std::vector<Real> coeff;
				
			coeff.resize(dims.size());
				
			begin=((proj_p._poly).constraints()).begin();
			end=((proj_p._poly).constraints()).end();
	
			this->_os << "C"<< this->_polyhedron_n << " = [ " ;
				
			bool first=true;
				
			for (j_cs=begin; j_cs!=end; j_cs++) {
				
				if (!first) {
					this->_os << "; ";
				} 
				
				bool not_zero=false;
				
				for (size_t i=0; i< coeff.size(); i++) {
					coeff[i]=j_cs->coefficient(Parma_Polyhedra_Library::Variable(dims[i]));
					
					not_zero=not_zero || (coeff[i]!=0);
				}
				
				if (not_zero) {
				
					first=false;
					
					for (size_t i=0; i< coeff.size(); i++) {
						this->_os << -coeff[i] << " ";
					}
				}
					
			}
	
			this->_os << " ]" <<std::endl;
			
			this->_os << "d"<< this->_polyhedron_n << " = [ ";
			
			first=true;
			
			for (j_cs=begin; j_cs!=end; j_cs++) {
		
				if (!first) {
					this->_os << " ; ";
				} 
				
				known_term=j_cs->inhomogeneous_term();
		
				bool not_zero=false;
				
				for (size_t i=0; i< coeff.size(); i++) {
					coeff[i]=j_cs->coefficient(Parma_Polyhedra_Library::Variable(dims[i]));
					
					not_zero=not_zero || (coeff[i]!=0);
				}
				
				if (not_zero) {
				
					first=false;
				
					this->_os << known_term;
				}
					
			}
	
			this->_os << " ]" <<std::endl;
			
			this->_os << "P" << this->_polyhedron_n << 
				" = polytope( C"<< this->_polyhedron_n <<  
				" , d"<< this->_polyhedron_n << " )"<<std::endl;
			
			this->_polyhedron_n++;
		}

		inline std::string &file_name() { return this->_f_name; }
		
		inline std::ofstream &stream() { return this->_os; }
		
		inline bool &os_given() { return this->_os_given; }
};


template <typename S>
class PolyhedronOneDimMatlabExporter {
		
	public:
		typedef S State;
		typedef typename S::Real Real;
		typedef Polyhedron<State> BasicSet;
		typedef Rectangle<State> Rectangle;
	
	private:
		
		bool _os_given;
	
		std::ofstream _os;
	
		std::string _f_name;
	
		bool _with_dims;
	
		std::vector<uint> _dims;
	
		uint _data_line;
		
		bool _stream_clean;
	
		void _open_stream_if_not_opened() {
		
			if (_os_given) return;
			
			this->_os_given=true;
			this->_stream_clean=true;
			
			open_OneDim_matlab_stream(this->_f_name,this->_os);
		}
	
	public:
		PolyhedronOneDimMatlabExporter(): _os_given(false), _with_dims(false), 
					_data_line(0) {}
		
		PolyhedronOneDimMatlabExporter(const std::string &file_name, 
						const Rectangle &bbox): _os_given(false), 
				 _f_name(file_name), _with_dims(false), _data_line(0) {}
					
		PolyhedronOneDimMatlabExporter(const std::string &file_name, 
					 const Rectangle &bbox, const std::vector<uint> &dims): 
				_os_given(false), 
				_f_name(file_name), _with_dims(true),
				_dims(dims), _data_line(0) {
					
			if (dims.size()!=1) {
				throw "The PolyhedronOneDimMatlabExporter export on one dimension only.";
			}				
		}
					
		PolyhedronOneDimMatlabExporter(const std::string &file_name, 
					const std::vector<uint> &dims): 
				_os_given(false), 
				_f_name(file_name), _with_dims(true),
				_dims(dims), _data_line(0) {
					
			if (dims.size()!=1) {
				throw "The PolyhedronOneDimMatlabExporter export on one dimension only.";
			}		
		}
					
		PolyhedronOneDimMatlabExporter(std::ofstream os, 
					const Rectangle &bbox): _os_given(true), 
				_os(os), _with_dims(false), _data_line(0){}
		
		PolyhedronOneDimMatlabExporter(std::ofstream os, const Rectangle &bbox, 
					const std::vector<uint> &dims): _os_given(true), 
				_os(os), _with_dims(true), _dims(dims),_data_line(0){
					
			if (dims.size()!=1) {
				throw "The PolyhedronOneDimMatlabExporter export on one dimension only.";
			}				
		}
			
		~PolyhedronOneDimMatlabExporter() {
			if (_os_given) {
				close_OneDim_mathlab_stream(this->_os);	
			}
		}
		
		void export_basicset(const BasicSet &p) {
			
			if (!this->_with_dims) {
				
				if (p.dim()!=1) 
					throw std::invalid_argument("No dimensions to project");
			
				std::vector<uint> dims;

				dims.resize(1);

				dims[0]=0;
			
				this->export_basicset(p,dims);
				return;	
			}
			
			this->export_basicset(p,_dims);
		}
		
		void export_basicset(const BasicSet &p,
			const std::vector<uint> &dims) {
				
			this->_open_stream_if_not_opened();

			if (!this->_stream_clean) {
				this->_os << " ; ";
			}
			
			BasicSet proj_p=p.project_on_dimentions(dims);
			
			Parma_Polyhedra_Library::Constraint_System::const_iterator j_cs, begin, end;
	
			Real coeff, known_term, error(1,1000);	
				
			proj_p.set_precision_to_upperapproximating(error);
				
			begin=((proj_p._poly).constraints()).begin();
			end=((proj_p._poly).constraints()).end();
	
			this->_os << this->_data_line << ", ";
			
			bool first=true;
					
			for (j_cs=begin; j_cs!=end; j_cs++) {
				
				coeff=j_cs->coefficient(Parma_Polyhedra_Library::Variable(dims[0]));
				known_term=j_cs->inhomogeneous_term();
				
				if (coeff!=0) {
				
					if (!first) {
						this->_os << ", ";
					} 
					
					first=false;
					
					this->_os << -(mpf_class)(known_term/coeff);
				}
					
			}
			
			this->_stream_clean=false;
			
			this->_data_line++;
		}

		inline std::string &file_name() { return this->_f_name; }
		
		inline std::ofstream &stream() { return this->_os; }
		
		inline bool &os_given() { return this->_os_given; }
};

}
}
}

#endif /* _POLYHEDRON_IO_H */
