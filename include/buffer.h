/***************************************************************************
 *            buffer.h
 *
 *  Mon May  3 11:57:51 2004
 *  Copyright  2004  Alberto Casagrande
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
 
#ifndef _BUFFER_H
#define _BUFFER_H

namespace Ariadne {
	
	
class Buffer{
		void *buf;
		unsigned int index;
		unsigned int buffer_dim;
	public:
		Buffer();
		Buffer(unsigned int new_buf_dim);
		Buffer(void *new_buf, unsigned int new_buf_dim);
		~Buffer();
		void set_index(unsigned int idx);
		unsigned int get_index();
		/* \brief Extracts a vector of byte from the buffer_dim
		 *
		 * \param byte_num is the wanted dimention of the output vector.
		 * After the execution of this method \a byte_num contains the
		 * real dimention of the output vector.
		 * \return A new \a void-type vector of lenght equal to the
		 * minimum between \a byte_num and \a buffer_num-index. The 
		 * returned vector contains exactly what is stored into the 
		 * buffer from \a index to the minimum between 
		 * \a index+byte_num and \a	buffer_dim.
		 */
		void *read_buffer(unsigned int *byte_num);
	
		/* \brief Returns a pointer to the buffer part pointed byte
		 * \a index.
		 */
		void *read_from_buffer();
};


}

#endif /* _BUFFER_H */
