/***************************************************************************
 *            io_error.h
 *
 *  Mon May  3 11:44:56 2004
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
 
#ifndef _IO_ERROR_H
#define _IO_ERROR_H

namespace Ariadne {
	
typedef enum {
	NO_READ_ERROR,
	PARSE_ERROR,
	LEXICAL_ERROR,
	NO_FILE_ERROR,
	READ_INTERNAL_ERROR
} _Read_Error_Type;
	
	
typedef struct {
		unsigned int line;
		unsigned int from_byte;
		_Read_Error_Type type;
} Read_Error;

typedef enum {
	NO_WRITE_ERROR,
	NO_SPACE_ERROR,
	NO_DIR_ERROR,
	WRITE_INTERNAL_ERROR
} _Write_Error_Type;

typedef struct {
		_Write_Error_Type type;
} Write_Error;


}

#endif /* _IO_ERROR_H */
