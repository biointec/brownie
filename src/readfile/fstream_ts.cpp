/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  <copyright holder> <email>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#include "fstream_ts.h"
#include <string.h> // memcpy
#include <iostream>
#include <assert.h>

ofstream_ts& ofstream_ts::write ( const char* s , std::streamsize n ) {
	omp_set_lock(&write_lock);
	if((buffer_size - current) < n) {
		ofs.write(buffer, current);
		current = 0;
	}
	memcpy(&buffer[current], s, n);
	current += n;
	
// 	// old way:
// 	ofs.write(s, n);
	omp_unset_lock(&write_lock);
	return *this;
}

ifstream_ts& ifstream_ts::read ( char* s , std::streamsize n ) {
	omp_set_lock(&read_lock);
	readn(s, n);
	omp_unset_lock(&read_lock);
	return *this;
}

ifstream_ts& ifstream_ts::read_c ( char* s , std::streamsize n ) {
	omp_set_lock(&read_lock);
	readn(s, n);
	return *this;
}

ifstream_ts& ifstream_ts::read_s ( char* s , std::streamsize n ) {
	readn(s, n);
	omp_unset_lock(&read_lock);
	return *this;
}

void ifstream_ts::readn ( char* s , std::streamsize n ) {
	if((buffer_size - current) < n) { // read more in, save last part to first
		memcpy(buffer, &buffer[current], buffer_size - current);
		current = buffer_size - current;
		ifs.read(&buffer[current], buffer_size - current);
		buffer_size = current + ifs.gcount(); // if buffer isnt full adjust size!
		current = 0;
	}
	memcpy(s, &buffer[current], n);
	current += n;
	
// // 	//Old way:
// 	ifs.read(s, n);
	
}