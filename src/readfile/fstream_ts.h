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


#ifndef FSTREAM_TS_H
#define FSTREAM_TS_H

#include <fstream>
#include <stdint.h>
#include <omp.h>
#include <iostream>

#define BUFFERSIZE 1048576

/**
 * Wrapper class for ofstream to make a threadsafe output filestream
 */
class ofstream_ts
{

public:
	ofstream_ts (int buffer_size_ = BUFFERSIZE) : ofs() {
		omp_init_lock(&write_lock);
		init_buffer(buffer_size_);
	}
	
	explicit ofstream_ts ( const char * filename, std::ios_base::openmode mode = std::ios_base::out, int buffer_size_ = BUFFERSIZE) : 
	ofs(filename, mode) {
		omp_init_lock(&write_lock);
		init_buffer(buffer_size_);
	}
	
	void open( const char * filename, std::ios_base::openmode mode = std::ios_base::out ) {
		ofs.open(filename, mode);
	}
	
	void close() {
		ofs.write(buffer, current);
		ofs.close();
	}
	
	bool is_open() {
		ofs.is_open();
	}
	
    ~ofstream_ts() {
		omp_destroy_lock(&write_lock);
		delete[] buffer;
	}
	
	/**
	 * Thread safe write operation
	 * stores data in buffer untill buffer is full, then writes all data at once
	 */
	ofstream_ts& write ( const char* s , std::streamsize n );
	
private:
	char *buffer;
	size_t buffer_size, current;
	omp_lock_t write_lock;
	std::ofstream ofs;
	
	void init_buffer(int buffer_size_) {
		buffer_size = buffer_size_;
		buffer = new char[buffer_size];
		current = 0;	
	}
};

/**
 * Wrapper class for ofstream to make a threadsafe output filestream
 */
class ifstream_ts
{

public:
	ifstream_ts (int buffer_size_ = BUFFERSIZE) : ifs() {
		omp_init_lock(&read_lock);
		init_buffer(buffer_size_);
	}
	
	explicit ifstream_ts ( const char * filename, std::ios_base::openmode mode = std::ios_base::in, int buffer_size_ = BUFFERSIZE) : 
	ifs(filename, mode) {
		omp_init_lock(&read_lock);
		init_buffer(buffer_size_);
	}
	
	void open( const char * filename, std::ios_base::openmode mode = std::ios_base::in ) {
		ifs.open(filename, mode);
	}
	
	void close() {
		ifs.close();
	}
	
	bool is_open() {
		ifs.is_open();
	}
	
	bool good() {
		return ifs.good() || (!ifs.good() && current < buffer_size);
	}
	
    ~ifstream_ts() {
		omp_destroy_lock(&read_lock);
		delete[] buffer;
	}
	
	/**
	 * Thread safe read operation
	 */
	ifstream_ts& read ( char* s , std::streamsize n );
	
	/**
	 * Thread safe read operation
	 * sets the lock but doesnt unset it yet cause more must be read sequentially by same thread
	 */
	ifstream_ts& read_c ( char* s , std::streamsize n );
	
	/**
	 * Thread safe read operation
	 * unsets the lock, assumes a read_c was done before read_s
	 */
	ifstream_ts& read_s ( char* s , std::streamsize n );
	
	/**
	 * If no more data can be ready, this function should be called to close the lock!
	 */
	void unset_readlock() {
		omp_unset_lock(&read_lock);
	}
private:
	char *buffer;
	size_t buffer_size, current;
	omp_lock_t read_lock;
	std::ifstream ifs;
	
	void init_buffer(int buffer_size_) {
		buffer_size = buffer_size_;
		buffer = new char[buffer_size];
		current = buffer_size;
	}
	
	/**
	 * Help function for the Thread safe read operation
	 */
	void readn ( char* s , std::streamsize n );
};

#endif // FSTREAM_TS_H
