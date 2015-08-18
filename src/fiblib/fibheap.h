/***************************************************************************
 *   Copyright (C) 2010 Jan Fostier (jan.fostier@intec.ugent.be)           *
 *   Original Velvet code by Daniel Zerbino (zerbino@ebi.ac.uk)            *
 *                                                                         *
 *   This file is part of Velvet 2.0                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef FIBHEAP_H
#define FIBHEAP_H

extern "C" {
        #include "fib.h"
}

// ============================================================================
// TYPEDEFS
// ============================================================================

typedef fibheap_el FibHeapEl;

// ============================================================================
// FIBONACCI HEAP WRAPPER
// ============================================================================

template <class T>
class FibHeap {

private:
        struct fibheap *fh;     // actual fibonacci heap

public:
        /**
         * Default constructor
         */
        FibHeap() {
                fh = fh_makekeyheap();
        }

        /**
         * Destuctor
         */
        ~FibHeap() {
                fh_deleteheap(fh);
        }

        /**
         * Insert a < key, value > pair
         * @param key Key
         * @param value Value
         * @return A pointer to the fibonacci heap element
         */
        FibHeapEl* insert(Time key, T value) {
                return fh_insertkey(fh, key, value);
        }

        /**
         * Return the minimal key in the heap
         * @return The minimal key in the heap
         */
        Time getMinKey() {
                return fh_minkey(fh);
        }

        /**
         * Removes the node with the smallest key and returns the value
         * @return The value of the node with the smallest key
         */
        T extractMinKey() {
                return fh_extractmin(fh);
        }

        /**
         * Replace the value in a given node with a new value
         * @param el Fibonacci heap element
         * @param T new value
         */
        void replaceValue(FibHeapEl* el, T newValue) {
                fh_replacedata(fh, el, newValue);
        }

        /**
         * Replace the key in a given node with a new value
         * @param el Fibonacci heap element
         * @param key new key
         */
        void replaceKey(FibHeapEl* el, Time newKey) {
                fh_replacekey(fh, el, newKey);
        }
};

#endif
