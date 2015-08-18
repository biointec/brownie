/*-
 * Copyright 1997, 1998-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *	$Id: fib.h,v 1.10 2003/01/14 10:11:30 jmg Exp $
 *
 */

#ifndef _FIB_H_
#define _FIB_H_

#ifdef _MSC_VER
        #include "msstdint.h"     /* uint64_t types, etc.  */
#else
        #include <stdint.h>
#endif

typedef float Time;    /* "Time", as defined by Daniel Zerbino */
typedef int32_t NodeID; /* NodeID definition */
typedef NodeID Data;    /* "Data", in this case a nodeID */

struct fibheap;
struct fibheap_el;
typedef int (*voidcmp)(Data, Data);

/* functions for key heaps */
struct fibheap *fh_makekeyheap(void);
struct fibheap_el *fh_insertkey(struct fibheap *, Time, Data);
Time fh_minkey(struct fibheap *);
Time fh_replacekey(struct fibheap *, struct fibheap_el *, Time);
Data fh_replacekeydata(struct fibheap *, struct fibheap_el *, Time, Data);

/* functions for void * heaps */
struct fibheap *fh_makeheap(void);
voidcmp fh_setcmp(struct fibheap *, voidcmp);
Data fh_setneginf(struct fibheap *, Data);
struct fibheap_el *fh_insert(struct fibheap *, Data);

/* shared functions */
Data fh_extractmin(struct fibheap *);
Data fh_min(struct fibheap *);
Data fh_replacedata(struct fibheap *, struct fibheap_el *, Data);
Data fh_delete(struct fibheap *, struct fibheap_el *);
void fh_deleteheap(struct fibheap *);
struct fibheap *fh_union(struct fibheap *, struct fibheap *);

#ifdef FH_STATS
int fh_maxn(struct fibheap *);
int fh_ninserts(struct fibheap *);
int fh_nextracts(struct fibheap *);
#endif

#endif /* _FIB_H_ */
