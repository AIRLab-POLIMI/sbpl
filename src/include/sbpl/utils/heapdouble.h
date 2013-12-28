/*
 * Copyright (c) 2008, Maxim Likhachev
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Pennsylvania nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __HEAPDOUBLE_H_
#define __HEAPDOUBLE_H_

#include <sbpl/planners/planner.h>
#include <sbpl/utils/keydouble.h>

//the maximum size of the heap
#define HEAPSIZE 20000000 
#define HEAPSIZE_INIT 5000
#define HEAPDOUBLE_PERFORMANCE_TEST GLOBAL_PERFORMANCE_TEST

struct HEAPELEMENT
{
    AbstractSearchState *heapstate;
    CKeyDouble key;
};

typedef struct HEAPELEMENT heapelement;

class CHeapDouble
{
    //data
public:
    int percolates; //for counting purposes
    heapelement* heap;
    int currentsize;
    int allocated;

    //constructors
public:
    CHeapDouble();
    ~CHeapDouble();

    //functions
public:
    bool emptyheap();
    bool fullheap();
    bool inheap(AbstractSearchState *AbstractSearchState);
    CKeyDouble getkeyheap(AbstractSearchState *AbstractSearchState);
    void makeemptyheap();
    void insertheap(AbstractSearchState *AbstractSearchState, CKeyDouble key);
    void deleteheap(AbstractSearchState *AbstractSearchState);
    void updateheap(AbstractSearchState *AbstractSearchState, CKeyDouble NewKey);
    AbstractSearchState *getminheap();
    AbstractSearchState *getminheap(CKeyDouble& ReturnKey);
    CKeyDouble getminkeyheap();
    AbstractSearchState *deleteminheap();
    void makeheap();
    void insert_unsafe(AbstractSearchState* state, CKeyDouble key);
    void updateheap_unsafe(AbstractSearchState* AbstractSearchState, CKeyDouble NewKey);
    void deleteheap_unsafe(AbstractSearchState* AbstractSearchState);

private:
    void percolatedown(int hole, heapelement tmp);
    void percolateup(int hole, heapelement tmp);
    void percolateupordown(int hole, heapelement tmp);

    void growheap();
    void sizecheck();
};

#endif

