/***************************************************************************
 *            test_expansion.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include <iostream>
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "expansion.h"

#include "test.h"
using namespace std;
using namespace Ariadne;


class TestExpansion
{
    typedef MultiIndex MI;
    typedef Expansion<Float> E;
  public:
    void test();
  private:
    void test_concept();
    void test_iterator_concept();
    void test_data_access();
    void test_cleanup();
    void test_constructors();
    void test_indexing();
    void test_find();
};


void TestExpansion::test()
{
    ARIADNE_TEST_CALL(test_data_access());
    ARIADNE_TEST_CALL(test_cleanup());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_indexing());
    ARIADNE_TEST_CALL(test_find());
}


void TestExpansion::test_concept()
{
    Float x=0;
    MI a(3);
    Expansion<Float> e(3);
    const Expansion<Float> ce(3);

    e=Expansion<Float>();
    e=Expansion<Float>(3);
    e=Expansion<Float>(ce);

    e=Expansion<Float>(3,1, 0.0, 0.0,0.0,0.0);
    e=Expansion<Float>(3,1, 1,2,3,5.0);

    //e=Expansion<Float>::variable(3u,0u);
    //e=Expansion<Float>::variables(3u)[0u];

    e=x;

    e.reserve(2u);
    e.insert(a,x);
    e.prepend(a,x);
    e.append(a,x);
    e.append(a,a,x);
    e.clear();

    x=ce[a];
    e[a]=1.0;

    ce.number_of_nonzeros();
    ce.argument_size();

    e.erase(e.begin());

    ce.check();
    e.cleanup();
}

void TestExpansion::test_iterator_concept()
{
    MI a(3);
    Expansion<Float> e(3);
    const Expansion<Float> cp(3);

    Expansion<Float>::iterator iter=e.begin(); iter=e.end(); iter=e.find(a);
    Expansion<Float>::const_iterator citer=e.begin(); citer=e.end(); citer=e.find(a);
    citer=e.begin(); citer=cp.end(); citer=cp.find(a);

    Expansion<Float>::value_type val=*iter;
    Expansion<Float>::reference ref=*iter;
    //Expansion<Float>::pointer ptr=iter.operator->();
    Expansion<Float>::const_reference ncref=*iter;

    // WARNING: Cannot convert non-constant pointer to constant pointer
    //Expansion<Float>::const_pointer ncptr=iter.operator->();

    Expansion<Float>::value_type cval=*citer;
    Expansion<Float>::const_reference cref=*citer;
    //Expansion<Float>::const_pointer cptr=citer.operator->();

    ++iter; --iter;
    ++citer; --citer;

    iter==iter; iter!=iter; citer==citer; citer!=citer;
    citer==iter; citer!=iter; iter==citer; iter!=citer;

}

// Test dereferencing of iterators
void TestExpansion::test_data_access()
{
    Expansion<Float> e(3);
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 0,0,0),2.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 1,0,0),3.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 0,1,0),5.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 1,0,1),7.0));

    // Test iterator difference
    ARIADNE_TEST_PRINT(e.begin());
    ARIADNE_TEST_PRINT(e.end());
    ARIADNE_TEST_EQUAL(e.begin()+e.number_of_nonzeros(),e.end());
    ARIADNE_TEST_EQUAL(e.end()-e.number_of_nonzeros(),e.begin());
    ARIADNE_TEST_EQUAL(e.end()-e.begin(),e.number_of_nonzeros());

    // Test derefencing of iterators
    ARIADNE_TEST_PRINT(e.begin());
    ARIADNE_TEST_PRINT(*e.begin());
    ARIADNE_TEST_EQUAL(e.begin()->key(),MI(3, 0,0,0));
    ARIADNE_TEST_EQUAL(e.begin()->data(),2.0);

    // Test hand-coded swap of values
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>::iterator,iter1,(e.begin()+1));
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>::iterator,iter2,(e.begin()+3));

    // Perform swap
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>::value_type,tmp,(*iter2));
    ARIADNE_TEST_ASSERT(tmp.key()==MI(3,1,0,1));
    ARIADNE_TEST_ASSERT(tmp.data()==7.0);
    ARIADNE_TEST_EXECUTE(*iter2=*iter1);
    ARIADNE_TEST_ASSERT(iter2->key()==MI(3,1,0,0));
    ARIADNE_TEST_ASSERT(iter2->data()==3.0);
    ARIADNE_TEST_EXECUTE(*iter1=tmp);
    ARIADNE_TEST_ASSERT(iter1->key()==MI(3,1,0,1));
    ARIADNE_TEST_ASSERT(iter1->data()==7.0);

}

void TestExpansion::test_cleanup()
{
    // Test to see if the cleanup/sort operations work.
    // Since these are used in the constructors, we can't use the main constructors to test this
    MI a(3);
    MI b(3); ++b;

    Expansion<Float> e(3);
    for(uint i=0; i!=2; ++i) {
        if(i%2) { e.append(a,1/(1.+i)); ++b; ++b; a=b; ++b; } else { e.append(b,1/(1.+i));}
    }

    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.sort());
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.remove_zeros());
    ARIADNE_TEST_PRINT(e);

}

void TestExpansion::test_constructors()
{
    // Empty polynomial
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,e1,(3));
    // Dense polynomial
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,e2,(3,2, 0.,0.,0.,0., 5.,2.,0.,0.,3.,0.));
    ARIADNE_TEST_EQUAL(e2[MI(3, 2,0,0)],5.0);
    // Sparse polynomial with unordered indiced
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,p3,(2,5, 1,2,5.0, 0,0,2.0, 1,0,3.0, 3,0,7.0, 0,1,11.0));
    ARIADNE_TEST_EQUAL(p3[MI(2, 1,2)],5.0);
    ARIADNE_TEST_EQUAL(p3[MI(2, 0,0)],2.0);

    // Unordered indices
    ARIADNE_TEST_EQUAL(E(2,5, 1,2,5.0, 0,0,2.0, 1,0,3.0, 3,0,7.0, 0,1,11.0),E(2,5, 0,0,2.0, 1,0,3.0, 0,1,11.0, 3,0,7.0, 1,2,5.0));
    // Repeated indices; do not sum in expansion class
    ARIADNE_TEST_ASSERT(E(2,3, 1,0,2.0, 0,2,7.0, 1,0,3.0)==E(2,3, 1,0,2.0, 1,0,3.0, 0,2,7.0)
        || E(2,3, 1,0,2.0, 0,2,7.0, 1,0,3.0)==E(2,3, 1,0,3.0, 1,0,2.0, 0,2,7.0));

    // Regression tests
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,pr2,(3,0, 0.0));
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,pr1,(3,1, 3,0,0, 0.0));
}




void TestExpansion::test_indexing()
{
    Expansion<Float> e(3,4, 0,0,0,2.0,  1,0,0,3.0, 1,0,1,5.0, 2,1,0,7.0);
    const Expansion<Float>& pc=e;
    ARIADNE_TEST_EQUAL(e[MI(3, 1,0,0)],3.0);

    e[MI(3, 1,0,0)]-=0.5;
    ARIADNE_TEST_EQUAL(e[MI(3, 1,0,0)],2.5);

    e[MI(3, 1,1,0)]=11.0;
    ARIADNE_TEST_EQUAL(e[MI(3, 1,1,0)],11.0);

    e.clear();
    e[MI(3, 0,0,0)]=2.0;
    e[MI(3, 0,1,0)]=3.0;
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EQUAL(e.number_of_nonzeros(),2);
    ARIADNE_TEST_EQUAL(e[MI(3, 0,0,0)],2.0);
    ARIADNE_TEST_EQUAL(e[MI(3, 0,1,0)],3.0);
    ARIADNE_TEST_EXECUTE(e[MI(3, 1,0,0)]=5.0);
    ARIADNE_TEST_EQUAL(e[MI(3, 1,0,0)],5.0);
    ARIADNE_TEST_EQUAL(e[MI(3, 0,1,0)],3.0);
    ARIADNE_TEST_EXECUTE(e[MI(3, 0,0,1)]=7.0);
    ARIADNE_TEST_EQUAL(e[MI(3, 0,0,1)],7.0);
    
    // Test insert at beginning
    e.clear();
    e[MI(3, 0,1,0)]=2.0;
    e[MI(3, 0,0,1)]=3.0;
    e[MI(3, 2,1,0)]=5.0;
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EQUAL(pc[MI(3, 0,0,0)],0.0);
    ARIADNE_TEST_EQUAL(e.number_of_nonzeros(),3);
    ARIADNE_TEST_EXECUTE(e[MI(3, 0,0,0)]=7);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EQUAL(e.number_of_nonzeros(),4);
    Expansion<Float>::const_iterator iter=e.begin();
    ARIADNE_TEST_EQUAL(iter->key(),MI(3, 0,0,0));
    ARIADNE_TEST_EQUAL(iter->data(),7.0);
    ++iter;
    if(iter->key()==MI(3, 2,1,0)) {
        std::cerr<<"Error in Expansion<X>::insert(MultiIndex,X) probably due to incorrect compiler optimization. Please inform the developers."; }
    ARIADNE_TEST_EQUAL(iter->key(),MI(3, 0,1,0));
    ARIADNE_TEST_EQUAL(iter->data(),2.0);
    ++iter;
    ARIADNE_TEST_EQUAL(iter->key(),MI(3, 0,0,1));
    ARIADNE_TEST_EQUAL(iter->data(),3.0);
    ++iter;
    ARIADNE_TEST_EQUAL(iter->key(),MI(3, 2,1,0));
    ARIADNE_TEST_EQUAL(iter->data(),5.0);



}

void TestExpansion::test_find()
{
    Expansion<Float> e(2,5, 1,2,5.0, 0,0,2.0, 1,0,3.0, 3,0,7.0, 0,1,11.0);
    MI a(2);
    a[0]=1; a[1]=2;
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_PRINT(e.find(a)-e.begin());
    ARIADNE_TEST_COMPARE(e.find(a),!=,e.end());
    ARIADNE_TEST_EQUAL(e.find(a)->key(),a);
    ARIADNE_TEST_EQUAL(e.find(a)->data(),5.0);
    a[1]=1;
    ARIADNE_TEST_EQUAL(e.find(a),e.end());
}

int main() {
    TestExpansion().test();
    return ARIADNE_TEST_FAILURES;
}
