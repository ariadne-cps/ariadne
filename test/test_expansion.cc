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
#include <vector>
#include "numeric.h"
#include "expansion.h"

#include "test.h"
using namespace std;
using namespace Ariadne;

template<class K, class V> struct MapValue {
    typedef K key_type;
    K first; V second;
    MapValue(const K& k, const V& v) : first(k), second(v) { }
    const K& key() const { return first; }
    const V& data() const { return second; }
};


class TestExpansion
{
    typedef MultiIndex MI;
    typedef Expansion<Float> E;
  public:
    void test();
  private:
    void test_working();
    void test_concept();
    void test_iterator_concept();
    void test_data_access();
    void test_cleanup();
    void test_constructors();
    void test_indexing();
    void test_find();
    void test_embed();
};


void TestExpansion::test()
{
    ARIADNE_TEST_CALL(test_working());
    ARIADNE_TEST_CALL(test_data_access());
    ARIADNE_TEST_CALL(test_cleanup());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_indexing());
    ARIADNE_TEST_CALL(test_find());
    ARIADNE_TEST_CALL(test_embed());
}


void TestExpansion::test_working()
{
    Expansion<Float> e(3);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_PRINT(e.size());
    // Append values
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 0,0,0),2.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 1,0,0),3.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 0,1,0),5.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 1,0,1),7.0));

    //ARIADNE_TEST_EQUAL(e.find(MI(3, 0,0,0)),e.begin());
    //assert(false);
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
    Expansion<Float>::pointer ptr=iter.operator->();
    Expansion<Float>::const_reference ncref=*iter;

    //WARNING: Cannot convert non-constant pointer to constant pointer
    Expansion<Float>::const_pointer ncptr=iter.operator->();
    ncptr=citer.operator->();

    Expansion<Float>::value_type cval=*citer;
    Expansion<Float>::const_reference cref=*citer;
    Expansion<Float>::const_pointer cptr=citer.operator->();

    ++iter; --iter;
    ++citer; --citer;

    iter==iter; iter!=iter; citer==citer; citer!=citer;
    citer==iter; citer!=iter; iter==citer; iter!=citer;

    ref=cref; cptr=ptr; ref=ncref; ncptr=ptr;
}

// Test dereferencing of iterators
void TestExpansion::test_data_access()
{
    // Test index and element size values for a large expansion
    Expansion<Float> e5(5);
    ARIADNE_TEST_EQUAL(e5._index_size(),2);
    ARIADNE_TEST_EQUAL(e5._element_size(),4);

    Expansion<Float> e(3);

    // Test index and element size values
    ARIADNE_TEST_EQUAL(e._index_size(),1);
    ARIADNE_TEST_EQUAL(e._element_size(),3);

    // Append values
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 0,0,0),2.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 1,0,0),3.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 0,1,0),5.0));
    ARIADNE_TEST_EXECUTE(e.append(MI(3, 1,0,1),7.0));

    // Test iterator difference
    ARIADNE_TEST_PRINT(e.begin());
    ARIADNE_TEST_PRINT(e.end());
    ARIADNE_TEST_EQUAL(e.begin()+e.number_of_nonzeros(),e.end());
    ARIADNE_TEST_EQUAL(e.end()-e.number_of_nonzeros(),e.begin());
    ARIADNE_TEST_EQUAL(e.end()-e.begin(),Expansion<Float>::difference_type(e.number_of_nonzeros()));

    // Test derefencing of iterators
    ARIADNE_TEST_PRINT(e.begin());
    ARIADNE_TEST_PRINT(*e.begin());
    ARIADNE_TEST_EQUAL(e.begin()->key(),MI(3, 0,0,0));
    ARIADNE_TEST_EQUAL(e.begin()->data(),2.0);



    // The behaviour of iterators is rather odd and not what might be expected
    // A MultiIndex reference assigned to by iter->key() changes its value
    // when the iterator is incremented, but a Float reference does not.
    // This behaviour should be changed in future versions if technologically
    // feasible.
    Expansion<Float>::iterator iter=e.begin();
    const MultiIndex& aref=iter->key();
    const Float& xref=iter->data();
    Float x1=iter->data();
    ++iter;
    MultiIndex a2=iter->key();
    ARIADNE_TEST_ASSERT(a2==aref);
    ARIADNE_TEST_ASSERT(x1==xref);

    // Test finding lower bound of multi index
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(iter=e.lower_bound(MI(3,1,0,0)));
    ARIADNE_TEST_PRINT(iter);
    ARIADNE_TEST_EQUAL(iter->key(),MI(3,1,0,0));
    ARIADNE_TEST_EQUAL(iter->data(),3.0);
    ARIADNE_TEST_EXECUTE(iter=e.lower_bound(MI(3,1,1,0)));
    ARIADNE_TEST_EQUAL(iter->key(),MI(3,1,0,1));
    ARIADNE_TEST_EQUAL(iter->data(),7.0);


    // Test finding of values of iterators
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(iter=e.find(MI(3,1,0,0)));
    ARIADNE_TEST_EQUAL(iter->key(),MI(3,1,0,0));
    ARIADNE_TEST_EQUAL(iter->data(),3.0);



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
    ARIADNE_TEST_PRINT(a);
    ARIADNE_TEST_PRINT(b);

    Expansion<Float> e(3);
    ARIADNE_TEST_PRINT(e);
    for(uint i=0; i!=2; ++i) {
        if(i%2) { e.append(a,1/(1.+i)); ++b; ++b; a=b; ++b; } else { e.append(b,1/(1.+i));}
        ARIADNE_TEST_PRINT(e);
    }

    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.sort());
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.remove_zeros());
    ARIADNE_TEST_PRINT(e);

}

void TestExpansion::test_constructors()
{
    // Empty expansion
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,e1,(3));
    // Expansion with all entries; useful for checking ordering of indices
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,e2,(3,4, 1., 2.,3.,4., 5.,6.,7.,8.,9.,10.,
        11.,12.,13.,14.,15.,16.,17.,18.,19.,20., 21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35));

    // Dense expansion
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,e3,(3,2, 0., 0.,0.,0., 5.,2.,0.,0.,3.,0.));
    ARIADNE_TEST_PRINT(e3);
    ARIADNE_TEST_COMPARE(e3.find(MI(3, 2,0,0)),!=,e3.end());
    ARIADNE_TEST_EQUAL(e3[MI(3, 2,0,0)],5.0);

    // Sparse expansion with unordered indiced
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,p3,(2,5, 1,2,5.0, 0,0,2.0, 1,0,3.0, 3,0,7.0, 0,1,11.0));
    ARIADNE_TEST_EQUAL(p3[MI(2, 1,2)],5.0);
    ARIADNE_TEST_EQUAL(p3[MI(2, 0,0)],2.0);

    // Unordered indices
    ARIADNE_TEST_EQUAL(E(2,5, 1,2,5.0, 0,0,2.0, 1,0,3.0, 3,0,7.0, 0,1,11.0),E(2,5, 0,0,2.0, 1,0,3.0, 0,1,11.0, 3,0,7.0, 1,2,5.0));
    // Repeated indices; do not sum in expansion class
    ARIADNE_TEST_ASSERT(E(2,3, 1,0,2.0, 0,2,7.0, 1,0,3.0)==E(2,3, 1,0,2.0, 1,0,3.0, 0,2,7.0)
        || E(2,3, 1,0,2.0, 0,2,7.0, 1,0,3.0)==E(2,3, 1,0,3.0, 1,0,2.0, 0,2,7.0));

    // Regression tests for expansions with only zeroes
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,pr2,(3,0, 0.0));
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,pr1,(3,1, 3,0,0, 0.0));

    // Regression tests for higher-order expansions
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,ho1,(5,4, 0,1,0,0,0,2.0, 0,1,0,0,1,3.0, 2,0,1,0,0,5.0, 0,0,0,0,0,7.0));
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

void TestExpansion::test_embed()
{
    ARIADNE_TEST_CONSTRUCT(Expansion<Float>,e,(2,5, 0,0,2.0, 1,0,3.0, 0,1,11.0, 3,0,7.0, 1,2,5.0));
    ARIADNE_TEST_EQUAL(embed(0,e,2),Expansion<Float>(4,5, 1,2,0,0,5.0, 0,0,0,0,2.0, 1,0,0,0,3.0, 3,0,0,0,7.0, 0,1,0,0,11.0));
    ARIADNE_TEST_EQUAL(embed(1,e,0),Expansion<Float>(3,5, 0,1,2,5.0, 0,0,0,2.0, 0,1,0,3.0, 0,3,0,7.0, 0,0,1,11.0));
    ARIADNE_TEST_EQUAL(embed(1,e,2),Expansion<Float>(5,5, 0,1,2,0,0,5.0, 0,0,0,0,0,2.0, 0,1,0,0,0,3.0, 0,3,0,0,0,7.0, 0,0,1,0,0,11.0));



Expansion<Float> x(5,17,
 0,0,0,0,0,0.0220457,1,0,0,0,0,0.00025196,0,1,0,0,0,0.000142431,0,0,1,0,0,-0.000391323,0,0,0,1,0,0.0483952,0,0,0,0,1,-0.028653,2,0,0,0,0,9.17162e-09,1,1,0,0,0,1.08229e-08,1,0,1,0,0,0.000248011,1,0,0,1,0,3.97819e-06,1,0,0,0,1,-1.9891e-06,0,2,0,0,0,3.19145e-09,0,1,1,0,0,0.000143275,0,1,0,1,0,-3.95255e-06,0,1,0,0,1,1.97627e-06,0,0,2,0,0,-2.31918e-08,0,0,1,1,0,6.27623e-06,0,0,1,0,1,-3.13812e-06,0,0,0,2,0,-0.0036531,0,0,0,1,1,0.0036531,0,0,0,0,2,-0.00269802,3,0,0,0,0,-1.36733e-12,2,1,0,0,0,-2.37889e-13,2,0,1,0,0,1.4827e-10,2,0,0,1,0,-2.96541e-10,2,0,0,0,1,1.4827e-10,1,2,0,0,0,1.12021e-12,1,1,1,0,0,-6.16868e-11,1,1,0,1,0,1.23374e-10,1,1,0,0,1,-6.16868e-11,1,0,2,0,0,2.94003e-08,1,0,1,1,0,-3.97819e-06,1,0,1,0,1,1.9891e-06,0,3,0,0,0,4.62908e-13,0,2,1,0,0,-8.50758e-11,0,2,0,1,0,1.70152e-10,0,2,0,0,1,-8.50758e-11,0,1,2,0,0,1.6979e-08,0,1,1,1,0,-2.29745e-06,0,1,1,0,1,1.14873e-06,4,0,0,0,0,2.54808e-16,3,1,0,0,0,-4.9344e-17,3,0,1,0,0,-2.21046e-14,3,0,0,1,0,4.42092e-14,3,0,0,0,1,-2.21046e-14,2,2,0,0,0,-2.24043e-16,2,1,1,0,0,1.37947e-14,2,1,0,1,0,-2.75894e-14,2,1,0,0,1,1.37947e-14,2,0,2,0,0,-9.31989e-09,2,0,1,1,0,2.96541e-10,2,0,1,0,1,-1.4827e-10,1,3,0,0,0,-9.94029e-18,1,2,1,0,0,1.07703e-14,1,2,0,1,0,-2.15406e-14,1,2,0,0,1,1.07703e-14,1,1,2,0,0,-1.07612e-08,1,1,1,1,0,-1.23374e-10,1,1,1,0,1,6.16868e-11,0,4,0,0,0,3.13236e-17,0,3,1,0,0,-2.63841e-15,0,3,0,1,0,5.27682e-15,0,3,0,0,1,-2.63841e-15,0,2,2,0,0,-3.10637e-09,0,2,1,1,0,-1.70152e-10,0,2,1,0,1,8.50758e-11,4,0,1,0,0,4.11928e-18,4,0,0,1,0,-8.23857e-18,4,0,0,0,1,4.11928e-18,3,1,1,0,0,-3.4276e-18,3,1,0,1,0,6.85519e-18,3,1,0,0,1,-3.4276e-18,3,0,2,0,0,1.38944e-12,3,0,1,1,0,-4.42092e-14,3,0,1,0,1,2.21046e-14,2,2,1,0,0,-1.9807e-18,2,2,0,1,0,3.9614e-18,2,2,0,0,1,-1.9807e-18,2,1,2,0,0,2.24094e-13,2,1,1,1,0,2.75894e-14,2,1,1,0,1,-1.37947e-14,1,3,1,0,0,1.1207e-18,1,3,0,1,0,-2.2414e-18,1,3,0,0,1,1.1207e-18,1,2,2,0,0,-1.13098e-12,1,2,1,1,0,2.15406e-14,1,2,1,0,1,-1.07703e-14,0,4,1,0,0,-0,0,3,2,0,0,-4.6027e-13,0,3,1,1,0,-5.27682e-15,0,3,1,0,1,2.63841e-15,4,0,2,0,0,-2.58927e-16,4,0,1,1,0,8.23857e-18,4,0,1,0,1,-4.11928e-18,3,1,2,0,0,5.27712e-17,3,1,1,1,0,-6.85519e-18,3,1,1,0,1,3.4276e-18,2,2,2,0,0,2.26023e-16,2,2,1,1,0,-3.9614e-18,2,2,1,0,1,1.9807e-18,1,3,2,0,0,8.81961e-18,1,3,1,1,0,2.2414e-18,1,3,1,0,1,-1.1207e-18,0,4,2,0,0,-3.13208e-17);

ARIADNE_TEST_PRINT(x);
ARIADNE_TEST_EXECUTE(embed(0,x,1));


}

int main() {
    TestExpansion().test();
    return ARIADNE_TEST_FAILURES;
}
