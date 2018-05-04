/*
 * Copyright 2018 <Benjamin Santos> <caos21@gmail.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define BOOST_TEST_MODULE enhancement_test

#define PTOL 1e-6

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <functional>
#include <cstdlib>

#include "../include/enhancement.h"
// #include "../include/eint.h"
// #include "../include/constants.h"

using namespace std;
using namespace boost::unit_test;
using namespace boost::test_tools;

using namespace eint;
using namespace enhancement;


//******************************************************************************
//**                              Test enhancement                            **
//******************************************************************************
BOOST_AUTO_TEST_CASE( enhancement_test ) {
  {
    darray rarray = {1, 10};
    darray qarray = {1, 20};
    double eps = 10.0;
    
    
    Enhancement enh(rarray, qarray, eps);

    enh.compute_reducedpairs();

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing particle pairs..." );
    for (auto pp: enh.particle_pairs) {
      pp.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing reduced particle pairs..." );
    for (auto rpp: enh.reduced_pairs) {
      rpp.print();
    }

    BOOST_TEST_MESSAGE( "\n\n[ii] Testing contact potentials..." );
    enh.compute_contact_potentials();
    for (auto cpot: enh.contact_potentials) {
      cpot.print();
    }
    
  }

}
