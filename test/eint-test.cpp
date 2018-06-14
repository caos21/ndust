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

#define BOOST_TEST_MODULE eint_test

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

#include "../include/array.h"
#include "../include/eint.h"
#include "../include/constants.h"

using namespace std;
using namespace boost::unit_test;
using namespace boost::test_tools;

using namespace eint;


//******************************************************************************
//**                              Test constants                              **
//******************************************************************************
BOOST_AUTO_TEST_CASE( constants ) {
  {
    BOOST_TEST_MESSAGE( "[ii] Testing electron charge..." );
    double true_val = eCharge;
    double test_val = 1.6021766208e-19;
    BOOST_TEST_MESSAGE( "[ii] e = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing Coulomb electron constant..." );
    double true_val = 8987551787.368177;
    double test_val = Kcoul;
    BOOST_TEST_MESSAGE( "[ii] K = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }  
}

//******************************************************************************
//**                           Test microparticles                            **
//******************************************************************************
BOOST_AUTO_TEST_CASE( microparticles ) {
  // parameters
  double r1o = 20e-6;
  double r2o = 10e-6;
  double q1o = 1e4*eCharge;
  double q2o = 2e4*eCharge;
  
  double r21 = r2o/r1o;
  double rt = 1 + r21;
  double q21 = q2o/q1o;
  double eps = 20.0;
  
  {
    BOOST_TEST_MESSAGE( "[ii] Testing r21..." );
    double true_val = 0.5;
    double test_val = r21;
    BOOST_TEST_MESSAGE( "[ii] r21 = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing rt..." );
    double true_val = 1.5;
    double test_val = rt;
    BOOST_TEST_MESSAGE( "[ii] r21 = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing q21..." );
    double true_val = 2.0;
    double test_val = q21;
    BOOST_TEST_MESSAGE( "[ii] q21 = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing potential_ipa_fact..." );
    double true_val = 6355969389.021359;
    double test_val = potential_ipa_fact(rt, r21, q21, eps);
    BOOST_TEST_MESSAGE( "[ii] Pipa = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }
  
  {
    BOOST_TEST_MESSAGE( "[ii] Testing force_ipa_fact..." );
    double true_val = -12923123549.208027;
    double test_val = force_ipa_fact(rt, r21, q21, eps);
    BOOST_TEST_MESSAGE( "[ii] Fipa = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing potential_mpc_funct..." );
    double true_val = 663443100.19583511;
    potential_mpc_funct pmpcfunct(r21, q21, eps, 50);
    double test_val = pmpcfunct(rt);
    BOOST_TEST_MESSAGE( "[ii] Pmpc = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
    true_val = 8.5151924232299881e-17;
    test_val *= q1o*q1o/r1o;
    BOOST_TEST_MESSAGE( "[ii] Potential_mpc = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);    
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing force_mpc_funct..." );
    double true_val = -80053707397.578506;
    force_mpc_funct fmpcfunct(r21, q21, eps, 40);
    double test_val = fmpcfunct(rt);
    BOOST_TEST_MESSAGE( "[ii] Fmpc = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
    true_val = -5.1373864803335422e-10;
    test_val *= q1o*q1o/(r1o*r1o);
    BOOST_TEST_MESSAGE( "[ii] Force_mpc = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
    
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing potential and force..." );
  
    darray rtarray = linear<double>(rt, 4.0*rt, 100);

    potential_mpc_funct pmpcfunct(r21, q21, eps, 40);

    potential_ipa_funct pipafunct(r21, q21, eps);    
    
    double potential_prefactor = q1o*q1o/r1o;
    
    vector<double> potential_mpc;
    vector<double> potential_ipa;

    for (auto ra: rtarray) {
      double pmpc = potential_prefactor * pmpcfunct(ra);
      double pipa = potential_prefactor * pipafunct(ra);
      potential_mpc.push_back(pmpc);
      potential_ipa.push_back(pipa);
    }
   
    vector<double> force_frompotmpc(99);
    vector<double> force_frompotipa(99);

    for (unsigned int i=0; i<99; ++i) {
      double dr = r1o*(rtarray[i+1]-rtarray[i]);
      force_frompotipa[i] = -1e12*(potential_ipa[i+1]-potential_ipa[i])/dr;
      force_frompotmpc[i] = -1e12*(potential_mpc[i+1]-potential_mpc[i])/dr;
    }
    
    force_mpc_funct fmpcfunct(r21, q21, eps, 40);
    force_ipa_funct fipafunct(r21, q21, eps);    
   
    vector<double> force_mpc;
    vector<double> force_ipa;

    double force_prefactor = 1e12*q1o*q1o/(r1o*r1o);
    for (auto ra: rtarray) {
      double fmpc = force_prefactor * fmpcfunct(ra);
      double fipa = force_prefactor * fipafunct(ra);      
      force_mpc.push_back(fmpc);
      force_ipa.push_back(fipa);
    }

    double ffipa[] =  {-70.40855498, -49.93330369, -34.75092439, -23.34306769,
		       -14.68081809,  -8.04956907,  -2.94242429,   1.0068095 ,
		       4.06693486,   6.43805338,   8.27098398,   9.68066072,
		       10.75553651,  11.56428607,  12.16065081,  12.58698517,
		       12.87688258,  13.05714044,  13.14924479,  13.17050199,
		       13.13490858,  13.05382481,  12.93649991,  12.79048445,
		       12.6219559 ,  12.43597725,  12.23670344,  12.02754699,
		       11.81131149,  11.59029967,  11.36640116,  11.14116418,
		       10.91585411,  10.69150168,  10.46894261,  10.24885039,
		       10.0317634 ,   9.81810743,   9.60821446,   9.40233827,
		       9.20066748,   9.0033365 ,   8.81043468,   8.622014  ,
		       8.43809551,   8.25867476,   8.08372637,   7.91320789,
		       7.74706301,   7.5852243 ,   7.42761554,   7.27415357,
		       7.12475002,   6.97931257,   6.8377462 ,   6.69995406,
		       6.56583832,   6.43530082,   6.30824362,   6.18456947,
		       6.06418216,   5.94698684,   5.83289025,   5.72180094,
		       5.61362939,   5.50828813,   5.40569184,   5.30575741,
		       5.20840396,   5.11355287,   5.02112778,   4.93105459,
		       4.84326143,   4.75767865,   4.67423876,   4.59287643,
		       4.51352841,   4.43613348,   4.36063245,   4.28696804,
		       4.21508489,   4.14492946,   4.07644999,   4.00959646,
		       3.94432053,   3.88057546,   3.8183161 ,   3.75749881,
		       3.69808141,   3.64002314,   3.5832846 ,   3.52782774,
		       3.47361573,   3.42061303,   3.36878524,   3.31809912,
		       3.26852254,   3.22002442,   3.1725747};

    //BOOST_TEST(ffipa == force_frompotipa, per_element() );
    for (unsigned int i=0; i < 99; ++i) {
      BOOST_TEST_CONTEXT("index " << i)
	{
	  BOOST_REQUIRE_CLOSE(ffipa[i], force_frompotipa[i], PTOL);
	}
    }
    
    
    double ffpmpc[] = { -6.78290394e+02,  -1.56344112e+02,  -7.64389141e+01,
			-4.26104968e+01,  -2.41794187e+01,  -1.28087725e+01,
			-5.26411327e+00,  -2.21937736e-02,   3.73028956e+00,
			6.46800328e+00,   8.48707216e+00,   9.98203321e+00,
			1.10860284e+01,   1.18932316e+01,   1.24720611e+01,
			1.28733250e+01,   1.31354424e+01,   1.32879026e+01,
			1.33536250e+01,   1.33506091e+01,   1.32931145e+01,
			1.31925204e+01,   1.30579628e+01,   1.28968134e+01,
			1.27150446e+01,   1.25175109e+01,   1.23081676e+01,
			1.20902428e+01,   1.18663736e+01,   1.16387146e+01,
			1.14090252e+01,   1.11787395e+01,   1.09490238e+01,
			1.07208229e+01,   1.04948982e+01,   1.02718590e+01,
			1.00521883e+01,   9.83626421e+00,   9.62437740e+00,
			9.41674612e+00,   9.21352828e+00,   9.01483176e+00,
			8.82072296e+00,   8.63123400e+00,   8.44636871e+00,
			8.26610770e+00,   8.09041262e+00,   7.91922973e+00,
			7.75249290e+00,   7.59012612e+00,   7.43204568e+00,
			7.27816191e+00,   7.12838067e+00,   6.98260466e+00,
			6.84073442e+00,   6.70266926e+00,   6.56830793e+00,
			6.43754928e+00,   6.31029271e+00,   6.18643863e+00,
			6.06588878e+00,   5.94854648e+00,   5.83431684e+00,
			5.72310696e+00,   5.61482606e+00,   5.50938552e+00,
			5.40669903e+00,   5.30668255e+00,   5.20925441e+00,
			5.11433526e+00,   5.02184811e+00,   4.93171827e+00,
			4.84387337e+00,   4.75824329e+00,   4.67476013e+00,
			4.59335817e+00,   4.51397383e+00,   4.43654560e+00,
			4.36101400e+00,   4.28732152e+00,   4.21541257e+00,
			4.14523341e+00,   4.07673210e+00,   4.00985846e+00,
			3.94456399e+00,   3.88080183e+00,   3.81852669e+00,
			3.75769483e+00,   3.69826397e+00,   3.64019326e+00,
			3.58344321e+00,   3.52797568e+00,   3.47375381e+00,
			3.42074196e+00,   3.36890569e+00,   3.31821170e+00,
			3.26862781e+00,   3.22012290e+00,   3.17266689e+00};

    //BOOST_TEST(ffpmpc == force_frompotmpc, per_element() );
    //BOOST_TEST(ffpmpc == force_frompotmpc);
    for (unsigned int i=0; i < 99; ++i) {
      BOOST_TEST_CONTEXT("index " << i)
	{
	  BOOST_REQUIRE_CLOSE(ffpmpc[i], force_frompotmpc[i], PTOL);
	}
    }
    
    double fmpc [] = {-513.73864803, -119.62859791,  -64.42540596,  -38.53108894,
		      -23.18272831,  -13.076627  ,   -6.0242446 ,   -0.92673472,
		      2.83958008,    5.65969312,    7.78612941,    9.39222037,
		      10.60115644,   11.50290964,   12.16462957,   12.63731554,
		      12.96025725,   13.16408175,   13.27289861,   13.30584343,
		      13.27820944,   13.20229056,   13.08801817,   12.94344802,
		      12.77513646,   12.58843391,   12.38771563,   12.17656463,
		      11.9579175 ,   11.73418144,   11.50732866,   11.27897295,
		      11.05043193,   10.82277808,   10.5968805 ,   10.37343938,
		      10.15301439,    9.93604831,    9.72288651,    9.51379331,
		      9.3089655 ,    9.10854374,    8.91262204,    8.72125573,
		      8.53446814,    8.35225626,    8.17459539,    8.00144319,
		      7.83274299,    7.66842663,    7.50841679,    7.35262906,
		      7.20097353,    7.05335627,    6.90968048,    6.76984749,
		      6.63375757,    6.50131064,    6.37240682,    6.24694691,
		      6.12483276,    6.0059676 ,    5.8902563 ,    5.77760553,
		      5.66792396,    5.56112235,    5.45711363,    5.35581302,
		      5.25713799,    5.16100834,    5.06734617,    4.97607588,
		      4.88712418,    4.80042002,    4.71589457,    4.6334812 ,
		      4.5531154 ,    4.47473478,    4.39827897,    4.32368958,
		      4.25091018,    4.17988618,    4.11056486,    4.04289522,
		      3.97682799,    3.91231558,    3.84931197,    3.7877727 ,
		      3.72765481,    3.66891678,    3.61151851,    3.55542121,
		      3.50058742,    3.44698092,    3.39456671,    3.34331096,
		      3.29318095,    3.24414506,    3.1961727 ,    3.14923431};

    //    BOOST_TEST(fmpc == force_mpc, per_element() );
    for (unsigned int i=0; i < 100; ++i) {
      BOOST_TEST_CONTEXT("index " << i)
	{
	  BOOST_REQUIRE_CLOSE(fmpc[i], force_mpc[i], PTOL);
	}
    }    
    
    ofstream outfile("force_test.dat");
    outfile << "#r\tforce_mpc\tforce_ipa";    
    
    for (unsigned int i=0; i<99; ++i) {
      outfile << '\n' << rtarray[i]*r1o
	      << '\t' << force_mpc[i] << '\t' << force_ipa[i]
	      << '\t' << force_frompotmpc[i] << '\t' << force_frompotipa[i];
    }
    outfile << '\n';
    outfile.close();
  }
}



//******************************************************************************
//**                           Test nanoparticles                             **
//******************************************************************************
BOOST_AUTO_TEST_CASE( nanoparticles ) {

  using namespace eint;
  
  double r2o = 0.5e-9;
  double r1o = 1.5e-9;
  double q2o = -2.0*eCharge;
  double q1o = eCharge;

  // WARNING
  // CRITICAL
  // r2o > r1o increase stability and fast convergence
  
  double r21 = r2o/r1o;
  double rt = 1 + r21;
  double q21 = q2o/q1o;
  double eps = 11.68;

  {
    BOOST_TEST_MESSAGE( "[ii] Testing r21..." );
    double true_val = 0.333333333333;
    double test_val = r21;
    BOOST_TEST_MESSAGE( "[ii] r21 = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing rt..." );
    double true_val = 1.333333333333;
    double test_val = rt;
    BOOST_TEST_MESSAGE( "[ii] r21 = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing q21..." );
    double true_val = -2.0;
    double test_val = q21;
    BOOST_TEST_MESSAGE( "[ii] r21 = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing potential_ipa_fact..." );
    double true_val = -23674188438.15027;
    double test_val = potential_ipa_fact(rt, r21, q21, eps);
    BOOST_TEST_MESSAGE( "[ii] Fipa = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }
  
  {
    BOOST_TEST_MESSAGE( "[ii] Testing force_ipa_fact..." );
    double true_val = -60267048378.50898;
    double test_val = force_ipa_fact(rt, r21, q21, eps);
    BOOST_TEST_MESSAGE( "[ii] Fipa = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing potential_mpc_funct..." );
    double true_val = -34282694806.041397;
    potential_mpc_funct pmpcfunct(r21, q21, eps, 50);
    double test_val = pmpcfunct(rt);
    BOOST_TEST_MESSAGE( "[ii] Pmpc = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
    true_val = -5.8668430992628157e-19;
    test_val *= q1o*q1o/r1o;
    BOOST_TEST_MESSAGE( "[ii] Potential_mpc = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);    
  }

  {
    BOOST_TEST_MESSAGE( "[ii] Testing force_mpc_funct..." );
    double true_val = -190887137462.55887;
    force_mpc_funct fmpcfunct(r21, q21, eps, 50);
    double test_val = fmpcfunct(rt);
    BOOST_TEST_MESSAGE( "[ii] Fmpc = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
    true_val = -2.1777846257346395e-09;
    test_val *= q1o*q1o/(r1o*r1o);
    BOOST_TEST_MESSAGE( "[ii] Force_mpc = " << true_val << " ? " << test_val);
    BOOST_REQUIRE_CLOSE(true_val, test_val, PTOL);
    
  }
}


