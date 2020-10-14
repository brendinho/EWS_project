/*
	Brendon Phillips
	PhD Candidate
	Bauch computational epidemiology research group
	Department of Applied Mathematics
	Faculty of Mathematics
	Universiity of Waterloo
*/

#ifndef PARAMETER_VALUES_HEADER_
#define PARAMETER_VALUES_HEADER_

#define NDEBUG

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <exception>
#include <experimental/filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <stdexcept>
#include <sys/times.h>
#include <unordered_set>
#include <unistd.h>
#include <utility>
#include <vector>
#include "prettyprint.hpp"
#include "stdafx.h"
#include "Snap.h"
#include "Timing.hpp"
#include "unistd.h"

namespace fs = std::experimental::filesystem;

#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wignored-attributes"

std::random_device h_rd;
std::mt19937_64 h_Mersenne(h_rd());
std::uniform_real_distribution<double> h_float(0,1);

int h_mod(int a,int b){ int temp = a%b; return (temp>=0)?temp:temp+b; }

int mod(int a,int b){ int temp = a%b; return (temp>=0)?temp:temp+b; }

int h_seed = getpid();

const float h_fastfloat(void)
{
	return h_float(h_Mersenne);
}

const int h_fastrand(int i) // i is the ceiling or the generated integer
{ //  a random number generator by Intel - using the bitshift operator

	// h_seed = (214013*h_seed+2531011);
	// return ((h_seed>>16)&0x7FFF)/32767.0*i;
	return h_mod( (int)(h_fastfloat()*i), i);
}

template <typename Type> void remove_element(std::vector<Type>& vec, Type target)
{ // loads faster than the find-erase idiom
	typename std::vector<Type>::iterator it = std::find(vec.begin(), vec.end(), target);
	if (it != vec.end()) { std::swap(*it, vec.back()); vec.pop_back(); }
	return;
}

template <typename Type> void random_insert(std::vector<Type>& vec, Type target)
{
	int offset = h_mod( h_fastrand(vec.size()), vec.size() );
	vec.insert( std::next(vec.begin(), offset), target );
	return;
}

template <typename Iter> Iter random_member(Iter head, Iter tail)
{ // returns a random element of a container
	int upper_bound = std::distance(head, tail), iterator_offset = 0;
	float temp = h_fastfloat()*upper_bound;
	int temp2 = (h_fastfloat()<0.5) ? floor(temp) : ceil(temp) ;
	iterator_offset = (temp2==upper_bound) ? 0 : temp2;
	std::advance(head, iterator_offset);
	return head;
}

template <typename Iter> double sum(Iter head, Iter tail)
{ // get the mean of the elements in a container
	return std::accumulate(head, tail, 0.0);
}

template <typename Iter> double mean(Iter head, Iter tail)
{ // get the mean of the elements in a container
	return sum(head, tail) / (float) distance(head, tail);
}

template <typename Iter> float unbiased_sd(const Iter head, const Iter tail)
{ // get the unbiased estimator for the standard deviation of the elements in a container
	const float temp_average = mean(head, tail);
	const float N = distance(head,tail);
	float squared = 0.0;
	for(Iter pos = head; pos < tail; ++pos){ squared += (*pos)*(*pos); }
	return sqrt( (squared - temp_average*temp_average*N) / (float)(N-1)  );
}

const std::vector<int> extract(std::map<int, int> const Container, const std::string vals)
{ // adapted from the site http://www.lonecpluspluscoder.com/2015/08/13/an-elegant-way-to-extract-keys-from-a-c-map/
	if( Container.empty() ){ return std::vector<int>(); }
	if( std::set<std::string>({"keys", "values"}).count(vals) == 0 ){ std::cerr << "improper type of values requested for extraction. options are 'keys' and 'values'" << std::endl; return std::vector<int>(); }
  	std::vector<int> the_vals;
 	if(vals == "keys") 			for(std::pair<int, int> element : Container){ the_vals.push_back( element.first ); }
	else if(vals == "values")	for(std::pair<int, int> element : Container){ the_vals.push_back( element.second ); }
  	return the_vals;
}

template<typename Type> std::vector<Type> range(Type length_from_zero)
{
	std::vector<Type> output(length_from_zero);
	std::iota(output.begin(), output.end(), 0);
	return output;
}

template<typename Type> std::vector<Type> range(Type start_val, Type max_val_plus_one)
{
	std::vector<Type> output(max_val_plus_one - start_val);
	std::iota(output.begin(), output.end(), start_val);
	return output;
}

/****************************** REPORTING CONSTANTS *******************************/


const int	global_infection_duration = 2; // -1 for unbounded infection

const int 	global_random_opinion_switch_interval = 1;

const int	global_case_importation_interval = 1;
const int	global_case_importation_delay = 2; // measured in weeks

const int	global_replenishment_interval = 1;

const float global_rewiring_probability = 0.001;

const int	global_neighbour_sampling_interval = 1;
const int 	global_correlation_origin = 0;

const int	global_infection_latency = 0;

const int	global_beta_parameter = 1.;

const int 	global_stored_history_length = 500;
const float	global_equilibrium_threshold = 0.005;

const int 	global_time_limit = 40000;
const int	global_time_minimum = 500;

const std::string global_physical_process = "SIRV\0";
const std::string global_social_process = "NHV\0";

const std::set<std::string> global_check_these_metrics =
{
	"num_S", "num_I", "num_R", "num_N", "num_H",
	"joint_dist_VI", "joint_dist_VR",
	"joint_dist_NS", "joint_dist_NI", "joint_dist_NR",
	"joint_dist_HS", "joint_dist_HI", "joint_dist_HR",
	"HH_joins", "NN_joins", "VV_joins", "HN_joins", "HV_joins", "NV_joins"
};

const std::set<std::string> global_legal_layers = {"physical", "social", "\0"};
const std::map<std::string, std::string> global_legal_states =
{
	{"physical", global_physical_process},
	{"social", global_social_process}
};

const bool check_layer(const std::string layer_name, const char state = '\0', const bool verbose = true)
{ // global function to check if layer names and intra-layer state names are correct

 	if( global_legal_layers.find(layer_name) == global_legal_layers.end() )
	{
		if(verbose){ std::cerr << "illegal layer name " << layer_name << std::endl; }
		exit(EXIT_FAILURE);
	}
 	if(state != '\0')
	{
		if( global_legal_states.at(layer_name).find(state) == std::string::npos)
		{
			if(verbose){ std::cerr << "illegal state " << state << " requested from layer " << layer_name << std::endl; }
			exit(EXIT_FAILURE);
		}
	}
 	return true;
}

const int num_files_matching(const std::string search_location, const std::string substring)
{
	// a non-recursive function to return the number of files in the requested location containing the substring given
	int count = 0;
	for(const auto & file : fs::directory_iterator(search_location))
	{
		if( file.path().string().find(substring) != std::string::npos )
		{
			count++;
		}
	}
	return count;
};

const int max_int(const std::vector<int> x)
{
	if( x.empty() ){ return std::numeric_limits<int>::quiet_NaN(); }
	return *std::max_element(std::cbegin(x), std::cend(x));
}
const int min_int(const std::vector<int> x)
{
	if( x.empty() ){ return std::numeric_limits<int>::quiet_NaN(); }
	return *std::min_element(std::cbegin(x), std::cend(x));
}
const int sum_int(const std::vector<int> x)
{
	if( x.empty() ){ return std::numeric_limits<int>::quiet_NaN(); }
	return std::accumulate(std::cbegin(x), std::cend(x), 0);
};
const double avg_int(const std::vector<int> x)
{
	if( x.empty() ){ return std::numeric_limits<int>::quiet_NaN(); }
	return std::accumulate(std::cbegin(x), std::cend(x), 0.0) / (double) x.size();
};

template <typename Type> std::vector<Type> set_to_vector(std::set<Type> input_set)
{
	std::set<Type> temp = input_set;
	std::vector<Type> output_vec;
	std::copy(temp.begin(), temp.end(), std::back_inserter(output_vec));
	return output_vec;
}

template <typename Type> std::set<Type> vector_to_set(std::set<Type> input_vec)
{
	std::vector<Type> temp = input_vec;
	std::set<Type> output_set;
	std::copy(temp.begin(), temp.end(), std::inserter(output_set, output_set.begin()));
	return output_set;
}

TVec<TInt> to_TVec(std::set<int> Container, int num = -1)
{
	std::set<int> temp = Container;
	TVec<TInt> output_vec;
	int stop = (num == -1) ? temp.size() : num;
	for(int i = 0; i < stop; ++i) output_vec.Add( *std::next(temp.begin(), i) );
	return output_vec;
}

TVec<TInt> to_TVec(std::vector<int> Container, int num = -1)
{
	std::vector<int> temp = Container;
	TVec<TInt> output_vec;
	int stop = (num == -1) ? temp.size() : num;
	for(int i = 0; i < stop; ++i) output_vec.Add(temp[i]);
	return output_vec;
}

float see_if_number(std::string input)
{
	// https://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
	// praises be to Ben Voigt
	// returns nan if input is not a number, so outpt value needs to be tested before assignment when function is used
	char* dummy;
	float converted = strtod(input.c_str(), &dummy);
	if(*dummy){ return std::numeric_limits<float>::quiet_NaN(); }
	return std::stof(input);
}

#endif
