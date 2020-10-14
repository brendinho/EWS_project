/*
    USE:

    . start the timint with the call TIC( <time_stop_number> )

    . finish the timing with try{ < var_to_store_time > = TOC( <time_stop_number> ) }

        the time stop number is an integer with unrestricted values
        the duration is given in milliseconds, and the variable used to store can be any numeric

    . =trying to retrieve time from an undefined stop number throws an exception, and returns a NaN value
*/

#ifndef TIMIMG_HPP_
#define TIMING_HPP_

// timing functions
#include <algorithm>
#include <chrono>
#include <exception>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

// get the keys of a map
// adapted from http://www.lonecpluspluscoder.com/2015/08/13/an-elegant-way-to-extract-keys-from-a-c-map/
template<typename Key_type, typename Value_type> std::vector<Key_type> timing_keys( std::map<Key_type,Value_type> const& input_map )
{
  	std::vector<Key_type> key_vector;
	for (auto const& element : input_map){ key_vector.push_back(element.first); }
  	return key_vector;
}

// get the values of a map - adapted from http://www.lonecpluspluscoder.com/2015/08/13/an-elegant-way-to-extract-keys-from-a-c-map/
template<typename TK, typename TV> std::vector<TV> timing_values(std::map<TK, TV> const& input_map)
{
	std::vector<TV> retval;
	for (auto const& element : input_map){ retval.push_back(element.second); }
	return retval;
}

std::map<int, time_point<Clock>> start;
std::vector<int> first_start_exception_thrown, first_finish_exception_thrown;

inline void TIC(int i)
{
	std::vector<int> time_stops = timing_keys(start);

	if( std::find(first_start_exception_thrown.begin(), first_start_exception_thrown.end(), i) != first_start_exception_thrown.end() )
	{
		return;
	}

	start[i] = Clock::now();
}

inline int TOC(int i, bool verbose = false)
{
	std::vector<int> time_stops = timing_keys(start);

	if( std::find( first_finish_exception_thrown.begin(), first_finish_exception_thrown.end(), i) != first_finish_exception_thrown.end() )
	{
		return std::numeric_limits<int>::quiet_NaN();
	}

	if( std::find( time_stops.begin(), time_stops.end(), i) == time_stops.end() )
	{
		first_finish_exception_thrown.push_back(i);
		throw std::invalid_argument( "TOC: time stop " + std::to_string(i) + " was not created." );
	}

	time_point<Clock> finish = Clock::now();

    auto milliseconds_taken = duration_cast<milliseconds>(finish - start[i]).count();

	int diff2 = milliseconds_taken;

	if(diff2==0) return 0;

	int diff = diff2/1000;

	int time_days = diff/86400;

	int time_hours = (diff - 86400*time_days)/3600;

	int time_minutes = (diff - 86400*time_days - 3600*time_hours)/60;

	int time_seconds = diff - 86400*time_days - 3600*time_hours - 60*time_minutes;

    if(verbose)
    {
        std::cout << "\t"
            << diff2 << " (total) milliseconds, "
        	<< time_days << " days, "
        	<< time_hours << " hours, "
        	<< time_minutes << " minutes, "
        	<< time_seconds << " seconds."
        	<< std::endl << std::flush;
    }

    return std::ceil(milliseconds_taken*1./1000);
}

#define TIME_CATCH catch(std::invalid_argument& e) { std::cout << e.what() << std::endl; }

#endif
