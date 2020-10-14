/*
Brendon Phillips
PhD candidate
Department of Applied Mathematics
University of Waterloo
*/

#ifndef SOCIAL_DATA_HEADER_
#define SOCIAL_DATA_HEADER_

#include "Parameters_hesitance.hpp"

class C_SocialData
{ // these objects do nothing for tehmselves, Jesus Christ
	private:

		// std::set<int> 	_like_minded_neighbours;
		char 			_state;
		int 			_opinion_changes;

	public:

		// void 				Agreeing_Neighbour_Add(int person) { _like_minded_neighbours.insert(person); }
		// const int 			Agreeing_Neighbour_Count() const { return _like_minded_neighbours.size(); }
		// void 				Agreeing_Neighbour_Delete(const int person) { _like_minded_neighbours.erase(person); }
		// void 				Agreeing_Neighbour_Reset() { _like_minded_neighbours.clear(); }
		// const std::set<int>	Agreeing_Neighbour_All(void) const { return _like_minded_neighbours; }

		bool State_Change_To(char new_state)
		{
			if(State_Get() == new_state) { return false; }
			// Agreeing_Neighbour_Reset();
			if(State_Get() !=  '\0') ++ _opinion_changes;
			State_Set(new_state);
			return true;
		}

		const char 	State_Get() const { return _state; }

		void 		State_Set(const char new_state) { _state = new_state; }

		const int 	Opinion_Changes_Count() const { return _opinion_changes; }
		void 		Opinion_Changes_Reset() { _opinion_changes = 0; }

		void 		Save(TSOut& SOut) const { ; }

		const int 	Spin_Get() const
		{
			if(State_Get() == char('V')) return +1;
			if(State_Get() == char('H')) return  0;
			if(State_Get() == char('N')) return -1;
			return std::numeric_limits<int>::quiet_NaN();
		}

		const bool	Is_Vaccinator() const { return (State_Get() == char('V')); }
		const bool 	Is_Hesitant()	const { return (State_Get() == char('H')); }
		const bool	Is_Evil()	const { return (State_Get() == char('N')); }

		C_SocialData& operator = (const C_SocialData& other)
		{
			if(this == &other) return *this;
			// _like_minded_neighbours = other.Agreeing_Neighbour_All();
			_state = other.State_Get();
			_opinion_changes = other.Opinion_Changes_Count();
			return *this;
		}

		C_SocialData()
		{
			_state = '\0';
			// _like_minded_neighbours = {};
			_opinion_changes = 0;
		}
		C_SocialData(char initial_state)
		{
			_state = initial_state;
			// _like_minded_neighbours = {};
			_opinion_changes = 0;
		}
		C_SocialData(const C_SocialData& other)
		{
			_state = other.State_Get();
			// _like_minded_neighbours = other.Agreeing_Neighbour_All();
			_opinion_changes = other.Opinion_Changes_Count();
		}
};


#endif
