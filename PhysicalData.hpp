/*
Brendon Phillips
PhD candidate
Department of Applied Mathematics
University of Waterloo
*/

#ifndef PHYSICAL_DATA_HEADER_
#define PHYSICAL_DATA_HEADER_

#include "Parameters_hesitance.hpp"

class C_PhysicalData
{
	private:

		std::set<int>	_infected_neighbours;
		char 			_state; // check the enum in the Parameters header file
		int 			_time_infected;

	public:

		const std::set<int> Infected_Neighbour_All(void) const { return _infected_neighbours; }
		void 				Infected_Neighbour_Add(const int person) { _infected_neighbours.insert(person); }
		const int 			Infected_Neighbour_Count() const { return _infected_neighbours.size(); }
		void 				Infected_Neighbour_Delete(const int person) { _infected_neighbours.erase(person); }
		void 				Infected_Neighbour_Reset() { _infected_neighbours.clear(); }

		const bool			Is_Infected() const { return (_state == 'I'); }

		void 				Save(TSOut& SOut) const { ; }

		const char 			State_Get() const { return _state; }

		bool State_Change_To(char new_state)
		{
			if(_state == new_state) { return false; }
			if(new_state == char('I')) { Time_Infected_Set(1); }
			else { Time_Infected_Set(0); }
			State_Set(new_state);
			return true;
		}
		void 				State_Set(char new_state) { _state = new_state; }

		const int 			Time_Infected_Get() const { return _time_infected; }
		void 				Time_Infected_Increment(void) { ++ _time_infected; }
		void 				Time_Infected_Set(const int time) { _time_infected = time; }

		C_PhysicalData& operator = (const C_PhysicalData& other)
		{
			if(this == &other) return *this;
			_infected_neighbours = other.Infected_Neighbour_All();
			_state = other.State_Get();
			_time_infected = other.Time_Infected_Get();
			return *this;
		}

		C_PhysicalData()
		{
			State_Set('\0');
			Time_Infected_Set(-1);
		}
		C_PhysicalData(char initial_state)
		{
			State_Set(initial_state);
			if(initial_state == char('I')) { Time_Infected_Set(1); }
		}
		C_PhysicalData(const C_PhysicalData& other)
		{
			_infected_neighbours = other.Infected_Neighbour_All();
			_state = other.State_Get();
			_time_infected = other.Time_Infected_Get();
		}
};

#endif
