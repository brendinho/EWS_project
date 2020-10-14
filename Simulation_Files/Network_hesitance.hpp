/*
	Brendon Phillips
	PhD Candidate
	Bauch computational epidemiology research group
	Department of Applied Mathematics
	Faculty of Mathematics
	Universiity of Waterloo
*/

#ifndef NETWORK_IMPLEMENTATION_HPP_
#define NETWORK_IMPLEMENTATION_HPP_

#include "Layer_hesitance.hpp"
#include "SocialLayer_hesitance.hpp"

std::random_device private_rd;
std::mt19937_64 private_Mersenne(private_rd());
std::uniform_real_distribution<double> random_double(0,1);

class C_Network // : private C_Layer<C_SocialData>, private C_Layer<C_PhysicalData>
{
    private:

		C_Connectivity_Layer<C_SocialData>	_social_layer;
	    C_Layer<C_PhysicalData>				_physical_layer;

		const double 		_birth_death_ratio;
	   	const double 		_case_importation_ratio;
		const double 		_social_norm;
		const double 		_initial_vaccinator_proportion;
		const double 		_random_opinion_switch_probability;
		const int	 		_instance;
		const double 		_perceived_vaccine_risk;
		const int 			_average_node_degree;
		const double 		_infection_prob;
		const std::string 	_network_structure;

	    std::map<std::multiset<char>, int> 	_joint_state_distribution;
	    std::vector<std::stringstream>    	_stats_file_buffer;

		// vectors to check for the equilibrium if it's reached
		// keys and spaces defined in the class constructor
		std::map<std::string, std::vector<float>> _the_check_vectors;

    public:

	    const bool Connect_Nodes(const std::string layer_name, const int person1, const int person2)
	    {
	        assert(check_layer(layer_name));
			if( person1 == person2){ return false; }
	        if(layer_name == "social"){ return _social_layer.Edge_Add(person1, person2); }
	        if(layer_name == "physical"){ return _physical_layer.Edge_Add(person1, person2); }
			return false;
	    }

		const int Num_Edges(const std::string layer_name) const
		{
		    assert(check_layer(layer_name));
		    if(layer_name == "social"){ return _social_layer.Edge_Count("real"); }
		    if(layer_name == "physical"){ return _physical_layer.Edge_Count(); }
		    return 0;
		}

		const int Num_Friendships(void) const
		{
			return _social_layer.Edge_Count("virtual");
		}

		const int Num_Neighbours(const std::string layer_name, const int node) const
		{
			assert(check_layer(layer_name));
			return Neighbours(layer_name, node).size();
		}

		const std::set<int> Neighbours(const std::string layer_name, const int person) const
		{
			assert(check_layer(layer_name));
			if(layer_name == "social"){ return _social_layer.Neighbours_All("real", person); }
			if(layer_name == "physical"){ return _physical_layer.Neighbours_All(person); }
			return {};
		}

		const int Neighbour_Random_ID(const std::string layer_name, const int person) const
		{
			assert(check_layer(layer_name));
			if(layer_name == "social"){ return _social_layer.Neighbour_Random_ID("real", person); }
			if(layer_name == "physical"){ return _physical_layer.Neighbour_Random_ID(person); }
			return std::numeric_limits<int>::quiet_NaN();
		}

		const std::set<int> Agreeing_Neighbours(const int person) const
		{
			return _social_layer.Neighbours_All("virtual", person);
		}

		const std::set<int> Nodes(const std::string layer_name) const
		{
			assert(check_layer(layer_name));
			if(layer_name == "social"){ return _social_layer.Nodes_All(); }
			if(layer_name == "physical"){ return _physical_layer.Nodes_All(); }
			return {};
		}

		const int Size(const std::string layer_name = "\0") const
		{
			assert(check_layer(layer_name));
			if(layer_name == "social"){ return _social_layer.Node_Count(); }
			if(layer_name == "physical"){ return _physical_layer.Node_Count(); }
			return std::max(_social_layer.Node_Count(), _physical_layer.Node_Count());
		}

		const char State(const std::string layer_name, const int person)
		{
			if(layer_name == "social"){ return _social_layer.Node_State_Get(person); }
			if(layer_name == "physical"){ return _physical_layer.Node_State_Get(person); }
			return '\0';
		}

		const bool Is_Vaccinator(const int person) const
		{
			assert(check_layer("social"));
			return (_social_layer.Node_State_Get(person) == 'V');
		}

		const bool Is_Evil(const int person) const
		{
			assert(check_layer("social"));
			return (_social_layer.Node_State_Get(person) == 'N');
		}

		const bool Is_Hesitant(const int person) const
		{
			assert(check_layer("social"));
			return (_social_layer.Node_State_Get(person) == 'H');
		}

		const bool Is_Infected(const int person) const
		{
			assert(check_layer("physical"));
			return (_physical_layer.Node_State_Get(person) == 'I');
		}

		const bool Node_Is_Infectious(const int person) const
		{
			return (_physical_layer.Node(person).Time_Infected_Get() > global_infection_latency);
		}

	    const bool Increment_Infected_Time(void)
	    {
		    if(global_infection_duration <= 0){ return false; }
		    for(int infected_node : _physical_layer.Nodes_In_State('I'))
		    {
			    if(_physical_layer.Node(infected_node).Time_Infected_Get() == global_infection_duration)
			    {
				    _physical_layer.Node_State_Change_To(infected_node, 'R');
				    -- _joint_state_distribution[{ _social_layer.Node_State_Get(infected_node), 'I' }];
				    ++ _joint_state_distribution[{ _social_layer.Node_State_Get(infected_node), 'R' }];
			    }
			    else
			    {
				    _physical_layer.Node(infected_node).Time_Infected_Increment();
			    }
		    }
		    return true;
	    }

		const std::map<std::multiset<char>, int> Joint_State_Distribution(void) const
		{
			return _joint_state_distribution;
		}

		const std::map<std::multiset<char>, int> Joint_State_Distribution(const char state) const
		{
			assert( check_layer("social", state) or check_layer("physical", state) );
			std::map<std::multiset<char>, int> temp_map;
			for(auto state_pair : _joint_state_distribution){ if( state_pair.first.count(state) ){ temp_map.insert(state_pair); } }
			return temp_map;
		}

		const int Joint_State_Distribution(const char state1, const char state2) const
		{
			// // the states must be in one state set of the other
			assert( check_layer("social", state1, false) or check_layer("physical", state1, false) );
			assert( check_layer("social", state2, false) or check_layer("physical", state2, false) );
			// both states can't be in the same set
			auto check = [=](void) -> bool
			{
				if( (state1 != 'V') and (state2 != 'V') )
				{
					if( check_layer("social", state1, false) and check_layer("physical", state1, false) ){ return false; }
					if( check_layer("social", state2, false) and check_layer("physical", state2, false) ){ return false; }
				}
				return true;
			};
			assert( check() );

			// if(_joint_state_distribution.find({state1, state2}) == _joint_state_distribution.end()) return 0; // at least one of the compartments requested has no nodes
			return _joint_state_distribution[{state1, state2}];
		}

		void Generate_Network_Random(const int number_of_edges)
		{
			const int number_of_nodes = std::min(_social_layer.Node_Count(), _physical_layer.Node_Count());
			const double ran_edge_prob = number_of_edges/(double)std::pow(number_of_nodes, 2);
	        // following the algorithm given in the paper "Efficent Generation of Large Random Networks"
	        int v = 1, w = -1;
	        while( v < number_of_nodes )
	        {
	            w += 1 + std::floor( std::log(1-random_double(private_Mersenne)) / std::log(1-ran_edge_prob) );
	            while( (w >= v) and (v < number_of_nodes) )
	            {
	                w -= v;
	                ++ v;
	            }
	            if(v < number_of_nodes)
	            {
					Connect_Nodes("social", v, w);
					Connect_Nodes("physical", v, w);
	            }
	        }
	        for(auto person : Nodes("social"))
	        { // can use eiyher the physical or social layers for the iteration of a single loop, since
	            if( Num_Neighbours("social", person) == 0 )
	            {
	                // randomly generate a new neighbour, going until we find one that doesn't create a loop
	                int new_neighbour;
	                do
	                { new_neighbour = _social_layer.Node_Random_ID(); }
	                while(new_neighbour == person);

	                Connect_Nodes("social", person, new_neighbour);
	                Connect_Nodes("physical", person, new_neighbour);
	            }
	        }
		    return;
		}

		void Generate_Network_Smallworld(const float rewiring_prob, const int avg_node_degree)
		{
			for(int person : range(Size()))
			{
				for(int neighbour : range(person+1, person+avg_node_degree+1))
				{
					int new_neighbour = neighbour%Size();

					if(random_double(private_Mersenne) < rewiring_prob)
					{
						int rand_node = -1;
						// rewire the link by never creating it in the first replacement
						// just finding them a random neighbour that ain't this dude
						do
						{
							rand_node = _social_layer.Node_Random_ID();
						} while( (rand_node == new_neighbour) or (rand_node==person) );

						new_neighbour = rand_node;
					}

					// std::cout << "connecting nodes " << person << " and " << neighbour << std::endl;
					Connect_Nodes("social", person, new_neighbour);
					Connect_Nodes("physical", person, new_neighbour);
				}
			}
			// std::cout << "stats of the social network" << std::endl; _social_layer.Statistics_Print();
			// std::cout << std::endl;
			std::cout << "stats of the physical network" << std::endl; _physical_layer.Statistics_Print();
		}

		const bool Change_State(const std::string layer_name, const int person, const char new_state, char killed_and_reborn = false)
		{
			/*
				uses evolutionary rules on the layer to determine the state changes of each node
			*/
			char old_social_state = State("social", person);
			char old_physical_state = State("physical", person);

			auto check_passed = [=](void) -> bool
			{
				if( killed_and_reborn ) return true;
				if(layer_name == "physical")
				{
					if( old_physical_state == 'R' ){ std::cout << "attempted illegal transition out of Recovered state." << std::endl; return false; }
					else if( old_physical_state == 'V' ){ std::cout << "attempted illegal transition out of Vaccinated state." << std::endl; return false; }
					else if( (old_physical_state != 'S') and (new_state == 'I') )
					{
						std::cout << "attempted illegal transition from state " << old_physical_state << " to physical state I." << std::endl;
						exit(EXIT_FAILURE);
					}
					else if( (old_physical_state != 'I') and (new_state == 'R') )
					{
						std::cout << "attempted illegal transition from state " << old_physical_state << " to physical state R." << std::endl;
						exit(EXIT_FAILURE);
					}
					else if( (old_physical_state != 'S') and (new_state == 'V') )
					{
						std::cout << "attempted illegal transition from state " << old_physical_state << " to physical state V." << std::endl;
						exit(EXIT_FAILURE);
					}
				}
				else if(layer_name == "social")
				{
					if( std::set<char>({old_social_state, new_state}) == std::set<char>({'N','V'}) )
					{ // person must become a straw man before fully switching opinions
						std::cout << "illegal social transition from " << old_social_state << " to " << new_state << " without passing through the Hesitant state." << std::endl;
						exit(EXIT_FAILURE);
					}
				}
				return true;
			};
			assert( check_passed() );

		    if(layer_name == "physical")
		    {
		        -- _joint_state_distribution[{old_social_state, old_physical_state}];
				++ _joint_state_distribution[{old_social_state, new_state}];
				_physical_layer.Node_State_Change_To(person, new_state);
		    }
		    else  // layer_name == social
		    {
				// printf("\trequested transition %c to %c. physical state is %c.\n", old_social_state, new_state, old_physical_state);
		        if((new_state == 'V') and (State("physical", person) == 'S')) // if they were susceptible before
	            {
					// std::cout << "\tgoing to be a susceptiple vaccinator. we can't let this happen." << std::endl;
					// changing from their old state to the new one, while becoming vaccinated
	                -- _joint_state_distribution[{old_social_state, 'S'}];
	                ++ _joint_state_distribution[{'V', 'V'}];
	                _physical_layer.Node_State_Change_To(person, 'V');
					_social_layer.Node_State_Change_To(person, new_state);
					// printf("\tchanged from (%c,%c) to (%c,%c).\n", old_social_state, old_physical_state, State("social", person), State("physical", person));
	            }
				else if(old_social_state != new_state) // the only special case to be considered was the above, otherwise SNAFU
				{
					// std::cout << "\told state isn't equal to the new state" << std::endl;
					_social_layer.Node_State_Change_To(person, new_state);
					// std::cout << "\tperson state changed to " << old_social_state << " to " << State("social", person) << std::endl;
					-- _joint_state_distribution[{old_physical_state, old_social_state}];
					++ _joint_state_distribution[{old_physical_state, new_state}];
					// std::cout << "\tchanged the joint state distribution" << std::endl;
				}
		    }
		    return true;
		}

		const std::set<int> Nodes_In_State(const std::string layer_name, const char state) const
		{
		    assert(check_layer(layer_name, state));
		    if(layer_name == "social"){ return _social_layer.Nodes_In_State(state); }
		    if(layer_name == "physical"){ return _physical_layer.Nodes_In_State(state); }
		    return {};
		}

		const int Nodes_In_State_Count(const std::string layer_name, const char state) const
		{
			assert(check_layer(layer_name, state));
		    if(layer_name == "social"){ return _social_layer.Nodes_In_State_Count(state); }
		    if(layer_name == "physical"){ return _physical_layer.Nodes_In_State_Count(state); }
		    return std::numeric_limits<int>::quiet_NaN();
		}

		const double Prob_Has_Been_Sick(const char social_state) const
		{ // returns the probability of a node in the given state having a history of illness, rather than
		    // never being infected
		    assert(check_layer("social", social_state));
		    return (
						Joint_State_Distribution('I', social_state) +
						Joint_State_Distribution('R', social_state)
					)
		        	/ (double) Nodes_In_State_Count("social", social_state);
		}

		const int Random_Node_In_State(const std::string layer_name, const char state)
		{ // returns an iterator to a random node in the requested state on a layer
		    assert(check_layer(layer_name, state));
		    if(layer_name == "social")
		    {
			    return _social_layer.Node_In_State_Random(state);
		    }
		    else if(layer_name == "physical")
		    {
			    return _physical_layer.Node_In_State_Random(state);
		    }
			return std::numeric_limits<int>::quiet_NaN();
		}

		const double Ratio_Infected_Neighbours(const char social_state) const
		{ // tells us how many of their social neighbours are infected with the disease
		    assert(check_layer("social", social_state));
			double sum_int_of_ratios = 0;
		    for(auto person : _social_layer.Nodes_All())
		    {
			    double num_infected_social_neighbours = 0.;
			    for(auto social_neighbour : _social_layer.Neighbours_All("real", person))
			    {
			    	if ( Is_Infected(social_neighbour) ){ ++ num_infected_social_neighbours; }
			    }
			    sum_int_of_ratios += num_infected_social_neighbours / (double) _social_layer.Neighbour_Count("real", person);
		    }
		    return sum_int_of_ratios / (double) _social_layer.Node_Count();
		}

		const std::vector<int> Shuffled_Nodes(const std::string layer_name) const
		{
			assert( check_layer(layer_name) );
		    if(layer_name == "social") return _social_layer.Nodes_Shuffled();
		    if(layer_name == "physical") return _physical_layer.Nodes_Shuffled();
			return std::vector<int>({});
		}

		const int Join_Count(const float proportion, const std::set<char> states) const
		{
			assert(check_layer("social"));
			return _social_layer.Join_Count(proportion, states);
		}

		const std::set<int> Node_Infected_Neighbours(const int person) const
		{
			std::set<int> sick_neighbours = {};
			for(auto social_neighbour : _social_layer.Neighbours_All("real", person))
			{
				if ( Is_Infected(social_neighbour) ){ sick_neighbours.insert(social_neighbour); }
			}
			return sick_neighbours;
		}

		const int Node_Infected_Neighbour_Count(const int person) const
		{
			int sick_neighbours = 0;
			for(auto social_neighbour : _social_layer.Neighbours_All("real", person))
			{
				if ( Is_Infected(social_neighbour) ){ ++ sick_neighbours; }
			}
			return sick_neighbours;
		}

		const int Square_Sum_Degrees(const std::string mode = "virtual") const
		{ // include a check to make sure that the states are not in the same state, except for the 'V'
		 	return _social_layer.Square_Sum_Degrees(mode);
		}

		const double Modularity_Score(const float proportion, const char state = '\0') const
		{
			assert( check_layer("social", state) );
			return _social_layer.Modularity_Score(proportion, state);
		}

		const double Watts_Strogatz_Avg_Clustering_Coeff(const float proportion, const char state = '\0') const
		{
			if(state != '\0') assert(check_layer("social", state));
			return _social_layer.Watts_Strogatz_Avg_Clustering_Coefficient(proportion, state);
		}

		const double Mutual_Information(const float proportion) const
		{ // calculates the Mutual Information statistic of the layer interaction of the network

			const std::vector<int> temp = Shuffled_Nodes("social");
			std::vector<int> all_nodes( temp.begin(), std::next(temp.begin(), std::ceil(proportion*Size())) );

			auto filtered = [=](const std::string the_layer, const char desired_state) -> float
			{
				return 1.*std::accumulate
				(
					all_nodes.begin(),
					all_nodes.end(),
					0.,
					[&](float accumulator, int person){ return accumulator + (State(the_layer, person) == desired_state); }
				) / all_nodes.size();
			};

			auto joint_prop = [=](const char desired_social_state, const char desired_physical_state) -> float
			{
				std::set<char> state_set { desired_social_state, desired_physical_state };
				return 1.*std::accumulate
				(
					all_nodes.begin(),
					all_nodes.end(),
					0.,
					[=](float accumulator, int person)
					{
						std::set<char> person_set = std::set<char>({State("social", person), State("physical", person)});
						return accumulator + ( person_set == state_set );
					}
				) / all_nodes.size();
			};

		    double network_size = Size();
		    double mutual_information_running_total = 0;

			// printf("\t\tMutual Information function\n");
		    for(char social_state : global_social_process)
		    {
				double p_soc = filtered("social", social_state);
		        for(char physical_state : global_physical_process)
		        {
					// printf("\t\t\t#%c (social), #%i; #%c (physical), %i; joint (%c,%c): %i.\n", social_state, Nodes_In_State_Count("social", social_state), physical_state, Nodes_In_State_Count("physical", physical_state), social_state, physical_state, Joint_State_Distribution(social_state, physical_state));
		            double p_phys = filtered("physical", physical_state);
		            double p_joint = joint_prop(social_state, physical_state);
		            if( p_joint > 0 ){ mutual_information_running_total += p_joint * std::log2( p_joint / (double)(p_soc*p_phys) ); }
					// printf("\t\t\trunning total: %f.\n", mutual_information_running_total);
		        }
		    }

		    return mutual_information_running_total;
		}

		const std::map<int, int> Connected_Component_Sizes(const float proportion, const char state = '\0') const
		{
			assert( check_layer("social", state) );
			return _social_layer.Connected_Component_Sizes(proportion, state);
		}

		const std::map<int, int> Echo_Chamber_Sizes(const float proportion, const char state = '\0') const
		{
			assert( check_layer("social", state) );
			return _social_layer.Echo_Chamber_Sizes(proportion, state);
		}

		const int Number_Triads(const float proportion, const char social_state = '\0')
		{
			return _social_layer.Number_Triads(proportion, social_state);
		}

		const float Opinion_Changes_Per_Person(const float proportion, const char social_state = '\0') const
		{
			return _social_layer.Opinion_Changes_Per_Person(proportion, social_state);
		}

		const float Diameter(const float proportion, const char state = '\0')
		{
			return _social_layer.Diameter(proportion, state);
		}

		void Push_Title_Lines(void)
		{ // write the title lines of the files to the respective buffers

		    // printing pure numbers, with the data being normalised during post-processing, so views can be changeed
		    assert( global_legal_layers.find("social") != global_legal_layers.end() );

			_stats_file_buffer[0]
				<< "# Brendon Phillips"
				<< "\n# PhD student"
				<< "\n# Computational epidemoilogy research group"
				<< "\n# Supervisor: Dr Chris Bauch"
				<< "\n# Department of Applied Mathematics"
				<< "\n# University of Waterloo"

			<< "\n\n"

				<< "N"                              		<< ","
				<< "instance"                       		<< ","
				<< "network_structure"                		<< ","
				<< "mean_degree"    						<< ","
				<< "perceived_vaccine_risk"                 << ","
				<< "infec_prob"                     		<< ","
				<< "social_norm"                    		<< ","
				<< "importation"                    		<< ","
				<< "initial_vacc_proportion"        		<< ","
				<< "random_opinion_switch"					<< ","
				<< "proportion_of_nodes"					<< ","
				<< "time"                           		<< ","

				<< "phys_S"                         		<< ","
				<< "phys_I"                         		<< ","
				<< "phys_R"                         		<< ","
				<< "phys_V"                         		<< ","
				<< "soc_H"									<< ","
				<< "soc_N"                          		<< ","
				<< "soc_V"                          		<< ","

				<< "mutual_info"                    		<< ","

				<< "watts_strogatz_hesitant"				<< ","
				<< "watts_strogatz_nonvacc" 				<< ","
				<< "watts_strogatz_vacc" 					<< ","
				<< "watts_strogatz_all" 					<< ","

				<< "prob_sick_hesitant"						<< ","
				<< "prob_sick_nonvacc"              		<< ","
				<< "prob_sick_vacc"                 		<< ","

				<< "HH_join_count"							<< ","
				<< "NN_join_count"                  		<< ","
				<< "VV_join_count"                  		<< ","
				<< "HN_join_count"							<< ","
				<< "HV_join_count"							<< ","
				<< "NV_join_count"                  		<< ","

				<< "hesitant_min_conn_comp_size"          	<< ","
				<< "hesitant_max_conn_comp_size"          	<< ","
				<< "hesitant_avg_conn_comp_size"         	<< ","
				<< "hesitant_number_conn_comps"          	<< ","

				<< "nonvacc_min_conn_comp_size"          	<< ","
				<< "nonvacc_max_conn_comp_size"          	<< ","
				<< "nonvacc_avg_conn_comp_size"         	<< ","
				<< "nonvacc_number_conn_comps"          	<< ","

				<< "vacc_min_conn_comp_size"            	<< ","
				<< "vacc_max_conn_comp_size"            	<< ","
				<< "vacc_avg_conn_comp_size"            	<< ","
				<< "vacc_number_conn_comps"             	<< ","

				<< "hesitant_min_chamber_size"            	<< ","
				<< "hesitant_max_chamber_size"            	<< ","
				<< "hesitant_avg_chamber_size"            	<< ","
				<< "hesitant_number_chambers"            	<< ","

				<< "nonvacc_min_chamber_size"            	<< ","
				<< "nonvacc_max_chamber_size"            	<< ","
				<< "nonvacc_avg_chamber_size"            	<< ","
				<< "nonvacc_number_chambers"            	<< ","

				<< "vacc_min_chamber_size"               	<< ","
				<< "vacc_max_chamber_size"               	<< ","
				<< "vacc_avg_chamber_size"              	<< ","
				<< "vacc_number_chambers"               	<< ","

				<< "hesitant_opinion_change"				<< ","
				<< "nonvacc_opinion_change"        		 	<< ","
				<< "vacc_opinion_change"            		<< ","
				<< "total_opinion_change"           		<< ","

				<< "ratio_hesitant_infected_neighbours"		<< ","
				<< "ratio_nonvacc_infected_neighbours"    	<< ","
				<< "ratio_vacc_infected_neighbours"			<< ","

				<< "num_hesitant_triads"					<< ","
				<< "num_nonvacc_triads"						<< ","
				<< "num_vacc_triads"						<< ","
				<< "total_num_triads"						<< ","

				<< "hesitant_diameter"						<< ","
				<< "nonvacc_diameter"						<< ","
				<< "vacc_diameter"

			<<  "\n";
		    // printing the header for the stats file
		    return;
		}

		const bool Record_Time_Step(const int time, const float proportion)
		{ // write the title lines of the files to the respective buffers

			bool we_should_stop = true;

			auto variance_ceiling_of_metric = [=](std::string name) -> float
			{
				if( std::set<std::string>({"num_S", "num_I", "num_R"}).count(name) )
				{ return Size("physical") * global_equilibrium_threshold; }

				if( std::set<std::string>({"num_N", "num_H"}).count(name) )
				{ return Size("social") * global_equilibrium_threshold; }

				if( std::set<std::string>({"HH_joins", "NN_joins", "VV_joins", "HN_joins", "HV_joins", "NV_joins"}).count(name) ) // , "HV_joins"
				{ return Num_Edges("social") * global_equilibrium_threshold; }

				if( std::set<std::string>({"joint_dist_VI", "joint_dist_VR", "joint_dist_NS", "joint_dist_NI", "joint_dist_NR", "joint_dist_HS", "joint_dist_HI", "joint_dist_HR"}).count(name) )
				{ return Size() * global_equilibrium_threshold; }

				return std::numeric_limits<float>::quiet_NaN();
			};

			if(time > global_stored_history_length)
			{
				for(std::string name : global_check_these_metrics)
				{
					// keep the length the same by popping the front and inserting at the back
					// at the time of coding, I found this to be better time than std::vector::pop_front(), which is O(n)
					_the_check_vectors.at(name).erase( _the_check_vectors.at(name).begin() );
				}
			}

			const double network_size = Size();
			assert( global_legal_layers.find("social") != global_legal_layers.end() );

			// vectors of the last few results of the important matrics to check for equilibrium
			// the lambdas for min, max, sum and avg of the vectors are below
			// if all the mmetrics are under the thresholds, then we stop

			// TIC(1000);
			_the_check_vectors.at("num_S").push_back( Nodes_In_State_Count("physical", 'S') );
			_the_check_vectors.at("num_I").push_back( Nodes_In_State_Count("physical", 'I') );
			_the_check_vectors.at("num_R").push_back( Nodes_In_State_Count("physical", 'R') );
			_the_check_vectors.at("num_N").push_back( Nodes_In_State_Count("social",   'N') );
			_the_check_vectors.at("num_H").push_back( Nodes_In_State_Count("social",   'H') );

			_the_check_vectors.at("joint_dist_VI").push_back( Joint_State_Distribution('V', 'I') );
			_the_check_vectors.at("joint_dist_VR").push_back( Joint_State_Distribution('V', 'R') );

			_the_check_vectors.at("joint_dist_NS").push_back( Joint_State_Distribution('N', 'S') );
			_the_check_vectors.at("joint_dist_NI").push_back( Joint_State_Distribution('N', 'I') );
			_the_check_vectors.at("joint_dist_NR").push_back( Joint_State_Distribution('N', 'R') );

			_the_check_vectors.at("joint_dist_HS").push_back( Joint_State_Distribution('H', 'S') );
			_the_check_vectors.at("joint_dist_HI").push_back( Joint_State_Distribution('H', 'I') );
			_the_check_vectors.at("joint_dist_HR").push_back( Joint_State_Distribution('H', 'R') );

			_the_check_vectors.at("HH_joins").push_back( Join_Count(1, {'H', 'H'}) );
			_the_check_vectors.at("NN_joins").push_back( Join_Count(1, {'N', 'N'}) );
			_the_check_vectors.at("VV_joins").push_back( Join_Count(1, {'V', 'V'}) );
			_the_check_vectors.at("HN_joins").push_back( Join_Count(1, {'H', 'N'}) );
			_the_check_vectors.at("HV_joins").push_back( Join_Count(1, {'H', 'V'}) );
			_the_check_vectors.at("NV_joins").push_back( Join_Count(1, {'N', 'V'}) );
			// try{ std::cout << "\t\tgrabbing the check vector metrics: " << TOC(1000) << " secs." << std::endl; } TIME_CATCH

			// TIC(1001);
			const std::map<int, int> hesitant_conn_comps = Connected_Component_Sizes(proportion, 'H');
			const std::map<int, int> nonvacc_conn_comps = Connected_Component_Sizes(proportion, 'N');
			const std::map<int, int> vacc_conn_comps = Connected_Component_Sizes(proportion, 'V');
			// try{ std::cout << "\t\tgrabbing the connected components: " << TOC(1001) << " secs." << std::endl; } TIME_CATCH

			// TIC(1002);
			const std::map<int, int> hesitant_echo_chambs = Echo_Chamber_Sizes(proportion, 'H');
			const std::map<int, int> nonvacc_echo_chambs = Echo_Chamber_Sizes(proportion, 'N');
			const std::map<int, int> vacc_echo_chambs = Echo_Chamber_Sizes(proportion, 'V');
			// try{ std::cout << "\t\tgrabbing the echo chambers: " << TOC(1002) << " secs." << std::endl; } TIME_CATCH

			// TIC(1005);
			const int num_hesitant_triads = Number_Triads(proportion, 'H');
			const int num_nonvacc_triads = Number_Triads(proportion, 'N');
			const int num_vacc_triads =  Number_Triads(proportion, 'V');
			const int num_total_triads = Number_Triads(proportion);
			// try{ std::cout << "\t\tgrabbing the numbers of triads: " << TOC(1005) << " secs." << std::endl; } TIME_CATCH

			// TIC(1008);
			const float hesitant_GCC = Watts_Strogatz_Avg_Clustering_Coeff(proportion, 'H');
			const float nonvacc_GCC = Watts_Strogatz_Avg_Clustering_Coeff(proportion, 'N');
			const float vacc_GCC = Watts_Strogatz_Avg_Clustering_Coeff(proportion, 'V');
			const float total_GCC = Watts_Strogatz_Avg_Clustering_Coeff(proportion);
			// try{ std::cout << "\t\tgrabbing the CCs: " << TOC(1008) << " secs." << std::endl; } TIME_CATCH

			// TIC(1009);
			const int square_sum_degrees = Square_Sum_Degrees();
			// try{ std::cout << "\t\tgetting the square sum of degrees: " << TOC(1009) << " secs."<< std::endl; } TIME_CATCH

			// TIC(1010);
			const int HH_join_count = Join_Count(proportion, {'H','H'});
			const int NN_join_count = Join_Count(proportion, {'N','N'});
			const int VV_join_count = Join_Count(proportion, {'V','V'});
			const int HN_join_count = Join_Count(proportion, {'H','N'});
			const int HV_join_count = Join_Count(proportion, {'H','V'});
			const int NV_join_count = Join_Count(proportion, {'N','V'});
			// try{ std::cout << "\t\tgetting the join counts: " << TOC(1010) << " secs." << std::endl; } TIME_CATCH

			// TIC(1011);
			const float mutual_information = Mutual_Information(proportion);
			// try{ std::cout << "\t\tcalculating the mutual information: " << TOC(1011) << " secs." << std::endl; } TIME_CATCH

			// TIC(1012);
			const float hesitant_diameter = Diameter(proportion, 'H');
			const float nonvacc_diameter = Diameter(proportion, 'N');
			const float vacc_diameter = Diameter(proportion, 'V');
			const float total_diameter = Diameter(proportion);
			// try{ std::cout << "\t\tcalculating the path lengths: " << TOC(1012) << " secs." << std::endl; } TIME_CATCH

			// TIC(1003);
			const std::vector<int> hesitant_conn_comp_sizes = extract( hesitant_conn_comps, "keys");
			const std::vector<int> nonvacc_conn_comp_sizes = extract( nonvacc_conn_comps, "keys" );
			const std::vector<int> vacc_conn_comp_sizes = extract( vacc_conn_comps, "keys" );

			const std::vector<int> hesitant_echo_chamber_sizes = extract( hesitant_echo_chambs, "keys" );
			const std::vector<int> nonvacc_echo_chamber_sizes = extract( nonvacc_echo_chambs, "keys" );
			const std::vector<int> vacc_echo_chamber_sizes = extract( vacc_echo_chambs, "keys" );

			const std::vector<int> hesitant_conn_comp_counts = extract( hesitant_conn_comps, "values");
			const std::vector<int> nonvacc_conn_comp_counts = extract( nonvacc_conn_comps, "values" );
			const std::vector<int> vacc_conn_comp_counts = extract( vacc_conn_comps, "values" );

			const std::vector<int> hesitant_echo_chamber_counts = extract( hesitant_echo_chambs, "values" );
			const std::vector<int> nonvacc_echo_chamber_counts = extract( nonvacc_echo_chambs, "values" );
			const std::vector<int> vacc_echo_chamber_counts = extract( vacc_echo_chambs, "values" );
			// try{ std::cout << "\t\tgrabbing the keys and values: " << TOC(1003) << " secs." << std::endl; } TIME_CATCH

			// TIC(1004);
			_stats_file_buffer.at(time+1)
				<< network_size                             				<< ","
				<< _instance                      							<< ","
				<< _network_structure                							<< ","
				<< _average_node_degree    									<< ","
				<< _perceived_vaccine_risk                      			<< ","
				<< _infection_prob                    						<< ","
				<< _social_norm                    							<< ","
				<< _case_importation_ratio	        		 				<< ","
				<< _initial_vaccinator_proportion        					<< ","
				<< _random_opinion_switch_probability 						<< ","
				<< proportion												<< ","
				<< time                           							<< ","

				<< Nodes_In_State_Count("physical", 'S') 					<< ","
				<< Nodes_In_State_Count("physical", 'I') 					<< ","
				<< Nodes_In_State_Count("physical", 'R')  					<< ","
				<< Nodes_In_State_Count("physical", 'V')   	 				<< ","
				<< Nodes_In_State_Count("social", 'H')						<< ","
				<< Nodes_In_State_Count("social", 'N')  					<< ","
				<< Nodes_In_State_Count("social", 'V')  					<< ","

				<< mutual_information										<< ","

				<< hesitant_GCC												<< ","
				<< nonvacc_GCC												<< ","
				<< vacc_GCC													<< ","
				<< total_GCC												<< ","

				<< Prob_Has_Been_Sick('H')									<< ","
				<< Prob_Has_Been_Sick('N')  								<< ","
				<< Prob_Has_Been_Sick('V')  								<< ","

				<< HH_join_count											<< ","
				<< NN_join_count                							<< ","
				<< VV_join_count											<< ","
				<< HN_join_count 											<< ","
				<< HV_join_count                							<< ","
				<< NV_join_count                							<< ","

				<< min_int(hesitant_conn_comp_sizes)             			<< ","
				<< max_int(hesitant_conn_comp_sizes)             			<< ","
				<< avg_int(hesitant_conn_comp_sizes)             			<< ","
				<< sum_int(hesitant_conn_comp_counts)           			<< ","

				<< min_int(nonvacc_conn_comp_sizes)             			<< ","
				<< max_int(nonvacc_conn_comp_sizes)             			<< ","
				<< avg_int(nonvacc_conn_comp_sizes)             			<< ","
				<< sum_int(nonvacc_conn_comp_counts)           				<< ","

				<< min_int(vacc_conn_comp_sizes)             				<< ","
				<< max_int(vacc_conn_comp_sizes)             				<< ","
				<< avg_int(vacc_conn_comp_sizes)             				<< ","
				<< sum_int(vacc_conn_comp_counts)             				<< ","

				<< min_int(hesitant_echo_chamber_sizes)         			<< ","
				<< max_int(hesitant_echo_chamber_sizes)          			<< ","
				<< avg_int(hesitant_echo_chamber_sizes)          			<< ","
				<< sum_int(hesitant_echo_chamber_counts)         			<< ","

				<< min_int(nonvacc_echo_chamber_sizes)         				<< ","
				<< max_int(nonvacc_echo_chamber_sizes)          			<< ","
				<< avg_int(nonvacc_echo_chamber_sizes)          			<< ","
				<< sum_int(nonvacc_echo_chamber_counts)         			<< ","

				<< min_int(vacc_echo_chamber_sizes)             			<< ","
				<< max_int(vacc_echo_chamber_sizes)             			<< ","
				<< avg_int(vacc_echo_chamber_sizes)             			<< ","
				<< sum_int(vacc_echo_chamber_counts)            			<< ","

				<< Opinion_Changes_Per_Person(proportion, 'H')				<< ","
				<< Opinion_Changes_Per_Person(proportion, 'N')        		<< ","
				<< Opinion_Changes_Per_Person(proportion, 'V')            	<< ","
				<< Opinion_Changes_Per_Person(proportion)          			<< ","

				<< Ratio_Infected_Neighbours('H')							<< ","
				<< Ratio_Infected_Neighbours('N') 							<< ","
				<< Ratio_Infected_Neighbours('V')  							<< ","

				<< num_hesitant_triads										<< ","
				<< num_nonvacc_triads 										<< ","
				<< num_vacc_triads 											<< ","
				<< num_total_triads 										<< ","

				<< hesitant_diameter										<< ","
				<< nonvacc_diameter											<< ","
				<< vacc_diameter

			<<  "\n";
			// try{ std::cout << "\t\twriting to the file buffer: " << TOC(1004) << " secs." << std::endl;; } TIME_CATCH

			if( (time <= global_stored_history_length) or (time <= global_time_minimum))
			{
				we_should_stop = false;
			}
			else for(auto name : global_check_these_metrics)
			{
				if( unbiased_sd(_the_check_vectors.at(name).cbegin(), _the_check_vectors.at(name).cend()) >= variance_ceiling_of_metric(name) )
				{
					we_should_stop = false;
					break;
				}
			}
			return we_should_stop;
		}

		void Print_Data_File(const std::string output_folder_name, const std::string file_name_parameter_stem, const bool reached_equilib) const
		{
		    Push_Title_Lines();
		    Print_Data(output_folder_name, file_name_parameter_stem, reached_equilib);
		    return;
		}

		void Print_Data(const std::string output_folder_name, const std::string file_name_global_stem, const bool reached_equilib) const
		{ // general function f0r p rinting all the desired output files from the simulation

			std::string data_file_name = output_folder_name + "/stats" + file_name_global_stem;
			if( reached_equilib ){ data_file_name += "_equilib"; }
			else { data_file_name += "_timed_out"; }
			data_file_name += ".bin";

			FILE* data_file = fopen(data_file_name.c_str(), "ab");
		    for(auto& line : _stats_file_buffer)
		    {
		        std::string temp(line.str());
		        temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
		        fwrite(temp.c_str(), sizeof temp.at(0), temp.size(), data_file);
			}
		    fclose(data_file);
		    std::cout << "Simulation data file name: " << data_file_name << std::endl << std::endl;
		}

		const int Spin(const int person) const
		{
			if( Is_Vaccinator(person) ) return +1;
			if( Is_Hesitant(person) ) return  0;
			if( Is_Evil(person) ) return -1;
		}

		const double Sentiment_Change_Probability(const int person, const int neighbour) const
		{
			assert(check_layer("social"));
			if( Is_Hesitant(neighbour) )
			{
				// we get nothing from talking to strawmen
				// no action from the social norm or anything else
				return 0.0;
			}

			const char current_opinion = State("social", person);
			std::set<int> all_social_neighbours = Neighbours("social", person);

			// float angel_neighbours = 0.; float evil_neighbours = 0.;
			// for(int friend : all_social_neighbours)
			// {
			// 	if( Is_Vaccinator(friend) ) ++angel_neighbours;
			// 	else if( Is_Evil(friend) ) ++evil_neighbours;
			// }

			const float angel_neighbours = std::accumulate(
				all_social_neighbours.begin(),
				all_social_neighbours.end(),
				0.0,
				[=](float accumulator, int neighbour) -> float { return (accumulator + Is_Vaccinator(neighbour)); }
			);
			const float evil_neighbours = std::accumulate(
				all_social_neighbours.begin(),
				all_social_neighbours.end(),
				0.0,
				[=](float accumulator, int neighbour) -> float { return (accumulator + Is_Evil(neighbour)); }
			);

			float big_delta  = 0.;

			if( Is_Vaccinator(neighbour) )
			{
				// if speaking to a vaccinator, return probability of * -> V
				const float ratio = (evil_neighbours - angel_neighbours)/all_social_neighbours.size();
				big_delta = -_social_norm*ratio - (_perceived_vaccine_risk - Node_Infected_Neighbour_Count(person));
			}
			else if( Is_Evil(neighbour) )
			{
				// if speaking to an idiot, return prob of * -> N
				const float ratio = (angel_neighbours - evil_neighbours)/all_social_neighbours.size();
				big_delta = -_social_norm*ratio + (_perceived_vaccine_risk - Node_Infected_Neighbour_Count(person));
			}

			return 1./(1. + std::exp(-global_beta_parameter*big_delta));
		}

		C_Network
		(
			const std::string structure,
			const int number_of_nodes,
			const int num_neighbours,
		    const double infection,
			const double perceived_vaccine_risk,
		    const double birth_death,
		    const double importation,
		    const double social_norm,
			const double random_switch,
			const double initial_vacc_proportion,
			const int this_instance
		) :
		    _birth_death_ratio(birth_death),
		    _case_importation_ratio(importation),
			_social_norm(social_norm),
			_perceived_vaccine_risk(perceived_vaccine_risk),
			_infection_prob(infection),
			_network_structure(structure),
			_random_opinion_switch_probability(random_switch),
		    _initial_vaccinator_proportion(initial_vacc_proportion),
			_average_node_degree(num_neighbours),
			_instance(this_instance)
		{
			int number_of_edges = number_of_nodes*num_neighbours;

			_social_layer	= C_Connectivity_Layer<C_SocialData>(number_of_nodes);
			_physical_layer = C_Layer<C_PhysicalData>(number_of_nodes);

			if(structure == "smallworld") Generate_Network_Smallworld(global_rewiring_probability, num_neighbours);
			else if(structure == "random") Generate_Network_Random(number_of_edges);

			std::vector<int> all_the_nodes = _social_layer.Nodes_Shuffled();

			int initial_vaccinator_quota = ceil( _initial_vaccinator_proportion * number_of_nodes );

			auto split_iterator = std::next(all_the_nodes.begin(), initial_vacc_proportion);

			std::vector<int> vaccers( all_the_nodes.begin(), all_the_nodes.begin()+initial_vaccinator_quota );
			std::vector<int> baby_murderers( all_the_nodes.begin() + initial_vaccinator_quota, all_the_nodes.end() );

			for(int good_guy : vaccers)
			{
				_social_layer.Node_State_Change_To(good_guy, 'V');
				_physical_layer.Node_State_Change_To(good_guy, 'V');
				++ _joint_state_distribution[{'V','V'}];
			}
			for(int motherfucker : baby_murderers)
			{
			   _social_layer.Node_State_Change_To(motherfucker, 'N');
			   _physical_layer.Node_State_Change_To(motherfucker, 'S');
			   ++ _joint_state_distribution[{'N','S'}];
			}

			for(auto name : global_check_these_metrics){ _the_check_vectors[name].reserve(global_stored_history_length); }

			// plus 2, because one for the title row, and one for time zero (the start of the simulation)
		    _stats_file_buffer.resize( global_time_limit + 2 );
		}

};

#endif
