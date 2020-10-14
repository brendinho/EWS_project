/*
	Brendon Phillips
	PhD Candidate
	Bauch computational epidemiology research group
	Department of Applied Mathematics
	Faculty of Mathematics
	Universiity of Waterloo
*/

#ifndef LAYER_DATA_HPP_
#define LAYER_DATA_HPP_

#include "SocialData_hesitance.hpp"
#include "PhysicalData.hpp"

// solution adapted from ildjarn's answer at
// https://stackoverflow.com/questions/41806062/c-conditional-templates-compilation-based-on-data-type
// template<bool B, typename T = void> using enable_if_t = typename std::enable_if<B, T>::type;

template <typename NodeType> class C_Layer
{
	private:
		std::random_device _layer_random_device;
		std::mt19937_64 _layer_Mersenne;
		std::uniform_int_distribution<int> _rand_int;

		// std::map<std::string, std::string> 	_layer_information;
		std::map<char, std::set<int>>		_nodes_in_state;
		TPt<TNodeNet<NodeType>>		 		_population; // argument: SocialData/PhysicalData, as it stands
		std::vector<int> 					_shuffled_ids;

		template<typename T = NodeType> typename std::enable_if<std::is_same<T, C_PhysicalData>::value, void>::type
			Edge_Add_Adjust_Node_Neighbour_Lists(const int person1, const int person2)
		{
			if(Node_State_Get(person2) == 'I'){ Node(person1).Infected_Neighbour_Add(person2); }
			if(Node_State_Get(person1) == 'I'){ Node(person2).Infected_Neighbour_Add(person1); }
			return;
		}

		template<typename T = NodeType> typename std::enable_if<std::is_same<T, C_SocialData>::value, void>::type
			Edge_Add_Adjust_Node_Neighbour_Lists(const int person1, const int person2)
		{
			return;
		}

		template<typename T = NodeType> typename std::enable_if<std::is_same<T, C_PhysicalData>::value, void>::type
			State_Change_Adjust_Node_Neighbour_Lists(int person, const char old_state, const char new_state)
		{
			if(old_state == 'I')
			{
				for(int neighbour : Neighbours_All(person)) { Node(neighbour).Infected_Neighbour_Delete(person); }
			}
			else if(new_state == 'I')
			{
				for(int neighbour : Neighbours_All(person)){ Node(neighbour).Infected_Neighbour_Add(person); }
			}
			return;
		}

		template<typename T = NodeType> typename std::enable_if<std::is_same<T, C_SocialData>::value, void>::type
			State_Change_Adjust_Node_Neighbour_Lists(int person, const char old_state, const char new_state)
		{
			return;
		}

		void Nodes_In_State_Add(const int person, const char new_state)
		{
			_nodes_in_state[new_state].insert(person);
		}

		void Nodes_In_State_Delete(const int person, const char old_state)
		{
			_nodes_in_state[old_state].erase(person);
		}

	public:

		const int random_int(void) { return _rand_int(_layer_Mersenne); }

		const TPt<TNodeNet<NodeType>> The_Network(void) const
		{
			return _population;
		}

		template<typename T = NodeType> typename std::enable_if<std::is_same<T, C_SocialData>::value, const float>::type Opinion_Changes_Per_Person(const float proportion, const char desired_social_state = '\0') const
		{
			// solution adapted from jpihl's answer at
			// https://stackoverflow.com/questions/6972368/stdenable-if-to-conditionally-compile-a-member-function

			std::vector<int> all_nodes;
			std::vector<int> temp_rando_vec;

			if(proportion == 1)
			{
				if(desired_social_state != '\0') all_nodes = set_to_vector(Nodes_In_State(desired_social_state));
				else all_nodes = set_to_vector(Nodes_All());
			}
			else
			{
				all_nodes = Nodes_Shuffled();
				std::vector<int> temp_rando_vec(all_nodes.begin(), std::next(all_nodes.begin(), std::ceil(proportion*Node_Count())) );

				if(desired_social_state != '\0')
				{
					temp_rando_vec.erase(
						std::remove_if(
							temp_rando_vec.begin(),
							temp_rando_vec.end(),
							[=](int person) -> bool { return (Node_State_Get(person) != desired_social_state); }
						),
						temp_rando_vec.end()
					);
				}
				all_nodes = temp_rando_vec;
			}
			const float pre = std::accumulate
			(
				all_nodes.begin(),
				all_nodes.end(),
				0.,
				[&](double sum_so_far, auto new_node){ return (sum_so_far +  Node(new_node).Opinion_Changes_Count()); }
			);
			if(desired_social_state == '\0'){ return pre; }
			return pre * 1. / all_nodes.size();
		}

		const bool Edge_Add(const int person1, const int person2)
		{  // the network is undirected right now
			if(person1 == person2){ return false; }
			if( _population->IsEdge(person1, person2, false) ){ return false; }
			_population->AddEdge(person1, person2);
			_population->AddEdge(person2, person1);
			Edge_Add_Adjust_Node_Neighbour_Lists(person1, person2);
			return true;
		}

		const int Edge_Count(void) const
		{
			return _population->GetEdges()/2;
		}

		const bool Edge_Delete(const int person1, const int person2)
		{ // FUNCTION NOT TESTED
			if( (not _population->IsEdge(person1, person2)) and (not _population->IsEdge(person2, person1)) ){ return false; }
			if( _population->IsEdge(person1, person2) ) _population->DelEdge(person1, person2);
			if( _population->IsEdge(person2, person1) ) _population->DelEdge(person2, person1);
			Edge_Add_Adjust_Node_Neighbour_Lists(person1, person2);
			return true;
		}

		const std::set<std::set<int>> Edges_All(void) const
		{
			const int num_edges = _population->GetEdges();
			std::set<std::set<int>> edge_set = {};

			for (auto EI = _population->BegEI(); EI < _population->EndEI(); EI++)
			{
				edge_set.insert({ EI.GetSrcNId(), EI.GetDstNId() });
			}

			return edge_set;
		}

		const std::set<int> Neighbours_All(const int person) const
		{
			int count = Neighbour_Count(person);
			if(count == 0){ return {}; }
			std::set<int> temp;
			for(int i = {}; i < count; ++i){ temp.insert( _population->GetNI(person).GetInNId(i) ); }
			return temp;
		}

		const std::set<int> Nodes_All(void) const
		{
			TIntV temp_vec_ids;
			_population->GetNIdV(temp_vec_ids);
			std::set<int> return_vec(temp_vec_ids.BegI(), temp_vec_ids.EndI());
			return return_vec;
		}

		const int Neighbour_At(const int person, const int position) const
		{
			return _population->GetNI(person).GetInNId(position);
		}

		const int Neighbour_Count(const int person) const
		{
			return _population->GetNI(person).GetInDeg();
		}

		NodeType& Neighbour_Random(const int person) const
		{
			return _population->GetNDat(Neighbour_Random_ID(person));
		}

		const int Neighbour_Random_ID(const int person) const
		{
			return  Neighbour_At(person, h_fastrand(Neighbour_Count(person)));
		}

		NodeType& Node(const int person) const
		{
			return _population->GetNDat(person);
		}

		NodeType& Node_Random(void) const
		{
			return _population->GetNI( Node_Random_ID() );
		}

		const int Node_Random_ID(void) const
		{
			return h_fastrand( Node_Count() );
		}

		const bool Node_State_Change_To(const int person, const char new_state)
		{
			const char old_state = Node(person).State_Get();
			if(old_state == new_state){ return false; }
			State_Change_Adjust_Node_Neighbour_Lists(person, old_state, new_state);
			Node(person).State_Change_To(new_state);
			Nodes_In_State_Delete(person, old_state);
			Nodes_In_State_Add(person, new_state);
			return true;
		}

		const char Node_State_Get(const int person) const
		{
			return Node(person).State_Get();
		}

		const std::set<int> Nodes_In_State(char state) const
		{
			if(not _nodes_in_state.count(state)) return {};
			// if(state == '\0') return Nodes_All();
			if(state == '\0') return std::set<int>({std::numeric_limits<int>::quiet_NaN()});
			return _nodes_in_state.at(state);
		}

		const int Nodes_In_State_Count(char state) const
		{
			if(not _nodes_in_state.count(state)) return 0;
			return _nodes_in_state.at(state).size();
		}

		const int Node_In_State_Random(const char state) const
		{
			if(not Nodes_In_State_Count(state)) return std::numeric_limits<int>::quiet_NaN();
			return *random_member(_nodes_in_state.at(state).begin(), _nodes_in_state.at(state).end());
		}

		const std::vector<int> Nodes_Shuffled(void)
		{
			std::shuffle(_shuffled_ids.begin(), _shuffled_ids.end(), _layer_Mersenne);
			return _shuffled_ids;
		}

		const int Node_Count(void) const
		{
			return _population->GetNodes();
		}

		const float Cluster_Coeff(void) const
		{
			return TSnap::GetClustCf(_population);
		}

		const float Diameter(void) const
		{
			return TSnap::GetBfsFullDiam(_population, std::log(Node_Count()));
		}

		// Print network statistics
		// using the prettyprint library for the STL containers
		const void Statistics_Print(void) const
		{
			std::cout << std::endl;
			std::cout << "number of nodes: " << _population->GetNodes() << std::endl;
			std::cout << "number of edges: " << _population->GetEdges() << std::endl;

			std::cout << std::endl;

			std::cout << "clustering coefficient: " << Cluster_Coeff() << std::endl;
			std::cout << "average path length: " << Diameter() << std::endl;

			return;
		}

		const std::map<int, int> Degree_Distribution(void) const
		{
			// degree is mentioned first
			TVec<TPair<TInt, TInt>> temp;
			TSnap::GetDegCnt(_population, temp);
			std::map<int, int> temp_map;
			for(int i {}; i < temp.Len(); ++i){ temp_map[ temp[i].GetVal1() ] = temp[i].GetVal2(); }
			return temp_map;
		}

		const int Stat_Square_Sum_Degrees(void) const
		{
			TIntV in_degree_sequence;
			TIntV out_degree_seqduence;
			TSnap::GetDegSeqV(_population, in_degree_sequence, out_degree_seqduence);
			return std::accumulate
			(
				in_degree_sequence.BegI(),
				in_degree_sequence.EndI(),
				0,
				[&](int sum_so_far, auto new_degree){ return (sum_so_far + new_degree*new_degree); }
			);
		}

		void Reset(void)
		{
			for(int person : Nodes_All())
			{
				Node_State_Change_To(person, '\0');
				for(int neighbour : Neighbours_All(person))
				{
					Edge_Delete(person, neighbour);
					Edge_Delete(neighbour, person);
				}
			}
		}

		void Generate_Network_Random(const int number_of_edges)
		{
			Reset();
			const int number_of_nodes = Node_Count();
			const double ran_edge_prob = number_of_edges/(double)std::pow(number_of_nodes, 2);
	        // following the algorithm given in the paper "Efficent Generation of Large Random Networks"
	        int v = 1, w = -1;
	        while( v < number_of_nodes )
	        {
	            w += 1 + std::floor( std::log(1-h_fastfloat()) / std::log(1-ran_edge_prob) );
	            while( (w >= v) and (v < number_of_nodes) )
	            {
	                w -= v;
	                ++ v;
	            }
	            if(v < number_of_nodes)
	            {
					Edge_Add(v, w);
	            }
	        }
	        for(auto person : Nodes_All())
	        { // can use either the physical or social layers for the iteration of a single loop, since
	            if( Neighbour_Count(person) == 0 )
	            {
	                // randomly generate a new neighbour, going until we find one that doesn't create a loop
	                int new_neighbour;
	                do
	                { new_neighbour = Node_Random_ID(); }
	                while(new_neighbour == person);

	                Edge_Add(person, new_neighbour);
	            }
	        }
		    return;
		}

		void Generate_Network_Smallworld(const float rewiring_prob, const int avg_node_degree)
		{
			Reset();
			const int number_of_nodes = Node_Count();

			for(int person : Nodes_All()){
			for(int neighbour : range(person+1, person+(avg_node_degree/2)+1))
			{
				Edge_Add(person,  mod(neighbour, number_of_nodes));
			}}

			for(int person : Nodes_All())
			{
				std::set<int> neighbour_set = Neighbours_All(person);
				for(int old_neighbour : neighbour_set)
				{
					if(h_fastfloat() < rewiring_prob)
					{
						int new_neighbour = -1;
						do {
							new_neighbour = h_fastrand(number_of_nodes);
						} while( (new_neighbour==person) or (new_neighbour==old_neighbour) or neighbour_set.count(new_neighbour) );

						Edge_Delete(person, old_neighbour);
						Edge_Delete(old_neighbour, person);

						Edge_Add(person, new_neighbour);
						Edge_Add(new_neighbour, person);
					}
				}
			}
		}

		~C_Layer() { ; }
		// C_Layer(void){ ; }
		// C_Layer(void) = delete;
		C_Layer(void) : _layer_Mersenne(_layer_random_device()), _rand_int(std::uniform_int_distribution<>(1, 1))
		{
			_population = _population = TNodeNet<NodeType>::New();
			_shuffled_ids = {};
		}
		C_Layer& operator = (const C_Layer& other)
		{
			if(this == &other) return *this;
			_layer_Mersenne = other._layer_Mersenne;
			_rand_int = std::uniform_int_distribution<>(1, other.Node_Count());
			_layer_Mersenne.seed(::time(NULL));
			_nodes_in_state =  other._nodes_in_state;
			_population = other._population;
			_shuffled_ids = other._shuffled_ids;
			return *this;
		}
		C_Layer(int number_of_nodes) : _layer_Mersenne(_layer_random_device()), _rand_int(std::uniform_int_distribution<>(1, number_of_nodes))
		{
			_layer_Mersenne.seed(::time(NULL));
			_population = TNodeNet<NodeType>::New();
			for(int i = 0; i < number_of_nodes; ++i)
			{
				_population->AddNode(i, NodeType('\0'));
				_shuffled_ids.push_back(i);
				_nodes_in_state['\0'].insert(i);
			}
			std::shuffle(_shuffled_ids.begin(), _shuffled_ids.end(), _layer_Mersenne);
		}
		C_Layer(const C_Layer& other) : _layer_Mersenne(_layer_random_device()), _rand_int(std::uniform_int_distribution<>(1, other._shuffled_ids.size()))
		{
			_layer_Mersenne.seed(::time(NULL));
			_population = other._population;
			_shuffled_ids = other._shuffled_ids;
			_nodes_in_state = other._nodes_in_state;
			return;
		}
		C_Layer(C_Layer&& other) : _layer_Mersenne(other._layer_Mersenne), _rand_int(std::uniform_int_distribution<>(1, other.Node_Count()))
		{
			_layer_Mersenne.seed(::time(NULL));
			_nodes_in_state =  other._nodes_in_state;
			_population = other._population;
			_shuffled_ids = other._shuffled_ids;
			return;
		}
};

#endif
