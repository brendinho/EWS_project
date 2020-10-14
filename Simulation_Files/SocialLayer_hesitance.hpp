/*
	Brendon Phillips
	PhD Candidate
	Bauch computational epidemiology research group
	Department of Applied Mathematics
	Faculty of Mathematics
	Universiity of Waterloo
*/

#ifndef SOCIAL_LAYER_HEADER_
#define SOCIAL_LAYER_HEADER_

template <typename NodeType> class C_Connectivity_Layer : public C_Layer<NodeType>
{
	private:

		// std::map<std::multiset<char>, int>	_join_count; // for the virtual layer
		std::map<std::set<char>, int>	_join_count; // for the virtual layer


		// the real layer is the original structure of the networl
		C_Layer<NodeType>	_real_layer;
		// the virtual layer is a network where nodes are connected only if they share opinions
		TPt<TUNGraph>		_virtual_layer;

	public:

		// void Information_Add(const std::string characteristic, const std::string description)
		// {
		// 	_real_layer.Information_Add(characteristic, description);
		// }

		const int random_int(void)
		{
			return _real_layer.random_int();
		}

		const int Edge_Count(const std::string layer_name) const
		{
			if(layer_name == "real") return _real_layer.Edge_Count();
			if(layer_name == "virtual") return _virtual_layer->GetEdges();
			else return std::numeric_limits<int>::quiet_NaN();
		}

		// const std::string Information_Get(const std::string characteristic) const
		// { // FUNCTION NOT CHECKED
		// 	return _real_layer.Information_Get(characteristic);
		// }

		const std::set<int> Neighbours_All(const std::string layer_name, const int person) const
		{
			if(layer_name == "real") return _real_layer.Neighbours_All(person);
			if(layer_name == "virtual")
			{
				std::set<int> temp;
				const int count = Neighbour_Count("virtual", person);
				if(count == 0){ return {}; }
				for(int i = {}; i < count; ++i){ temp.insert( _virtual_layer->GetNI(person).GetInNId(i) ); }
				return temp;
			}
			return std::set<int>({});
		}

		const std::set<int> Nodes_All(void) const
		{
			return _real_layer.Nodes_All();
		}

		const int Neighbour_At(const std::string layer_name, const int person, const int position) const
		{
			if(layer_name == "real") return _real_layer.Neighbour_At(person, position);
			if(layer_name == "virtual") return _virtual_layer->GetNI(person).GetInNId(position);
			return std::numeric_limits<int>::quiet_NaN();
		}

		const int Neighbour_Count(const std::string layer_name, const int person) const
		{
			if(layer_name == "real") return _real_layer.Neighbour_Count(person);
			if(layer_name == "virtual") return _virtual_layer->GetNI(person).GetDeg();
			return std::numeric_limits<int>::quiet_NaN();
		}

		NodeType& Neighbour_Random_Node_Info(const int person) const
		{
			return _real_layer.Neighbour_Random(person);
		}

		const int Neighbour_Random_ID(const std::string layer_name, const int person) const
		{
			if(layer_name == "real") return _real_layer.Neighbour_Random_ID(person);
			if(layer_name == "virtual") return Neighbour_At(layer_name, person, h_fastrand(Neighbour_Count("virtual", person)));
			return std::numeric_limits<int>::quiet_NaN();
		}

		NodeType& Node(const int person) const
		{
			return _real_layer.Node(person);
		}

		NodeType& Node_Random(void) const
		{
			return _real_layer.Node_Random();
		}

		const int Node_Random_ID(void) const
		{
			return _real_layer.Node_Random_ID();
		}

		const char Node_State_Get(const int person)
		{
			return _real_layer.Node_State_Get(person);
		}

		const std::set<int> Nodes_In_State(char state) const
		{
			return _real_layer.Nodes_In_State(state);
		}

		const int Nodes_In_State_Count(char state) const
		{
			return _real_layer.Nodes_In_State_Count(state);
		}

		const int Node_In_State_Random(const char state) const
		{
			return _real_layer.Node_In_State_Random(state);
		}

		const std::vector<int> Nodes_Shuffled(void)
		{
			return _real_layer.Nodes_Shuffled();
		}

		const int Node_Count(void) const
		{
			return _real_layer.Node_Count();
		}

		const int Square_Sum_Degrees(const std::string layer_name) const
		{
			if(layer_name == "real") {return _real_layer.Stat_Square_Sum_Degrees(); }
			TIntV in_degree_sequence;
			TIntV out_degree_seqduence;
			TSnap::GetDegSeqV(_virtual_layer, in_degree_sequence, out_degree_seqduence);
			return std::accumulate
			(
				in_degree_sequence.BegI(),
				in_degree_sequence.EndI(),
				0,
				[=](int sum_so_far, auto new_degree){ return (sum_so_far + new_degree*new_degree); }
			);
		}

		const bool Edge_Add(const int person1, const int person2)
		{ // FIX
			if(person1 == person2){ return false; }
			if( not _real_layer.Edge_Add(person1, person2) ){ return false; }
			++ _join_count[{ Node_State_Get(person1), Node_State_Get(person2) }];
			if( Node_State_Get(person1) == Node_State_Get(person2) )
			{
				_virtual_layer->AddEdge(person1, person2);
				_virtual_layer->AddEdge(person2, person1);
			}
			return true;
		}

		const bool Edge_Delete(const int person1, const int person2)
		{ // FUNCTION NOT TESTED
			if(person1 == person2){ return false; }
			if( not _real_layer.Delete_Edge(person1, person2) ){ return false; }
			_virtual_layer->DelEdge(person1, person2);
			_virtual_layer->DelEdge(person2, person1);
			-- _join_count[{ Node_State_Get(person1), Node_State_Get(person2) }];
			return true;
		}

		const bool Node_State_Change_To(const int person, const char new_state)
		{
			auto old_state = Node_State_Get(person);
			if( not _real_layer.Node_State_Change_To(person, new_state) ){ return false; }
			for(auto nbr : Neighbours_All("real", person))
			{
				--  _join_count[{old_state, Node_State_Get(nbr)}];
				++  _join_count[{new_state, Node_State_Get(nbr)}];
				if( Node_State_Get(person) != Node_State_Get(nbr) )
				{ // then it was previously a neighbour
					_virtual_layer->DelEdge(person, nbr);
				}
				else
				{
					_virtual_layer->AddEdge(person, nbr);
				}
			}
			return true;
		}

		const int Join_Count(const float proportion, const std::set<char> state_pair) const
		{
			assert( check_layer("social",*state_pair.begin()) );
			assert( check_layer("social",*next(state_pair.begin(),1)) );
			if( not _join_count.count(state_pair) ){ return 0; }
			if(proportion == 1) return _join_count.at(state_pair);

			/*
				simple check. get the quota of ranmly selected nodes. if we're looking for a singe-state join count, then filter for only nodes in that state. then use this list of random nodes to retrieve a subgraph of the virual sentiment network and count the edges there.
			*/
			std::vector<int> temp = Nodes_Shuffled();
			std::vector<int> randos;

			if( state_pair.size() == 1)
			{
				// if looking for the NN or VV join count, we only need nodes in that one state. so, filler out the nodes not in that state
				std::copy_if(
					temp.begin(),
					std::next(temp.begin(), std::ceil(proportion*Node_Count())),
					std::back_inserter(randos),
					[=](int person){ return (Node_State_Get(person) == *state_pair.begin()); }
				);
			}
			else{ randos = temp; }

			// much faster than a vaive loop, thanks be to the gods of Stanford
			TPt<TUNGraph> subnet = TSnap::GetSubGraph(_virtual_layer, to_TVec(randos));
			return subnet->GetEdges();
		}

		const float Opinion_Changes_Per_Person(const float proportion, const char social_state = '\0') const
		{
			return _real_layer.Opinion_Changes_Per_Person(proportion, social_state);
		}

		const std::map<int, int> Degree_Distribution(const std::string layer_name) const
		{
			// degree is mentioned first
			if(layer_name == "real") { return _real_layer.Degree_Distribution(); }
			TVec<TPair<TInt, TInt>> temp;
			TSnap::GetDegCnt(_virtual_layer, temp);
			std::map<int, int> temp_map;
			for(int i {}; i < temp.Len(); ++i){ temp_map[ temp[i].GetVal1() ] = temp[i].GetVal2(); }
			return temp_map;
		}

		const float Diameter(const float proportion, const char desired_social_state = '\0') const
		{
			std::set<int> nodes_for_subgraph;
			// if we're using all the nodes, then just extract the components from the virtual network and be done with it
			if(proportion == 1)
			{
				if(desired_social_state == '\0') TSnap::GetBfsEffDiam(_virtual_layer, std::log(Node_Count()));
				else nodes_for_subgraph = Nodes_In_State(desired_social_state);
			}
			else
			{
				std::vector<int> temp = Nodes_Shuffled();
				if(desired_social_state == '\0')
				{
					std::copy(
						temp.begin(),
						std::next( temp.begin(), std::ceil(proportion*Node_Count()) ),
						std::inserter(nodes_for_subgraph, nodes_for_subgraph.begin())
					);
				}
				else
				{
					std::copy_if(
						temp.begin(),
						std::next( temp.begin(), std::ceil(proportion*Node_Count()) ),
						std::inserter(nodes_for_subgraph, nodes_for_subgraph.begin()),
						[=](int person) -> bool { return(Node_State_Get(person) == desired_social_state); }
					);
				}
			}

			TPt<TUNGraph> the_subgraph = TSnap::GetSubGraph(_virtual_layer, to_TVec(nodes_for_subgraph));
			if(nodes_for_subgraph.size() <= 2) return std::numeric_limits<float>::quiet_NaN();
			return TSnap::GetBfsFullDiam(the_subgraph, std::max(std::log(nodes_for_subgraph.size()), 1.));
		}

		const std::map<int, int> Connected_Component_Sizes(const float proportion, const char desired_state = '\0') const
		{
			// return variable to return the size and numbers of the components
			TVec<TPair<TInt, TInt>> component_size_and_count;
			TPt<TUNGraph> the_subgraph;

			// if we're using all the nodes, then just extract the components from the virtual network and be done with it
			if(proportion == 1)
			{
				// if looking at the total numbers, we don't even need a subgraph
				if(desired_state == '\0') the_subgraph = _virtual_layer;
				// create a subgraph of only nodes in the desired state
				else the_subgraph = TSnap::GetSubGraph(_virtual_layer, to_TVec(Nodes_In_State(desired_state)));
			}
			else
			{
				/*
					there are a couple tasks here. first get a number of random nodes (from the real graph irrespective of opinion). then filter them to find which of them for the opinion we want. form a subgraph with these compliant randos and then find the connected components in this subgraph. report them.
				*/

				// get the shuffled nodes and select the proportion we're looking for
				std::vector<int> temp = Nodes_Shuffled();
				TVec<TInt> nodes_for_subgraph;

				for(int i = 0; i <= std::ceil(proportion*Node_Count()); ++i)
				{
					// if state is not given, we'll take it
					if(desired_state == '\0') nodes_for_subgraph.Add(temp[i]);
					// if a state is given and the person agrees with the required state, then we'll take it
					else if(Node_State_Get(temp[i]) == desired_state) nodes_for_subgraph.Add(temp[i]);
					// in all other situations, we don't know her
				}
				// make a subgraph with these compliant nodes
				the_subgraph = TSnap::GetSubGraph(_virtual_layer, nodes_for_subgraph);
			}

			// get the connected components of this compliant subgraph
			TSnap::GetWccSzCnt(the_subgraph, component_size_and_count);
			std::map<int, int> component_sizes;
			for(int i = 0; i < component_size_and_count.Len(); ++i)
			{
				component_sizes[ component_size_and_count[i].GetVal1() ] = component_size_and_count[i].GetVal2();
			}
			return component_sizes;
		}

		const std::map<int, int> Echo_Chamber_Sizes(const float proportion, const char desired_state = '\0') const
		{
			auto chambered = [=](std::set<int> people) -> std::set<int>
			{
				/* get the chambered nodes in the given list, which are nodes that have no disagreeing neighbours */
				std::set<int> accepted_nodes;
				for(int person : people)
				{
					bool in_a_chamber = true;
					/*
						NO NOT USE THE VIRTUAL NETWORK
						If you only consider agreeing neighbours, you've intentionally lost all information about how many disagreeing neighbours they had and the problem becomes trivial. Use the real network and check all their neighbours.
					*/
					for(int neighbour : Neighbours_All("real", person))
					{
						if(Node_State_Get(person) != Node_State_Get(neighbour))
						{
							in_a_chamber = false;
							break;
						}
					}
					if( in_a_chamber )
					{
						accepted_nodes.insert(person);
					}
				}
				return accepted_nodes;
			};

			// variable to store the sizes and numbers of the chambers
			TVec<TPair<TInt, TInt>> chamber_size_and_count;

			// store the nodes that are in a chamber
			std::set<int> chambered_nodes;

			if(proportion == 1)
			{
				/*
					nothing special to be done here if using all the nodes. get the nodes in the desired state, filter them to get the chambered nodes. get the subgraph with these nodes, and then find the conencted components.
				*/

				// if we're interested in the whole network, then just get all the nodes
				if(desired_state == '\0') chambered_nodes = chambered(Nodes_All());
				// get the nodes in the state, and filter them to get the chambered nodes
				else chambered_nodes = chambered(Nodes_In_State(desired_state));
			}
			else
			{
				/*
					Same deal as before, but usiong a random subset. CHoose some random nodes from the entire set, find the ones in the desired state, filter them to find the chambered ones, create a subnoetwork with those nodes, and then find the connected components. piece of piss.
				*/

				// get some random nodes
				std::vector<int> temp = Nodes_Shuffled();
				// container for the random variables
				std::set<int> randos;

				// if we're just looking for a total, just get the first n nodes
				if(desired_state == '\0') std::copy_n(temp.begin(), std::ceil(proportion*Node_Count()), std::inserter(randos, randos.begin()));
				// else get the random nodes with opinion specified in the function call
				else std::copy_if(
						temp.begin(),
						std::next( temp.begin(), std::ceil(proportion*Node_Count()) ),
						std::inserter(randos, randos.begin()),
						[=](int person){return (Node_State_Get(person) == desired_state);}
					);

				// get the chambered nodes from the randos
				chambered_nodes = chambered(randos);
			}

			/*
				get the subgraph of the chambered nodes formed from the VIRTUAL GRAPHS. Since we know which nodes are chambered, the only connections we need now are those between asgreeing friends
			*/
			// TVec<TInt> snap_chambereds; for(int trapped : chambered_nodes) snap_chambereds.Add(trapped);
			TPt<TUNGraph> subgraph_for_this_state = TSnap::GetSubGraph(_virtual_layer, to_TVec(chambered_nodes));

			// find the connected components of this chambered subgraph
			TSnap::GetWccSzCnt(subgraph_for_this_state, chamber_size_and_count);

			// self-explanatory
			std::map<int, int> component_sizes;
			for(int i {}; i < chamber_size_and_count.Len(); ++i)
			{
				component_sizes[ chamber_size_and_count[i].GetVal1() ] = chamber_size_and_count[i].GetVal2();
			}
			return component_sizes;
		}

		const double Modularity_Score(const float proportion, const char desired_state = '\0') const
		{
			std::set<int> nodes_for_subgraph;
			if(proportion == 1)
			{
				if(desired_state == '\0') nodes_for_subgraph = Nodes_All();
				else nodes_for_subgraph = Nodes_In_State(desired_state);
			}
			else
			{
				std::vector<int> temp = Nodes_Shuffled();

				if(desired_state == '\0') std::copy_n(temp.begin(), std::ceil(proportion*Node_Count()), std::inserter(nodes_for_subgraph, nodes_for_subgraph.begin()));
				else std::copy_if(
					temp.begin(),
					std::next(temp.begin(), std::ceil(proportion*Node_Count())),
					std::inserter(nodes_for_subgraph, nodes_for_subgraph.begin()),
					[=](int person){ return(Node_State_Get(person) == desired_state); }
				);
			}
			return TSnap::GetModularity(_virtual_layer, to_TVec(nodes_for_subgraph), Edge_Count("virtual"));
		}

		const int Number_Triads(const float proportion, const char desired_state = '\0') const
		{

			std::vector<int> nodes_for_subgraph;
			if(proportion == 1)
			{
				if(desired_state == '\0') return TSnap::GetTriads(_virtual_layer);
				else nodes_for_subgraph = set_to_vector(Nodes_In_State(desired_state));
			}
			else
			{
				const std::vector<int> temp = Nodes_Shuffled();

				if(desired_state == '\0') nodes_for_subgraph = temp;
				else std::copy_if(
					temp.begin(),
					std::next(temp.begin(), std::ceil(proportion*Node_Count())),
					std::back_inserter(nodes_for_subgraph),
					[=](int person){ return(Node_State_Get(person) == desired_state); }
				);
			}

			TPt<TUNGraph> subgraph_for_this_state = TSnap::GetSubGraph(_virtual_layer, to_TVec(nodes_for_subgraph));
			return TSnap::GetTriads(subgraph_for_this_state);
		}

		const double Watts_Strogatz_Avg_Clustering_Coefficient(const float proportion, const char desired_state = '\0') const
		{

			std::set<int> nodes_for_subgraph;
			if(proportion == 1)
			{
				if(desired_state == '\0') return TSnap::GetClustCf(_virtual_layer);
				else nodes_for_subgraph = Nodes_In_State(desired_state);
			}
			else
			{
				std::vector<int> temp = Nodes_Shuffled();

				if(desired_state == '\0') std::copy_n(temp.begin(), std::ceil(proportion*Node_Count()), std::inserter(nodes_for_subgraph, nodes_for_subgraph.begin()));
				else std::copy_if(
					temp.begin(),
					std::next(temp.begin(), std::ceil(proportion*Node_Count())),
					std::inserter(nodes_for_subgraph, nodes_for_subgraph.begin()),
					[=](int person){ return(Node_State_Get(person) == desired_state); }
				);
			}

			TPt<TUNGraph> the_subgraph = TSnap::GetSubGraph(_virtual_layer, to_TVec(nodes_for_subgraph));
			return TSnap::GetClustCf(the_subgraph);
		}

		C_Connectivity_Layer& operator = (const C_Connectivity_Layer& other)
		{
			_join_count = other._join_count;
			_virtual_layer = other._virtual_layer;
			_real_layer = other._real_layer;
		}
		// C_Connectivity_Layer(void) = delete;
		C_Connectivity_Layer(void){ ; }
		C_Connectivity_Layer(int number_of_nodes) // : _real_layer(C_Layer<NodeType>(number_of_nodes, default_state))
		{
			_virtual_layer = TUNGraph::New();
			for(int i = 0; i < number_of_nodes; ++i)
			{
				_virtual_layer->AddNode(i);
			}
			_real_layer = C_Layer<NodeType>(number_of_nodes);
		}
		C_Connectivity_Layer(const C_Connectivity_Layer& other)
		{
			_join_count = other._join_count;
			_virtual_layer = other._virtual_layer;
			_real_layer = other._real_layer;
		}

};

#endif
