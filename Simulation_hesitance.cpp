/*
	Brendon Phillips
	PhD Candidate
	Bauch computational epidemiology research group
	Department of Applied Mathematics
	Faculty of Mathematics
	Universiity of Waterloo
*/

#include "Network_hesitance.hpp"

int main(int argc, char* argv[])
{
	// trial values for testing the code without command line arguments
	/*
	const std::string local_network_structure = "smallworld";
	const int local_network_size = 10000;

	const int local_average_node_degree = 50;
	const int instance = 4;
	const double local_infection_probability = 0.2;

	const double local_perceived_vaccine_risk = -0.0125;
	const double local_replenishment_ratio = 2.4e-4;
	const double local_case_importation_ratio = 2.5e-5;
	const double local_social_norm = 0.125;

	const double local_initial_vaccinator_proportion = 0.75; // 0.05;
	const double local_random_opinion_switch_probability = 0.0001;

	const double local_proportion_to_be_measured = 0.6;

	const int max_hours = 0;
	const int max_minutes = 5;
	*/

	// /*
	const std::string local_network_structure = argv[1];
	std::cout << "network structure: " << local_network_structure << std::endl;

	const int local_network_size = atoi(argv[2]);
	std::cout << "network size: " << local_network_size << std::endl;

	const int local_average_node_degree = atoi(argv[3]);
	std::cout << "average neighbourhood size: " << local_average_node_degree << std::endl;

	const int instance = atoi(argv[4]);
	std::cout << "instance: " << instance << std::endl;

	const double local_infection_probability = atof(argv[5]);
	std::cout << "infection probability: " << local_infection_probability << std::endl;

    const double local_perceived_vaccine_risk = atof(argv[6]);
	std::cout << "perceived vaccine risk: " << local_perceived_vaccine_risk << std::endl;

	const double local_replenishment_ratio = atof(argv[7]);
	std::cout << "replenishment rate: " << local_replenishment_ratio << std::endl;

	const double local_case_importation_ratio = atof(argv[8]);
	std::cout << "case importation ratio: " << local_case_importation_ratio << std::endl;

	const double local_social_norm = atof(argv[9]);
	std::cout << "social norm: " << local_social_norm << std::endl;

	const double local_initial_vaccinator_proportion = atof(argv[10]);
	std::cout << "initial vaccinator proportion: " << local_initial_vaccinator_proportion << std::endl;

	const double local_random_opinion_switch_probability = atof(argv[11]);
	std::cout << "random opinion switch prob: " << local_random_opinion_switch_probability << std::endl;

	const double local_proportion_to_be_measured = atof(argv[12]);
	std::cout << "proportion of nodes: " << local_proportion_to_be_measured << std::endl;

	const int max_hours = atoi(argv[13]);
	const int max_minutes = atoi(argv[14]);
	std::cout << "time limit: " << max_hours << " hours, " << max_minutes << " minutes" << std::endl;

	char szPath[128] = "";
    gethostname(szPath, sizeof(szPath));

	const std::string host = szPath;
	std::string stem = "";

	std::cerr << "host name is " << host << std::endl;

	if(host.find("gra") != std::string::npos)
	{
		stem = "home/b2philli/scratch";
	}
	else if((host.find("hagrid") != std::string::npos) or (host.find("hpc") != std::string::npos) or (host.find("pr1") != std::string::npos))
	{
		stem = "work/b2philli/project/b2philli";
	}
	else
	{
		std::cerr << "no folder name found, host name is " << host << "."<< std::endl;
		exit(EXIT_FAILURE);
	}

	const std::string local_output_folder_name = ("/" + stem);

	// std::cout << "local output " << local_output_folder_name << std::endl; return 0;

	// */
	if(
		(instance < 0) or
		(local_infection_probability < 0) or
		(local_infection_probability > 1) or
		(local_replenishment_ratio < 0) or
		(local_replenishment_ratio > 1) or
		(local_case_importation_ratio < 0) or
		(local_case_importation_ratio > 1) or
		(local_initial_vaccinator_proportion < 0) or
		(local_initial_vaccinator_proportion > 1) or
		(local_random_opinion_switch_probability < 0) or
		(local_random_opinion_switch_probability > 1) or
		(local_proportion_to_be_measured < 0) or
		(local_proportion_to_be_measured > 1) or
		((max_hours == 0) and (max_minutes == 0)) or
		(max_hours < 0) or
		(max_minutes < 0)
	)
	{
		std::cout << "Parameter value error..." << std::endl;
		exit(EXIT_FAILURE);
	}

    std::random_device rd;
    std::mt19937_64 Mersenne(rd());
    std::uniform_real_distribution<double> random(0,1);

    std::stringstream local_perceived_vaccine_risk_stream;
    std::stringstream infection_probability_stream;
    std::stringstream replenishment_stream;
    std::stringstream importation_stream;
    std::stringstream social_norm_stream;
    std::stringstream init_vacc_prop_stream;
    std::stringstream random_switch_prob_stream;

    local_perceived_vaccine_risk_stream	<< std::fixed << std::setprecision(4) 	<< local_perceived_vaccine_risk;
    infection_probability_stream   		<< std::setprecision(5) 				<< local_infection_probability;
    importation_stream             		<< std::setprecision(5) 				<< local_case_importation_ratio;
    replenishment_stream           		<< std::setprecision(5) 				<< local_replenishment_ratio;
    social_norm_stream             		<< std::setprecision(5) 				<< local_social_norm;
    init_vacc_prop_stream          		<< std::fixed << std::setprecision(2) 	<< local_initial_vaccinator_proportion;
    random_switch_prob_stream      		<< std::fixed << std::setprecision(4) 	<< local_random_opinion_switch_probability;

    std::string name_addition = "";

    std::string file_name_local_stem =
        name_addition   +
        "_N_"           + std::to_string(local_network_size) +
        "_dur_"         + std::to_string(global_infection_duration) +
        "_vaccprop_"    + init_vacc_prop_stream.str() +
        "_pay_"         + local_perceived_vaccine_risk_stream.str() +
        "_inf_"         + infection_probability_stream.str() +
        "_imp_"         + importation_stream.str() +
        "_rep_"         + replenishment_stream.str() +
        "_struct_"      + local_network_structure +
        "_deg_"        	+ std::to_string( local_average_node_degree ) +
        "_norm_"        + social_norm_stream.str() +
        "_switch_"		+ random_switch_prob_stream.str() +
		"_proportion_"	+ std::to_string(local_proportion_to_be_measured) +
        "_inst_"        + std::to_string(instance)
    ;

    std::string data_file_name_stub = local_output_folder_name + "/hes" + file_name_local_stem;

	std::cout << "Preliminary data file name: " << data_file_name_stub << std::endl << std::endl;

	std::ofstream flag_file;
	if( num_files_matching(local_output_folder_name, file_name_local_stem) )
	{
    	std::cout << "These parameter values have already been run." << std::endl;
    	exit(EXIT_SUCCESS);
	}
	flag_file.open(data_file_name_stub);
	if( not flag_file.is_open() )
	{
		std::cout << "Flag file could not be created." << std::endl;
    	exit(EXIT_SUCCESS);
	}
	flag_file.close();

	const int max_seconds = 3600*max_hours + 60*max_minutes;

	TIC(-1); // timing the construction of the network
	C_Network City
    (
		local_network_structure,
        local_network_size,
		local_average_node_degree,
        local_infection_probability,
		local_perceived_vaccine_risk,
        local_replenishment_ratio,
        local_case_importation_ratio,
        local_social_norm,
		local_random_opinion_switch_probability,
		local_initial_vaccinator_proportion,
		instance
    );
	try{ std::cout << "Network generation time: " << TOC(-1) << " secs." << std::endl; } TIME_CATCH

	TIC(0); // time the entire program

	bool reached_equilib = false;

    auto patient_zero = City.Random_Node_In_State("physical", 'S');

    City.Change_State("physical", patient_zero, 'I');

    City.Record_Time_Step(0, local_proportion_to_be_measured); // since patient zero is infected

    for(int t = 1; t <= global_time_limit; ++t)
    {
        std::cout << "time step " << t << std::endl;

        City.Increment_Infected_Time();

		std:: cout << "\timporting people" << std::endl;

        // imports a percentage of cases by infecting susceptibles
        if(global_case_importation_interval != 0){
		if(mod(t, global_case_importation_interval) == 0){
		if(t > global_case_importation_delay){
		for(auto person : City.Nodes_In_State("physical", 'S'))
		{
			if(random(Mersenne) < local_case_importation_ratio)
			{
			    City.Change_State("physical", person, 'I');
			}
		}}}}

		std:: cout << "\tinfecting people" << std::endl;

        // shuffles the agents so that the infection spread appears random
        std::set<int> temp_susceptibles = City.Nodes_In_State("physical", 'S');

        std::vector<int> unprocessed_susceptibles(temp_susceptibles.begin(), temp_susceptibles.end());
        temp_susceptibles.clear();
        std::shuffle(unprocessed_susceptibles.begin(), unprocessed_susceptibles.end(), Mersenne);

        // read the shuffled list in order, as an alternative to true random sampling without replacement
        for(int susceptible : unprocessed_susceptibles){
        for(int infected_neighbour : City.Node_Infected_Neighbours(susceptible))
        {
            if(City.Node_Is_Infectious(infected_neighbour)){
            if(random(Mersenne) < local_infection_probability)
            {
                City.Change_State("physical", susceptible, 'I');
                break;
            }}
        }}

		unprocessed_susceptibles.clear();

		std:: cout << "\topinion sampling and change" << std::endl;

	    // nothing fancy, just so it's clear what the fuck is going on with this shit
	    for(auto person : City.Shuffled_Nodes("social"))
		{
			const int neighbour = City.Neighbour_Random_ID("social", person);

		    if( City.State("social", person) != City.State("social", neighbour) )
	        {
				// if talking to a straw man, nothing happens so move on without doing anything
				if( City.Is_Hesitant(neighbour) ) continue;

				// if the magical probability function says that we have to reconsider
				if(random(Mersenne) < City.Sentiment_Change_Probability(person, neighbour))
				{
					// a clash of two extreme views
						if( 	City.Is_Vaccinator(person) 	and City.Is_Evil(neighbour) ) 		City.Change_State("social", person, 'H');
					else if( 	City.Is_Evil(person) 		and City.Is_Vaccinator(neighbour) ) City.Change_State("social", person, 'H');
					// straw man vulnerable to any sentiment expressed by neighbours
					else if(	City.Is_Hesitant(person)	and City.Is_Evil(neighbour) )		City.Change_State("social", person, 'N');
					else if(	City.Is_Hesitant(person)	and City.Is_Vaccinator(neighbour) )	City.Change_State("social", person, 'V');
					// no rules written for hesitant neighbour, since it doesn't do anything, and the continue statement should take care of that
	            }
	        }
		}

		std:: cout << "\topinion random switch" << std::endl;

        // people switch their opinions randomly without stimulus
        if( mod(t, global_random_opinion_switch_interval) == 0 ){
		for(auto person : City.Nodes("social")){
        if( random(Mersenne) < local_random_opinion_switch_probability )
        {
				 if( City.Is_Vaccinator(person) )	City.Change_State("social", person, 'H');
			else if( City.Is_Evil(person) ) 		City.Change_State("social", person, 'H');
			else if( City.Is_Hesitant(person) )
			{
				const char new_opinion = (random(Mersenne) < local_initial_vaccinator_proportion) ? 'V' : 'N';
				City.Change_State("social", person, new_opinion);
			}
        }}}

		std:: cout << "\tkilling and birthing people" << std::endl;

        // we kill and birth people at the end of the week
		// it's possible for infected people to die straightway
        if( mod(t, global_replenishment_interval) == 0)
        { // working with expectation here
            const std::vector<int> alive_people  = City.Shuffled_Nodes("physical");
			const float kill_quota = std::ceil( local_replenishment_ratio * City.Size() );

			std::for_each(
				alive_people.begin(),
				std::next(alive_people.begin(), std::floor(kill_quota*local_initial_vaccinator_proportion)),
				[&](int sensible){
					// reborn as vaccinated pro-vaccine citizen
					City.Change_State("social",   sensible, 'V', true);
					City.Change_State("physical", sensible, 'V', true);
				}
			);

			std::for_each(
				std::next(alive_people.begin(), std::floor(kill_quota*local_initial_vaccinator_proportion)),
				std::next(alive_people.begin(), kill_quota),
				[&](int asshole){
					// Jenny McCarthy, Wakefield, religious exemption or just some other run_of_the_mill fuckery messing with their head
					City.Change_State("social",   asshole, 'H', true);
					City.Change_State("physical", asshole, 'S', true);
				}
			);
        }

		reached_equilib = City.Record_Time_Step(t, local_proportion_to_be_measured);
		try{ std::cout << "total time elapsed" << std::endl; TOC(0, true); } TIME_CATCH

		if( reached_equilib or (TOC(0) >= max_seconds) ) { break; }

    } // END THE TIME LOOP

	try
	{
	    std::cout << std::endl << "Duration of the simulation." << std::endl;
		TOC(0, true);
	    std::cout << std::endl << std::endl;
	}
	TIME_CATCH

    remove(data_file_name_stub.c_str());
    City.Print_Data_File(local_output_folder_name, file_name_local_stem, reached_equilib);

    exit(EXIT_SUCCESS);
}
