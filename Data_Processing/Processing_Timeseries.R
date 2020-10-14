# Brendon Philips

# Time series generator

rm(list=ls());

user_name <- Sys.info()[8][[1]];
source(sprintf("/media/%s/Simulations/Processing/Parameter_Values.R", user_name));

# source("Processing_Parameter_Table.R");
Params <- fread(parameter_file, colClasses="character");

Parameters <- Params[num(beta)==1 & num(size)==40000 & num(infec)==0.2]

start_time <-  Sys.time();

# cl <- parallel::makeCluster( detectCores()[1]-1 ); #not to overload your computer
# doParallel::registerDoParallel(cl);

foreach( line_index=1:nrow(Parameters), .packages=c("data.table")) %dopar%
{
	this_tuple <- Parameters[line_index];

	beta <- this_tuple[, beta];
	initial_vacc_prop <- this_tuple[, init_prop];
	size <- this_tuple[, size];
	duration <- this_tuple[, duration];
	risk <- this_tuple[, risk];
	infection_prob <- this_tuple[, infec];
	importation_rate <- this_tuple[, import];
	replenishment_rate <-  this_tuple[, birth_death];
	physical_topology <- this_tuple[, phys_top];
	physical_degree <- this_tuple[, phys_deg];
	social_topology <- this_tuple[, soc_top];
	social_degree <- this_tuple[, soc_deg];
	random_opinion_switch <- this_tuple[, random_switch];
	social_norm <- this_tuple[, norm];

	file_end_line <- "_inst*bin";

	output_file_name <- paste(
		CSV_file_path,
		"Time_Series",
		"_N_", size,
		"_dur_", duration,
		"_beta_", beta,
		"_vaccprop_", initial_vacc_prop,
		"_risk_", risk,
		"_inf_", infection_prob,
		"_imp_", importation_rate,
		"_rep_", replenishment_rate,
		"_ptop_", physical_topology,
		"_pdeg_", physical_degree,
		"_stop_", social_topology,
		"_sdeg_", social_degree,
		"_norm_", social_norm,
		"_switch_", random_opinion_switch,
		".csv",
		sep=""
	);

	# all the files for this trial, across all the payoffs
	data_file_list <- Sys.glob(
		paste(
			# data_file_path,
			# "stats_cleaned",
			"/home/b2philli/Desktop/Files/",
			"stats*",
			"_N_", size,
			"_dur_", duration,
			"_beta_", beta,
			"_vaccprop_", initial_vacc_prop,
			"_pay_", risk,
			"_inf_", infection_prob,
			"_imp_", importation_rate,
			"_rep_", replenishment_rate,
			"_ptop_", physical_topology,
			"_pdeg_", physical_degree,
			"_stop_", social_topology,
			"_sdeg_", social_degree,
			"_norm_", social_norm,
			"_switch_", random_opinion_switch,
			"_inst*bin",
			sep=""
		)
	);

	if(length(data_file_list) == 0) stop();
	# if(length(data_file_list) == 0) next;

	# read first file in glob vector, get the column titles, and process it like a regular file
	if(file.size(data_file_list[1]) == 0){ next; }
	columns <- names(fread(data_file_list[1]));

	averaged_column_names <-  c(setdiff(columns, c(fixed_values, "V1", "time", "imported_cases", "secondary_infections", "sum_squares_degrees")), "ratio_sec_to_all_infs_full_sim", "ratio_sec_to_all_infs_at_the_end");

	mean_names <- as.vector(sapply(
		averaged_column_names,
		function(x) paste(x, "_mean", sep="")
	));

	sd_names <- as.vector(sapply(
		averaged_column_names,
		function(x) paste(x, "_sd", sep="")
	));

	DT <- data.table(matrix(
		0.0,
		nrow=0,
		ncol=(length(columns))
	));
	names(DT) <- columns;

	for(file_index in seq_along(data_file_list))
	{
		Data <- fread(data_file_list[file_index]);
		Data <- Data[, V1:=NULL];

		if( ("secondary_infections" %in% names(Data)) & ("imported_cases" %in% names(Data)) )
		{
			# this piece of code only calculates the ratio for the last n time steps of the simulation
			# we've done this so that the code below can be used
			number_of_zero_padding_rows <- nrow(Data)-use_this_many_entries_for_ratio_sec_new;

			secondaries <- c( rep(0, number_of_zero_padding_rows), tail(Data[, secondary_infections], use_this_many_entries_for_ratio_sec_new));
			importeds <-   c( rep(0, number_of_zero_padding_rows), tail(Data[, imported_cases      ], use_this_many_entries_for_ratio_sec_new));

			# can also be used for the entire table, with a "by=instance" argument, and [[2]] for the cumulative_*
			# Data[, "cases_imported"	:=	ceiling(importation*head( c(0, phys_S), -1)) ];
			Data[, "ratio_sec_to_all_infs_at_the_end" := cumsum(secondaries)/cumsum(secondaries + importeds) ];
			Data[, "ratio_sec_to_all_infs_full_sim"   := cumsum(secondary_infections)/cumsum(secondary_infections + imported_cases) ];
		}
		else
		{
			averaged_column_names <- averaged_column_names[ ! averaged_column_names %in% c("ratio_sec_to_all_infs_full_sim", "ratio_sec_to_all_infs_at_the_end")]
		}

		DT <- rbind(DT, Data, fill=TRUE);

		rm(Data, data_file_list);
	}

	Time_Series <- data.table(matrix(0.0, nrow=(DT[, max(time)]+1), ncol=length(c(mean_names, sd_names)) ));
	names(Time_Series) <- c(mean_names, sd_names);

	for(col in averaged_column_names)
	{
		# writeLines(paste("\t\t\tcolumn:", col));
		Time_Series[, paste(col, "_mean", sep="") := DT[, mean(get(col), na.rm=TRUE), by=time][, 2] ];
		Time_Series[, paste(col, "_sd",   sep="") := DT[, sd(get(col), na.rm=TRUE),   by=time][, 2] ];
	}

	Time_Series[, "time":=list(0:max(DT$time)) ];
	Time_Series[, (reinserted_values) := DT[, unique(.SD), .SDcols=reinserted_values]];

	write.csv( Time_Series, file = output_file_name);

	rm(Time_Series, DT);
}

print(Sys.time() - start_time)

# parallel::stopCluster(cl);
