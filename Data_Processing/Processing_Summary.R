# Brendon Phillips

# TIme series generator

rm(list=ls());

source("/media/b2philli/Simulations/Processing/Parameter_Values.R");

time_start <- Sys.time();

# source(sprintf("%s/Processing/Processing_Parameter_Table.R", home_folder));
Params <- fread(parameter_file, colClasses="character");
Params[, c("V1", "time"):=NULL];

DO_SD <- TRUE;

if( exists("cl") ) parallel::stopCluster(cl);
# cl <- parallel::makeCluster( detectCores()[1]-1 ); #not to overload your computer
cl <- parallel::makeCluster( 5 ); #not to overload your computer
doParallel::registerDoParallel(cl);

for(the_size in c("10000")) # "40000"))
{
	data_file_path <- sprintf("%s/N_%s_cleaned/", home_folder, the_size);

	Parameters <- Params[num(beta)==1 & num(init_prop)==0.05 & num(size)==num(the_size) & num(norm)>=0 & num(norm)<=3 & num(risk)>=-1 & num(risk)<=1];
	Parameters <- unique(Parameters);

	# monitoring_file <- sprintf("processing_log_size_%s.txt", the_size);
	# writeLines(c(""), monitoring_file);
	# skipped_parameter_tuples <- c();

	# foreach(line_index=1400:nrow(Parameters), .packages=c("data.table")) %dopar%
	for(line_index in 6025:nrow(Parameters))
	{
		# sink(monitoring_file, append=TRUE);
		# cat(sprintf("iteration number: %i\n", line_index));

		writeLines(sprintf("iteration number: %i", line_index));

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

		file_end_line <- "_inst*";
		file_end_line <- sprintf("%s%s", file_end_line, ".bin");

		output_file_name <- paste(
			CSV_file_path,
			"Summary",
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

		data_file_list <- Sys.glob(
			paste(
				data_file_path,
				# "stats_cleaned",
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
				file_end_line,
				sep=""
			)
		);

		# # use this one to stop the slave if running in parallel
		# if(length(data_file_list) == 0) stop();
		# use this one of running the loop in serial
		if(length(data_file_list) == 0) next;

		Data <- data.table();

		# read first file in glob vector, get the column titles, and process it like a regular file
		columns <-names(fread(data_file_list[1]));

		averaged_column_names =  c(setdiff(columns, c(fixed_values, "V1", "time", "imported_cases", "secondary_infections", "sum_squares_degrees")), "ratio_sec_to_all_infs_full_sim", "ratio_sec_to_all_infs_at_the_end");

		Mean_DT <- data.table()

		mean_names <- as.vector(sapply(
			averaged_column_names,
			function(x) paste(x, "_mean", sep="")
		));

		if(DO_SD)
		{
			sd_names <- as.vector(sapply(
				averaged_column_names,
				function(x) paste(x, "_sd", sep="")
			));

			Mean_DT <- data.table(matrix(
				0.0,
				nrow=0,
				ncol=(length(c(mean_names, sd_names, "instance")))
			));
			names(Mean_DT) <- c(mean_names, sd_names, "instance");
		}
		else
		{
			Mean_DT <- data.table(matrix(
				0.0,
				nrow=0,
				ncol=(length(c(mean_names, "instance")))
			));
			names(Mean_DT) <- c(mean_names, "instance");
		}


		for( file_index in seq_along(data_file_list) )
		{
			if(file.size(data_file_list[file_index]) == 0){ next; }

			# Data is for the mean and sd calculations based on the values in the temp table
			Data <- data.table()

			if(DO_SD){ Data <- data.table(matrix( nrow=1, ncol=length(c(mean_names, sd_names)) )); }
			else{ Data <- data.table(matrix( nrow=1, ncol=length(mean_names) )); }

			file_name <- data_file_list[file_index];

			temp <- tryCatch(
				{ fread(data_file_list[file_index]); },
				error = function(err){ return(sprintf("warning on reading file %s"), file_name); } ,
				warning = function(war){ return(sprintf("error in reading the file %s", file_name)); }
			);

			if(! exists("temp")) next;
			if(! is.data.table(temp)) next;

			if("V1" %in% names(temp)) { temp <- temp[, V1:=NULL]; }

			# the most important line
			temp <- temp[time >=  max(time)-number_of_time_steps_averaged];

			# some of the topology names have errors in them
			temp[,  social_topology:=tolower(social_topology)];
			temp[, physical_topology:=tolower(physical_topology)];

			if( ("secondary_infections" %in% names(temp)) & ("imported_cases" %in% names(temp)) )
			{
				temp[, "ratio_sec_to_all_infs_at_the_end" := temp[, cumsum(secondary_infections)/cumsum(secondary_infections + imported_cases)] ];
			}
			# else
			# {
			# 	averaged_column_names <- averaged_column_names[ ! averaged_column_names %in% c("ratio_sec_to_all_infs_full_sim", "ratio_sec_to_all_infs_at_the_end")]
			# }

			if(DO_SD){ names(Data) <- c(mean_names, sd_names); }
			else{ names(Data) <- mean_names; }

			instance_number <- temp[, unique(instance)];

			for(col in averaged_column_names){ if(col %in% names(temp))
			{
				Data[, paste(col, "_mean", sep="") := temp[, mean(get(col), na.rm=TRUE)] ];
				if(DO_SD){ Data[, paste(col, "_sd",   sep="") := temp[, sd(get(col), na.rm=TRUE)] ]; }
			}}

			Data[, (reinserted_values) :=  unique(temp[, .SD, .SDcols=reinserted_values])];

			Mean_DT <- rbind(Mean_DT, Data, fill=TRUE);

			rm(temp, Data, instance_number);
		}

		rm(data_file_list);

		# since the "change_files" file may change the instance number of the file name, it may not match the actual
		# instance number of the data.table inside, and it may also be multiply defined.
		# choosing and assigning a unique one here avoids the hassle
		Mean_DT[, instance:=1:length(instance)];

		if(DO_SD)
		{
			mean_of_all <- Mean_DT[, lapply(.SD, function(inp){ mean(x=inp, na.rm=TRUE) }), .SDcols=c(mean_names, sd_names)];
		}
		else
		{
			mean_of_all <- Mean_DT[, lapply(.SD, function(inp){ mean(x=inp, na.rm=TRUE) }), .SDcols=mean_names];
		}

		mean_of_all[, "instance":=NaN];
		mean_of_all[, (reinserted_values) :=  Mean_DT[, unique(.SD), .SDcols=reinserted_values] ];
		Mean_DT <- rbind(Mean_DT, mean_of_all);

		write.csv( Mean_DT, file = output_file_name);

		rm(Mean_DT);
	}
}

if( exists("cl") ) parallel::stopCluster(cl);
print("done processing the files")
print(Sys.time()-time_start)
