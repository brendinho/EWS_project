# Brendon Phillips
# PhD candidate
# Department of Applied Mathematics
# University of Waterloo

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

time_start <- Sys.time();

# source(sprintf("%s/Processing/Processing_Parameter_Table.R", home_folder));
Parameters <- fread(parameter_file, colClasses="character"); Parameters[, V1:=NULL];
if("time" %in% names(Parameters)) Parameters[, time:=NULL];
Parameters <- unique(Parameters);

if( exists("cl") ) parallel::stopCluster(cl);
# cl <- parallel::makeCluster( detectCores()[1]-1 ); #not to overload your computer
cl <- parallel::makeCluster( 9 ); #not to overload your computer
doParallel::registerDoParallel(cl);

data_file_path <- sprintf("%s/", stored_data_path);

interspersal_order <- c(1, 11, 2, 12, 3, 13, 4, 14, 5, 15, 6, 16, 7, 17, 8, 18, 9, 19, 10, 20);

# pb <- txtProgressBar(min=0, max=nrow(Parameters), style=3);
# for(line_index in 1:nrow(Parameters))
foreach(line_index=1:nrow(Parameters), .packages=c("data.table", "comprehenr")) %dopar%
{
	# setTxtProgressBar(pb, line_index);
	this_tuple <- Parameters[line_index];

	structure <- this_tuple[, network_structure];
	initial_vacc_prop <- this_tuple[, initial_vacc_proportion];
	size <- this_tuple[, N];
	risk <- this_tuple[, perceived_vaccine_risk];
	infection_prob <- this_tuple[, infec_prob];
	importation_rate <- this_tuple[, importation];
	degree <- this_tuple[, mean_degree];
	random_opinion_switch <- this_tuple[, random_opinion_switch];
	social_norm <- this_tuple[, social_norm];
	proportion <- this_tuple[, proportion_of_nodes];

	the_values <- to_list(for(col in names(Parameters)) convert_if_number(this_tuple[[col]]) );
	new_stem <- sprintf("%sHes_%s_inst_*.csv", data_file_path, paste(unlist(c(fragments, the_values))[interspersal_order], collapse='_'));
	data_file_list <- Sys.glob(new_stem);

	# if(length(data_file_list) == 0) next;
	if(length(data_file_list) == 0) stop();

	output_file_name <- sprintf("%sHes_%s_summary.csv", CSV_file_path, paste(unlist(c(fragments, the_values))[interspersal_order], collapse='_'));

	# read first file in glob vector, get the column titles, and process it like a regular file
	columns <-names(fread(data_file_list[1]));
	averaged_column_names =  setdiff(columns, c(hes_fixed_values, "V1", "instance", "beta_parameter", "time"));

	mean_names	<- as.vector(sapply( averaged_column_names, function(x) paste(x, "_mean", sep="") ));
	sd_names 	<- as.vector(sapply( averaged_column_names, function(x) paste(x, "_sd",   sep="") ));

	Mean_DT <- data.table(matrix( 0.0, nrow=0, ncol=(length(c(mean_names, sd_names, "instance"))) ));
	names(Mean_DT) <- c("instance", mean_names, sd_names);

	# I'll use this for calculating the mean of all the last 500's of each trial, instead of finding a mean of means
	All_of_Them <- data.table()

	for( file_index in seq_along(data_file_list) )
	{
		if(file.size(data_file_list[file_index]) == 0){ stop(); }
		# if(file.size(data_file_list[file_index]) == 0){ next; }

		Data <- data.table(matrix( nrow=1, ncol=length(c(mean_names, sd_names)) ));
		names(Data) <- c(mean_names, sd_names);

		file_name <- data_file_list[file_index];

		temp <- tryCatch(
			{ fread(data_file_list[file_index]); },
			error = function(err){ return(sprintf("warning on reading file %s"), file_name); } ,
			warning = function(war){ return(sprintf("error in reading the file %s", file_name)); }
		);
		if(! exists("temp")) next;
		if(! is.data.table(temp)) next;
		if("V1" %in% names(temp)) temp[, V1:=NULL];

		temp <- temp[time >=  max(time)-number_of_time_steps_averaged];
		temp[,  network_structure:=tolower(network_structure)];

		All_of_Them <- rbind(All_of_Them, temp, fill=TRUE)

		for(col in averaged_column_names){ if(col %in% names(temp))
		{
			Data[, paste(col, "_mean", sep="") := temp[, mean(get(col), na.rm=TRUE)] ];
			Data[, paste(col, "_sd",   sep="") := temp[, sd(get(col), na.rm=TRUE)] ];
		}}

		# Data[, (hes_fixed_values) :=  unique(temp[, .SD, .SDcols=hes_fixed_values])];
		Mean_DT <- rbind(Mean_DT, Data, fill=TRUE);

		rm(temp, Data)
	}

	Mean_DT[, instance:=1:nrow(Mean_DT)];
	Mean_DT <- rbind(Mean_DT, as.list(rep(NA, ncol(Mean_DT))));

	for(col in averaged_column_names){ if(col %in% names(All_of_Them))
	{
		Mean_DT[nrow(Mean_DT), paste(col, "_mean", sep="") := All_of_Them[, mean(get(col), na.rm=TRUE)] ];
		Mean_DT[nrow(Mean_DT), paste(col, "_sd",   sep="") := All_of_Them[, sd(get(col), na.rm=TRUE)] ];
	}}

	Mean_DT[, (hes_fixed_values) :=  this_tuple];

	write.csv( Mean_DT, file = output_file_name);

	rm(Mean_DT);
}

writeLines("\n")

if( exists("cl") ) parallel::stopCluster(cl);
print("done processing the files")
print(Sys.time()-time_start)
