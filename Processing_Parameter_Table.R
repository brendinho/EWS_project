# Brendon Phillips
# Ph.D. candidate
# Department of Applied Mathematics
# University of Waterloo
# extract the parameters from the dataset, with limits

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

time_start <- Sys.time();

the_sizes <- c("10000", "40000", "562500")

Params <- data.table(matrix(nrow=0, ncol=15));
names(Params) <- c("size", "duration", "beta", "init_prop", "risk", "infec", "import", "birth_death", "phys_top", "phys_deg", "soc_top", "soc_deg", "norm", "random_switch", "time");

for(the_size in the_sizes)
{
	data_file_path <- sprintf("%s/N_%s_cleaned", stored_data_path, the_size);
	all_raw_files <- Sys.glob(sprintf("%s/stats*bin", data_file_path));

	for(file_index in seq_along(all_raw_files))
	{
		file_name_pieces <- strsplit(all_raw_files[file_index], "_")[[1]];

		size <- toString(file_name_pieces[which(file_name_pieces == "N")+1]);
		duration <- toString(file_name_pieces[which(file_name_pieces == "dur")+1]);
		beta <- toString(file_name_pieces[which(file_name_pieces == "beta")+1]);
		init_prop <- toString(file_name_pieces[which(file_name_pieces == "vaccprop")+1]);
		risk <- toString(file_name_pieces[which(file_name_pieces == "pay")+1]);
		infection <- toString(file_name_pieces[which(file_name_pieces == "inf")+1]);
		importation <- toString(file_name_pieces[which(file_name_pieces == "imp")+1]);
		birth_death <- toString(file_name_pieces[which(file_name_pieces == "rep")+1]);
		physical_top <- file_name_pieces[which(file_name_pieces == "ptop")+1];
		physical_deg <- toString(file_name_pieces[which(file_name_pieces == "pdeg")+1]);
		social_top <- file_name_pieces[which(file_name_pieces == "stop")+1];
		social_deg <- toString(file_name_pieces[which(file_name_pieces == "sdeg")+1]);
		norm <- toString(file_name_pieces[which(file_name_pieces == "norm")+1]);
		random_switch <- toString(file_name_pieces[which(file_name_pieces == "switch")+1]);
		time <- if("equilib.bin" %in% file_name_pieces ) "equilib" else "unconverged";

		Params <- rbind(Params, as.list(c(size, duration, beta, init_prop, risk, infection, importation, birth_death, physical_top, physical_deg, social_top, social_deg, norm, random_switch, time)));
	}
}

write.csv(unique(Params), file=parameter_file);
print(Sys.time() - time_start);
