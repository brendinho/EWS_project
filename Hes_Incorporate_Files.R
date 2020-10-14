# Brendon Phillips
# PhD candidate
# Department of Applied Mathematics
# University of Waterloo

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R")

start_time <- Sys.time();

# # takes 2 hours
# Processed_Files <- Sys.glob(sprintf("%s/*csv", stored_data_path));
# Processed_Hashes <- data.table(matrix("", ncol=2, nrow=length(Processed_Files)));
# names(Processed_Hashes) <- c("hash", "done"); Processed_Hashes[, done:=as.numeric(done)];
# for(file_index in 1:length(Processed_Files))
# {
# 	set(Processed_Hashes, file_index, names(Processed_Hashes), list(md5sum(Processed_Files[file_index]), 1));
# }
# Processed_Hashes <- unique(Processed_Hashes);
# fwrite(Processed_Hashes, "processed_hashes.csv");

Processed_Hashes <- data.table(fread( "processed_hashes.csv"));

Processed_Hashes <- unique(Processed_Hashes);
setkey(Processed_Hashes, hash);

print(Sys.time() - start_time);

data_folder_names <- c("/home/b2philli/MFCF"); # , "/home/b2philli/Sharcnet_project""/home/b2philli/Sharcnet_scratch",

number_of_files <- 0;
for(the_folder in data_folder_names){ number_of_files <- number_of_files + length(Sys.glob(sprintf("%s/*bin", the_folder))); }

interspersal_order <- c(1, 11, 2, 12, 3, 13, 4, 14, 5, 15, 6, 16, 7, 17, 8, 18, 9, 19, 10, 20);
for(the_folder in data_folder_names)
{
	writeLines(sprintf("\nTHE FOLDER: %s", the_folder));
	filelist <- Sys.glob(sprintf("%s/*bin", the_folder));

	# pb <- txtProgressBar(min=0, max=length(filelist), style=3);
	pb <- txtProgressBar(min=0, max=length(filelist), style=3);
	for(file_index in 1:length(filelist))
	{
		# print(sprintf("%s out of %s", file_index, number_of_files));
		setTxtProgressBar(pb, file_index);
		old_filename <- filelist[file_index];
		if(file.info(old_filename)$size == 0){ next; }
		old_file_hash <- md5sum(old_filename)[[1]];
		if(!is.na( Processed_Hashes[.(old_file_hash), unique(done)] )){ next; }
		DT <- data.table(fread(old_filename));
		name_parameters <- DT[, unique(.SD), .SDcols=hes_fixed_values];
		the_values <- to_list(for(col in hes_fixed_values) convert_if_number(name_parameters[[col]]) )
		new_stem <- sprintf("Hes_%s_inst", paste(unlist(c(fragments, the_values))[interspersal_order], collapse='_'));

		instance <- 0;
		while(instance <- instance + 1)
		{
			new_file_name <- sprintf("%s/%s_%i.csv", stored_data_path, new_stem, instance);
			if(!file.exists(new_file_name))
			{
				write.csv(DT, new_file_name);
				Processed_Hashes[old_file_hash] <- 1;
				break;
			}
		}
	}

	fwrite(Processed_Hashes, "processed_hashes.csv");

	# stop(Sys.time() - start_time)
}

fwrite(Processed_Hashes, "processed_hashes.csv");

print(Sys.time() - start_time);

# if(old_file_hash %in% Processed_Hashes){ next; }
# if(!is.null(Processed_Hashes[old_file_hash])){ next; }

# write.csv(unique(na.omit(parameters_of_the_new_files)), "/home/b2philli/Desktop/new_files.csv");
# write.csv(unique(na.omit(The_Parameters)), "/home/b2philli/Dropbox/CSV/Hes_Parameters.csv");

# Processed_Hashes <- scan("processed_hashes.txt", character());
# Processed_Hashes <- character(length(Processed_Files));
# for(index in seq_along(Processed_Files))
# {
# 	Processed_Hashes[index] <- md5sum(Processed_Files[index])[[1]];
# }

# gather all the bins from the three folders, process clean and remname them, dump tham into the hesitance folder

# print("Hashing files...");
# start_hash <- Sys.time();
# Processed_Hashes <- as.vector(md5sum(Processed_Files));
# Processed_Hashes <- character(length(Processed_Files));
# Processed_Hashes <- Processed_Hashes[Processed_Hashes!=""];
# Processed_Hashes <- data.table(read.csv("processed_hashes.txt"));

# Processed_Hashes <- list();
# for(file_index in 1:length(Processed_Files[1:50]))
# {
# 	print(file_index)
# 	Processed_Hashes[ md5sum(Processed_Files[file_index]) ] <- 1;
# }
# save(Processed_Hashes, file="processed_hashes.RData");
# Processed_Hashes <- load("processed_hashes.RData");

# The_Parameters <- data.table(matrix(nrow=number_of_files, ncol=length(hes_fixed_values)))
# names(The_Parameters) <- hes_fixed_values;
# number_cols <- setdiff(hes_fixed_values, "network_structure")
# The_Parameters[, (number_cols):=lapply(.SD, as.numeric), .SDcols=number_cols];
# The_Parameters[, ("network_structure"):=lapply(.SD, as.character), .SDcols="network_structure"];

# set(The_Parameters, file_counter, names(The_Parameters), the_values);
# The_Parameters[file_counter, (colnames(The_Parameters)):=the_values];

# if(!(new_file_name %in% Processed_Files))

# Processed_Files <- c(Processed_Files, new_file_name);
# Processed_Hashes <- c(Processed_Hashes, old_file_hash);

# parameters_of_the_new_files <- rbind(parameters_of_the_new_files,  DT[, unique(.SD), .SDcols=hes_fixed_values]);

# file_counter <- file_counter + 1;
