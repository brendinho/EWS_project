# Brendon Phillips
# Ph.D. candidate
# Department of Applied Mathematics
# University of Waterloo
# extract the parameters from the dataset, with limits

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

time_start <- Sys.time();

raw_files <- Sys.glob(sprintf("%s/*csv", stored_data_path));

Params <- data.table(matrix(nrow=length(raw_files), ncol=length(hes_fixed_values)+1, 0));
names(Params) <- c(hes_fixed_values, "instance");
Params[, network_structure:=as.character(network_structure)];

pb <- txtProgressBar(min=0, max=length(raw_files), style=3);
for(file_index in seq_along(raw_files))
{
	setTxtProgressBar(pb, file_index);
	the_split <- strsplit(strsplit(raw_files[file_index], ".csv")[[1]], "_")[[1]];
	getval <- function(x) return( convert_if_number(the_split[which(the_split == x)+1]) )

	# parallel short list to hes_fixed_values, since the file name has shorter keys than the table in the file
	hes_fixed_values_short <- c("N", "struct", "deg", "risk", "infec", "norm", "imp", "init", "switch", "prop", "inst");
	the_values <- to_list(for(val in hes_fixed_values_short) getval(val));

	set(Params, file_index, names(Params), the_values);
}
Params <- unique(Params[N!=0]);

write.csv(unique(Params), file=sprintf("%s_%s.csv", strsplit(parameter_file, '.csv')[[1]], "instances"));
write.csv(unique(Params[, -c("instance")]), file=parameter_file);
print(Sys.time() - time_start);
