# THIS FILE MUST BE RUN IN SERIAL _ MUST MUST MUST - BECAUSE OF THE FILE NAME CHECKING AND RACE CONDITIONS

start <- Sys.time();

get_from <- "/home/b2philli/Desktop/Holding";
send_to <- "/home/b2philli/Desktop/Files";

for(filename in Sys.glob(sprintf("%s/*bin", get_from)))
{
	# print(file.exists(filename));
	# print(file.exists(sprintf("%s/%s", send_to, tail(str_split(filename, "/")[[1]], n=1) )))

	if( file.exists(sprintf("%s/%s", send_to, tail(str_split(filename, "/")[[1]], n=1) )) )
	{
		first_part_of_the_name <- tail(str_split(spl[1], "/")[[1]], 1);

		spl <- str_split(filename, "inst_")[[1]];
		second_part_of_the_name <- paste(str_split(spl[2], "_")[[1]][-1], collapse="_");

		trial_number <- 0;
		repeat
		{
			if( file.exists(new_name) )
			{
				trial_number <- trial_number+1;
				new_name <- sprintf("%s/%sinst_%s_%s", send_to, first_part_of_the_name, trial_number, second_part_of_the_name);
			}
			else
			{
				temp <- fread(filename);
				temp[, instance:=trial_number];
				write.csv( fread(filename), new_name );
				break;
			}
		}
	}
	else
	{
		system(sprintf("cp %s %s", filename, send_to));
	}
	system(sprintf("rm %s", filename));
}

system("python3 remove_duplicates.py");

temp <- Sys.time();
print(end-start);
