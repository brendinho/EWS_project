# Brendon Phillips
# PhD candidate
# Department of Mathematics
# University of Waterloo

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Parameter_Values.R");

folder_name <- paste("Time_Series", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character")
Params[, "V1":=NULL]

# Parameters <- Params[num(size)==562500 & num(random_switch)==0.0001 & num(init_prop)==0.05 & num(norm)==0 & phys_top=="random" & num(risk)%in%c(0,0.01,-0.01)]
# Parameters <- Params[num(size)==562500 & num(random_switch)==0.001 & num(init_prop)==0.05 & num(norm)==0 & phys_top=="random" & num(risk)%in%c(0,0.005,-0.005)]
# Parameters <- Params[num(size)==562500 & num(random_switch)==0.01 & num(init_prop)==0.05 & num(norm)==0 & phys_top=="random" & num(risk)%in%c(0,0.02,-0.02)]
Parameters <- Params[num(size)==10000 & num(beta)==1 & num(init_prop)==0.05 & num(norm)==0 & num(risk)%in%c(0.0312, -0.0312, 0) & num(infec)%in%c(0.8)]; #
# Parameters <- Params[num(size)==40000 & num(beta)==1 & num(norm)%in%c(0) & num(infec)%in%c(0.2, 0.8) & num(init_prop)==0.05 & num(risk)%in%c(-0.0312, 0.0312, 0)]; #

# soc_V_EVO_N_40000_dur_2_beta_1.000000_vaccprop_0.05_risk_0.0312_in_0.8_im_2.5e-05_rep_0.00024_pt_random_pd_30_st_random_sd_30_sn_0_sw_1e-04

Parameters <- unique(Parameters);
quantities_to_plot <- c("soc_V", "phys_V"); #, "phys_R", "phys_S",

# "vacc_modularity", "nonvacc_modularity", "total_modularity", "prob_sick_nonvacc", "prob_sick_vacc", "vacc_min_conn_comp_size", "vacc_max_conn_comp_size", "vacc_avg_conn_comp_size", "nonvacc_min_conn_comp_size", "nonvacc_max_conn_comp_size", "nonvacc_avg_conn_comp_size");
# ,"mutual_info", "NV_join_count") #, "moran_i", "NN_join_counnt", "VV_join_count", "basics", "neighbours", "joins", "opinion", "prob_sick", "warning_signals", "opinion"); # "modularity", "nonvacc_conn_comps", "vacc_conn_comps",  "vacc_echo_chambers", "nonvacc_echo_chambers", "number", "watts_strogatz"

time_start <- Sys.time();

text_size <- 35;
title_size <-45;

for(line_index in 1:nrow(Parameters))
{
	this_tuple <- Parameters[line_index];

	beta <- this_tuple[, beta];
	initial_vacc_prop <- this_tuple[, init_prop];
	the_size <- this_tuple[, size];
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

	data_file_path <- sprintf("/media/%s/Simulations/N_%s_cleaned/", user_name, the_size);

	file_end_line <- "_inst*";
	file_end_line <- sprintf("%s%s", file_end_line, ".bin");

	data_file_list <- Sys.glob(
		paste(
			data_file_path,
			# "stats_cleaned",
			"stats*",
			"_N_", the_size,
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

	# painfully slow
	# 8 mins with processor running full tilt
	# preread all the files, figure out how much storage we need, make that size data.table, and then fill it

	# line_counter <- 0; column_counter <- 0;
	# for(file in data_file_list)
	# {
	# 	if(file.size(file) == 0) { next; }
	# 	temp <- fread(file);
	# 	line_counter <- line_counter + nrow(temp);
	# 	column_counter <- max(column_counter, ncol(temp));
	# }
	# just makesure that the columns are all there and lined up properly before blindly inserting int o the table
	# stop();

	if(length(data_file_list) == 0) next;

	Summary_Data <- data.table();

	for(file_index in 1:length(data_file_list))
	{
		if(file.size(data_file_list[file_index]) == 0) { next; }
		temp <- fread(data_file_list[file_index]);
		if(nrow(temp) == 0){ next; }
		if(nrow(Summary_Data) != 0)
		{
			# just in case the trial number is not unique, make it so
			if(temp[, unique(instance)] %in% unique(Summary_Data[, instance]))
			{
				temp_trial <- 0;
				while(temp_trial %in% unique(Summary_Data[, instance])) { temp_trial = temp_trial + 1; }
				temp[, ("instance"):=rep(temp_trial, nrow(temp))]
			}
		}
		Summary_Data <- rbind(temp, Summary_Data, fill=TRUE);
	}

	Summary_Data[, "V1":=0]
	props <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_N", "soc_V");
	Summary_Data[, (props):=lapply(.SD, function(x) x/num(the_size)), .SDcols=props];

	for(metric in quantities_to_plot)
	{
		plot_file_name <- paste(
			graph_file_path,
			folder_name,
			"/", metric,
			"_EVO",
			"_N_", the_size,
			"_dur_", duration,
			"_beta_", beta,
			"_vaccprop_", initial_vacc_prop,
			"_risk_", risk,
			"_in_", infection_prob,
			"_im_", importation_rate,
			"_rep_", replenishment_rate,
			"_pt_", physical_topology,
			"_pd_", physical_degree,
			"_st_", social_topology,
			"_sd_", social_degree,
			"_sn_", social_norm,
			"_sw_", random_opinion_switch,
			".png",
			sep=""
		);

		pl <- ggplot(Summary_Data, aes(time, get(metric), group=instance)) +
			geom_line() +
			labs(x=expression(tau), y=axis_label(metric)) +
			theme(
				panel.background=element_blank(),
				axis.line=element_line(colour = "black"),
				axis.text=element_text(size=text_size),
				axis.title=element_text(size=title_size),
				legend.position="none"
			);

		ggsave(
			plot_file_name,
			plot = pl, width=10, height=4, limitsize=FALSE, dpi=50
		);
		dev.off();

		system(sprintf("convert %s -trim %s", plot_file_name, plot_file_name));
	}
}

graphics.off()

print("done")
print(Sys.time()-time_start)
