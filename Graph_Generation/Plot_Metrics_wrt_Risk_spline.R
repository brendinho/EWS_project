# This is the file for the first diagram tableau in the figure - the oen for multiple transitions is in another file
rm(list=ls())

options(scipen=10000)

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Change_wrt_perceived_vaccine_risk_individual/Splines", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# source("Processing_Parameter_Table.R");
Params <- fread(parameter_file, colClasses="character")

metrics_to_plot <- c("plain", "watts_strogatz"); # "modularity", "watts_strogatz", "opinion_change", "echo_chamber_size", "conn_comp_size", "conn_comp_number", "echo_chamber_number", "watts_strogatz", "prob_sick"

# c("plain", "geary_c", "mutual_info", "NV_join_count", "moran_i", "join_counts"); # , "getis_ord""echo_chamber_number"); #, "echo_chamber_number", "conn_comp_number", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "vacc_conn_comp_size", "nonvacc_conn_comp_size"); # , "forced_new_ratio", "join_count",

time_start <- Sys.time();

for(the_size in c("10000")) # "40000",
{
	Parameters <- Params[num(size)==the_size & num(beta)==1 & soc_top=="random" & phys_top=="random"];

	if(num(the_size) == 10000) Parameters <- Parameters[num(init_prop)==0.05 & num(norm)%in%c(0, 0.5)];
	# if(num(the_size) == 10000)  Parameters <- Parameters[num(init_prop)%in%c(0.05, 0.95) & num(norm)%in%c(0, 0.5)];


	# if(num(the_size) == 10000)  Parameters <- Parameters[num(init_prop)==0.05 & num(norm)%in%c(0,0.25)];
	# if(num(the_size) == 40000)  Parameters <- Parameters[num(init_prop)==0.05 & num(norm)%in%c(0,0.25)];
	# if(num(the_size) == 562500) Parameters <- Parameters[num(init_prop)==0.05 & num(norm)%in%c(0,0.25) & num(random_switch)==0.1];

	Parameters[, c("V1", "risk", "time"):=NULL];
	Parameters <- unique(Parameters);

	# foreach( line_index=1:nrow(Parameters), .packages=c("data.table")) %dopar%
	for(line_index in 1:nrow(Parameters))
	{
		this_tuple <- Parameters[line_index];

		beta <- this_tuple[, beta];
		initial_vacc_prop <- this_tuple[, init_prop];
		size <- this_tuple[, size];
		duration <- this_tuple[, duration];
		infection_prob <- this_tuple[, infec];
		importation_rate <- this_tuple[, import];
		replenishment_rate <-  this_tuple[, birth_death];
		physical_topology <- this_tuple[, phys_top];
		physical_degree <- this_tuple[, phys_deg];
		social_topology <- this_tuple[, soc_top];
		social_degree <- this_tuple[, soc_deg];
		random_opinion_switch <- this_tuple[, random_switch];
		social_norm <- this_tuple[, norm];
		num_size <- num(size);

		filenames <- Sys.glob(paste(
			CSV_file_path,
			"Summary",
			"_N_", size,
			"_dur_", duration,
			"_beta_", beta,
			"_vaccprop_", initial_vacc_prop,
			"_risk_", "*",
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
		));

		if(length(filenames) == 0){ next; }

		Proportion = data.table();
		for(file in filenames)
		{
			temp <- data.table(read.csv(file));
			Proportion <- rbind(Proportion, temp[!is.na(instance), ], fill=TRUE);
		}
		Proportion[, X:=NULL];

		if(num(the_size) == 10000)  Proportion <- Proportion[perceived_vaccine_risk<=0.25 & perceived_vaccine_risk>=-0.25];
		# if(num(the_size) == 40000)  Proportion <- Proportion[perceived_vaccine_risk<=0.25 & perceived_vaccine_risk>=-0.25];
		# if(num(the_size) == 562500) Proportion <- Proportion[perceived_vaccine_risk<=0.25 & perceived_vaccine_risk>=-0.25];

		Trans_Table <- get_transitions(Proportion, initial_vacc_prop);

		file_name_ending <- paste(
			"_N_", size,
			"_dur_", duration,
			"_beta_", beta,
			"_vaccprop_", initial_vacc_prop,
			"_inf_", infection_prob,
			"_imp_", importation_rate,
			"_rep_", replenishment_rate,
			"_ptop_", physical_topology,
			"_pdeg_", physical_degree,
			"_stop_", social_topology,
			"_sdeg_", social_degree,
			"_norm_", social_norm,
			"_switch_", random_opinion_switch,
			".png",
			sep=""
		);

		the_size <- num(size);

		# if there are less than 5 points to be plotted, then cancel this plot and skip to the next parameter set
		if(nrow(Proportion)<5){ next; }

		for(metric in metrics_to_plot)
		{
			# # plot_these_here <- c("phys_R_mean")# , "phys_V_mean");
			plot_these_here <- sprintf("%s_mean", plot_these(metric));

			Proportion_Here <- melt(Proportion[, .SD, .SDcols=c("perceived_vaccine_risk", plot_these_here)], id="perceived_vaccine_risk");

			the_plot <- ggplot(Proportion_Here) +
				geom_spline(data=Proportion_Here, aes(x=perceived_vaccine_risk, y=value, colour=variable), na.rm=TRUE) # +
				# geom_ribbon(data=DT, aes(x=dep, ymin=lwr, ymax=upr, colour=grp), alpha=.3, linetype=0);

			if(isTRUE( Trans_Table$is_soc_trans ))
			{
				the_plot <-  the_plot + geom_vline(xintercept=Trans_Table[, soc_intersec_risk], size=1);
			}
			if(isTRUE( Trans_Table$is_phys_trans ))
			{
				the_plot <-  the_plot + geom_vline(xintercept=Trans_Table[, phys_intersec_risk], size=1);
			}
		 	the_plot<- the_plot +
				scale_y_continuous(label=scientific_format(digits=1)) +
				theme(axis.text=element_text(size=50), axis.title=element_text(size=50), legend.position="none") + #
				xlab(expression(kappa)) + ylab(axis_label(metric));

			pl <- the_plot # + geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.3, linetype=0);

			ggsave(
				paste(graph_file_path, folder_name, "/indi_spline_", metric, file_name_ending, sep=""),
				plot=pl, width=20, height=5, units="in", limitsize=FALSE #  height=5,
			)
		}
	}
}

print(Sys.time()-time_start);
