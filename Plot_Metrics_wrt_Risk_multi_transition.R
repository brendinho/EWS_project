# This is the file for the first diagram tableau in the figure - the oen for multiple transitions is in another file
rm(list=ls())

options(scipen=10000)

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Change_wrt_risk_multi", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# source("Processing_Parameter_Table.R");
Params <- fread(parameter_file, colClasses="character")
# Parameters <- Params[num(beta)==1 & num(init_prop)==0.05 & num(size)==40000 & num(infec)==0.2 & num(norm)%in%c(1.125, 1.5, 1.625, 1.75, 1.875, 2.125)];
Parameters <- Params[num(init_prop)==0.05 & num(size)==10000 & num(norm)%in%c(1.96875, 1.59375, 1.8125, 2.375, 2.03125, 2.125) & num(infec)==0.8];
Parameters[, c("V1", "risk", "time"):=NULL];

Parameters <- unique(Parameters);

metrics_to_plot <- c("plain", "NV_join_count", "mutual_info"); # "moran_i", "join_counts" , "geary_c", "getis_ord", "geary_c", "getis_ord"); # "echo_chamber_number"); #, "echo_chamber_number", "conn_comp_number", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "vacc_conn_comp_size", "nonvacc_conn_comp_size"); # , "forced_new_ratio", "join_count",

time_start <- Sys.time();

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

    Prop_vs_Risk = data.table();

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

	for(file in filenames)
	{
		temp <- data.table(read.csv(file));
		Prop_vs_Risk <- rbind(Prop_vs_Risk, temp[is.na(instance), ], fill=TRUE);
	}
	Prop_vs_Risk[, X:=NULL]

	Trans_Table <- get_transitions(Prop_vs_Risk, initial_vacc_prop);

	prediction <- change_point_test(Prop_vs_Risk, "moran_i_mean", test_to_use="buishand_r");

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
		"_multitransition.png",
		sep=""
	);

	the_size <- num(size);

	colours <- list(c("red", "lightpink"), c("blue", "lightskyblue1"), c("green", "palegreen"));

	Proportion 	<- Prop_vs_Risk[initial_vacc_proportion==initial_vacc_prop];
	Proportion	<- Proportion[perceived_vaccine_risk<=-0.05 & perceived_vaccine_risk>=-0.65];

	mains <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_V", "soc_N");
	mains_mean <- sapply(mains, function(x) sprintf("%s_mean",x));
	mains_sd <- sapply(mains, function(x) sprintf("%s_sd",x));

	Proportion[, (mains_sd)  :=lapply(.SD, "/", the_size), .SDcols=mains_sd  ];
	Proportion[, (mains_mean):=lapply(.SD, "/", the_size), .SDcols=mains_mean];

	joins <- c("NN_join_count", "NV_join_count", "VV_join_count");
	joins_mean <- sapply(joins, function(x) sprintf("%s_mean",x));
	joins_sd <- sapply(joins, function(x) sprintf("%s_sd",x));

	# do the sd first before normalising the means, since normalising the means changes to total 1
	# nirmalising the sd's by one doesn't do fuck all, then
	Proportion[, (joins_sd)  :=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean), .SDcols=joins_sd  ];
	Proportion[, (joins_mean):=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean), .SDcols=joins_mean];

	Proportion <- Proportion[order(perceived_vaccine_risk)];

	# if there are less than 5 points to be plotted, then cancel this plot and skip to the next parameter set
	if(nrow(Proportion)<5){ next; }

	for(metric in metrics_to_plot)
	{
		pl <- ggplot(Proportion, aes(x=perceived_vaccine_risk), aes_string(x="perceived_vaccine_risk")) + 
		  labs(x=expression(kappa), y= (if(metric!="NV_join_count") axis_label(metric) else expression("<"*scriptstyle(N)*","*scriptstyle(V[s])*">"))) + 
		  theme(
		    panel.background=element_blank(),
		    axis.title=element_text(size=80), 
		    axis.text=element_text(size=60)
		  );

		pl <- pl + theme();

		if(metric == "plain")
		{
			pl <- pl +
				geom_ribbon(aes(ymin = (soc_V_mean-soc_V_sd), ymax = (soc_V_mean+soc_V_sd)),fill = "grey70") +
				geom_ribbon(aes(ymin = (phys_R_mean-phys_R_sd), ymax = (phys_R_mean+phys_R_sd)),fill = "grey70") +
				geom_ribbon(aes(ymin = (phys_V_mean-phys_V_sd), ymax = (phys_V_mean+phys_V_sd)),fill = "grey70") +
				geom_ribbon(aes(ymin = (soc_N_mean-soc_N_sd), ymax = (soc_N_mean+soc_N_sd)), fill = "grey70") +
				geom_line(aes(y=soc_V_mean), colour="green", size=1.5) +
				geom_line(aes(y=phys_R_mean), colour="black", linetype="dotdash", size=1.5) +
				geom_line(aes(y=phys_V_mean), colour="darkblue", linetype="dotdash", size=1.5) +
				geom_line(aes(y=soc_N_mean), colour="red", size=1.5);
		}
		else
		{
			for(index in seq_along(plot_these(metric)))
			{
				prefix <- plot_these(metric)[index];
				main_col <- colours[[index]][1];
				ribb_col <- colours[[index]][2];

				tmean <- sprintf("%s_mean", prefix);
				tsd <- sprintf("%s_sd", prefix);

				if(!( tmean %in% names(Proportion) )){ next; }

				pl <- pl +
					geom_ribbon(aes_string(ymin=Proportion[, get(tmean)-get(tsd)], ymax=(Proportion[, get(tmean)+get(tsd)])), fill=ribb_col) +
					geom_line(aes_string(x="perceived_vaccine_risk", y=tmean), colour=main_col, size=1.5);
			}
		}

		if(isTRUE( Trans_Table$is_soc_trans ))
		{
			pl <-  pl + geom_vline(xintercept=Trans_Table[, soc_intersec_risk], size=1.5);
		}
		if(isTRUE( Trans_Table$is_phys_trans ))
		{
			pl <-  pl + geom_vline(xintercept=Trans_Table[, phys_intersec_risk], size=1.5);
		}
		# pl <- pl + geom_vline(xintercept=prediction$kappa, size=3, linetype="dotted", colour="black")

		ggsave(
			paste(graph_file_path, folder_name, "/indi_", metric, file_name_ending, sep=""),
			plot=pl, width=20, height=7, units="in", limitsize=FALSE #  height=5,
		)
	}
}

print(Sys.time()-time_start);
