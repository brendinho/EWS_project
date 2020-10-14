# This is the file for the first diagram tableau in the figure - the oen for multiple transitions is in another file
rm(list=ls())

options(scipen=10000)

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Kappa_Series", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# source("Processing_Parameter_Table.R");
Params <- fread(parameter_file, colClasses="character");

image_height <- 5.5;
image_width <- 22;

metrics_to_plot <- c("conn_comp_size", "echo_chamber_size", "plain", "number", "modularity", "watts_strogatz", "opinion_change", "prob_sick", "echo_chamber_number", "conn_comp_number", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "vacc_conn_comp_size", "nonvacc_conn_comp_size");
# , "geary_c", "mutual_info", "NV_join_count", "moran_i", "join_counts", "join_count");

time_start <- Sys.time();

for(the_size in c("10000", "40000"))
{
	Parameters <- Params[num(size)==the_size & num(beta)==1 & soc_top=="random" & phys_top=="random"];

	# if(num(the_size) == 10000) Parameters <- Parameters[num(init_prop)==0.05 & num(norm)%in%c(0, 0.5)];
	# if(num(the_size) == 10000)  Parameters <- Parameters[num(init_prop)%in%c(0.05, 0.95) & num(norm)%in%c(0, 0.5)];

	if(num(the_size) == 10000)  Parameters <- Parameters[num(init_prop)==0.05];
	# if(num(the_size) == 40000)  Parameters <- Parameters[num(init_prop)==0.05 & num(norm)%in%c(0,0.25)];
	# if(num(the_size) == 562500) Parameters <- Parameters[num(init_prop)==0.05 & num(norm)%in%c(0,0.25) & num(random_switch)==0.1];

	Parameters[, c("V1", "risk", "time", "norm"):=NULL];
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
			Proportion <- rbind(Proportion, temp[is.na(instance), ], fill=TRUE);
		}
		Proportion[, X:=NULL];

		if(num(the_size) == 10000)  Proportion <- Proportion[perceived_vaccine_risk<=0.25 & perceived_vaccine_risk>=-0.25];
		if(num(the_size) == 40000)  Proportion <- Proportion[perceived_vaccine_risk<=0.25 & perceived_vaccine_risk>=-0.25];
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

		colours <- list(c("red", "lightpink"), c("blue", "lightskyblue1"), c("green", "palegreen"));

		mains <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_V", "soc_N");
		mains_mean <- sapply(mains, function(x) sprintf("%s_mean",x));
		mains_sd <- sapply(mains, function(x) sprintf("%s_sd",x));

		Proportion[, (mains_sd)  :=lapply(.SD, "/", the_size), .SDcols=mains_sd  ];
		Proportion[, (mains_mean):=lapply(.SD, "/", the_size), .SDcols=mains_mean];

		joins <- c("NN_join_count", "NV_join_count", "VV_join_count");
		joins_mean <- sapply(joins, function(x) sprintf("%s_mean",x));
		joins_sd <- sapply(joins, function(x) sprintf("%s_sd",x));

		# do the sd first before normalising the means, since normalising the means changes to total 1
		# normalising the sd's by one doesn't do fuck all, then
		Proportion[, (joins_sd)  :=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean), .SDcols=joins_sd  ];
		Proportion[, (joins_mean):=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean), .SDcols=joins_mean];

		Proportion <- Proportion[order(perceived_vaccine_risk)];

		# if there are less than 5 points to be plotted, then cancel this plot and skip to the next parameter set
		if(nrow(Proportion)<5){ next; }

		for(metric in metrics_to_plot)
		{
			if(metric != "plain"){
			for(the_thing in plot_these(metric))
			{
				if(!( the_thing %in% names(metric) )) next;
			}}

			print(metric)

			DT <- data.table(grp=character(), dep=numeric(), val=numeric(), lwr=numeric(), upr=numeric()); # dep means dependent

			for(index in seq_along(plot_these(metric)))
			{
				prefix <- plot_these(metric)[index];
				tmean <- sprintf("%s_mean", prefix);
				tsd <- sprintf("%s_sd", prefix);

				  temp_table <- data.table(
				    grp = rep(prefix, nrow(Proportion)),
				    dep = Proportion[, perceived_vaccine_risk],
				    val = Proportion[, get(tmean)],
				    lwr = Proportion[, get(tmean)-get(tsd)],
				    upr = Proportion[, get(tmean)+get(tsd)]
				  );

				DT <- do.call(smartbind, list(DT, temp_table));
			}

			pl_wo_ribbon <- ggplot(data=DT, aes(x=dep, y=val, group=grp, colour=grp, fill=grp)) +
				geom_point() +
				geom_line(size=2) +
				# theme_minimal() +
				theme(
				  plot.margin = unit(c(0.5,0,0,0), "cm"),
			    axis.title.y = element_text(margin = margin(t=0, r=20, b=0, l=0)),
				  axis.title.x = element_text(margin = margin(t=10, r=0, b=0, l=0))
				);
				# scale_y_continuous(breaks = scales::pretty_breaks(n = 4));
				# scale_y_continuous(breaks=pretty(DT$val, 4))
				# ylim(c(-Inf, Inf+0.2))

			if(isTRUE( Trans_Table$is_soc_trans ))
			{
				pl_wo_ribbon <-  pl_wo_ribbon + geom_vline(xintercept=Trans_Table[, soc_intersec_risk], size=1);
			}
			if(isTRUE( Trans_Table$is_phys_trans ))
			{
				pl_wo_ribbon <-  pl_wo_ribbon + geom_vline(xintercept=Trans_Table[, phys_intersec_risk], size=1);
			}
			pl_wo_ribbon <- pl_wo_ribbon +
				theme(
				  axis.text=element_text(size=70),
				  axis.title=element_text(size=80),
				  legend.position="none"
				) + #
				xlab(expression(kappa)) + ylab(axis_label(metric));

			pl <- pl_wo_ribbon + geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.3, linetype=0);

			ggsave(
				paste(graph_file_path, folder_name, "/indi_", metric, file_name_ending, sep=""),
				plot=pl, width=image_width, height=image_height, units="in", limitsize=FALSE #  height=5,
			)
			dev.off();

			if(metric == "watts_strogatz")
			{
				ggsave(
					paste(graph_file_path, folder_name, "/indi_", metric, "_wo_ribbon", file_name_ending, sep=""),
					plot=pl_wo_ribbon, width=image_width, height=image_height, units="in", limitsize=FALSE
				);
			  dev.off();
			}

			# stop()
		}
	}
}

print(Sys.time()-time_start);
