# Summary Plot for the Topology
rm(list=ls())

source("~/Dropbox/Small_Data_Processing/Parameter_Values.R");

folder_name <- paste("Heatmaps", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# metrics_to_plot <- c(,, "geary_c_mean", "mutual_info_mean", "vacc_min_conn_comp_size_mean", "vacc_max_conn_comp_size_mean" "nonvacc_echo_chamber_size", "vacc_echo_chamber_size", "conn_comp_number", "nonvacc_conn_comp_size", "vacc_conn_comp_size", "join_count", "moran_i", "geary_c", "getis_ord", "mutual_info", "watts_strogatz", "modularity", "forced_new_ratio");

metrics_to_plot <- c("moran_i", "join_counts", "nonvacc_echo_chambers", "vacc_echo_chambers", "mutual_info", "conn_chamb_numbers", "nonvacc_conn_comps", "vacc_conn_comps"); #, "modularity", "watts_strogatz",  "geary_c", "getis_ord", "infections", "opinion");

start <- Sys.time();

wide <- 15000;
high <- 10000;
paragraph <- c(4.25,5.5,3,1);
# resolut <- 300;

for(line_index in 1:nrow(Parameters))
# foreach( line_index=1:nrow(Parameters), .packages=c("data.table")) %dopar%
{
	this_tuple <- Parameters[line_index];
#
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
    social_topology <- physical_topology;
    social_degree <- physical_degree;

    Proportion_vs_Risk_Table <- data.table();

	temp_start <- Sys.time();

	for(risk in vaccine_risk_vector){
	for(social_norm in social_norm_vector)
	{
		filename <- paste(
			CSV_file_path,
			"Summary",
			"_N_", size,
			"_dur_", duration,
			"_vaccprop_", init_vacc,
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
			if(vacc_allowed) "" else "_no_vaccination",
			".csv",
			sep=""
		);

		if( file.exists(filename) )
		{
			temp <- data.table(read.csv(filename));
			Proportion_vs_Risk_Table <- rbind(Proportion_vs_Risk_Table, temp[is.na(instance), ], fill=TRUE);
		}
	}}

	temp_end <- Sys.time();
	print(temp_end-temp_start);

	writeLines("\n\n");
	writeLines(paste("\treplenishment ratio: ", replenishment_rate, sep=""));
	writeLines(paste("\timportation rate: ", importation_rate, sep=""));
	writeLines(paste("\tnetwork size: ", size, sep=""));
	writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
	writeLines(paste("\tphysical topology: ", physical_topology, sep=""));
	writeLines(paste("\tphysical degree: ", physical_degree, sep=""));
	writeLines(paste("\tsocial topology: ", social_topology, sep=""));
	writeLines(paste("\tsocial degree: ", social_degree, sep=""));
	writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));

	# since we're plotting the sequence with lines, we need to make sure that the risks are in numerical order
	Proportion_vs_Risk_Table <- Proportion_vs_Risk_Table[order(perceived_vaccine_risk)];

	file_name_ending <- paste(
		"_N_", size,
		"_vacc_", init_vacc,
		"_dur_", duration,
		"_inf_", infection_prob,
		"_imp_", importation_rate,
		"_rep_", replenishment_rate,
		"_ptop_", physical_topology,
		"_pdeg_", physical_degree,
		"_stop_", social_topology,
		"_sdeg_", social_degree,
		"_norm_", sprintf("%1.2f", social_norm),
		"_switch_", random_opinion_switch,
		".png",
		sep=""
	);

	for(metric in metrics_to_plot)
	{
		plot_list <- list();

		if(metric == "moran_i")
		{
			DF <- Proportion_vs_Risk_Table[, .SD, .SDcol=c("social_norm", "perceived_vaccine_risk", "soc_V_mean", "moran_i_mean", "N")];

			local({

				plot_list[[1]] <<- ggplot(
						DF,
						aes(x=social_norm, y=perceived_vaccine_risk)
					) +
					geom_raster(aes(fill=soc_V_mean/N))  +
					scale_colour_gradientn(colours=terrain.colors(20)) +
					xlab(expression(kappa)) + ylab(expression(sigma)) +
					guides(fill=guide_legend(title="Vacc"));

				plot_list[[2]] <<- ggplot(
						DF,
						aes(x=social_norm, y=perceived_vaccine_risk)
					) +
					geom_raster(aes(fill=moran_i_mean))  +
					# scale_fill_gradient2(low="navy", high="red", mid="white") +
					scale_colour_gradientn(colours=terrain.colors(20)) +
					xlab(expression(kappa)) + ylab(expression(sigma)) +
					guides(fill=guide_legend(title="Moran's I"));

			});

			row_length <- 2;
		}
		else if(metric == "join_counts")
		{
			DF <- Proportion_vs_Risk_Table[, .SD, .SDcol=c("social_norm", "perceived_vaccine_risk", "soc_V_mean", "NN_join_count_mean", "NV_join_count_mean", "VV_join_count_mean")];

			local(
			{
				plot_list[[1]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=soc_V_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Vacc"));

				plot_list[[2]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=NN_join_count_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="[N,N]"));

				plot_list[[3]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=NV_join_count_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="[N,V]"));

				plot_list[[4]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=VV_join_count_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="[V,V]"));

			});

			row_length <- 4;
		}
		else if(metric == "nonvacc_echo_chambers")
		{
			DF <- Proportion_vs_Risk_Table[, .SD, .SDcol=c("social_norm", "perceived_vaccine_risk", "soc_V_mean", "nonvacc_min_chamber_size_mean", "nonvacc_max_chamber_size_mean", "nonvacc_avg_chamber_size_mean")];

			local(
			{
				plot_list[[1]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=soc_V_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Vacc"));

				plot_list[[2]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=nonvacc_min_chamber_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Min"));

				plot_list[[3]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=nonvacc_max_chamber_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Max"));

				plot_list[[4]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=nonvacc_avg_chamber_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Avg"));
			});

			row_length <- 4;
		}
		else if(metric == "vacc_echo_chambers")
		{
			DF <- Proportion_vs_Risk_Table[, .SD, .SDcol=c("social_norm", "perceived_vaccine_risk", "soc_V_mean", "vacc_min_chamber_size_mean", "vacc_max_chamber_size_mean", "vacc_avg_chamber_size_mean")];

			local(
			{
				plot_list[[1]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=soc_V_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Vacc"));

				plot_list[[2]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=vacc_min_chamber_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Min"));

				plot_list[[3]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=vacc_max_chamber_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Max"));

				plot_list[[4]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=vacc_avg_chamber_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Avg"));
			});

			row_length <- 4;
		}
		else if(metric == "mutual_info")
		{
			DF <- Proportion_vs_Risk_Table[, .SD, .SDcol=c("social_norm", "perceived_vaccine_risk", "soc_V_mean", "mutual_info_mean")];

			local(
			{
				plot_list[[1]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=soc_V_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Vacc"));

				plot_list[[2]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=mutual_info_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="M"));
			});

			row_length <- 2;
		}
		else if(metric == "conn_chamb_numbers")
		{
			DF <- Proportion_vs_Risk_Table[, .SD, .SDcol=c("social_norm", "perceived_vaccine_risk", "soc_V_mean", "vacc_number_conn_comps_mean", "nonvacc_number_conn_comps_mean", "vacc_number_chambers_mean", "nonvacc_number_chambers_mean")];

			local(
			{
				plot_list[[1]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=soc_V_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Vacc"));

				plot_list[[2]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=vacc_number_conn_comps_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="V Conn"));

				plot_list[[3]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=nonvacc_number_conn_comps_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="NVacc Conn"));

				plot_list[[4]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=vacc_number_chambers_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="V Chamb"));

				plot_list[[5]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=nonvacc_number_chambers_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="NVacc Chamb"));
			});

			row_length <- 5;
		}
		else if(metric == "nonvacc_conn_comps")
		{
			DF <- Proportion_vs_Risk_Table[, .SD, .SDcol=c("social_norm", "perceived_vaccine_risk", "soc_V_mean", "nonvacc_min_conn_comp_size_mean", "nonvacc_max_conn_comp_size_mean", "nonvacc_avg_conn_comp_size_mean")];

			local(
			{
				plot_list[[1]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=soc_V_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Vacc"));

				plot_list[[2]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=nonvacc_min_conn_comp_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Min"));

				plot_list[[3]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=nonvacc_max_conn_comp_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Max"));

				plot_list[[4]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=nonvacc_avg_conn_comp_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Avg"));

			});

			row_length <- 4;
		}
		else if(metric == "vacc_conn_comps")
		{
			DF <- Proportion_vs_Risk_Table[, .SD, .SDcol=c("social_norm", "perceived_vaccine_risk", "soc_V_mean", "vacc_min_conn_comp_size_mean", "vacc_max_conn_comp_size_mean", "vacc_avg_conn_comp_size_mean")];

			local(
			{
				plot_list[[1]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=soc_V_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Vacc"));

				plot_list[[2]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=vacc_min_conn_comp_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Min"));

				plot_list[[3]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=vacc_max_conn_comp_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Max"));

				plot_list[[4]] <<- ggplot(
							DF,
							aes(x=social_norm, y=perceived_vaccine_risk)
						) +
						geom_raster(aes(fill=vacc_avg_conn_comp_size_mean))  +
						scale_colour_gradientn(colours=terrain.colors(20)) +
						xlab(expression(kappa)) + ylab(expression(sigma)) +
						guides(fill=guide_legend(title="Avg"));

			});

			row_length <- 4;
		}
	}


	theme_set(theme_cowplot(font_size=20))

	ggsave(
		filename = paste(graph_file_path, folder_name, "/", metric, file_name_ending, sep=""),
		plot = plot_grid(
			plotlist=plot_list,
			ncol=row_length, align = 'h'
		),
		width=28*row_length, height=28, units="cm",
		limitsize=FALSE
	);

	dev.off()

	rm(plot_list, row_length);

}

end <- Sys.time();

print(end-start);
