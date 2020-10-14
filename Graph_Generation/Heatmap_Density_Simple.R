# Summary Plot for the Topology
rm(list=ls())

user_name <- Sys.info()[8][[1]];
source(sprintf("/media/%s/Simulations/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Heatmaps", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

metrics_to_plot <- c("phys_V_mean", "soc_V_mean") # "geary_c_mean") #, , "NN_join_count_mean", "VV_join_count_mean", "NV_join_count_mean", "mutual_info_mean", "moran_i_mean", "ratio_sec_to_all_infs_at_the_end_mean") # ) # ) # , ); #, "join_counts", "nonvacc_echo_chambers", "vacc_echo_chambers", "mutual_info", "conn_chamb_numbers", "nonvacc_conn_comps", "vacc_conn_comps"); #, "modularity", "watts_strogatz",  "geary_c", "getis_ord", "infections", "opinion");

start <- Sys.time();

Params <- fread(parameter_file, colClasses="character")
Parameters <- Params[num(beta)==1 & num(size)==40000 & num(init_prop)==0.05 & num(infec)==0.2, ];
Parameters[, c("V1", "risk", "norm", "time"):=NULL];
Parameters <- unique(Parameters);

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
	random_opinion_switch <- this_tuple[ , random_switch];
    social_topology <- physical_topology;
    social_degree <- physical_degree;

    DT <- data.table()

	file_list <- Sys.glob(paste(
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
		"_norm_", "*",
		"_switch_", random_opinion_switch,
		".csv",
		sep=""
	));

	for(filename in file_list)
	{
		temp <- data.table(read.csv(filename));
		DT <- rbind(DT, temp[is.na(instance), ], fill=TRUE);
	}

	if(nrow(DT)==0){ next; }

	DT <- unique(DT[
		num(perceived_vaccine_risk)<=0.2 &
		num(perceived_vaccine_risk)>=-1 &
		num(social_norm)>=0 &
		num(social_norm)<=2.4]
	);

	DT <- unique(DT[num(social_norm)%in%seq(0, 2.4, by=0.125) & num(perceived_vaccine_risk)%in%seq(-1, 1, by=0.03125)])

	for(metric in metrics_to_plot)
	{
		file_name_ending <- paste(
			"_N_", size,
			"_vacc_", initial_vacc_prop,
			"_dur_", duration,
			"_inf_", infection_prob,
			"_imp_", importation_rate,
			"_rep_", replenishment_rate,
			"_ptop_", physical_topology,
			"_pdeg_", physical_degree,
			"_stop_", social_topology,
			"_sdeg_", social_degree,
			"_switch_", random_opinion_switch,
			".png",
			sep=""
		);

		p <- ggplot(DT, aes(x=perceived_vaccine_risk, y=social_norm)) + geom_tile(aes(fill=get(metric)))

		stop()
	}
}

end <- Sys.time();

print(end-start);

# a <- list(
#   autotick = FALSE,
#   ticks = "outside",
#   tick0 = 0,
#   dtick = 0.25,
#   ticklen = 5,
#   tickwidth = 2,
#   tickcolor = toRGB("blue")
# )
# norm_range <- seq(0, 2, by=0.25)
# risk_range <- seq(-1, 1, by=0.25)
# # p <- plot_ly(z=DT_mat, type="contour") %>% layout(xaxis = a, yaxis = a)
# p <- plot_ly(x=~norm_range, z=~risk_range, type="contour") %>% add_markers() %>%
#   add_markers(y = ~rev(s)) %>%
#   layout(xaxis = a, yaxis = a)
#
# contour.mat <- ifelse(DT_mat < 20000, 0, DT_mat)
# contour.mat2 <- ifelse(DT_mat < 30000, 0, DT_mat)
#
# filled.contour(DT_mat, color = terrain.colors,
#            plot.axes = contour(contour.mat, levels = 1,
#                                drawlabels = FALSE, axes = FALSE,
#                                frame.plot = FFALSE, add = TRUE) +
# 							   contour(contour.mat2, levels = 1,
# 				                                   drawlabels = FALSE, axes = FALSE,
# 				                                   frame.plot = FFALSE, add = TRUE))
