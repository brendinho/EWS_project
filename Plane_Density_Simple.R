# Summary Plot for the Topology
rm(list=ls())

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Heatmaps", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

metrics_to_plot <- c("moran_i"); #"moran_i", "join_counts", "nonvacc_echo_chambers", "vacc_echo_chambers", "mutual_info", "conn_chamb_numbers", , "vacc_conn_comps"); #, "modularity", "watts_strogatz",  "geary_c", "getis_ord", "infections", "opinion");

start <- Sys.time();

Params <- fread(parameter_file, colClasses="character")
Parameters_all <- Params[num(beta)==1 & num(size)==40000 & num(init_prop)==0.05 & num(infec)==0.2, ];
Parameters <- Parameters_all;
Parameters[, c("V1", "risk", "norm", "time"):=NULL]; # , "connectivity"
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

	# for(risk in sprintf("%.4f", seq(-1, 1, by=0.03125))){
	# for(social_norm in seq(0, 3, by=0.125))
	# {
	# 	filename <- paste(
	# 		CSV_file_path,
	# 		"Summary",
	# 		"_N_", size,
	# 		"_dur_", duration,
	# 		"_beta_", beta,
	# 		"_vaccprop_", initial_vacc_prop,
	# 		"_risk_", risk,
	# 		"_inf_", infection_prob,
	# 		"_imp_", importation_rate,
	# 		"_rep_", replenishment_rate,
	# 		"_ptop_", physical_topology,
	# 		"_pdeg_", physical_degree,
	# 		"_stop_", social_topology,
	# 		"_sdeg_", social_degree,
	# 		"_norm_", social_norm,
	# 		"_switch_", random_opinion_switch,
	# 		".csv",
	# 		sep=""
	# 	);
	#
	# 	if( file.exists(filename) )
	# 	{
	# 		temp <- data.table(read.csv(filename));
	# 		DT <- rbind(DT, temp[is.na(instance), ], fill=TRUE);
	# 	}
	# }}

	# # # # # # write.csv(DT, file=sprintf("%sN_%s_table.csv", CSV_file_path, size));

	# if(nrow(DT)==0){ next; }

	DT <- fread(sprintf("%sN_40000_table.csv", CSV_file_path));
	DT[, "V1":=NULL];
	DT[, soc_V_mean:=soc_V_mean/N];
	DT[, phys_V_mean:=phys_V_mean/N];

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

	for (metric in c("soc_V_mean", "phys_V_mean")) #, "phys_I_mean"))
	{
		# ATTRIBUTION: https://colinbousige.github.io/post/3dplotr/
		DT.interp <- with(DT,
			interp(x = social_norm, y = perceived_vaccine_risk, z = get(metric),
				duplicate="median",
				xo=seq(min(DT$social_norm), max(DT$social_norm), length = 100),
				yo=seq(min(DT$perceived_vaccine_risk), max(DT$perceived_vaccine_risk), length = 100),
				extrap=FALSE, linear=FALSE)
			);
		# Regrouping this list to a 3-columns data.frame
		melt_x <- rep(DT.interp$x, times=length(DT.interp$y));
		melt_y <- rep(DT.interp$y, each=length(DT.interp$x));
		melt_z <- as.vector(DT.interp$z);
		DT.smooth <- data.table(na.omit(data.frame(social_norm=melt_x, perceived_vaccine_risk=melt_y, Z=melt_z)));
		DT.smooth <- DT.smooth[perceived_vaccine_risk<=0.25 & perceived_vaccine_risk>=-0.8];
		DT.smooth <- DT.smooth[social_norm>head(DT.smooth$social_norm, 1) & social_norm<tail(DT.smooth$social_norm, 1)];

		break_vector <- c(0, 0.5, 1)
		if(metric == "phys_I_mean"){ break_vector <- c(minn(DT.smooth$Z), maxx(DT.smooth$Z), mean(minn(DT.smooth$Z), maxx(DT.smooth$Z))); }

		pl <- ggplot(
				data = DT.smooth,
				aes(x=social_norm, y=perceived_vaccine_risk)
			) +
			geom_raster(aes(fill = Z)) +
			scale_fill_gradientn(colours=colorRampPalette(c("seagreen","brown"))(500), breaks=break_vector) +
			# scale_fill_gradientn(colours=colorRampPalette(c("red", "yellow", "blue"))(500), breaks=c(0, 0.5, 1)) +
			# scale_fill_gradientn(colours=colorBrewerPalette("RdYlBu")) +
			# scale_fill_distiller(palette="Spectral", direction=1) +
			# scale_fill_gradientn(colours = terrain.colors(10)) +
			# ggtitle("Number of agents with pro-vaccine sentiment") +
			labs(x=expression(sigma), y=expression(kappa), fill=axis_label(metric))
			theme_bw();

		theme_set(theme_cowplot(font_size=20))

		ggsave(
			filename = paste(graph_file_path, folder_name, "/", metric, "_DENSITY", file_name_ending, sep=""),
			plot = pl,
			width=28, height=20, units="cm",
			limitsize=FALSE
		);

		dev.off();
	}

	# for(metric in c("ratio_sec_to_all_infs_at_the_end_mean", "NV_join_count_mean", "moran_i_mean", "mutual_info_mean", "VV_join_count_mean", "NN_join_count_mean", "nonvacc_opinion_change_mean", "vacc_opinion_change_mean"))
 	# {
	# 	if(metric == "ratio_sec_to_all_infs_at_the_end_mean")
	# 	{ # the cmoothing process interp() doesn't take NAs, so set them all to zero
	# 		set(
	# 			DT,
	# 			which(is.na( DT[[metric]] )),
	# 			which(names(DT)==metric),
	# 			0
	# 		);
	# 	}
	#
	# 	# ATTRIBUTION: https://colinbousige.github.io/post/3dplotr/
	# 	DT.interp <- with(DT,
	# 		interp(x = social_norm, y = perceived_vaccine_risk, z = get(metric),
	# 			duplicate="median",
	# 			xo=seq(min(DT$social_norm), max(DT$social_norm), length = 100),
	# 			yo=seq(min(DT$perceived_vaccine_risk), max(DT$perceived_vaccine_risk), length = 100),
	# 			extrap=FALSE, linear=FALSE)
	# 		);
	# 	# Regrouping this list to a 3-columns data.frame
	# 	melt_x <- rep(DT.interp$x, times=length(DT.interp$y));
	# 	melt_y <- rep(DT.interp$y, each=length(DT.interp$x));
	# 	melt_z <- as.vector(DT.interp$z);
	# 	DT.smooth <- data.table(na.omit(data.frame(social_norm=melt_x, perceived_vaccine_risk=melt_y, Z=melt_z)));
	#
	# 	data_table <- DT.smooth;
	# 	if(metric == "nonvacc_opinion_change_mean"){ data_table <- DT.smooth[Z>0] }
	#
	# 	# DT.smooth <- DT[, .SD, .SDcols=c("social_norm", "perceived_vaccine_risk", metric)];
	# 	# names(DT.smooth) <- c("social_norm", "perceived_vaccine_risk", "Z");
	# 	# data_table <- DT.smooth
	#
	# 	data_table <- data_table[perceived_vaccine_risk<=0.25 & perceived_vaccine_risk>=-0.8];
	# 	data_table <- data_table[social_norm>head(data_table$social_norm, 1) & social_norm<tail(data_table$social_norm, 1)];
	#
	# 	if(metric == "ratio_sec_to_all_infs_at_the_end_mean")
	# 	{ # scrub the data to get rid if small negative values
	# 		set(data_table, which(data_table$Z<0), which(names(data_table)=="Z"), 0);
	# 	}
	#
	# 	pl <- ggplot(
	# 			data = data_table,
	# 			aes(x=social_norm, y=perceived_vaccine_risk)
	# 		) +
	# 		geom_raster(aes(fill = Z)) +
	# 		ggtitle(axis_label(metric)) +
	# 		labs(x=expression(sigma), y=expression(kappa), fill=axis_label(metric)) +
	# 		theme_bw();
	#
	# 	if(metric == "NV_join_count_mean")
	# 	{
	# 		pl <- pl + scale_fill_gradientn(
	# 				# colours = colorRampPalette(c("orange", "blue"))(500),
	# 				# colours = brewer.pal(11, "Spectral"),
	# 				colours = terrain.colors(500),
	# 				breaks = c(0, 500, DT.smooth[, max(Z)]),
	# 				values = scales::rescale(c(1, 500, DT.smooth[, max(Z)]))
	# 			);
	# 	}
	# 	else if(metric %in% c("nonvacc_opinion_change_mean", "vacc_opinion_change_mean"))
	# 	{
	# 		pl <- pl + scale_fill_gradientn(
	# 				# colours = colorRampPalette(c("orange", "blue"))(500),
	# 				colours = brewer.pal(11, "Spectral"),
	# 				# colours = terrain.colors(500),
	# 				# breaks = c(0, 51000, DT.smooth[, max(Z)]),
	# 				values = scales::rescale(c(0, 10000, DT.smooth[, max(Z)]))
	# 			);
	# 	}
	# 	else if(metric == "ratio_sec_to_all_infs_at_the_end_mean")
	# 	{
	# 		pl <- pl + scale_fill_gradientn(
	# 				# colours = colorRampPalette(c("orange", "blue"))(500),
	# 				colours = terrain.colors(500),
	# 				breaks=c(min(DT.smooth[, minn(Z)]), min(DT.smooth[, maxx(Z)]))
	# 			);
	# 	}
	# 	else
	# 	{
	# 		pl <- pl + scale_fill_gradientn(
	# 				# colours = colorRampPalette(c("orange", "blue"))(500),
	# 				colours = terrain.colors(500),
	# 				breaks=c(min(DT.smooth[, minn(Z)]), min(DT.smooth[, maxx(Z)]))
	# 			);
	# 	}
	#
	# 	theme_set(theme_cowplot(font_size=20));
	#
	# 	ggsave(
	# 		filename = paste(graph_file_path, folder_name, "/", metric, "_DENSITY", file_name_ending, sep=""),
	# 		plot = pl,
	# 		width=28, height=20, units="cm",
	# 		limitsize=FALSE
	# 	);
	#
	# 	dev.off();
	# }
}

end <- Sys.time();

print(end-start);
