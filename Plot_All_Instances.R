# Brendon Phillips
# PhD candidate
# Bahc lab cumputational epidemiology group
# Department of Mathematics
# University of Waterloo

rm(list=ls());

start <- Sys.time();

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Cloud_Plots", sep="");
dir.create(sprintf("%s%s", graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character"); # as.numeric(size)!=10000 &
Params[, c("V1", "risk"):=NULL];
Parameters <- Params[num(size)==10000 &num(beta)==1 & phys_top!="lattice" & soc_top!="lattice" & num(init_prop)==0.05 & num(norm)%in%c(0,0.5)];
Parameters <- unique(Parameters);

# metrics_to_plot <- c("soc_V_mean", "moran_i_mean", "moran_i_calc_mean", "mutual_info_mean", "phys_R_mean", "phys_V_mean", "phys_S_mean", "phys_I_mean", "NN_join_count_mean", "NV_join_count_mean", "VV_join_count_mean");

metrics_to_plot <- c("watts_strogatz"); # "watts_strogatz_nonvacc_mean", "watts_strogatz_vacc_mean", "watts_strogatz_all_mean", "vacc_number_conn_comps_mean", "nonvacc_number_conn_comps_mean", "vacc_number_chambers_mean", "nonvacc_number_chambers_mean");

for(line_index in 1:nrow(Parameters))
{
	this_tuple <- Parameters[line_index];

	steepness <- this_tuple[, beta];
	social_norm <- this_tuple[, norm];
	replenishment_rate <-  this_tuple[, birth_death];
	initial_vacc_prop <- this_tuple[, init_prop];
	size <- this_tuple[, size];
	duration <- this_tuple[, duration];
	infection_prob <- this_tuple[, infec];
	importation_rate <- this_tuple[, import];
	physical_topology <- this_tuple[, phys_top];
	physical_degree <- this_tuple[, phys_deg];
	social_topology <- this_tuple[, soc_top];
	social_degree <- this_tuple[, soc_deg];
	random_switch <- this_tuple[, random_switch];

	filenames <- Sys.glob(paste(
		CSV_file_path,
		"Summary",
		"_N_", size,
		"_dur_", duration,
		"_beta_", steepness,
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
		"_switch_", random_switch,
		".csv",
		sep=""
	));

	if(length(filenames) == 0){ next; }

	Summary_Data <- data.table();
	for(file in unique(filenames))
	{
		temp <- data.table(read.csv(file));
		Summary_Data <- rbind(Summary_Data, temp[!is.na(instance), ], fill=TRUE);
	}
	if(nrow(Summary_Data)==0){ next; }

	writeLines(paste("\treplenishment ratio: ", replenishment_rate, sep=""));
	writeLines(paste("\timportation rate: ", importation_rate, sep=""));
	writeLines(paste("\tsocial norm: ", social_norm, sep=""));
	writeLines(paste("\tnetwork size: ", size, sep=""));
	writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
	writeLines(paste("\tphysical topology: ", physical_topology, sep=""));
	writeLines(paste("\tphysical degree: ", physical_degree, sep=""));
	writeLines(paste("\tinitial vaccinator proportion: ", initial_vacc_prop, sep=""));
	writeLines(paste("\tsocial topology: ", social_topology, sep=""));
	writeLines(paste("\tsocial degree: ", social_degree, sep=""));
	writeLines(paste("\tbeta: ", steepness, sep=""));
	writeLines(paste("\topinion switching rate: ", random_switch, sep=""));

	# Summary_Data <- Summary_Data[, X:=NULL];
	Summary_Data[, "colour":=""];

	Parameter_Tuples <- Summary_Data[, unique(.SD), .SDcols=setdiff(fixed_values, c("instance"))];

	# colour_vector <- colorRampPalette(brewer.pal(11, "Paired"))(nrow(Parameter_Tuples));
	Summary_Data <- Summary_Data[-1<=perceived_vaccine_risk & perceived_vaccine_risk<=1, ];

	file_name_ending <- paste(
		"_N_", size,
		"_risk_", social_norm,
		"_dur_", duration,
		"_beta_", steepness,
		"_vaccprop_", initial_vacc_prop,
		"_inf_", infection_prob,
		"_imp_", importation_rate,
		"_rep_", replenishment_rate,
		"_top_", physical_topology,
		"_deg_", physical_degree,
		"_switch_", random_switch,
		".png",
		sep=""
	);

	risk_vector <- Summary_Data[, sort(unique(perceived_vaccine_risk))];
	suspected_ends <- numeric();

	for(risk in risk_vector)
	{
		high_risks <- risk_vector[which(Summary_Data[perceived_vaccine_risk==risk, soc_V_mean>0.8*as.numeric(size)])];
		low_risks <- risk_vector[which(Summary_Data[perceived_vaccine_risk==risk, soc_V_mean<0.2*as.numeric(size)])];

		if(isTRUE( (length(high_risks)>0) & (length(low_risks)>0) )){ suspected_ends <- c(suspected_ends, risk); }
	}

	for(metric in metrics_to_plot)
	{
		Table_Here <- melt(Summary_Data[, .SD, .SDcols=c("perceived_vaccine_risk", sprintf("%s_mean", plot_these(metric)))], id="perceived_vaccine_risk");

		print(metric);
		denom <- 1;
		if(grepl("soc", metric)) { denom <- as.numeric(size); }
		else if(grepl("phys", metric)) { denom <- as.numeric(size); }
		else if(grepl("join", metric)) { denom <- Summary_Data[, maxx(NN_join_count_mean+NV_join_count_mean+VV_join_count_mean)]; }

		the_plot <- ggplot() +
			geom_rect(aes(
						ymin=-Inf,
						ymax=Inf,
						xmin=minn(suspected_ends),
						xmax=maxx(suspected_ends)),
						alpha=0.6, fill='gray'
				) +
			geom_point(data=Table_Here, aes(x=perceived_vaccine_risk, y=value/denom, colour=variable), size=2) +
			geom_spline(aes(x=perceived_vaccine_risk, y=value, colour=variable), data=Table_Here) +
			theme(legend.position = "none", panel.background=element_blank(), axis.line=element_line(colour = "black"), axis.text=element_text(size=50), axis.title=element_text(size=70), axis.title.y = element_blank()) + labs(x=expression(kappa), y=axis_label(metric))

		metric_name <- metric;
		ggsave(paste(graph_file_path, folder_name, "/", metric_name, file_name_ending, sep=""), plot=the_plot, height=7, width=20);

		dev.off()
	}

	graphics.off()
}

end <- Sys.time();

print(end-start);

# melt(Summary_Data[, .SD, .SDcols=c("perceived_vaccine_risk", filter_these(names(Summary_Data), yes=c("watts"), no=c("_sd")))], id="perceived_vaccine_risk")
