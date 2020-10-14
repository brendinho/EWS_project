# Brendon Phillips
# PhD candidate
# Bahc lab cumputational epidemiology group
# Department of Mathematics
# University of Waterloo

rm(list=ls());

start <- Sys.time();

user_name <- Sys.info()[8][[1]];
source(sprintf("/media/%s/Simulations/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Cloud_Plots", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# axis_title2 <- function(x)
# {
# 	if(x == "phys_S_mean"){ return("Susceptibles"); }
# 	if(x == "phys_I_mean"){ return("Infecteds"); }
# 	if(x == "phys_R_mean"){ return("Recovereds"); }
# 	if(x == "phys_V_mean"){ return("Vaccinateds"); }
#
# 	if(x == "soc_V_mean"){  return("Vaxxers"); }
#
# 	if(x == "mutual_info_mean"){ return("Mutual Information"); }
# 	if(x %in% c("moran_i_mean", "morran_i_calc_mean")){ return("Moran's I"); }
#
# 	if(x == "NN_join_count_mean"){ return("[N,N]"); }
# 	if(x == "NV_join_count_mean"){ return("[N,V]"); }
# 	if(x == "VV_join_count_mean"){ return("[V,V]"); }
# }
#
# max_value <- function(x, size)
# {
# 	if( grepl("phys", x) ){ return(size); }
# 	if( grepl("soc", x) ){ return(size); }
# 	if( grepl("join", x) ){ return(size*(size)); }
# 	return(1);
# }


Params <- fread(parameter_file, colClasses="character") # as.numeric(size)!=10000 &
Parameters <- Params[num(beta)==1 & phys_top!="lattice" & soc_top!="lattice" & num(init_prop)==0.05, ];
Parameters[, c("V1", "risk", "random_switch", "norm", "time", "infec", "import", "phys_deg", "soc_deg"):=NULL]; # "init_prop"
Parameters <- unique(Parameters);

metrics_to_plot <- c("soc_V_mean", "moran_i_mean", "moran_i_calc_mean", "mutual_info_mean", "phys_R_mean", "phys_V_mean", "phys_S_mean", "phys_I_mean", "NN_join_count_mean", "NV_join_count_mean", "VV_join_count_mean");

for(line_index in 1:nrow(Parameters))
# foreach( line_index=1:nrow(Parameters), .packages=c("data.table")) %dopar%
{
	this_tuple <- Parameters[line_index];

	initial_vacc_prop <- this_tuple[, init_prop];
	size <- this_tuple[, size];
	duration <- this_tuple[, duration];
	beta <- this_tuple[, beta];
	# infection_prob <- this_tuple[, infec];
	# importation_rate <- this_tuple[, import];
	physical_topology <- this_tuple[, phys_top];
	# physical_degree <- this_tuple[, phys_deg];
	social_topology <- this_tuple[, soc_top];
	# social_degree <- this_tuple[, soc_deg];

	filenames <- Sys.glob(paste(
		CSV_file_path,
		"Summary",
		"_N_", size,
		"_dur_", duration,
		"_beta_", beta,
		"_vaccprop_", initial_vacc_prop,
		"_risk_", "*",
		# "_inf_", infection_prob,
		"_inf_", "*",
		# "_imp_", importation_rate,
		"_imp_", "*",
		"_rep_", "*",
		"_ptop_", physical_topology,
		# "_pdeg_", physical_degree,
		"_pdeg_", "*",
		"_stop_", social_topology,
		# "_sdeg_", social_degree,
		"_sdeg_", "*",
		"_norm_", "*",
		"_switch_", "*",
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

	for(metric in metrics_to_plot)
	{
		if(!isTRUE( metric %in% names(Summary_Data) )){ metrics_to_plot <- metrics_to_plot[!(metrics_to_plot==metric)] }
	}

	# writeLines(paste("\treplenishment ratio: ", replenishment_rate, sep=""));
	# writeLines(paste("\timportation rate: ", importation_rate, sep=""));
	# writeLines(paste("\tsocial norm: ", social_norm, sep=""));
	writeLines(paste("\tnetwork size: ", size, sep=""));
	# writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
	writeLines(paste("\tphysical topology: ", physical_topology, sep=""));
	# writeLines(paste("\tphysical degree: ", physical_degree, sep=""));
	writeLines(paste("\tinitial vaccinator proportion: ", initial_vacc_prop, sep=""));
	writeLines(paste("\tsocial topology: ", social_topology, sep=""));
	# writeLines(paste("\tsocial degree: ", social_degree, sep=""));
	writeLines(paste("\tbeta: ", beta, sep=""));
	# writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));

	Summary_Data <- Summary_Data[, .SD, .SDcols=c(fixed_values, metrics_to_plot)];

	# Summary_Data <- Summary_Data[, X:=NULL];
	Summary_Data[, "colour":=""];

	Parameter_Tuples <- Summary_Data[, unique(.SD), .SDcols=setdiff(fixed_values, c("instance"))];

	# colour_vector <- colorRampPalette(brewer.pal(11, "Paired"))(nrow(Parameter_Tuples));
	Summary_Data <- Summary_Data[-1<=perceived_vaccine_risk & perceived_vaccine_risk<=1, ];

	norm_vector <- Summary_Data[, sort(unique(social_norm))];
	colour_vector <- green2red(length(norm_vector));
	colour_row <- which(names(Summary_Data) == "colour");

	for(i in 1:nrow(Summary_Data))
	{
		set(
			Summary_Data,
			i,
			colour_row,
			colour_vector[Summary_Data[i, which(social_norm == norm_vector)]]
		);
	}

	file_name_ending <- paste(
		"_N_", size,
		"_dur_", duration,
		"_beta_", beta,
		"_vaccprop_", initial_vacc_prop,
		"_ptop_", physical_topology,
		"_stop_", social_topology,
		".png",
		sep=""
	);

	# find the end of the zone of bistability - we define it by the social transition
	# suspected_end <- Summary_Data[soc_V_mean>0.8*max(soc_V_mean), .SD, .SDcols=c("soc_V_mean", "perceived_vaccine_risk")][, max(perceived_vaccine_risk)];
	# sd_threshold <- Summary_Data[, sd(c(0.2*as.numeric(size), 0.8*as.numeric(size)))]
	#
	# standard_deviations <- numeric();
	# for(index in seq_along(risk_vector))
	# {
	# 	standard_deviations[index] <- Summary_Data[perceived_vaccine_risk==risk_vector[index], sd(soc_V_mean)];
	# }
	# suspected_ends <- risk_vector[which(standard_deviations > sd_threshold)];

	risk_vector <- Summary_Data[, sort(unique(perceived_vaccine_risk))];
	suspected_ends <- numeric();

	for(risk in risk_vector)
	{
		high_risks <- risk_vector[which(Summary_Data[perceived_vaccine_risk==risk, soc_V_mean>0.8*as.numeric(size)])];
		low_risks <- risk_vector[which(Summary_Data[perceived_vaccine_risk==risk, soc_V_mean<0.2*as.numeric(size)])];

		if(isTRUE( (length(high_risks)>0) & (length(low_risks)>0) )){ suspected_ends <- c(suspected_ends, risk); }
	}

	# calculate the lagged difference of the standard deviations to find the greatest difference
	# # have to use sd - can't use ordinary difference between values, since each risk valdue corresponds to many values in this clooud plot
	# # use to set the width of the shaded ribbon in the cloud plota below
	# differences <- abs(diff(standard_deviations, lag=1));
	# suspected_ends <- risk_vector[which(match(differences, sort(differences, TRUE)[1:4])!=0)];

	for(metric in metrics_to_plot)
	{
		print(metric);
		denom <- 1;
		if(grepl("soc", metric)) { denom <- as.numeric(size); }
		else if(grepl("phys", metric)) { denom <- as.numeric(size); }
		else if(grepl("join", metric)) { denom <- Summary_Data[, maxx(NN_join_count_mean+NV_join_count_mean+VV_join_count_mean)]; }

		the_plot <- ggplot() +
			geom_rect(aes(
						ymin=Summary_Data[, minn(get(metric)/denom)],
						ymax=Summary_Data[, maxx(get(metric)/denom)],
						xmin=minn(suspected_ends),
						xmax=maxx(suspected_ends)),
						alpha=0.6, fill='gray'
					) +
				geom_point(data=Summary_Data, aes(x=perceived_vaccine_risk, y=get(metric)/denom, col=colour), size=2) +
				theme(legend.position = "none", panel.background=element_blank(), axis.line=element_line(colour = "black"), axis.text=element_text(size=50), axis.title=element_text(size=70), axis.title.y = element_blank()) +
				labs(x=expression(kappa))


		metric_name <- metric;
		if(metric == "moran_i_calc_mean"){ metric_name <- "moran_i_mean"; }
		ggsave(paste(graph_file_path, folder_name, "/", metric_name, file_name_ending, sep=""), plot=the_plot, height=10, width=20)
	}

}

end <- Sys.time();

print(end-start);
