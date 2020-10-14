# This is the file for the first diagram tableau in the figure - the oen for multiple transitions is in another file
rm(list=ls())

options(scipen=10000)

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Skew_vs", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# source("Processing_Parameter_Table.R");
Params <- fread(parameter_file, colClasses="character")

graph_width <- 10;
graph_height <- 3;
axis_text_size <- 35;
axis_title_size <- 35;
point_size <- 30;

metrics_to_plot <- c("geary_c", "mutual_info", "NV_join_count", "moran_i", "NN_join_count", "VV_join_count");
metric_names <- c("<G>", "<M>", "<N,Vs>", "<I>", "<N,N>", "<Vs,Vs>");
# , "getis_ord""echo_chamber_number"); #, "echo_chamber_number", "conn_comp_number", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "vacc_conn_comp_size", "nonvacc_conn_comp_size"); # , "forced_new_ratio", "join_count",

time_start <- Sys.time();

for(the_size in c("10000", "40000")) #
{
	Parameters <- Params[num(size)==the_size & num(beta)==1 & soc_top=="random" & phys_top=="random" & num(init_prop)==0.05];

	if(num(the_size) == "40000") Parameters <- Parameters[num(infec)==0.2];
	if(num(the_size) == "10000") Parameters <- Parameters[num(import)==0.00025]

	norm_vector <- Parameters[, unique(norm)];
	infec_vector <- Parameters[, unique(infec)];

	Parameters[, c("V1", "risk", "time", "norm"):=NULL];
	Parameters <- unique(Parameters);

	for(line_index in 1:nrow(Parameters))
	{
		this_tuple <- Parameters[line_index];

		beta <- this_tuple[, beta];
		initial_vacc_prop <- this_tuple[, init_prop];
		size <- this_tuple[, size];
		duration <- this_tuple[, duration];
		importation_rate <- this_tuple[, import];
		replenishment_rate <-  this_tuple[, birth_death];
		physical_topology <- this_tuple[, phys_top];
		physical_degree <- this_tuple[, phys_deg];
		social_topology <- this_tuple[, soc_top];
		social_degree <- this_tuple[, soc_deg];
		random_opinion_switch <- this_tuple[, random_switch];
		infection_prob <- this_tuple[, infec];

		num_size <- num(size);

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
			"_switch_", random_opinion_switch,
			".png",
			sep=""
			);

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
			"_norm_", "*",
			"_switch_", random_opinion_switch,
			".csv",
			sep=""
		));

		if(length(filenames) == 0){ next; }

		Proportion <- data.table();
		for(file in filenames)
		{
			temp <- data.table(read.csv(file));
			Proportion <- rbind(Proportion, temp[is.na(instance), ], fill=TRUE);
		}
		Proportion[, X:=NULL];
		Proportion <- Proportion[social_norm<=3];

		for(metric_index in seq_along(metrics_to_plot))
		{
			metric <- metrics_to_plot[metric_index];

			Skew_vs <- data.table(
				norm = numeric(),
				gap = numeric(),
				skew = numeric()
			);

			Kurtosis_vs <- data.table(
				norm = numeric(),
				gap = numeric(),
				skew = numeric()
			)

			for(norm in unique(Proportion$social_norm))
			{
				DT <- Proportion[social_norm==norm][order(perceived_vaccine_risk)];
				Trans_Table <- get_transitions(DT, initial_vacc_prop);

				trans_gap <- Trans_Table[, phys_intersec_risk-soc_intersec_risk];

				if(length(DT) == 0) next;

				skew <- skewness(DT[[ sprintf("%s_mean", metric) ]]);
				Skew_vs <- rbind(Skew_vs, list(norm, trans_gap, skew));

				kurt <- kurtosis(DT[[ sprintf("%s_mean", metric) ]]);
				Kurtosis_vs <- rbind(Kurtosis_vs, list(norm, trans_gap, kurt));
			}

			skew_norm <- ggplot(Skew_vs, aes(x=norm, y=skew)) +
				geom_point(aes(size=point_size)) +
				geom_line(size=2, colour="blue") +
				# geom_smooth(method="loess", se=F, size=2) +
				xlab(expression(sigma)) +
				ylab(bquote(gamma['1']*"("*.(metric_names[metric_index])*")")) +
				theme(
					# panel.background=element_blank(),
					axis.text=element_text(size=axis_text_size),
					axis.title=element_text(size=axis_title_size),
					legend.position = "none",
					axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
					axis.line.y = element_line(colour="black", size=0.5, linetype="solid")
				);
			norm_plot_name <- paste(graph_file_path, folder_name, "/skew_v_norm_", metric, file_name_ending, sep="");

			ggsave(
				norm_plot_name,
				plot=skew_norm, width=graph_width, height=graph_height, dpi=50, limitsize=FALSE
			);
			dev.off()
			system(sprintf("convert %s -trim %s", norm_plot_name, norm_plot_name));

			skew_gap_outliers <- boxplot(Skew_vs$gap, plot=F)$out;
			Skew_vs <- Skew_vs[!is.na(gap) & !(num(gap)%in%skew_gap_outliers)][order(gap)];
			skew_gap <- ggplot(Skew_vs, aes(x=gap, y=skew)) +
				geom_point(aes(size=point_size)) +
				xlab(expression('K'['p']*'-K'['s'])) +
				ylab(bquote(gamma['1']*"("*.(metric_names[metric_index])*")")) +
				geom_smooth(method="loess", se=F, size=2, colour="blue") +
				theme(
					# panel.background=element_blank(),
					axis.text=element_text(size=axis_text_size),
					axis.title.y=element_text(size=axis_title_size),
					axis.title.x=element_text(size=axis_title_size-10),
					legend.position = "none",
					axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
					axis.line.y = element_line(colour="black", size=0.5, linetype="solid")
				);
			gap_plot_name <- paste(graph_file_path, folder_name, "/skew_v_gap_", metric, file_name_ending, sep="");

			ggsave(
				gap_plot_name,
				plot=skew_gap, width=graph_width, height=graph_height, limitsize=FALSE, dpi=50
			);
			dev.off()
			system(sprintf("convert %s -trim %s", gap_plot_name, gap_plot_name));

			# kurt_norm <- ggplot(Kurtosis_vs, aes(x=norm, y=skew)) +
			# 	geom_point(aes(size=point_size)) +
			# 	geom_smooth(method="loess", se=F) +
			# 	xlab(expression(sigma)) +
			# 	ylab(bquote(gamma['2'])) +
			# 	theme(
			# 		# panel.background=element_blank(),
			# 		axis.text=element_text(size=axis_text_size),
			# 		axis.title=element_text(size=axis_title_size),
			# 		legend.position = "none",
			# 		axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
			# 		axis.line.y = element_line(colour="black", size=0.5, linetype="solid")
			# 	);
			#
			# ggsave(
			# 	paste(graph_file_path, folder_name, "/kurt_v_norm_", metric, file_name_ending, sep=""),
			# 	plot=kurt_norm, width=graph_width, height=graph_height, units="in", limitsize=FALSE
			# );
			#
			# kurt_gap_outliers <- boxplot(Kurtosis_vs$gap, plot=F)$out;
			# Kurtosis_vs <- Kurtosis_vs[!is.na(gap) & !(num(gap)%in%kurt_gap_outliers)][order(gap)];
			#
			# kurt_gap <- ggplot(Kurtosis_vs, aes(x=gap, y=skew)) +
			# 	geom_point(aes(size=point_size)) +
			# 	xlab(expression('K'['p']*'-K'['s'])) +
			# 	ylab(bquote(gamma['2'])) +
			# 	geom_smooth(method="loess", se=F) +
			# 	theme(
			# 		# panel.background=element_blank(),
			# 		axis.text=element_text(size=axis_text_size),
			# 		axis.title=element_text(size=axis_title_size),
			# 		legend.position = "none",
			# 		axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
			# 		axis.line.y = element_line(colour="black", size=0.5, linetype="solid")
			# 	);
			#
			# ggsave(
			# 	paste(graph_file_path, folder_name, "/kurt_v_gap_", metric, file_name_ending, sep=""),
			# 	plot=kurt_gap, width=graph_width, height=graph_height, units="in", limitsize=FALSE
			# );
		}
	}
}

print(Sys.time()-time_start);
