# Brendon Phillips
# PhD candidate
# Bahc lab cumputational epidemiology group
# Department of Mathematics
# University of Waterloo

rm(list=ls());

start_time <- Sys.time();

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Cloud_Plots", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character"); Params[, V1:=NULL];
Parameters <- unique(Params[network_structure=="smallworld"]);
Parameters[, c("perceived_vaccine_risk", "random_opinion_switch", "social_norm", "infec_prob", "importation", "mean_degree", "initial_vacc_proportion", "proportion_of_nodes"):=NULL]; #
Parameters <- unique(Parameters);

metrics_to_plot <- c('phys_S_mean', 'phys_I_mean', 'phys_R_mean', 'phys_V_mean', 'soc_H_mean', 'soc_N_mean', 'soc_V_mean'); #, 'mutual_info_mean', 'watts_strogatz_hesitant_mean', 'watts_strogatz_nonvacc_mean', 'watts_strogatz_vacc_mean', 'watts_strogatz_all_mean', 'prob_sick_hesitant_mean', 'prob_sick_nonvacc_mean', 'prob_sick_vacc_mean', 'HH_join_count_mean', 'NN_join_count_mean', 'VV_join_count_mean', 'HN_join_count_mean', 'HV_join_count_mean', 'NV_join_count_mean', 'hesitant_min_conn_comp_size_mean', 'hesitant_max_conn_comp_size_mean', 'hesitant_avg_conn_comp_size_mean', 'hesitant_number_conn_comps_mean', 'nonvacc_min_conn_comp_size_mean', 'nonvacc_max_conn_comp_size_mean', 'nonvacc_avg_conn_comp_size_mean', 'nonvacc_number_conn_comps_mean', 'vacc_min_conn_comp_size_mean', 'vacc_max_conn_comp_size_mean', 'vacc_avg_conn_comp_size_mean', 'vacc_number_conn_comps_mean', 'hesitant_min_chamber_size_mean', 'hesitant_max_chamber_size_mean', 'hesitant_avg_chamber_size_mean', 'hesitant_number_chambers_mean', 'nonvacc_min_chamber_size_mean', 'nonvacc_max_chamber_size_mean', 'nonvacc_avg_chamber_size_mean', 'nonvacc_number_chambers_mean', 'vacc_min_chamber_size_mean', 'vacc_max_chamber_size_mean', 'vacc_avg_chamber_size_mean', 'vacc_number_chambers_mean', 'hesitant_opinion_change_mean', 'nonvacc_opinion_change_mean', 'vacc_opinion_change_mean', 'total_opinion_change_mean', 'ratio_hesitant_infected_neighbours_mean', 'ratio_nonvacc_infected_neighbours_mean', 'ratio_vacc_infected_neighbours_mean', 'num_hesitant_triads_mean', 'num_nonvacc_triads_mean', 'num_vacc_triads_mean', 'total_num_triads_mean', 'hesitant_diameter_mean', 'nonvacc_diameter_mean', 'vacc_diameter_mean');

stop();

for(line_index in 1:nrow(Parameters))
# foreach( line_index=1:nrow(Parameters), .packages=c("data.table")) %dopar%
{
	this_tuple <- Parameters[line_index];

	# initial_vacc_prop <- this_tuple$initial_vacc_proportion;
	size <- this_tuple$N;
	structure <- this_tuple$network_structure;
	filenames <- Sys.glob(paste(
		CSV_file_path,
		"Hes_N_", size,
		"_struct_", structure,
		"_deg_", "*",
		"_risk_", "*",
		"_infec_", "*",
		"_norm_", "*",
		"_imp_", "*",
		"_init_", "*",
		"_switch_", "*",
		"_prop_", "*",
		"_summary.csv",
		sep=""
	));

	if(length(filenames) == 0){ next; }

	col_names <- c(hes_fixed_values, "instance", metrics_to_plot);
	Summary_Data <- data.table(matrix(NaN, nrow=length(filenames)*100, ncol=length(col_names))); # assuming 100 instances per parameter set
	names(Summary_Data) <- col_names;
	Summary_Data[, network_structure:=as.character(network_structure)];

	row_number <- 1;
	for(file in filenames)
	{
		temp <- data.table(read.csv(file))[!is.na(instance), .SD, .SDcols=names(Summary_Data)];
		for(the_inst in 1:nrow(temp))
		{
			set(Summary_Data, row_number, names(Summary_Data), temp[the_inst, .SD, .SDcols=names(Summary_Data)]);
			row_number <- row_number + 1;
		}
	}
	Summary_Data <- na.omit(Summary_Data);
	if(nrow(Summary_Data)==0){ next; }

	for(metric in metrics_to_plot)
	{
		plot_name <- sprintf("%s%s/Hes_Cloud_%s_N_%s_struct_%s.png", graph_file_path, folder_name, metric, size, structure); # _networkprop_%s

		pl <- ggplot(data=Summary_Data, mapping=aes(x=perceived_vaccine_risk, y=get(metric), colour=factor(initial_vacc_proportion))) +
			#, group=factor(initial_vacc_proportion)
			geom_point(size=3) +
			geom_smooth(method="gam", size=2, aes(fill=get(metric))) +
			scale_colour_manual(values=c("blue", "red")) +
			scale_fill_discrete(labels=unique(Summary_Data$initial_vacc_proportion)) +
			xlab(expression(kappa)) + ylab(axis_label(metric)) +
			theme(
					panel.background=element_blank(),
					axis.text=element_text(size=45),
					axis.title=element_text(size=55),
					legend.text=element_text(size=40),
					legend.title=element_text(size=50),
					legend.key = element_rect(size = 2),
					legend.key.width=unit(2, 'cm'),
					legend.justification='center',
					legend.direction="vertical",
					# legend.position="right",
					# legend.position=c(-1, max(Summary_Data[[metric]])),
					legend.position=c(0.07,0.7),
					legend.title.align=0.5,
					legend.spacing.y = unit(3, "mm"),
			        panel.border = element_rect(colour = "black", fill=NA),
			        legend.box.background = element_rect(colour = "black")
		) +
		guides(colour=guide_legend(title=expression(alpha), size=30));

		ggsave(
			plot_name,
			plot = pl, width=25, height=5, dpi=50
		);
		dev.off()

		system(sprintf("convert %s -trim %s", plot_name, plot_name));
	}

	stop(Sys.time() - start_time);
}

end <- Sys.time();

print(end-start);
