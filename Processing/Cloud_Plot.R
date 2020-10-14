# Brendon Phillips
# PhD candidate
# Bauch lab cumputational epidemiology group
# Department of Mathematics
# University of Waterloo

rm(list=ls());

start_time <- Sys.time();

source("/home/b2philli/Dropbox/Processing/Parameter_Values.R");

folder_name <- paste("Cloud_Plots", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character"); Params[, V1:=NULL];
Parameters <- unique(Params[phys_top=="random", .SD, .SDcols=c("size", "phys_top")]);

metrics_to_plot <- c('phys_S_mean', 'phys_I_mean', 'phys_R_mean', 'phys_V_mean', 'soc_N_mean', 'soc_V_mean');

Parameters <- Parameters[size==10000];

for(line_index in 1:nrow(Parameters))
{
	this_tuple <- Parameters[line_index];

	filenames <- Sys.glob(paste(
		CSV_file_path,
		"Summary",
		"_N_", this_tuple$size,
		"_dur_", "*",
		"_beta_", "*",
		"_vaccprop_", "*",
		"_risk_", "*",
		"_inf_", "*",
		"_imp_", "*",
		"_rep_", "*",
		"_ptop_", this_tuple$phys_top,
		"_pdeg_", "*",
		"_stop_", this_tuple$phys_top,
		"_sdeg_", "*",
		"_norm_", "*",
		"_switch_", "*",
		".csv",
		sep=""
	));

	if(length(filenames) == 0){ next; }

	col_names <- c(fixed_values, metrics_to_plot);
	Summary_Data <- data.table(matrix(NaN, nrow=length(filenames)*100, ncol=length(col_names))); # assuming 100 instances per parameter set
	names(Summary_Data) <- col_names;
	Summary_Data[, physical_topology:=as.character(physical_topology)];
	Summary_Data[, social_topology:=as.character(social_topology)];

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

	initial_vacc_prop_count <- length(unique(Summary_Data$initial_vacc_proportion))
	# palette <- rainbow(initial_vacc_prop_count);
	pal <- colorRampPalette(c("blue", "green", "red"))
	palette <- pal(initial_vacc_prop_count);

	for(metric in metrics_to_plot)
	{
		plot_name <- sprintf("%s%s/Cloud_%s_N_%s_struct_%s.png", graph_file_path, folder_name, metric, this_tuple$size, this_tuple$phys_top); # _networkprop_%s

		pl <- ggplot(data=Summary_Data, mapping=aes(x=perceived_vaccine_risk, y=get(metric), colour=factor(initial_vacc_proportion))) +
			#, group=factor(initial_vacc_proportion)
			geom_point(size=3) +
			# geom_smooth(method="gam", size=2) +
			scale_colour_manual(values=palette) +
			scale_fill_discrete(labels=unique(Summary_Data$initial_vacc_proportion)) +
			xlab(expression(kappa)) + ylab(axis_label(metric)) +
			theme(
					panel.background=element_blank(),
					axis.text=element_text(size=45),
					axis.title=element_text(size=55),
					legend.text=element_text(size=40),
					legend.title=element_text(size=50),
					# legend.key = element_rect(size = 2),
					legend.key.width=unit(2, 'cm'),
					legend.justification='center',
					legend.direction="horizontal",
					legend.position="bottom",
					# legend.position=c(0.07,0.7),
					# legend.title.align=0.5,
					legend.spacing.y = unit(3, "mm"),
			        panel.border = element_rect(colour = "black", fill=NA),
			        legend.box.background = element_rect(colour = "black")
		) +
		guides(colour=guide_legend(title=expression(alpha), nrow=1, override.aes = list(size=15)));

		ggsave(
			plot_name,
			plot = pl, width=25, height=6, dpi=50
		);
		dev.off()

		system(sprintf("convert %s -trim %s", plot_name, plot_name));
	}
}

end <- Sys.time();

print(end-start);
