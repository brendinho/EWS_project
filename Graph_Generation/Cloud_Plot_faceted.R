# # Brendon Phillips
# # PhD candidate
# # Bauch lab cumputational epidemiology group
# # Department of Mathematics
# # University of Waterloo
#
# rm(list=ls());
#
# start_time <- Sys.time();
#
# source("/home/b2philli/Dropbox/Processing/Parameter_Values.R");
#
# folder_name <- paste("Cloud_Plots", sep="");
# dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);
#
# Params <- fread(parameter_file, colClasses="character"); Params[, V1:=NULL];
# Parameters <- unique(Params[phys_top=="random", .SD, .SDcols=c("size", "phys_top")]);
#
# metrics_to_plot <- c('phys_S_mean', 'phys_V_mean');
#
# Parameters <- Parameters[size==562500];
#
# for(line_index in 1:nrow(Parameters))
# {
# 	this_tuple <- Parameters[line_index];
#
# 	filenames <- Sys.glob(paste(
# 		CSV_file_path,
# 		"Summary",
# 		"_N_", this_tuple$size,
# 		"_dur_", "*",
# 		"_beta_", "*",
# 		"_vaccprop_", "*",
# 		"_risk_", "*",
# 		"_inf_", "*",
# 		"_imp_", "*",
# 		"_rep_", "*",
# 		"_ptop_", this_tuple$phys_top,
# 		"_pdeg_", "*",
# 		"_stop_", this_tuple$phys_top,
# 		"_sdeg_", "*",
# 		"_norm_", "*",
# 		"_switch_", "*",
# 		".csv",
# 		sep=""
# 	));
#
# 	if(length(filenames) == 0){ next; }
#
# 	col_names <- c("perceived_vaccine_risk", "initial_vacc_proportion", metrics_to_plot);
# 	Summary_Data <- data.table(matrix(NaN, nrow=length(filenames)*100, ncol=length(col_names))); # assuming 100 instances per parameter set
# 	names(Summary_Data) <- col_names;
#
# 	row_number <- 1;
# 	for(file in filenames)
# 	{
# 		temp <- data.table(read.csv(file))[!is.na(instance), .SD, .SDcols=names(Summary_Data)];
# 		for(the_inst in 1:nrow(temp))
# 		{
# 			set(Summary_Data, row_number, names(Summary_Data), as.list(temp[the_inst]));
# 			row_number <- row_number + 1;
# 		}
# 	}
# 	Summary_Data <- na.omit(Summary_Data[perceived_vaccine_risk<=1 & perceived_vaccine_risk>=-1]);
# 	if(nrow(Summary_Data)==0){ next; }
# 	Summary_Data <- melt(Summary_Data, id.vars=c("perceived_vaccine_risk", "initial_vacc_proportion"))
#
# 	write.csv(Summary_Data, "temp_cloud_plot_fata.csv");

	Summary_Data <- data.table(read.csv("temp_cloud_plot_fata.csv")); Summary_Data[, X:=NULL];

	initial_vacc_prop_count <- length(unique(Summary_Data$initial_vacc_proportion))
	pal <- colorRampPalette(c("blue", "green", "red"))
	palette <- pal(initial_vacc_prop_count);

	plot_name <- sprintf(
		"%s%s/Cloud_%s_N_%s_struct_%s.png",
		graph_file_path,
		folder_name,
		paste(metrics_to_plot, sep="", collapse="_"),
		this_tuple$size,
		this_tuple$phys_top
	); # _networkprop_%

	levels(Summary_Data$variable) <- c(expression("<S"*">"), expression("<V"['p']*">"));

	# the_label <- variable=c(phys_S_mean="<S>", phys_V_mean=bquote("<V"['p']*">"))

	pl <- ggplot(data=Summary_Data, mapping=aes(x=perceived_vaccine_risk, y=value, colour=factor(initial_vacc_proportion))) +
		#, group=factor(initial_vacc_proportion)
		geom_point(size=3) +
		scale_colour_manual(values=palette) +
		xlab(expression(kappa)) +
		# scale_fill_discrete(labels=unique(Summary_Data$initial_vacc_proportion)) +
		theme(
			panel.background=element_blank(),
			axis.title.y=element_blank(),
			axis.text=element_text(size=45),
			axis.title=element_text(size=55),
			legend.text=element_text(size=40),
			legend.title=element_text(size=50),
			# legend.key = element_rect(size = 2),
			strip.text=element_text(size=40),
			panel.spacing = unit(2, "lines"),
			strip.text.y.right=element_text(angle=0),
			legend.key.width=unit(2, 'cm'),
			legend.justification='center',
			legend.direction="vertical",
			legend.position="right",
			legend.title.align=0.5,
			legend.spacing.y = unit(3, "mm"),
			panel.border = element_rect(colour = "black", fill=NA),
			legend.box.background = element_rect(colour = "black"),
			axis.line=element_line(),
		) +
		facet_grid(variable~., scales="free_y", labeller=label_parsed) +
		annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
		guides(colour=guide_legend(title=expression(alpha), override.aes = list(size=15)));

		# geom_smooth(method="gam", size=2) +

	ggsave(
		plot_name,
		plot = pl, width=25, height=9, dpi=50
	);
	dev.off()

# 	system(sprintf("convert %s -trim %s", plot_name, plot_name));
#
# }
#
# end <- Sys.time();
#
# print(end-start);
