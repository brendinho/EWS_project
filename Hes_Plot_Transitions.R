# Brendon Phillips
# PhD candidate
# Department of Applied Mathematics
# University of Waterloo

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Location", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character")
Params[, c("V1", "perceived_vaccine_risk", "social_norm"):=NULL];
Params <- Params[network_structure=="smallworld" & num(initial_vacc_proportion)==0.25];
Params <- unique(Params[, proportion_of_nodes:=NULL]);

start_time <- Sys.time();

trend_width <- 25;
trend_height <- 7;
text_size <- 40;
title_size <- 40;

for(line_index in 1:nrow(Params))
{
    this_tuple <- Params[line_index];

    the_size <- convert_if_number(this_tuple$N);
    structure <- this_tuple$network_structure;
    degree <- convert_if_number(this_tuple$mean_degree);
    infection_prob <- convert_if_number(this_tuple$infec_prob);
    importation_rate <- "2.5e-05"; # convert_if_number(this_tuple$importation);
    initial_vacc_prop <- convert_if_number(this_tuple$initial_vacc_proportion);
    random_opinion_switch <- "1e-04"; convert_if_number(this_tuple$random_opinion_switch);

	writeLines("\n");
	writeLines(paste("\timportation rate: ", importation_rate, sep=""));
	writeLines(paste("\tnetwork size: ", the_size, sep=""));
	writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
	writeLines(paste("\tstructure: ", structure, sep=""));
	writeLines(paste("\tdegree: ", degree, sep=""));
	writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));

	Real_Transitions <- data.table(matrix(ncol=4));
	names(Real_Transitions) <- c("norm", "proportion", "phys_trans", "soc_trans");
	for(test_proportion in c("0.6", "0.8", "1"))
	{
		real_filenames <- Sys.glob(paste(CSV_file_path, "Hes_N_", the_size, "_struct_", structure, "_deg_", degree,  "_risk_*",
				"_infec_", infection_prob, "_norm_*",  "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_",
				random_opinion_switch, "_prop_", test_proportion, "_summary.csv", sep=""));

		if(length(real_filenames) < 10){ print("Not enough real shit to do the fakes. moving on..."); next; }

		Real_Prop <- data.table();
		for(file in real_filenames)
		{
			temp <- data.table(read.csv(file));
			Real_Prop <- rbind(Real_Prop, temp[is.na(instance)], fill=TRUE);
		}
		Real_Prop[, ("X"):=NULL];

		for(norm_index in seq_along(Real_Prop[, unique(social_norm)]))
		{
			soc_norm <- Real_Prop[, sort(unique(social_norm))][norm_index];
			transes <- get_transitions(Real_Prop[social_norm==soc_norm][order(perceived_vaccine_risk)], initial_vacc_prop);
			Real_Transitions <- rbind( Real_Transitions, list(soc_norm, test_proportion, transes$phys_intersec_risk, transes$soc_intersec_risk) );
		}
		Real_Transitions <- Real_Transitions[rowSums(is.na(Real_Transitions)) != ncol(Real_Transitions)][norm<2.5];
	}
	write.csv(Real_Transitions, "temp_transitions_plot.csv");

	# Real_Transitions <- data.table(read.csv("temp_transitions_plot.csv")); Real_Transitions[, X:=NULL];

	Real_Transitions <-  melt(Real_Transitions, id.vars=c("norm", "proportion"), measure.vars=c('phys_trans', 'soc_trans'));

	# the_labels <- c("Social", "Physical");
	the_labels <- c("Social Transition", "Physical Transition");
	names(the_labels) <- c("soc_trans", "phys_trans")

	pl_Trans <- ggplot(Real_Transitions, aes(x=norm, y=value, colour=proportion)) +
		geom_line(size=1.5) +
		labs(x=expression(sigma), y=expression('K'['*']^beta), colour=expression(beta)) +
		facet_wrap(variable~., dir='v', labeller=labeller(variable=the_labels), scales='free') + # strip.position="right",
		# guides(colour = guide_legend(nrow=1)) +
		theme(
			axis.line=element_line(colour = "black"),
			axis.text=element_text(size = text_size-5),
			axis.title=element_text(size = title_size),
			legend.text=element_text(size = text_size-5),
			legend.title=element_text(size = title_size),
			legend.position="right",
			legend.key.width=unit(1, 'cm'),
			legend.spacing.x=unit(1, 'cm'),
			legend.direction="vertical",
			legend.justification="center",
			panel.background=element_rect(fill='grey95'),
			legend.title.align=0.5,
			legend.box.background = element_rect(colour = "black"),
			# legend.box.margin=unit(c(0,1,0,0), "cm"),
			strip.text.x = element_text(size = text_size-5),
		);

	graph_name <- paste(graph_file_path, folder_name, "/Hes_Transitions_N_", the_size, "_struct_", structure, "_deg_", degree, "_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop, "_switch_", random_opinion_switch, "_all_props.png", sep="");

	ggsave(graph_name, plot=pl_Trans, width=trend_width, height=trend_height, limitsize=FALSE, dpi=50);
	dev.off();
}

print(Sys.time() - start_time);

graphics.off()
