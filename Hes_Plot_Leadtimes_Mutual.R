# Brendon Phillips
# PhD candidate
# Department of Applied Mathematics
# University of Waterloo

# Summary Plot for the Topology
rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Lead_Time", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character")
Params[, c("V1", "perceived_vaccine_risk", "social_norm"):=NULL];
Params <- Params[network_structure=="smallworld" & num(initial_vacc_proportion)==0.25];
Params <- unique(Params);

start_time <- Sys.time();

text_size <- 45;
title_size <- 45;

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
    proportion <- convert_if_number(this_tuple$proportion_of_nodes);

    filenames <- Sys.glob(paste(CSV_file_path, "Hes_N_", the_size, "_struct_", structure, "_deg_", degree,  "_risk_*", "_infec_", infection_prob,
                                "_norm_*",  "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_", proportion,
                                "_summary.csv", sep=""));

    if(length(filenames) < 10) next;

	writeLines("\n");
    writeLines(paste("\timportation rate: ", importation_rate, sep=""));
    writeLines(paste("\tnetwork size: ", the_size, sep=""));
    writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
    writeLines(paste("\tstructure: ", structure, sep=""));
    writeLines(paste("\tdegree: ", degree, sep=""));
    writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));
    writeLines(sprintf("\tproportion tested: %f", proportion));

    Prop_vs_Risk <- data.table();
    for(file in filenames)
	{
		temp <- data.table(read.csv(file));
    	Prop_vs_Risk <- rbind(Prop_vs_Risk, temp[is.na(instance)], fill=TRUE);
	}
    Prop_vs_Risk[, ("X"):=NULL];
    write.csv(Prop_vs_Risk, "mutual_risk.csv");

	# Prop_vs_Risk <- data.table(read.csv("mutual_risk.csv")); Prop_vs_Risk[, X:=NULL];

    col_names <- c("test", "norm", "phys_trans", "soc_trans", "mutual");
	Leads <- data.table(matrix(ncol=length(col_names), nrow=0));

    for(CHANGE_TEST in c("lanzante", "pettitt", "buishand_r", "snh"))
    {
        writeLines(sprintf("\ttest: %s", CHANGE_TEST));

		file_stem <- paste("_N_", the_size, "_struct_", structure, "_deg_", degree, "_infec_", infection_prob,
							"_imp_", importation_rate, "_init_", initial_vacc_prop,
							"_switch_", random_opinion_switch, "_prop_", proportion, sep=""
						);

        names(Leads) <- col_names;
        for(norm_index in seq_along(Prop_vs_Risk[, unique(social_norm)]))
        {
            soc_norm <- Prop_vs_Risk[, sort(unique(social_norm))][norm_index];
            Prop_Table_Here <- Prop_vs_Risk[social_norm==soc_norm][order(perceived_vaccine_risk)];
            transitions <- get_transitions(Prop_Table_Here, initial_vacc_prop);
            change_point <- change_point_test(Prop_Table_Here, "mutual_info_mean", CHANGE_TEST)$kappa;
			Leads <- rbind(Leads, list(CHANGE_TEST, soc_norm, transitions$phys_intersec_risk, transitions$soc_intersec_risk, change_point));
        }
	}
    write.csv(Leads, "temp_mutual_predictions.csv");

	Leads <- data.table(read.csv("temp_mutual_predictions.csv")); Leads[, X:=NULL];

	Leads$the_lead <- Leads[, min(phys_trans, soc_trans)-mutual, by=c("norm", "test")]$V1

    pl_EWS <- ggplot(Leads, aes(x=norm, y=the_lead, colour=test)) +
        geom_line(size=1.5) +
        # scale_colour_manual(values=colours_here) +
        # labs(x=expression(sigma), y=bquote('min(K'['s']*',K'['p']*')-B{'*.(axis_label("mutual_info"))*'}')) +
		labs(x=expression(sigma), y=bquote('K'['M']*'-B{'*.(axis_label("mutual_info"))*'}')) +
        theme(
			# panel.background=element_blank(),
			# panel.grid.major=element_line(size = 0.50, linetype = 'solid', colour = 'grey'),
			# panel.grid.minor=element_line(size = 0.25, linetype = 'solid', colour = 'grey'),
            axis.line=element_line(colour = "black"),
            axis.text=element_text(size=text_size-5),
            axis.title=element_text(size=title_size-7),
			legend.text=element_text(size=text_size-10),
			legend.title=element_blank(),
			legend.key.width=unit(3, 'cm'),
			legend.spacing.x=unit(1, 'cm'),
			# legend.position="bottom",
			# legend.direction="horizontal",
			legend.position="right",
			legend.direction="vertical",
			legend.justification="center"
        ) +
		scale_colour_discrete(breaks=unique(Leads$test), labels=as.vector(sapply(unique(Leads$test), test_name))) # +
		# guides(colour = guide_legend(nrow=1));

	graph_name <- paste(graph_file_path, folder_name, "/Hes_Mutuals_N_", the_size, "_struct_", structure, "_deg_", degree, "_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop, "_switch_", random_opinion_switch, "_prop_", proportion, "_all_tests.png", sep="");

    ggsave(graph_name, plot = pl_EWS, width=25, height=4, units="in", limitsize=FALSE);
    dev.off();

}
