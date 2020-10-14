# Summary Plot for the Topology
rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Location", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character");
Params <- Params[num(initial_vacc_proportion)==0.25 & network_structure=="smallworld"];
Params[, c("V1","perceived_vaccine_risk", "social_norm"):=NULL];
Params <- unique(Params);

start_time <- Sys.time();

text_width <- 40;

for(line_index in 1:nrow(Params))
{
	this_tuple <- Params[line_index];

	initial_vacc_prop <- convert_if_number(this_tuple$initial_vacc_prop);
    the_size <- convert_if_number(this_tuple$N);
    infection_prob <- convert_if_number(this_tuple$infec_prob);
    importation_rate <- "2.5e-05"; # convert_if_number(this_tuple$importation);
    structure <- convert_if_number(this_tuple$network_structure);
    degree <- convert_if_number(this_tuple$mean_degree);
    random_opinion_switch <- "1e-04"; # convert_if_number(this_tuple$random_opinion_switch);
    proportion <- convert_if_number(this_tuple$proportion_of_nodes);

    filenames <- Sys.glob(paste(CSV_file_path, "Hes_N_", the_size, "_struct_", structure, "_deg_", degree,  "_risk_*", "_infec_", infection_prob,
                                "_norm_*",  "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_", proportion,
                                "_summary.csv", sep=""));

	if(length(filenames) < 10){ next; }

	Table <- data.table();
	for(file in filenames)
	{
		temp <- data.table(read.csv(file));
		temp <- unique(temp[
					,
					lapply(.SD, num),
					.SDcols=c("perceived_vaccine_risk", "phys_V_mean", "phys_R_mean", "initial_vacc_proportion", "soc_V_mean", "soc_N_mean", "soc_H_mean", "instance", "social_norm")
				]
				[order(perceived_vaccine_risk)]);

		Table <- rbind(Table, temp[is.na(instance)], fill=TRUE);
	}
	write.csv(Table, "temp_intertransition_table.csv")

	# Table <- data.table(read.csv("temp_intertransition_table.csv"));
	# Table[, c("X", "instance"):=NULL];

	Locations <- data.table(matrix(ncol=3));
	names(Locations) <- c("norm", "phys_trans", "soc_trans");

	for(soc_norm in Table[, num(unique(social_norm))])
	{
		transitions <- get_transitions(Table[social_norm==soc_norm], initial_vacc_prop);

		Locations <- rbind(
			Locations,
			list(soc_norm, transitions[, phys_intersec_risk], transitions[, soc_intersec_risk])
		);
	}
	Locations <- Locations[order(norm)];

	if(num(proportion)==1)
	{
		# Locations <- Locations[num(norm)<=2.5]
		min_x_inset <- -0.13;
		max_x_inset <- 1.35;
		min_y_inset <- 0.0015;
		max_y_inset <- 0.008;
	} else if(num(proportion) == 0.8)
	{
		min_x_inset <- -0.13;
		max_x_inset <- 2;
		min_y_inset <- 0.005;
		max_y_inset <- 0.05;
	} else if(num(proportion) == 0.6)
	{
		min_x_inset <- 1.5;
		max_x_inset <- 2.5;
		min_y_inset <- -0.001;
		max_y_inset <- 0.0015;
	}

	# keep these lines in order!
	Spline_Table <- data.table(data.frame(spline(Locations[, .SD, .SDcols=c("norm", "phys_trans")], n=nrow(Locations)*5)));
	Spline_Table <- cbind(Spline_Table, data.table(data.frame(spline(Locations[, .SD, .SDcols=c("norm", "soc_trans")], n=nrow(Locations)*5)))[, y]);
	Spline_Table <- cbind(Spline_Table, Spline_Table[, y-V2]);
	names(Spline_Table) <- c("norm", "phys_trans", "soc_trans", "differences");

	pl_trans_vals <- ggplot(Spline_Table, aes(x=norm)) +
		geom_point(aes(y=phys_trans, colour=expression('K'['p'])), colour='red', size=1.5) +
		geom_line(data=Spline_Table, aes(norm, phys_trans), colour='red', size=0.75) +
		geom_point(aes(y=soc_trans, colour=expression('K'['s'])), colour='blue', size=1.5) +
		geom_line(data=Spline_Table, aes(norm, soc_trans), colour='blue', size=0.75) +
		geom_ribbon(Spline_Table, mapping=aes(ymin=phys_trans, ymax=soc_trans), fill='grey', alpha=0.5) +
		xlab(expression(sigma)) + ylab(expression(kappa)) +
		theme(
			panel.background=element_blank(),
			axis.line=element_line(colour = "black"),
			axis.text=element_text(size=text_width-5),
			axis.title=element_text(size=text_width+5)
		);

	trans_vals_plot_name <- paste(graph_file_path, folder_name, "/Hes_TransVals_N_", the_size, "_struct_", structure, "_deg_", degree,
		"_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop, "_switch_", random_opinion_switch, "_prop_", proportion, ".png", sep="");

	ggsave(
		trans_vals_plot_name,
		plot = pl_trans_vals, width=22, height=5, limitsize=FALSE, dpi=50
	);
	dev.off()

	pl_intertrans_dist <- ggplot(Spline_Table, aes(x=norm)) +
		# labs(x=expression(sigma), y=expression('K'['*']^beta), colour=expression(beta)) +
		xlab(expression(sigma)) + ylab(bquote('K'['p']^.(proportion)*'-'*'K'['s']^.(proportion))) +
		theme(
			panel.background=element_blank(),
			axis.line=element_line(colour = "black"),
			axis.text=element_text(size=text_width),
			axis.title=element_text(size=text_width+10)
		) +
		geom_point(aes(y=differences), colour='purple', size=2) +
		geom_line(data=Spline_Table, aes(norm, differences), colour='purple', size=1.5) +
		geom_hline(yintercept=0, size=1, colour="black", linetype="dashed");

	intertrans_dist_plot_name <- paste(graph_file_path, folder_name, "/Hes_Intertrans_Plain_N_", the_size, "_struct_", structure, "_deg_", degree,
		"_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop, "_switch_", random_opinion_switch, "_prop_", proportion, ".png", sep="");

	ggsave(
		intertrans_dist_plot_name,
		plot = pl_intertrans_dist, width=22, height=5, limitsize=FALSE, dpi=50
	);
	dev.off()

	# Locations[, ("differences"):=phys_trans-soc_trans];

	total_pl <- pl_intertrans_dist +
		ylim(c(min(Spline_Table$differences)*1.05, 1.05*max(Spline_Table$differences))) +
		annotation_custom(
			ggplotGrob(pl_trans_vals),
			xmin=min_x_inset, xmax=max_x_inset,
			ymin=min_y_inset, ymax=max_y_inset
		)

	total_plot_file_name <- paste(graph_file_path, folder_name, "/Hes_IntertransitionDist_N_", the_size, "_struct_", structure, "_deg_", degree,
		"_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop, "_switch_", random_opinion_switch, "_prop_", proportion, ".png", sep="");

	ggsave(
		total_plot_file_name,
		plot = total_pl,
		width=22, height=5, limitsize=FALSE, dpi=50
	); # height=4
	dev.off()

	# stop()
}

graphics.off();
print(Sys.time()-start_time);
