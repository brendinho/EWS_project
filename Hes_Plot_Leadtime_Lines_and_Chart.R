# Brendon Phillips
# PhD candidate
# Department of Applied Mathematics
# University of Waterloo

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Lead_Time", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character")
Params[, c("V1", "perceived_vaccine_risk", "social_norm"):=NULL];
Params <- Params[network_structure=="smallworld" & num(initial_vacc_proportion)==0.25];
Params <- unique(Params);

start_time <- Sys.time();

text_size <- 35;
title_size <- 35;

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

	writeLines("\n");
	writeLines(paste("\timportation rate: ", importation_rate, sep=""));
	writeLines(paste("\tnetwork size: ", the_size, sep=""));
	writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
	writeLines(paste("\tstructure: ", structure, sep=""));
	writeLines(paste("\tdegree: ", degree, sep=""));
	writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));
	writeLines(sprintf("\tproportion tested: %f", proportion));

    filenames <- Sys.glob(paste(CSV_file_path, "Hes_N_", the_size, "_struct_", structure, "_deg_", degree,  "_risk_*", "_infec_", infection_prob,
                                "_norm_*",  "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_", proportion,
                                "_summary.csv", sep=""));

    if(length(filenames) < 10) next;

	# if the tested proportion isn't 1, then we want to get ht real estimates of the transitions from the prop_1 trials
	if(proportion != 1)
	{
		real_filenames <- Sys.glob(paste(CSV_file_path, "Hes_N_", the_size, "_struct_", structure, "_deg_", degree,  "_risk_*", "_infec_", infection_prob,
									"_norm_*",  "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_1_summary.csv", sep=""));


		if(length(real_filenames) < 10){ print("Not enough real shit to do the fakes. moving on..."); next; }

		Real_Prop <- data.table();
		for(file in real_filenames)
		{
			temp <- data.table(read.csv(file));
			Real_Prop <- rbind(Real_Prop, temp[is.na(instance)], fill=TRUE);
		}
		Real_Prop[, ("X"):=NULL];

		Real_Transitions <- data.table(matrix(ncol=3));
		names(Real_Transitions) <- c("norm", "phys_trans", "soc_trans");
		for(norm_index in seq_along(Real_Prop[, unique(social_norm)]))
		{
			soc_norm <- Real_Prop[, sort(unique(social_norm))][norm_index];
			transes <- get_transitions(Real_Prop[social_norm==soc_norm][order(perceived_vaccine_risk)], initial_vacc_prop);
			Real_Transitions <- rbind( Real_Transitions, list(soc_norm, transes$phys_intersec_risk, transes$soc_intersec_risk) );
		}
		Real_Transitions <- Real_Transitions[rowSums(is.na(Real_Transitions)) != ncol(Real_Transitions)][order(norm)];

		rm(real_filenames, Real_Prop);
	}

	Prop_vs_Risk <- data.table();
	for(file in filenames)
	{
		temp <- data.table(read.csv(file));
		Prop_vs_Risk <- rbind(Prop_vs_Risk, temp[is.na(instance)], fill=TRUE);
	}
	Prop_vs_Risk[, ("X"):=NULL];

	if(proportion != 1)
	{  # only take the norm values for which we can compare with the real estimates
		Prop_vs_Risk <- Prop_vs_Risk[social_norm %in% intersect(Prop_vs_Risk$social_norm, Real_Transitions$norm)];
		# make sure that they now have values only for the same norms, so there's no checking later on
		Real_Transitions <- Real_Transitions[norm %in% intersect(Prop_vs_Risk$social_norm, Real_Transitions$norm)];
	}
	write.csv(Prop_vs_Risk, "temp_risk.csv");

	# Prop_vs_Risk <- data.table(read.csv("temp_risk.csv")); Prop_vs_Risk[, X:=NULL];

	metrics <- setdiff(names(Prop_vs_Risk), c(hes_fixed_values, "instance"));
	metrics <- filter_these(metrics, no=c("_sd", "ratio_", "soc", "phys"));

    col_names <- c("norm", "phys_trans", "soc_trans", names.change(metrics));

	the_names_to_group_together <- c("watts", "mutual", "diameter", "triads", "join.count", "opinion.change", "prob.sick",
					"avg.conn.comp", "min.conn.comp", "max.conn.comp", "number.conn", "avg.chamber.size", "min.chamber.size", "max.chamber.size", "number.chambers"
	); # "phys", "soc",#

    for(CHANGE_TEST in c("snh")) # , "lanzante", "pettitt", "buishand_r")) # "binseg")) #
    {
		file_stem <- paste("_N_", the_size, "_struct_", structure, "_deg_", degree, "_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop,
							"_switch_", random_opinion_switch, "_prop_", proportion, sep=""
						);

		writeLines(sprintf("\ttest: %s", CHANGE_TEST));
        Change_Predictions <- data.table(matrix(ncol=length(col_names)));
        names(Change_Predictions) <- col_names;
        for(norm_index in seq_along(Prop_vs_Risk[, unique(social_norm)]))
        {
            soc_norm <- Prop_vs_Risk[, sort(unique(social_norm))][norm_index];
            Prop_Table_Here <- Prop_vs_Risk[social_norm==soc_norm][order(perceived_vaccine_risk)];
            transitions <- get_transitions(Prop_Table_Here, initial_vacc_prop);
            change_points <- c();
            for(metric in metrics) change_points[[metric]] <- change_point_test(Prop_Table_Here, metric, CHANGE_TEST)$kappa;
            Change_Predictions <- rbind(
                Change_Predictions,
                as.list(c(
                    soc_norm,
					transitions$phys_intersec_risk,
					transitions$soc_intersec_risk,
                    data.table(change_points)[, change_points]
                ))
            );
        }
        Change_Predictions <- Change_Predictions[rowSums(is.na(Change_Predictions)) != ncol(Change_Predictions)][order(norm)];
		if(proportion !=1) # if testing only a subset of the nodes, we're getting the earlier predicted transition for either dynamic
		{	# reorder the norms for the real_transitions vector so that the two tables have the same order of norms
			Real_Transitions <- Real_Transitions[order(norm)];
			# get the first predicted physical and social transitions
			Change_Predictions[, phys_trans:=pmin(Real_Transitions$phys_trans, Change_Predictions$phys_trans)];
			Change_Predictions[,  soc_trans:=pmin(Real_Transitions$soc_trans,  Change_Predictions$soc_trans )];
		}
        write.csv(Change_Predictions, "temp_change_predictions.csv");

        # Change_Predictions <- data.table(read.csv("temp_change_predictions.csv")); Change_Predictions[, X:=NULL];

        All_Lead_Times <- Change_Predictions[, min(soc_trans, phys_trans)] - Change_Predictions[, .SD, .SDcols=names.change(metrics)];
        All_Lead_Times <- cbind(Change_Predictions$norm, All_Lead_Times); names(All_Lead_Times)[1] <- "norm";

        Max_or_Min <- data.table(matrix(nrow=nrow(All_Lead_Times), ncol=ncol(All_Lead_Times)-1, "neither"));
        Max_or_Min <- cbind(Change_Predictions$norm, Max_or_Min);
        names(Max_or_Min) <- names(All_Lead_Times);

		for(num_norm in unique(All_Lead_Times$norm)) { for(name in setdiff(names(All_Lead_Times), "norm")) {
            if(is.na( All_Lead_Times[norm==num_norm, get(name)] )){ Max_or_Min[which(norm==num_norm), (name):="none"]; }
            else if(All_Lead_Times[norm==num_norm, get(name)]<0)
			{
				Max_or_Min[which(norm==num_norm), (name):="fail"];
				All_Lead_Times[which(norm==num_norm), (name):=NA];
			}
        }}

		for(i in 1:nrow(All_Lead_Times))
        {
            rows_with_max <- which(All_Lead_Times[i, ] == maxx(All_Lead_Times[i, -c("norm")]));
            for(here in rows_with_max)
            {
                if(Max_or_Min[i, ..here][[1]] == "neither"){ Max_or_Min[i, here] <- "max"; }
            }
            rows_with_min <- which(All_Lead_Times[i, ] == minn(All_Lead_Times[i, -c("norm")]));
            for(here in rows_with_min)
            {
                if(Max_or_Min[i, ..here][[1]] == "max"){ Max_or_Min[i, here] <- "both"; }
                else if(Max_or_Min[i, ..here][[1]] == "neither"){ Max_or_Min[i, here] <- "min"; }
            }
        }

		Max_Min_Grid <- Max_or_Min %>% gather(col_name, value, -norm) %>% ggplot(aes(factor(norm), col_name, fill=value)) +
			geom_tile(colour = 'black') +
			labs(x=expression(sigma), y=test_name(CHANGE_TEST, "new")) +
			scale_x_discrete(breaks=seq(min(Max_or_Min$norm), max(Max_or_Min$norm), by=.25)) +
			scale_y_discrete(labels=as.vector(sapply(setdiff(names(Max_or_Min), "norm"), axis_label))) +
			scale_fill_manual(values = c('max'='green', 'neither'='grey', 'min'='red', 'both'='yellow', 'none'='white', 'fail'='black')) +
			theme(
				legend.position="none",
				axis.text.x=element_text(size=text_size+5),
				axis.text.y=element_text(size=text_size+5),
				axis.title=element_text(size=title_size+15)
			);

		ggsave(
			paste(graph_file_path, folder_name, "/", "Hes_LeadTime_grid_", CHANGE_TEST, file_stem, "_total.png", sep=''),
			plot = Max_Min_Grid, width=ncol(Max_or_Min)/2, height=ncol(Max_or_Min)/1.5, limitsize=FALSE, dpi=50
		);
		dev.off();

        Joint <- data.table(EWS=character(), Maxima=numeric(), Minima=numeric());
        for(row_name in setdiff(names(Max_or_Min), "norm"))
        {
            Joint <- rbind(Joint,
                           list(
                               row_name,
                               length(which(Max_or_Min[, get(row_name)] == "max"))/nrow(Max_or_Min),
                               length(which(Max_or_Min[, get(row_name)] %in% c("min", "fail", "none")))/nrow(Max_or_Min)
                           )
            );
        }

		ordered_rows <- c();
		for(i in the_names_to_group_together)
        {
            ordered_rows <- c(ordered_rows, which(grepl(i, Joint$EWS))); # group all the like metrics together
        }
        Joint <- Joint[ordered_rows]; # put the rows in order
        Joint[, Position:=1:nrow(Joint)]; # put in a row that preserves this order, so we can use the reorder() command

		Joint <- melt(Joint, id.vars=c("EWS", "Position"), measure.vars=c("Maxima", "Minima"));

		pl_max <- ggplot(Joint[variable=="Maxima"], aes(x=reorder(EWS, Position), y=value)) +
			geom_bar(stat='identity', fill='green') +
			labs(title="maximal lead times") +
			scale_x_discrete(labels=as.vector(sapply(Joint$EWS, axis_label))) +
			geom_text(aes(label=sprintf("%i%%", round(100*abs(value), 0))), size=10, colour="black") +
			coord_flip() +
			theme(
				legend.position='none',
				plot.title=element_text(size=title_size),
				axis.text.y=element_text(size=text_size-5, hjust=0.5),
				axis.text.x=element_text(size=text_size-5),
				axis.ticks.y=element_blank(),
				axis.line.y=element_blank(),
				axis.title=element_blank(),
				panel.grid.major=element_line(size = 0.50, linetype = 'solid', colour = 'grey'),
				panel.grid.minor=element_line(size = 0.25, linetype = 'solid', colour = 'grey'),
			);

		pl_min <- ggplot(Joint[variable=="Minima"], aes(x=reorder(EWS, Position), y=value)) +
			geom_bar(stat='identity', fill='red') +
			labs(title="minimal, undefined, negative lead times") +
			scale_y_continuous('', trans='reverse') +
			geom_text(aes(label=sprintf("%i%%", round(100*abs(value), 0))), size=10, colour="black") +
			coord_flip() +
			theme(
				legend.position='none',
				plot.title=element_text(size=title_size),
				axis.text.x=element_text(size=text_size-5),
				axis.text.y=element_blank(),
				axis.ticks.y=element_blank(),
				axis.line.y=element_blank(),
				axis.title=element_blank(),
				panel.grid.major=element_line(size = 0.50, linetype = 'solid', colour = 'grey'),
				panel.grid.minor=element_line(size = 0.25, linetype = 'solid', colour = 'grey'),
			);

		pl_performance <- grid.arrange(
			pl_min,
			pl_max,
			widths=c(0.46, 0.54),
			ncol=2,
			bottom=expression("Proportion of social("*sigma*") norm values")
		);

		graph_name <- paste(graph_file_path, folder_name, "/", "Hes_LeadTime_bars_", CHANGE_TEST, file_stem, "_total.png", sep='');
		ggsave(
			graph_name,
			plot = pl_performance, width=25, height=length(unique(Joint$EWS))/1.8, limitsize=FALSE, dpi=50
		);

        warning_set <- c("conn_comp_size", "num_triads", "number", "mutual", "prob_sick", "join_counts", "echo_chamber_size", "watts_strogatz", "opinion_change", "net_diameter"); # , "outputs"

		double_row_legend <- c(); # c("conn_comp_size")

		Joint <- Joint[variable=="Maxima"];
		Joint[, variable:=NULL];

        for(warn in warning_set)
        {
			writeLines(sprintf("\t\twarnings: %s", warn));
            lets_find_the_maxes <- plot_these(warn);

            getPalette <- colorRampPalette(brewer.pal(8, "Dark2"));
            colours_here <- getPalette(length(lets_find_the_maxes));
            names(colours_here) <- names.change(lets_find_the_maxes);

            Plot_From_This_Table <- Joint[EWS %in% names.change(lets_find_the_maxes)];

            pl_best <- ggplot(
	                Plot_From_This_Table,
	                aes(x=EWS, y=value, fill=EWS) # , colour=EWS
	            ) +
                geom_bar(stat="identity") +
                geom_text(aes(label=sprintf("%i%%", round(100*value, 0))), vjust=-0.1, size=8, colour="black") +
                scale_y_continuous(limits=c(0, max(Plot_From_This_Table$value, na.rm=TRUE)+0.015)) + # eyJ1IjoiYjJwaGlsbGkiLCJhIjoiY2s1NjZ5bzJlMDI4dzNucGtkdjA2NjNydiJ9
                scale_fill_manual(values=colours_here, name="", breaks=Plot_From_This_Table$EWS, labels=as.vector(sapply(Plot_From_This_Table$EWS, axis_label))) +
                ylab("# Maxima (LLT)") +
                theme(
                    panel.background=element_blank(),
                    axis.line=element_line(colour = "black"),
                    axis.text=element_blank(),
                    axis.title=element_text(size=title_size),
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
					# legend.position="none",
					legend.position="right",
					legend.direction="vertical",
					legend.justification="center",
					legend.key.size = unit(2,"line"),
					legend.key.width=unit(2, "line"),
					legend.spacing.x = unit(0, 'cm'),
					legend.spacing.y = unit(0, 'cm'),
					legend.text=element_text(size=text_size),
					legend.title=element_blank(),
					legend.box.background = element_rect(colour = "black"),
					legend.box.margin=unit(c(0,1,0,0), "cm"),

                );
				# if(warn %in% double_row_legend){ pl_best <- pl_best + guides(fill=guide_legend(nrow=2,byrow=TRUE)); }
				# else { pl_best <- pl_best + guides(fill = guide_legend(nrow=1)); }

            pl_legend <- g_legend(pl_best);
            dev.off();

            Lead_Times_Here <- All_Lead_Times[, .SD, .SDcols=intersect(names.change(lets_find_the_maxes), names(All_Lead_Times))];
            Lead_Times_Here <- cbind(Change_Predictions$norm, Lead_Times_Here);
            names(Lead_Times_Here)[1] <- "norm";
            Lead_Times_Here <- melt(Lead_Times_Here, id="norm");

            pl_EWS <- ggplot(Lead_Times_Here, aes(norm, value)) +
                geom_line(aes(colour=variable), size=1) +
				scale_colour_manual(values=colours_here, labels=as.vector(sapply(Plot_From_This_Table$EWS, axis_label))) +
                ylim(c(min(Lead_Times_Here$value, na.rm=TRUE), max(Lead_Times_Here$value, na.rm=TRUE)+0.04)) +
                # labs(x=expression(sigma), y=bquote('min(K'['*']*')-'*.(axis_label(CHANGE_TEST))*'{'*.(axis_label(warn))*'}')) +
				labs(x=expression(sigma), y=bquote('K'['M']*'-'*.(axis_label(CHANGE_TEST))*'{'*.(axis_label(warn))*'}')) +
                theme(
                    panel.background=element_blank(),
                    axis.line=element_line(colour = "black"),
                    axis.text=element_text(size=text_size),
                    axis.title=element_text(size=title_size),
					axis.title.y=element_text(hjust=1.2),
					# legend.position="bottom",
					# legend.direction="horizontal",
					legend.position="right",
					legend.direction="vertical",
					legend.justification="center",
					legend.key.size = unit(2,"line"),
					legend.key.width=unit(2, "line"),
					legend.spacing.x = unit(0.3, 'cm'),
					legend.spacing.y = unit(0, 'cm'),
					legend.text=element_text(size=text_size),
					legend.title=element_blank(),
					legend.box.background = element_rect(colour = "black"),
					legend.box.margin=unit(c(0,1,0,0), "cm"),
                ) +
                geom_hline(yintercept=0, size=1, colour="black", linetype="dashed") +
				guides(fill = guide_legend(override.aes = list(shape = 5, size=9))) # +

			# pl_legend <- g_legend(pl_EWS);
			# dev.off();

			# pair_heights <- {if(warn %in% double_row_legend) c(10,3) else c(10,1)};
			# left_pl <- grid.arrange(
			#     pl_EWS + theme(legend.position="none"),
			#     pl_legend,
			#     nrow=2,
			#     heights=pair_heights
			# );
			# left_pl <- pl_EWS;
            # right_pl <- pl_best + theme(legend.position="none");
            # pl_truth <- plot_grid( left_pl, right_pl, rel_widths=c(7.5,2.5) );
			pl_truth <- plot_grid(
				pl_EWS+theme(legend.position="none"),
				pl_best+theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm")), # top, right, bottom, left),
				pl_legend,
				rel_widths=c(7,2,1),
				nrow=1
			) #
            pl_truth <- cowplot::ggdraw(pl_truth) + theme(plot.background=element_blank());

			graph_name <- paste(graph_file_path, folder_name, "/Hes_Leadtime_", warn, "_", CHANGE_TEST, "_N_", the_size, "_struct_", structure, "_deg_", degree, "_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_", proportion, ".png", sep="");

            ggsave(
                graph_name,
                plot = pl_truth, width=25, height=4, limitsize=FALSE, dpi=50
            );
            dev.off();
			# system(sprintf("convert %s -trim %s", graph_name, graph_name));
			# stop()

        }
    }
}

print(Sys.time() - start_time);

graphics.off()
