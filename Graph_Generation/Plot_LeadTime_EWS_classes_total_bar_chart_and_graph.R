# Dear God, why have I retooled this code so many fucking times???

# Summary Plot for the Topology
rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Parameter_Values.R");

folder_name <- paste("Lead_Time", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character")
Params[, c("V1", "risk", "norm"):=NULL];

start_time <- Sys.time();

main_or_supplement <- "main";

if(main_or_supplement == "main")
{
  image_width <- 25;
  image_height <- 5.5;
  text_size <- 40;
  title_size <- 40;
} else if(main_or_supplement == "supplement")
{
  image_width <- 10;
  image_height <- 10;
  text_size <- 30;
  title_size <- 40;
}

for(the_pair in c( list(c(10000, "new")) )) # , list(c(10000, "new")),, list(c(10000, "old"))
{
    the_size <- num(the_pair[[1]])
    old_or_new <- the_pair[[2]]

	Parameters <- unique(Params[phys_top=="random" & num(init_prop)==0.05 & num(beta)==1 & num(size)==the_size]);

	if(num(the_size)==10000){ Parameters <- Parameters[num(import)==0.00025]; }
	else if(num(the_size)==40000){ Parameters <- Parameters[num(infec)==0.2]; }
	else if(num(the_size)==562500){ Parameters <- Parameters[]; }

	Parameters <- unique(Parameters);

	for(line_index in 1:nrow(Parameters))
	{
		this_tuple <- Parameters[line_index];

		initial_vacc_prop <- this_tuple[, init_prop]; beta <- this_tuple[, beta]; duration <- this_tuple[, duration];
		infection_prob <- this_tuple[, infec]; importation_rate <- this_tuple[, import]; replenishment_rate <-  this_tuple[, birth_death];
		random_opinion_switch <- this_tuple[, random_switch]; physical_topology <- this_tuple[, phys_top]; physical_degree <- this_tuple[, phys_deg];
		social_topology <- this_tuple[, soc_top]; social_degree <- this_tuple[, soc_deg];

		filenames <- Sys.glob(paste(CSV_file_path, "Summary", "_N_", the_size, "_dur_", duration, "_beta_", beta, "_vaccprop_", initial_vacc_prop, "_risk_", "*", "_inf_",
		                            infection_prob, "_imp_", importation_rate, "_rep_", replenishment_rate, "_ptop_", physical_topology, "_pdeg_", physical_degree, "_stop_",
		                            social_topology, "_sdeg_", social_degree, "_norm_", "*", "_switch_", random_opinion_switch, ".csv", sep=""));

		if(length(filenames) < 10) next;

		writeLines("\n"); writeLines(paste("\treplenishment ratio: ", replenishment_rate, sep="")); writeLines(paste("\timportation rate: ", importation_rate, sep=""));
		writeLines(paste("\tnetwork size: ", the_size, sep="")); writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
		writeLines(paste("\tphysical topology: ", physical_topology, sep=""));  writeLines(paste("\tphysical degree: ", physical_degree, sep=""));
		writeLines(paste("\tsocial topology: ", social_topology, sep="")); writeLines(paste("\tsocial degree: ", social_degree, sep=""));
		writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));

		Prop_vs_Risk <- data.table();
		for(file in filenames){ temp <- data.table(read.csv(file));
		Prop_vs_Risk <- rbind(Prop_vs_Risk, temp[is.na(instance)], fill=TRUE); }
		Prop_vs_Risk[, ("X"):=NULL];
		write.csv(Prop_vs_Risk, "temp_risk.csv");

		# Prop_vs_Risk <- data.table(read.csv("temp_risk.csv")); Prop_vs_Risk[, X:=NULL];

		metrics <- c();

		if(the_size == 40000)
		{
		    metrics <- filter_these(names(Prop_vs_Risk), no=c("sd", fixed_values[fixed_values!="N"], "dist", "chamber", "conn", "mod", "watts", "ratio", "prob",
		                                                      "opinion", "getis"));
		}
		else if(the_size == 10000)
		{
		    if(old_or_new == "old"){
		        metrics <- filter_these(names(Prop_vs_Risk), no=c("sd", fixed_values[fixed_values!="N"], "dist", "chamber", "conn", "mod",
		                                                          "watts", "ratio", "prob", "opinion", "getis"));
		    } else {
		        metrics <- filter_these(names(Prop_vs_Risk), no=c(fixed_values, "sd", "ratio", "moran", "getis", "join", "geary", "mutual",
		                                                          "phys", "soc"));
		    }
		}
		metrics <- metrics[metrics!="N"];

		the_names_to_group_together <- c();
		if(the_size == 10000){
		    if(old_or_new == "old"){
		        the_names_to_group_together <- c("join.count", "phys", "soc", "moran", "geary", "mutual");
		    }
		    else if(old_or_new == "new"){
		        the_names_to_group_together <- c("phys", "soc", "watts", "number.chambers", "avg.conn.comp", "min.conn.comp", "max.conn.comp",
		                                         "number.conn", "prob", "modularity", "opinion", "avg.chamber.size", "min.chamber.size",
		                                         "max.chamber.size");
		    }
		} else if(the_size == 40000){
		    the_names_to_group_together <- c("join.count", "phys", "soc", "moran", "geary", "mutual");
		}

		col_names <- c("norm", "phys_trans", "soc_trans", names.change(metrics));

		for(CHANGE_TEST in c("lanzante", "pettitt", "buishand_r", "snh")) #
		{
		    writeLines(sprintf("\t\ttest: %s", CHANGE_TEST));
			file_name_ending <- paste("_N_", the_size, "_dur_", duration, "_beta_", beta, "_vaccprop_", initial_vacc_prop, "_inf_", infection_prob,
			                          "_imp_", importation_rate, "_rep_", replenishment_rate, "_top_", physical_topology, "_deg_", physical_degree,
			                          "_switch_", random_opinion_switch, "_mode_", main_or_supplement, ".png", sep=""
			);

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
						soc_norm, transitions$phys_intersec_risk, transitions$soc_intersec_risk,
						data.table(change_points)[, change_points]
					))
				);
			}
			Change_Predictions <- Change_Predictions[rowSums(is.na(Change_Predictions)) != ncol(Change_Predictions)][order(norm)];
			write.csv(Change_Predictions, "temp_change_predictions.csv");

			# Change_Predictions <- data.table(read.csv("temp_change_predictions.csv")); Change_Predictions[, X:=NULL];

			if(the_size == 40000){ Change_Predictions <- Change_Predictions[norm<=2.625]; }

			All_Lead_Times <- Change_Predictions$soc_trans - Change_Predictions[, .SD, .SDcols=names.change(metrics)];
			All_Lead_Times <- cbind(Change_Predictions$norm, All_Lead_Times); names(All_Lead_Times)[1] <- "norm";
			# Filter(function(All_Leadx)!all(is.na(x)), All_Lead_Times);

			Max_or_Min <- data.table(matrix(nrow=nrow(All_Lead_Times), ncol=ncol(All_Lead_Times)-1, "neither"));
			Max_or_Min <- cbind(Change_Predictions$norm, Max_or_Min);
			names(Max_or_Min) <- names(All_Lead_Times);

			for(num_norm in unique(All_Lead_Times$norm)) { for(name in setdiff(names(All_Lead_Times), "norm")) {
			    if(is.na( All_Lead_Times[norm==num_norm, get(name)] )){ Max_or_Min[which(norm==num_norm), (name):="empty"]; }
			}}

			for(i in 1:nrow(All_Lead_Times))
			{
			    rows_with_max <- which(All_Lead_Times[i, ] == max(All_Lead_Times[i, -c("norm")], na.rm=TRUE));
			    rows_with_min <- which(All_Lead_Times[i, ] == min(All_Lead_Times[i, -c("norm")], na.rm=TRUE));
			    for(here in rows_with_max){ Max_or_Min[i, here] <- "max"; }
			    for(here in rows_with_min)
			    {
			        if(Max_or_Min[i, ..here][[1]] != "neither"){ Max_or_Min[i, here] <- "both"; next; }
			        Max_or_Min[i, here] <- "min";
			    }
			}

			Joint <- data.table(EWS=character(), Maxima=numeric(), Minima=numeric());
			for(row_name in setdiff(names(Max_or_Min), "norm"))
			{
			    Joint <- rbind(Joint,
			                   list(
			                       row_name,
			                       length(which(Max_or_Min[, get(row_name)] == "max"))/nrow(Max_or_Min),
			                       length(which(Max_or_Min[, get(row_name)] == "min"))/nrow(Max_or_Min)
			                   )
			    );
			}

			# Joint <- Joint[, difference:=count_max-count_min][order(difference)];

			ordered_rows <- numeric(); # order the table so that the bar chart places related metrics together instead of spreading them around the graph
			for(i in the_names_to_group_together)
			{
			  ordered_rows <- c(ordered_rows, which(grepl(i, Joint$EWS))); # group all the like metrics together
			}
			Joint <- Joint[ordered_rows]; # put the rows in order
			Joint[, Position:=1:nrow(Joint)]; # put in a row that preserves this order, so we can use the reorder() command

			Comparison <- data.table(norm=All_Lead_Times$norm);

			warning_set <- c();
			if(the_size == 10000){
			    if(old_or_new=="old"){ warning_set <- c("Dyn", "WS"); }
			    if(old_or_new=="new"){
			        warning_set <- c("Conn_Echo_Count", "Conn_Size", "Echo_Size", "Prob_Sick", "WattsGCC", "Modu", "Senti");
			        # "V_Conn_Size", "NV_Conn_Size", ", "V_Echo_Size", "NV_Echo_Size");
                }
			}
			else if(the_size == 40000){ warning_set <- c("WS", "Dyn"); }

			for(warnings in warning_set)
			{
			    lets_find_the_maxes <- c();
			    writeLines(sprintf("\t\t\twarnings: %s", warnings));

			    if(warnings == "V_Conn_Size"){ lets_find_the_maxes <- c("vacc_min_conn_comp_size", "vacc_max_conn_comp_size", "vacc_avg_conn_comp_size"); }
				if(warnings == "NV_Conn_Size"){ lets_find_the_maxes <- c("nonvacc_min_conn_comp_size", "nonvacc_max_conn_comp_size", "nonvacc_avg_conn_comp_size"); }
			    if(warnings == "Conn_Size"){ lets_find_the_maxes <- c("vacc_min_conn_comp_size", "vacc_max_conn_comp_size", "vacc_avg_conn_comp_size",
			                                                          "nonvacc_min_conn_comp_size", "nonvacc_max_conn_comp_size", "nonvacc_avg_conn_comp_size"); }

			    if(warnings == "V_Echo_Size"){ lets_find_the_maxes <- c("vacc_min_chamber_size", "vacc_max_chamber_size", "vacc_avg_chamber_size"); }
			    if(warnings == "NV_Echo_Size"){ lets_find_the_maxes <- c("nonvacc_min_chamber_size", "nonvacc_max_chamber_size", "nonvacc_avg_chamber_size"); }
			    if(warnings == "Echo_Size"){ lets_find_the_maxes <- c("vacc_min_chamber_size", "vacc_max_chamber_size", "vacc_avg_chamber_size", "nonvacc_min_chamber_size",
			                                                          "nonvacc_max_chamber_size", "nonvacc_avg_chamber_size") }

			    if(warnings == "Prob_Sick") lets_find_the_maxes <- c("prob_sick_nonvacc_mean", "prob_sick_vacc_mean");
			    if(warnings == "WattsGCC") lets_find_the_maxes <- c("watts_strogatz_nonvacc_mean", "watts_strogatz_vacc_mean", "watts_strogatz_all_mean");
			    if(warnings == "Modu") lets_find_the_maxes <- c("vacc_modularity_mean", "nonvacc_modularity_mean", "total_modularity_mean");
			    if(warnings == "Senti") lets_find_the_maxes <- c("nonvacc_opinion_change_mean", "vacc_opinion_change_mean", "total_opinion_change_mean");

				if(warnings == "WS") lets_find_the_maxes <- c("NN_join_count", "NV_join_count", "VV_join_count", "geary_c", "moran_i", "mutual_info");
				if(warnings == "Dyn") lets_find_the_maxes <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_N", "soc_V");

				if(warnings == "Conn_Echo_Count") lets_find_the_maxes <- c("vacc_number_conn_comps_mean", "nonvacc_number_conn_comps_mean",
				                                                           "nonvacc_number_chambers_mean", "vacc_number_chambers_mean");
				getPalette <- colorRampPalette(brewer.pal(11, "Dark2"));
				colours_here <- getPalette(length(lets_find_the_maxes));
				names(colours_here) <- names.change(lets_find_the_maxes);

				Plot_From_This_Table <- Joint[EWS %in% names.change(lets_find_the_maxes)];

				mode <- "Maxima";

				pl_best <- ggplot(
						Plot_From_This_Table,
						aes(x=EWS, y=get(mode), fill=EWS)
					) +
					geom_bar(stat="identity") +
					geom_text(aes(label=sprintf("%i%%", round(100*get(mode), 0))), vjust=-0.1, size=10, colour="black") +
					scale_y_continuous(limits=c(0, Plot_From_This_Table[, max(get(mode))+0.015])) +
					scale_fill_manual(values=colours_here, name="", breaks=Plot_From_This_Table$EWS, labels=as.vector(sapply(Plot_From_This_Table$EWS, axis_label))) +
					ylab("# Maxima") +
					theme(
					  panel.background=element_blank(),
					  axis.line=element_line(colour = "black"),
					  axis.text=element_blank(),
					  axis.title=element_text(size=title_size),
					  axis.title.x=element_blank(),
					  legend.text=element_text(size=text_size),
					  legend.title=element_blank(),
					  axis.text.x=element_blank(),
					  axis.ticks.x=element_blank(),
					  legend.position=c(0.03, 0.5)
					) +
					guides(fill = guide_legend(nrow=1));

				if(warnings == "Dyn"){
				  pl_best <- pl_best + theme(legend.key.size = unit(1,"line"), legend.key.width=unit(5, "line"), legend.spacing.x = unit(1, 'cm'));
				}	else{
				  pl_best <- pl_best + theme(legend.key.size = unit(2,"line"), legend.key.width=unit(3.5, "line"), legend.spacing.x = unit(0.6, 'cm'));
				}

				pl_legend <- g_legend(pl_best);
				dev.off();

				Lead_Times_Here <- All_Lead_Times[, .SD, .SDcols=intersect(names.change(lets_find_the_maxes), names(All_Lead_Times))];
				Lead_Times_Here <- cbind(Change_Predictions$norm, Lead_Times_Here);
				names(Lead_Times_Here)[1] <- "norm";
				Lead_Times_Here <- melt(Lead_Times_Here, id="norm");

				pl_EWS <- ggplot(Lead_Times_Here, aes(norm, value)) +
					geom_line(aes(colour=variable), size=1) +
					scale_colour_manual(values=colours_here) +
                    ylim(c(min(Lead_Times_Here$value, na.rm=TRUE), max(Lead_Times_Here$value, na.rm=TRUE)+0.04)) +
				    labs(x=expression(sigma), y=bquote('K'['s']*'-'*.(axis_label(CHANGE_TEST))*'{'*.(axis_label(warnings))*'}')) +
				    theme(
						panel.background=element_blank(),
						axis.line=element_line(colour = "black"),
						axis.text=element_text(size=text_size),
						axis.title=element_text(size=title_size),
						legend.position="none"
					) +
					geom_hline(yintercept=0, size=1, colour="black", linetype="dashed");

				left_pl <- grid.arrange(
					pl_EWS + theme(legend.position="none"),
					pl_legend,
					nrow=2,
					heights=c(10,1)
				);
				right_pl <- pl_best + theme(legend.position="none");
				pl_truth <- plot_grid( left_pl, right_pl, rel_widths=c(7.5,2.5) )
				pl_truth <- cowplot::ggdraw(pl_truth) + theme(plot.background=element_blank());

				graph_name <- paste(graph_file_path, folder_name, "/", "LeadTime", sep='');
				# if(the_size == 10000){ graph_name <- paste(graph_name, '_', old_or_new, sep=''); }
				graph_name <- paste(graph_name, '_', warnings, "_", CHANGE_TEST, file_name_ending, sep="");

				ggsave(
					graph_name,
					plot = pl_truth, width=image_width, height=image_height, limitsize=FALSE, dpi=50
				);
				dev.off();
			}
		}
	}
}
print(Sys.time() - start_time);

graphics.off()
