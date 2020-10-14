# Dear God, why have I retooled this code so many fucking times???

# Summary Plot for the Topology
rm(list=ls());

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Lead_Time", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character")
Params[, c("V1", "risk", "norm"):=NULL];

start_time <- Sys.time();

# cl <- parallel::makeCluster(6)
# doParallel::registerDoParallel(cl)

main_or_supplement <- "main";

if(main_or_supplement == "main")
{
	image_width <- 25;
	image_height <- 5;
	text_size <- 35;
	title_size <-45;
} else if(main_or_supplement == "supplement")
{
  image_width <- 25;
  image_height <- 5;
  text_size <- 45;
  title_size <-60;
}

for(the_size in c("40000")) # "10000",
{
	Parameters <- unique(Params[phys_top=="random" & num(init_prop)==0.05 & num(beta)==1 & num(size)==the_size]);

	if(num(the_size)==10000){ Parameters <- Parameters[num(import)==0.00025]; }
	else if(num(the_size)==40000){ Parameters <- Parameters[num(infec)==0.2]; }
	else if(num(the_size)==562500){ Parameters <- Parameters[]; }

	Parameters <- unique(Parameters);

	for(line_index in 1:nrow(Parameters))
	{
		this_tuple <- Parameters[line_index];

		initial_vacc_prop <- this_tuple[, init_prop];
		size <- this_tuple[, size];
		beta <- this_tuple[, beta];
		duration <- this_tuple[, duration];
		infection_prob <- this_tuple[, infec];
		importation_rate <- this_tuple[, import];
		replenishment_rate <-  this_tuple[, birth_death];
		random_opinion_switch <- this_tuple[, random_switch];
		physical_topology <- this_tuple[, phys_top];
		physical_degree <- this_tuple[, phys_deg];
		social_topology <- this_tuple[, soc_top];
		social_degree <- this_tuple[, soc_deg];

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

		if(length(filenames)==0){ next;	}

		writeLines("\n");
		writeLines(paste("\treplenishment ratio: ", replenishment_rate, sep=""));
		writeLines(paste("\timportation rate: ", importation_rate, sep=""));
		writeLines(paste("\tnetwork size: ", size, sep=""));
		writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
		writeLines(paste("\tphysical topology: ", physical_topology, sep=""));
		writeLines(paste("\tphysical degree: ", physical_degree, sep=""));
		writeLines(paste("\tsocial topology: ", social_topology, sep=""));
		writeLines(paste("\tsocial degree: ", social_degree, sep=""));
		writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));

		if(length(filenames) < 10) next;

		Prop_vs_Risk <- data.table();
		for(file in filenames)
		{
			temp <- data.table(read.csv(file));
			Prop_vs_Risk <- rbind(Prop_vs_Risk, temp[is.na(instance)], fill=TRUE);
		}
		Prop_vs_Risk[, ("X"):=NULL];
		write.csv(Prop_vs_Risk, "temp_risk.csv");

		# Prop_vs_Risk <- data.table(read.csv("temp_risk.csv"));
		# Prop_vs_Risk[, X:=NULL];

		metrics <- c();

		if(size == "40000")
		{
			metrics <- filter_these(names(Prop_vs_Risk), no=c("sd", fixed_values[fixed_values!="N"], "dist", "chamber", "conn", "mod", "watts", "ratio", "prob", "opinion"));
			metrics <- metrics[metrics!="N"];
		}
		else if(size == "10000")
		{
			# metrics <- filter_these(names(Prop_vs_Risk), no=c(fixed_values, "sd", "ratio", "moran", "getis", "join", "geary", "mutual"));
		  metrics <- filter_these(names(Prop_vs_Risk), no=c("sd", fixed_values[fixed_values!="N"], "dist", "chamber", "conn", "mod", "watts", "ratio", "prob", "opinion"));

			metrics <- metrics[metrics!="N"];
		}

		# setting the colour scheme for the plotted warning signals
		# this comes in handy later, when we do the side-by-side of the graph and the bar chart
		# we can make sure that the colours of the variables remain the same

		colours <- rainbow(length(metrics));
		names(colours) <- names.change(metrics);

		col_names <- c("norm", "phys_trans", "soc_trans", names.change(metrics));

		for(CHANGE_TEST in c("lanzante")) # , "pettitt", "snh", "buishand_r"
		{
			file_name_ending <- paste(
				"_N_", size,
				"_dur_", duration,
				"_beta_", beta,
				"_vaccprop_", initial_vacc_prop,
				"_inf_", infection_prob,
				"_imp_", importation_rate,
				"_rep_", replenishment_rate,
				"_top_", physical_topology,
				"_deg_", physical_degree,
				"_switch_", random_opinion_switch,
				"_mode_", main_or_supplement,
				".png",
				sep=""
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

			# Change_Predictions <- data.table(read.csv("temp_change_predictions.csv"));
			# Change_Predictions[, X:=NULL];

			Change_Predictions <- Change_Predictions[norm<=2.625];

			All_Lead_Times <- Change_Predictions$soc_trans - Change_Predictions[, .SD, .SDcols=names.change(metrics)];

			All_Maxima <- data.table(table(names(unlist(apply(All_Lead_Times, 1, function(x) which(x==max(x, na.rm=TRUE)))))));
			All_Maxima[, N:=round(N/nrow(Change_Predictions),2)]
			names(All_Maxima) <- c("EWS", "count_max");

			All_Minima <- data.table(table(names(unlist(apply(All_Lead_Times, 1, function(x) which(x==min(x, na.rm=TRUE)))))));
			All_Minima[, N:=round(N/nrow(Change_Predictions),2)]
			names(All_Minima) <- c("EWS", "count_min");

			All_Lead_Times <- cbind(Change_Predictions$norm, All_Lead_Times);
			names(All_Lead_Times)[1] <- "norm";

			Joint <- merge(All_Maxima, All_Minima, by="EWS", all=TRUE);
			Joint[is.na(Joint)] <- 0;
			Joint <- Joint[, difference:=count_max-count_min][order(difference)];

			ordered_rows <- numeric(); # order the table so that the bar chart places related metrics together instead of spreading them around the graph
			the_names_to_group_together <- c();
			# if(size == "10000"){ the_names_to_group_together <- c("phys", "soc", "watts", "avg.chamber.size", "min.chamber.size", "max.chamber.size", "number.chambers", "avg.conn.comp", "min.conn.comp", "max.conn.comp", "number.conn", "prob", "modularity", "opinion"); }
			if(size == "10000"){ the_names_to_group_together <- c("phys", "soc", "geary", "moran", "mutual", "join.count"); }

			else if(size == "40000"){ the_names_to_group_together <- c("join.count", "phys", "soc", "moran", "geary", "mutual"); }
			for(i in the_names_to_group_together)
			{
				ordered_rows <- c(ordered_rows, which(grepl(i, Joint$EWS))); # group all the like metrics together
			}
			Joint <- Joint[ordered_rows]; # put the rows in order
			Joint[, Position:=1:nrow(Joint)]; # put in a row that preserves this order, so we can use the reorder() command

			################################################################################################################################################

			for(mode in c("max", "min"))
			{
				pl_performance <- ggplot(
						Joint,
						aes(x=reorder(EWS, Position), y=get(sprintf("count_%s", mode)), fill=EWS)
					) +
					geom_bar(stat="identity") +
					geom_hline(yintercept=0.25, linetype="dashed", colour="black") +
					geom_text(aes(label=sprintf("%i%%", round(100*get(sprintf("count_%s", mode)), 0))), vjust=-0.1, size=(if(size=="40000") 13 else 8), colour="black") +
					scale_x_discrete(labels=as.vector(sapply(Joint$EWS, axis_label))) +
					scale_y_continuous(limits=c(0, Joint[, max(get(sprintf("count_%s", mode)))+0.04])) +
					scale_fill_manual("legend", values=colours) +
					ylab(sprintf( "Ratio of %s lead times", (if(mode=="max") "best" else "worst") )) +
					theme(
						axis.text.y=element_blank(),
						axis.text=element_text(size=text_size),
						axis.title.y=element_text(size=title_size),
						axis.title.x=element_blank(),
						panel.background=element_blank(),
						legend.position="none",
						axis.text.x=(
							if(size=="10000") element_text(angle=45, hjust=1) else if(size=="40000") element_text(angle=30, hjust=1)
						)
					);

				# ggsave(
				# 	paste(graph_file_path, folder_name, "/", "LeadTime_", mode, "_performance_", CHANGE_TEST, file_name_ending, sep=""),
				# 	plot = pl_performance, width=image_width, height=image_height+2, units="in", limitsize=FALSE
				# );
				# dev.off();
			}

			################################################################################################################################################

			differences <- ggplot(Joint, aes(x=EWS, y=difference, fill=EWS)) + geom_bar(stat="identity") + theme(legend.position="none");

			Comparison <- data.table(norm=All_Lead_Times$norm);

			warning_set <- c();
			# if(size == "10000"){ warning_set <- c("V_Conn_Size", "NV_Conn_Size", "Conn_Echo_Count", "V_Echo", "NV_Echo", "Prob_Sick", "WattsGCC", "Modu", "Senti"); }
			# else if(size == "40000"){ warning_set <- c("WS"); }
			warning_set <- c("WS");
			warning_set <- c("Dyn", warning_set);

			for(warnings in warning_set)
			{
				if(warnings == "WS") lets_find_the_maxes <- c("NN_join_count", "NV_join_count", "VV_join_count", "geary_c", "moran_i", "mutual_info");
				if(warnings == "Dyn") lets_find_the_maxes <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_N", "soc_V");
				if(warnings == "V_Conn_Size") lets_find_the_maxes <- c("vacc_min_conn_comp_size", "vacc_max_conn_comp_size", "vacc_avg_conn_comp_size");
				if(warnings == "NV_Conn_Size") lets_find_the_maxes <- c("nonvacc_min_conn_comp_size", "nonvacc_max_conn_comp_size", "nonvacc_avg_conn_comp_size")
				if(warnings == "Conn_Echo_Count") lets_find_the_maxes <- c("vacc_number_conn_comps_mean", "nonvacc_number_conn_comps_mean", "nonvacc_number_chambers_mean", "vacc_number_chambers_mean");
				if(warnings == "V_Echo") lets_find_the_maxes <- c("vacc_min_chamber_size", "vacc_max_chamber_size", "vacc_avg_chamber_size");
				if(warnings == "NV_Echo") lets_find_the_maxes <- c("nonvacc_min_chamber_size", "nonvacc_max_chamber_size", "nonvacc_avg_chamber_size");
				if(warnings == "Prob_Sick") lets_find_the_maxes <- c("prob_sick_nonvacc_mean", "prob_sick_vacc_mean");
				if(warnings == "WattsGCC") lets_find_the_maxes <- c("watts_strogatz_nonvacc_mean", "watts_strogatz_vacc_mean", "watts_strogatz_all_mean");
				if(warnings == "Modu") lets_find_the_maxes <- c("vacc_modularity_mean", "nonvacc_modularity_mean", "total_modularity_mean");
				if(warnings == "Senti") lets_find_the_maxes <- c("nonvacc_opinion_change_mean", "vacc_opinion_change_mean", "total_opinion_change_mean");

				colours_here <- rainbow(length(lets_find_the_maxes));
				names(colours_here) <- names.change(lets_find_the_maxes);

				Plot_From_This_Table <- Joint[EWS %in% names.change(lets_find_the_maxes)];

				mode <- "max";

				pl_best <- ggplot(
						Plot_From_This_Table,
						aes(x=EWS, y=count_max, fill=EWS)
					) +
					geom_bar(stat="identity") +
					geom_text(aes(label=sprintf("%i%%", round(100*get(sprintf("count_%s", mode)), 0))), vjust=-0.1, size=9, colour="black") +
					scale_y_continuous(limits=c(0, Plot_From_This_Table[, max(get(sprintf("count_%s", mode)))+0.01])) +
					scale_fill_manual(values=colours_here, name="", breaks=Plot_From_This_Table$EWS, labels=as.vector(sapply(Plot_From_This_Table$EWS, axis_label))) +
					ylab("Ratio of best lead times") +
					theme(
						panel.background=element_blank(),
						axis.line=element_line(colour = "black"),
						axis.text=element_blank(),
						axis.title=element_text(size=text_size),
						axis.title.x=element_blank(),
						legend.text=element_text(size=text_size),
						legend.title=element_blank(),
						axis.text.x=element_blank(),
						axis.ticks.x=element_blank(),
						legend.position=c(0.03, 0.5)
					) +
					guides(fill = guide_legend(nrow=1));

				if(warnings == "Dyn"){ pl_best <- pl_best + theme(legend.key.size = unit(1.5,"line"), legend.key.width=unit(2, "line"), legend.spacing.x = unit(1, 'cm')); }
				else{ pl_best <- pl_best + theme(legend.key.size = unit(2,"line"), legend.key.width=unit(2, "line"), legend.spacing.x = unit(0.6, 'cm')); }

				pl_legend <- g_legend(pl_best); dev.off();

				Lead_Times_Here <- All_Lead_Times[, .SD, .SDcols=intersect(names.change(lets_find_the_maxes), names(All_Lead_Times))];
				Lead_Times_Here <- cbind(Change_Predictions$norm, Lead_Times_Here);
				names(Lead_Times_Here)[1] <- "norm";
				Lead_Times_Here <- melt(Lead_Times_Here, id="norm");

				pl_EWS <- ggplot(Lead_Times_Here, aes(norm, value)) +
					geom_line(aes(colour=variable), size=1) +
					scale_fill_manual("legend", values=colours_here) +
					labs(x=expression(sigma), y=bquote('K'['s']*'-B'[.(axis_label(warnings))]^.(axis_label(CHANGE_TEST)))) +
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

				# ggsave(
				# 	paste(graph_file_path, folder_name, "/", "LeadTime_", warnings, "_", CHANGE_TEST, file_name_ending, sep=""),
				# 	plot = pl_truth, width=image_width, height=image_height, units="in", limitsize=FALSE
				# );
				# dev.off();

				# R list comprehension from the "comprehenr" package
				Max_min_lead_times_from_this_EWS_class <- data.table(
					m = to_vec(for(i in 1:nrow(All_Lead_Times)) min(All_Lead_Times[i, .SD, .SDcols=intersect(names.change(lets_find_the_maxes), names(All_Lead_Times))], na.rm=T)),
					M = to_vec(for(i in 1:nrow(All_Lead_Times)) max(All_Lead_Times[i, .SD, .SDcols=intersect(names.change(lets_find_the_maxes), names(All_Lead_Times))], na.rm=T))
				);
				names(Max_min_lead_times_from_this_EWS_class) <- sprintf("%s_%s", warnings, c("min", "max"));

				Comparison <- cbind(Comparison, Max_min_lead_times_from_this_EWS_class);

				if(warnings != "Dyn")
				{
					Chi_DT <- data.table(
						norm = Comparison$norm,
						chi_min = Comparison[, get(sprintf("%s_min", warnings))-Dyn_min],
						chi_max = Comparison[, get(sprintf("%s_max", warnings))-Dyn_max]
					);

					epsilon_min = abs(max(Chi_DT$chi_min)-min(Chi_DT$chi_min))/100;
					epsilon_max = abs(max(Chi_DT$chi_max)-min(Chi_DT$chi_max))/100;

					warnings_min_bigger <- 	round(sum(Chi_DT$chi_min > epsilon_min)/nrow(Comparison)*100, 2);
					warnings_min_equal <-	round(sum(abs(Chi_DT$chi_min) < epsilon_min)/nrow(Comparison)*100, 2);
					warnings_min_smaller <-	round(sum(Chi_DT$chi_min < -epsilon_min)/nrow(Comparison)*100, 2);

					warnings_max_bigger <- 	round(sum(Chi_DT$chi_max > epsilon_max)/nrow(Comparison)*100, 2);
					warnings_max_equal <-	round(sum(abs(Chi_DT$chi_max) < epsilon_max)/nrow(Comparison)*100, 2);
					warnings_max_smaller <- round(sum(Chi_DT$chi_max < -epsilon_max)/nrow(Comparison)*100, 2);

					chi_min_neg <- round(sum(Chi_DT$chi_min < -epsilon_min)/nrow(Comparison)*100, 1);
					chi_min_zer <- round(sum(abs(Chi_DT$chi_min) <= epsilon_min)/nrow(Comparison)*100, 1);
					chi_min_pos <- 100-chi_min_neg-chi_min_zer;

					chi_max_neg <- round(sum(Chi_DT$chi_max < -epsilon_max)/nrow(Comparison)*100, 1);
					chi_max_zer <- round(sum(abs(Chi_DT$chi_max) <= epsilon_max)/nrow(Comparison)*100, 1);
					chi_max_pos <- 100-chi_max_neg-chi_max_zer;

					mytable <- data.table();

					if(main_or_supplement == "main"){
					  mytable <- data.table(
  						" "=c("Xmin", 	"Xmax"),
  						neg=c(sprintf("%g%%", chi_min_neg),		sprintf("%g%%", chi_max_neg)),
  						zero=c(sprintf("%g%%", chi_min_zer),	sprintf("%g%%", chi_max_zer)),
  						pos=c(sprintf("%g%%", chi_min_pos),		sprintf("%g%%", chi_max_pos))
					  );
					} else if(main_or_supplement == "supplement"){
					  mytable <- data.table(
					    " "=c("neg", 	"zero", "pos"),
					    Xmin=c(sprintf("%g%%", chi_min_neg), sprintf("%g%%", chi_min_zer), sprintf("%g%%", chi_min_pos)),
					    Xmax=c(sprintf("%g%%", chi_max_neg), sprintf("%g%%", chi_max_zer), sprintf("%g%%", chi_max_pos))
					  );
					}

					tt <- gridExtra::ttheme_default(base_size=text_size, core = list(padding=unit(c(4, 9), "mm")));

					pl <- ggplot(Chi_DT);

				 pl_leadtime_comparison <- ggplot(Chi_DT, aes(x=norm)) +
					  ylim(c(Chi_DT[, min(chi_min)], Chi_DT[, max(chi_min)+0.1])) +
						geom_hline(yintercept=0, "dashed") +
						geom_rect(data=Chi_DT, mapping=aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf), fill='green', alpha=0.008) +
						geom_rect(data=Chi_DT, mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), fill='red', alpha=0.008) +
						labs(x=expression(sigma), y=bquote(chi[''['*']]^.(axis_label(CHANGE_TEST)))) +
					  geom_point(aes(y=chi_min), colour="blue", size=3) + geom_line(aes(y=chi_min), colour="blue", size=1.5) +
					  geom_point(aes(y=chi_max), colour="red", size=3)  + geom_line(aes(y=chi_max), colour="red", size=1.5)   +
						theme(
							panel.background=element_blank(),
							axis.line=element_line(colour = "black"),
							axis.text=element_text(size=text_size+15),
							axis.title=element_text(size=title_size+15),
							legend.position="none"
						) +
						geom_hline(yintercept=0, size=0.5, colour="black", linetype="dashed");

          if(main_or_supplement == "main")
          {
            pl <- pl_leadtime_comparison + annotation_custom(tableGrob(mytable, rows=NULL ,theme=tt), xmin=1.7, xmax=Inf, ymin=0.1, ymax=Inf); # ymax=0.8*max(Comparison[, -c("norm")])
          } else if(main_or_supplement == "supplement"){
            tbl <- tableGrob(mytable, rows=NULL, theme=tt) # , widths=unit(1, "null"), heights=unit(1/nrow(mytable)),"npc");
            pl <- plot_grid(pl_leadtime_comparison, tbl, rel_widths=c(7,2), ncol=2, nrow=1);
            pl <- cowplot::ggdraw(pl) + theme(plot.background=element_blank());
          }

					ggsave(
						paste(graph_file_path, folder_name, "/", "LeadTime_comparison_", warnings, "_", CHANGE_TEST, file_name_ending, sep=""),
						plot = pl, width=image_width, height=image_height, units="in", limitsize=FALSE
					);
					dev.off();

					# stop(Sys.time() - start_time);
				}
			}
			# stop()
		}
	}
	# stop()
}
print(Sys.time() - start_time);

graphics.off()
