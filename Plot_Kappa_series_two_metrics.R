# This is the file for the first diagram tableau in the figure - the oen for multiple transitions is in another file
rm(list=ls())

options(scipen=10000)

source("/home/b2philli/Dropbox/Processing/Parameter_Values.R");

folder_name <- paste("Kappa_Series", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# source("Processing_Parameter_Table.R");
Params <- fread(parameter_file, colClasses="character");

image_height <- 4;
image_width <- 25;

# metrics_to_plot <- c("conn_comp_size", "echo_chamber_size", "plain", "number", "modularity", "watts_strogatz", "opinion_change", "prob_sick", "echo_chamber_number", "conn_comp_number", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "vacc_conn_comp_size", "nonvacc_conn_comp_size");

time_start <- Sys.time();

mains <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_V", "soc_N");
mains_mean <- unname(sapply(mains, function(x) sprintf("%s_mean",x)));
mains_sd <- unname(sapply(mains, function(x) sprintf("%s_sd",x)));

joins <- c("NN_join_count", "NV_join_count", "VV_join_count");
joins_mean <- unname(sapply(joins, function(x) sprintf("%s_mean",x)));
joins_sd <- unname(sapply(joins, function(x) sprintf("%s_sd",x)));

for(the_pair in c( list(c(size=40000, mode="old")) )) # "10000",list(c(size=10000, mode="old")),
{
	the_size <- the_pair[["size"]];
	old_or_new <- the_pair[["mode"]];
	classes_to_plot <- c();

	Parameters <- Params[num(size)==the_size & num(beta)==1 & soc_top=="random" & phys_top=="random" & num(init_prop)==0.05];

	if(the_size == 40000){
		Parameters <- Parameters[num(beta)==1 & num(init_prop==0.05) & num(birth_death)==0.00024 & num(risk)>=-0.6 & num(risk)<=-0.3 & soc_top=="random"];
		classes_to_plot <- c("plain", "mutual_info");
		norms_to_use <- c(1.125,1.5,1.625,1.75,2.125,2.25);
	}
	else if(the_size == 10000){
		if(old_or_new == "old"){ classes_to_plot <- c("geary_c", "mutual_info", "NV_join_count", "moran_i", "join_counts", "plain"); }
		if(old_or_new == "new"){ Parameters <- Parameters[num(init_prop)==0.05 & birth_death==0.00024]; }
	}
	# if(num(the_size) == 562500) Parameters <- Parameters[num(init_prop)==0.05 & num(norm)%in%c(0,0.25) & num(random_switch)==0.1];

	Parameters[, c("V1", "risk", "time"):=NULL];
	Parameters <- unique(Parameters);

	for(line_index in 1:nrow(Parameters))
	{
		this_tuple <- Parameters[line_index];

		beta <- this_tuple[, beta];
		initial_vacc_prop <- this_tuple[, init_prop];
		size <- this_tuple[, size];
		duration <- this_tuple[, duration];
		infection_prob <- this_tuple[, infec];
		importation_rate <- this_tuple[, import];
		replenishment_rate <-  this_tuple[, birth_death];
		physical_topology <- this_tuple[, phys_top];
		physical_degree <- this_tuple[, phys_deg];
		social_topology <- this_tuple[, soc_top];
		social_degree <- this_tuple[, soc_deg];
		random_opinion_switch <- this_tuple[, random_switch];
		num_size <- num(size);
		the_norm <- this_tuple[, norm];

		Proportion <- data.table();
		Transitions <-  data.table(norm=numeric(), soc=numeric(), phys=numeric());

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
			"_norm_", the_norm,
			"_switch_", random_opinion_switch,
			".csv",
			sep=""
		));

		for(file in filenames)
		{
			temp <- data.table(read.csv(file));
			Proportion <- rbind(Proportion, temp[is.na(instance)], fill=TRUE);
		}

		Trans_Table <- get_transitions(Proportion, initial_vacc_prop);

		# Proportion[, (mains_sd)  :=lapply(.SD, "/", num(the_size)), .SDcols=mains_sd  ];
		# Proportion[, (mains_mean):=lapply(.SD, "/", num(the_size)), .SDcols=mains_mean];
		#
		# # do the sd first before normalising the means, since normalising the means changes to total 1
		# # normalising the sd's by one doesn't do fuck all, then
		# Proportion[, (joins_sd)  :=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean), .SDcols=joins_sd ];
		# Proportion[, (joins_mean):=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean), .SDcols=joins_mean];

		# if there are less than 5 points to be plotted, then cancel this plot and skip to the next parameter set
		if(nrow(Proportion)<5){ next; }

		text_size <- 30;

		The_Double <- data.table(risk=numeric(), quant=numeric(), variable=numeric(), value=numeric());

		for(index in 1:nrow(Proportion))
		{
			the_row <- Proportion[index];
			for(the_class in classes_to_plot)
			{
				for(metric in plot_these(the_class))
				{
					The_Double <- rbind(The_Double, list(the_row$perceived_vaccine_risk, metric, the_class, the_row[[sprintf("%s_mean", metric)]]));
				}
			}
		}

		pl <- ggplot(data=The_Double, aes(x=risk, y=value, colour=quant)) + # , ymin=lwr, ymax=upr
		geom_line(size=2) +
		xlab(expression(kappa)) + #ylab(axis_label(variable)) +
		theme(
			# legend.position="bottom",
			# legend.direction="horizontal",
			# legend.key.width=unit(4, 'cm'),
			panel.background=element_blank(),
			axis.text=element_text(size=text_size),
			legend.text=element_text(size=text_size+10),
			# strip.text.x = element_text(size = text_size),
			axis.title=element_text(size=text_size+10),
			panel.spacing = unit(3, "lines"),
			legend.title=element_blank(),
			legend.key.width=unit(1, 'cm'),
			# legend.spacing.x=unit(1, 'cm'),
			legend.direction="vertical",
			legend.justification='center',
			legend.position="right",
		) +
		facet_wrap(~variable, nrow=1, scales='free', labeller=label_bquote(sigma*' = '*.(norm)));

		# scale_colour_discrete(labels=as.vector(sapply(plots_to_make, axis_label))) +
		# scale_fill_discrete(labels=as.vector(sapply(plots_to_make, axis_label))) +
		# geom_vline(aes(xintercept=phys), data=Transitions, size=1, linetype="dashed") +
		# geom_vline(aes(xintercept=soc),  data=Transitions, size=1, linetype="dotted") +

		plot_file_name <- paste(graph_file_path, folder_name, "/indi_twin_", metric, "_N_", size, "_dur_", duration, "_beta_", beta,
			"_vaccprop_", initial_vacc_prop, "_inf_", infection_prob, "_imp_", importation_rate, "_rep_", replenishment_rate,
			"_ptop_", physical_topology, "_pdeg_", physical_degree, "_stop_", social_topology, "_sdeg_", social_degree,
			"_norm_", min(norms_to_use), "_", max(norms_to_use), "_switch_", random_opinion_switch, ".png", sep=""
		);


		stop();


		stop()

		# for(metric in metrics_to_plot)
		# {
		# 	# if(metric != "plain"){
		# 	# for(the_thing in plot_these(metric))
		# 	# {
		# 	# 	if(!( the_thing %in% names(metric) )) next;
		# 	# }}
		#
		# 	plots_to_make <- plot_these(metric);
		# 	DT <- data.table(grp=character(), norm=numeric(), risk=numeric(), val=numeric(), lwr=numeric(), upr=numeric());
		#
		# 	for(the_norm in norms_to_use){ for(index in seq_along(plots_to_make))
		# 	{
		# 		Proportion_Here <- Proportion[social_norm==the_norm];
		#
		# 		prefix <- plot_these(metric)[index];
		# 		tmean <- sprintf("%s_mean", prefix);
		# 		tsd <- sprintf("%s_sd", prefix);
		#
		# 		temp_table <- data.table(
		# 			grp = rep(prefix, nrow(Proportion_Here)),
		# 			norm = rep(the_norm, nrow(Proportion_Here)),
		# 			risk = Proportion_Here[, perceived_vaccine_risk],
		# 			val = Proportion_Here[, get(tmean)],
		# 			lwr = Proportion_Here[, get(tmean)-get(tsd)],
		# 			upr = Proportion_Here[, get(tmean)+get(tsd)]
		# 		);
		#
		# 		DT <- data.table(do.call(smartbind, list(DT, temp_table)));
		# 	}}
		#
		# 	text_size <- 30;
		#
		# 	scientific_10 <- function(x) {
		# 		parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
		# 	}
		#
		# 	if(metric %in% c("vacc_conn_comp_size", "nonvacc_conn_comp_size", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "num_triads")){ pl <- pl + scale_y_continuous(labels=scientific_10); }

			# plot_file_name <- paste(graph_file_path, folder_name, "/indi_twin_", metric, "_N_", size, "_dur_", duration, "_beta_", beta, "_vaccprop_", initial_vacc_prop,
			# 	"_inf_", infection_prob, "_imp_", importation_rate, "_rep_", replenishment_rate, "_ptop_", physical_topology, "_pdeg_", physical_degree,
			# 	"_stop_", social_topology, "_sdeg_", social_degree, "_norm_", min(norms_to_use), "_", max(norms_to_use), "_switch_", random_opinion_switch, ".png", sep=""
			# );

			ggsave(
				plot_file_name,
				plot=pl, width=image_width, height=image_height, limitsize=FALSE, dpi=50 #  height=5,
			)
			dev.off();

			system(sprintf("convert %s -trim %s", plot_file_name, plot_file_name));
		}
	}
}

print(Sys.time()-time_start);
