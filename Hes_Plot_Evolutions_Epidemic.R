# Brendon Phillips
# PhD candidate
# Department of Mathematics
# University of Waterloo

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Time_Series", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Parameters <- fread(parameter_file, colClasses="character")
Parameters[, "V1":=NULL];
Parameters[network_structure=="smallworld"];

Parameters <- Parameters[num(initial_vacc_proportion)==0.25 & num(social_norm)==0 & num(perceived_vaccine_risk)%in%c(-0.0125, 0, 0.00625)]; # num(proportion_of_nodes)==1 &
Parameters <- unique(Parameters);

quantities_to_plot <- c("phys_I");

text_size <- 35;
title_size <-45;

time_start <- Sys.time();

for(line_index in 1:nrow(Parameters))
{
	this_tuple <- Parameters[line_index];

	initial_vacc_prop <- (this_tuple$initial_vacc_proportion);
	size <- convert_if_number(this_tuple$N);
	risk <- (this_tuple$perceived_vaccine_risk);
	infection_prob <- (this_tuple$infec_prob);
	importation_rate <- (this_tuple$importation);
	structure <- this_tuple$network_structure;
	degree <- (this_tuple$mean_degree);
	random_opinion_switch <- this_tuple$random_opinion_switch;
	social_norm <- (this_tuple$social_norm);
	proportion <- (this_tuple$proportion_of_nodes);

	data_file_list <- Sys.glob(
		paste(
			stored_data_path,
			"/Hes",
			"_N_", size,
			"_struct_", structure,
			"_deg_", degree,
			"_risk_", risk,
			"_infec_", infection_prob,
			"_norm_", social_norm,
			"_imp_", importation_rate,
			"_init_", initial_vacc_prop,
			"_switch_", random_opinion_switch,
			"_prop_", proportion,
			"_inst_*.csv",
			sep=""
		)
	);

    if(length(data_file_list) == 0){ next; }

	Summary_Data <- data.table();

	for(file_index in 1:length(data_file_list))
	{
		if(file.size(data_file_list[file_index]) == 0) { next; }
		temp <- fread(data_file_list[file_index]);
		if(nrow(temp) == 0){ next; }
		if(nrow(Summary_Data) != 0)
		{
			# just in case the trial number is not unique, make it so
			if(temp[, unique(instance)] %in% unique(Summary_Data[, instance]))
			{
				temp_trial <- 0;
				while(temp_trial %in% unique(Summary_Data[, instance])) { temp_trial = temp_trial + 1; }
				temp[, ("instance"):=rep(temp_trial, nrow(temp))]
			}
		}
		Summary_Data <- rbind(temp[time<=6], Summary_Data, fill=TRUE);
  	}

	Summary_Data[, "V1":=0]
	props <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_N", "soc_H", "soc_V");
	Summary_Data[, (props):=lapply(.SD, function(x) x/size), .SDcols=props];

	for(metric in quantities_to_plot){ if(metric %in% names(Summary_Data))
	{
		if(length(Summary_Data[[metric]]) == 0) next;

		plot_file_name <-paste(
		    graph_file_path,
		    folder_name, "/",
		    "Hes_", metric,
		    "_N_", size,
		    "_struct_", structure,
		    "_deg_", degree,
		    "_risk_", risk,
		    "_infec_", infection_prob,
		    "_norm_", social_norm,
		    "_imp_", importation_rate,
		    "_init_", initial_vacc_prop,
		    "_switch_", random_opinion_switch,
		    "_prop_", proportion,
		    "_Epidemic.png",
		    sep=""
		);

		pl <- ggplot(Summary_Data, aes(time, get(metric), group=instance)) +
		    geom_line() +
		    labs(x=expression(tau), y=axis_label(metric)) +
		    theme(
		        panel.background=element_blank(),
		        axis.line=element_line(colour = "black"),
		        axis.text=element_text(size=text_size),
		        axis.title=element_text(size=title_size),
		        legend.position="none"
		    )

		ggsave(
		    plot_file_name,
		    plot = pl, width=10, height=4, limitsize=FALSE, dpi=50
		);
		dev.off();
	}}
}

# if( exists("cl") ) parallel::stopCluster(cl);
# print("done")
# print(Sys.time()-time_start)
