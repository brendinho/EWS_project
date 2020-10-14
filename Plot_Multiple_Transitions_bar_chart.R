# This is the file for the first diagram tableau in the figure - the oen for multiple transitions is in another file
rm(list=ls())

options(scipen=10000)

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Number_of_Transitions", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# source("Processing_Parameter_Table.R");
Params <- fread(parameter_file, colClasses="character")
# Parameters <- Params[num(beta)==1 & num(init_prop)==0.05 & num(size)==40000 & num(infec)==0.2];
Parameters <- Params[num(beta)==1 & num(init_prop)==0.05 & num(size)==10000 & num(infec)==0.2 & num(import)==0.00025];
Parameters[, c("V1", "risk", "time", "norm"):=NULL];

Parameters <- unique(Parameters);

time_start <- Sys.time();

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

    Prop_vs_Risk = data.table();

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

	file_name_ending <- paste(
		"_N_", size,
		"_dur_", duration,
		"_beta_", beta,
		"_vaccprop_", initial_vacc_prop,
		"_inf_", infection_prob,
		"_imp_", importation_rate,
		"_rep_", replenishment_rate,
		"_ptop_", physical_topology,
		"_pdeg_", physical_degree,
		"_stop_", social_topology,
		"_sdeg_", social_degree,
		"_switch_", random_opinion_switch,
		".png",
		sep=""
	);

	if(length(filenames) == 0){ next; }

	for(file in filenames)
	{
		temp <- data.table(read.csv(file));
		Prop_vs_Risk <- rbind(Prop_vs_Risk, temp[is.na(instance)], fill=TRUE);
	}
	Prop_vs_Risk[, X:=NULL];

	Multiple_Transition <- data.table(norms=numeric(), type=character());

	for(norm in unique(Prop_vs_Risk$social_norm))
	{
		if((norm>3.1) || (norm<0)) next;
		DT <- Prop_vs_Risk[social_norm==norm][order(perceived_vaccine_risk)];
		Trans_Table <- get_transitions(DT, initial_vacc_prop);

		social_kappas <- Trans_Table[, unlist(soc_trans_risks)];
		for(i in unlist(social_kappas)){ Multiple_Transition <- rbind(Multiple_Transition, list(norm, "social")) }

		physical_kappas <- Trans_Table[, unlist(phys_trans_risks)];
		for(k in unlist(physical_kappas)){ Multiple_Transition <- rbind(Multiple_Transition, list(norm, "physical")) }
	}

	BINS <- 0;
	if(size==40000) BINS <- 0.125;
	if(size==10000) BINS <- 0.03125;

	colours <- c("#6666FF", "#FF6699");
	names(colours) <- unique(Multiple_Transition$type);

	# types <- c("social", "physical");
	types <- colours;
	names(types) <- c("#K_s", "#K_p");

	pl <- ggplot(Multiple_Transition, aes(norms, fill=type, colour=group)) +
			geom_histogram(binwidth=BINS, colour="black", alpha=0.8) +
			labs(x=expression(sigma), y="Count") +
  	  		scale_fill_discrete(labels=names(types)) +
			theme(
			  axis.line = element_line(size=1, colour = "black"),
				panel.grid.major = element_line(colour = "#d3d3d3"),
				panel.grid.minor = element_line(colour = "grey"),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.text=element_text(size=40),
				axis.title=element_text(size=40),
				legend.text=element_text(size=40),
				legend.position="right",
				legend.title=element_blank(),
				legend.direction="vertical"
			);

	ggsave(
		paste(graph_file_path, folder_name, "/", "number_of_transitions", file_name_ending, sep=""),
		plot=pl, width=20, height=4, limitsize=FALSE, dpi=50 #  height=5,
	)
}

print(Sys.time()-time_start);
