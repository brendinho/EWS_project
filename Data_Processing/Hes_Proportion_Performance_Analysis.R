# Brendon Phillips
# PhD candidate
# Department of Applied Mathematics
# University of Waterloo

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Lead_Time", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character");
Params[, c("V1", "perceived_vaccine_risk", "social_norm", "proportion_of_nodes"):=NULL];
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

	writeLines("\n");
	writeLines(paste("\timportation rate: ", importation_rate, sep=""));
	writeLines(paste("\tnetwork size: ", the_size, sep=""));
	writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
	writeLines(paste("\tstructure: ", structure, sep=""));
	writeLines(paste("\tdegree: ", degree, sep=""));
	writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));

	# ALL_DATA <- data.table();
	# metrics <- character();
	#
	# Real_Transitions <- data.table(matrix(ncol=3));
	# for(the_proportion in c(1, 0.8, 0.6))
	# {
	# 	DT <- data.table();
	# 	filenames <- Sys.glob(paste(CSV_file_path, "Hes_N_", the_size, "_struct_", structure, "_deg_", degree,  "_risk_*",
	# 						"_infec_", infection_prob, "_norm_*",  "_imp_", importation_rate, "_init_", initial_vacc_prop,
	# 						"_switch_", random_opinion_switch, "_prop_", the_proportion, "_summary.csv", sep=""));
	#
	# 	for(file in filenames)
	# 	{
	# 		temp <- data.table(read.csv(file));
	# 		DT <- rbind(DT, temp[is.na(instance)], fill=TRUE);
	# 	}
	# 	DT[, c("X", "instance"):=NULL];
	#
	# 	metrics <- setdiff(names(DT), c(hes_fixed_values, "instance"));
	# 	metrics <- filter_these(metrics, no=c("_sd", "ratio_", "soc", "phys"));
	#
	# 	if(the_proportion == 1)
	# 	{
	# 		names(Real_Transitions) <- c("norm", "phys_trans", "soc_trans");
	# 		for(norm_index in seq_along(DT[, unique(social_norm)]))
	# 		{
	# 			soc_norm <- DT[, sort(unique(social_norm))][norm_index];
	# 			transes <- get_transitions(DT[social_norm==soc_norm][order(perceived_vaccine_risk)], initial_vacc_prop);
	# 			Real_Transitions <- rbind( Real_Transitions, list(soc_norm, transes$phys_intersec_risk, transes$soc_intersec_risk) );
	# 		}
	# 		Real_Transitions <- Real_Transitions[rowSums(is.na(Real_Transitions)) != ncol(Real_Transitions)][order(norm)];
	# 		Real_Transitions[, proportion:=the_proportion];
	# 	}

		CHANGE_TEST <- "snh";

	# 	col_names <- c("norm", "phys_trans", "soc_trans", "proportion", names.change(metrics));
	#
	# 	# get all the rest of the change times
	# 	writeLines(sprintf("\ttest: %s, proportion: %f", CHANGE_TEST, the_proportion));
	# 	Change_Predictions <- data.table(matrix(ncol=length(col_names)));
	# 	names(Change_Predictions) <- col_names;
	# 	for(norm_index in seq_along(DT[, unique(social_norm)]))
	# 	{
	# 		soc_norm <- DT[, sort(unique(social_norm))][norm_index];
	# 		Prop_Table_Here <- DT[social_norm==soc_norm][order(perceived_vaccine_risk)];
	# 		transitions <- get_transitions(Prop_Table_Here, initial_vacc_prop);
	# 		change_points <- c();
	# 		for(metric in metrics) change_points[[metric]] <- change_point_test(Prop_Table_Here, metric, CHANGE_TEST)$kappa;
	# 		Change_Predictions <- rbind(
	# 			Change_Predictions,
	# 			as.list(c(
	# 				soc_norm,
	# 				transitions$phys_intersec_risk,
	# 				transitions$soc_intersec_risk,
	# 				the_proportion,
	# 				data.table(change_points)[, change_points]
	# 			))
	# 		);
	# 	}
	# 	Change_Predictions <- Change_Predictions[!is.na(norm)];
	#
	# 	if(the_proportion != 1)
	# 	{
	# 		for(the_norm in DT[, unique(social_norm)])
	# 		{
	# 			Change_Predictions[
	# 				norm==the_norm,
	# 				soc_trans:=suppressWarnings(minn(c(Change_Predictions[norm==the_norm, soc_trans], Real_Transitions[norm==the_norm, soc_trans])))
	# 			];
	# 			Change_Predictions[
	# 				norm==the_norm,
	# 				phys_trans:=suppressWarnings(minn(c(Change_Predictions[norm==the_norm, phys_trans], Real_Transitions[norm==the_norm, phys_trans])))
	# 			];
	# 		}
	# 	}
	# 	All_Lead_Times <- Change_Predictions[, min(soc_trans, phys_trans)] - Change_Predictions[, .SD, .SDcols=names.change(metrics)];
	# 	All_Lead_Times <- cbind(Change_Predictions[, .SD, .SDcols=c("norm", "soc_trans", "phys_trans", "proportion")], All_Lead_Times);
	# 	names(All_Lead_Times)[1:3] <- c("norm", "soc_trans", "phys_trans", "proportion");
	#
	# 	ALL_DATA <- rbind(ALL_DATA, All_Lead_Times, fill=TRUE);
	# }
	# write.csv(ALL_DATA, "all_the_data_for_proportion_analysis.csv");

	ALL_DATA <- data.table(read.csv("all_the_data_for_proportion_analysis.csv"));
	ALL_DATA[, c("X"):=NULL];

	metrics <- filter_these(names(ALL_DATA), yes=c(".change"))

	Performance <- data.table(matrix(nrow=length(metrics)*length(unique(ALL_DATA$proportion)), ncol=9, 0.));
	names(Performance) <- c("test", "proportion", "metric", "min_diff", "mean_diff", "max_diff", "min_lead", "mean_lead", "max_lead");
	Performance[, metric:=as.character(metric)]; Performance[, test:=as.character(test)];

	file_index <- 1;
	for(the_proportion in c(0.6, 0.8, 1)){
	for(the_metric in unique(metrics)){

		lead_times_right_now <- ALL_DATA[proportion==the_proportion, .SD, .SDcols=c("norm", the_metric)][order(norm)];
		the_merger <- merge(lead_times_right_now, ALL_DATA[proportion==1, .SD, .SDcols=c("norm", the_metric)], by=c("norm"));
		names(the_merger) <- c("norm", "sampled", "real");

		set(
			Performance, file_index, names(Performance),
			list(
				CHANGE_TEST, the_proportion, the_metric,
				if(the_proportion != 1) 100*the_merger[, minn(abs(real-sampled))/mean(real)] else NA,
				if(the_proportion != 1) 100*the_merger[, mean(abs(real-sampled))/mean(real)] else NA,
				if(the_proportion != 1) 100*the_merger[, maxx(abs(real-sampled))/mean(real)] else NA,
				minn(lead_times_right_now[[the_metric]]),
				mean(lead_times_right_now[[the_metric]], na.rm=TRUE),
				maxx(lead_times_right_now[[the_metric]])
			)
		);
		file_index <- file_index + 1;
		# if(the_proportion != 1) 100*the_merger[, minn(abs(real-sampled))]/mean_lead else NA,
		# if(the_proportion != 1) 100*the_merger[, mean(abs(real-sampled))]/mean_lead else NA,
		# if(the_proportion != 1) 100*the_merger[, maxx(abs(real-sampled))]/mean_lead else NA,

		# stop()
	}}
	write.csv(Performance, "Proportion_Performance_Comparison.csv");

	# Performance[proportion==1 & metric%in%filter_these(Performance$metric, yes=c("watts")),]

	stop(Sys.time() - start_time)

}

print(Sys.time() - start_time);

graphics.off()
