# Summary Plot for the Topology
rm(list=ls());

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Location", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# source("Processing_Parameter_Table.R");
Params <- fread(parameter_file, colClasses="character")
Params[, c("V1", "risk", "norm"):=NULL];

text_width <- 40;

start_time <- Sys.time();

for(the_size in c("10000", "562500", "40000")) #
{
	DT <- unique(Params[num(beta)==1 & num(init_prop)==0.05 & num(size)==num(the_size), ])

	if(num(the_size)==10000){  DT <- DT[num(import)==0.00025]; }
	if(num(the_size)==40000){  DT <- DT[num(infec)!=0.8]; }
	if(num(the_size)==562500){ DT <- DT[num(random_switch)%in%c(0.1, 0.01,0.001, 0.0001) & phys_top=="random"] };

	DT <- unique(DT);

	for(line_index in 1:nrow(DT))
	{
		this_tuple <- DT[line_index];

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

		if(length(filenames) < 10)
		{
			# print("less than 50 observations. continue (yes/No)?");
			# if(isTRUE( readline() != "yes" )){ next; }
			next;
		}

		Table <- data.table();
		for(file in filenames)
		{
			temp <- data.table(read.csv(file));
			temp <- unique(temp[
					,
					lapply(.SD, num),
					.SDcols=c("perceived_vaccine_risk", "phys_V_mean", "phys_R_mean", "initial_vacc_proportion", "soc_V_mean", "soc_N_mean", "instance", "social_norm")
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

		if(num(the_size)==10000)
		{
			if(num(infection_prob)==0.8)
			{
				Locations <- Locations[num(norm)<=2.5]
				min_x_inset <- 0.06;
				max_x_inset <- 2;
				min_y_inset <- 0;
				max_y_inset <- 0.043;
			}
			if(num(infection_prob)==0.2)
			{
				Locations <- Locations[num(norm)<=2.5]
				min_x_inset <- 0.06;
				max_x_inset <- 2;
				min_y_inset <- 0.004;
				max_y_inset <- 0.20;
			}
		}
		if(num(the_size)==40000)
		{
			Locations <- Locations[num(norm)<=2.5]
			min_x_inset <- 0.06;
			max_x_inset <- 2.5;
			min_y_inset <- 0.005;
			max_y_inset <- 0.09;
		}
		if(num(the_size)==562500)
		{
			if(num(random_opinion_switch)==0.1)
			{
				Locations <- Locations[num(norm)<=2.5]
				min_x_inset <- 0.1;
				max_x_inset <- 2.5;
				min_y_inset <- 0.05;
				max_y_inset <- 0.45;
			}
			if(num(random_opinion_switch)==0.01)
			{
				Locations <- Locations[num(norm)<=2]
				min_x_inset <- 0.05;
				max_x_inset <- 2;
				min_y_inset <- 0.005;
				max_y_inset <- 0.051;
			}
			if(num(random_opinion_switch)==0.001)
			{
				Locations <- Locations[num(norm)<=2];
				min_x_inset <- 0.1;
				max_x_inset <- 2;
				min_y_inset <- 0.002;
				max_y_inset <- 0.019;
			}
			if(num(random_opinion_switch)==1e-4)
			{
				Locations <- Locations[num(norm)<=2];
				min_x_inset <- 0.07;
				max_x_inset <- 2;
				min_y_inset <- 0.003;
				max_y_inset <- 0.019;
			}
		}

		# keep these lines in order!
		Spline_Table <- data.table(data.frame(spline(Locations[, .SD, .SDcols=c("norm", "phys_trans")], n=nrow(Locations)*5)));
		Spline_Table <- cbind(Spline_Table, data.table(data.frame(spline(Locations[, .SD, .SDcols=c("norm", "soc_trans")], n=nrow(Locations)*5)))[, y]);
		Spline_Table <- cbind(Spline_Table, Spline_Table[, y-V2]);
		names(Spline_Table) <- c("norm", "phys_trans", "soc_trans", "differences");

		pl_inset <- ggplot(Spline_Table, aes(x=norm)) +
			geom_point(aes(y=phys_trans, colour=expression('K'['p'])), colour='red', size=1.5) +
			geom_line(data=Spline_Table, aes(norm, phys_trans), colour='red', size=0.75) +
			geom_point(aes(y=soc_trans, colour=expression('K'['s'])), colour='blue', size=1.5) +
			geom_line(data=Spline_Table, aes(norm, soc_trans), colour='blue', size=0.75) +
			geom_ribbon(Spline_Table, mapping=aes(ymin=phys_trans, ymax=soc_trans), fill='grey', alpha=0.5) +
			xlab(expression(sigma)) + ylab(expression(kappa)) +
			theme(
				panel.background=element_blank(),
				axis.line=element_line(colour = "black"),
				axis.text=element_text(size=text_width),
				axis.title=element_text(size=text_width+15)
			);

		Locations[, ("differences"):=phys_trans-soc_trans];

		pl <- ggplot(Spline_Table, aes(x=norm)) +
			xlab(expression(sigma)) + ylab(expression('K'['p']*'-'*'K'['s'])) +
			theme(
			  panel.background=element_blank(),
			  axis.line=element_line(colour = "black"),
			  axis.text=element_text(size=text_width),
			  axis.title=element_text(size=text_width+10)
			 )+
			ylim(c(min(Spline_Table$differences)-0.002, 1.05*max(Spline_Table$differences))) +
			geom_point(aes(y=differences), colour='purple', size=2) +
			geom_line(data=Spline_Table, aes(norm, differences), colour='purple', size=1.5) +
		  annotation_custom(
		    ggplotGrob(pl_inset),
		    xmin=min_x_inset, xmax=max_x_inset,
		    ymin=min_y_inset, ymax=max_y_inset
		  );

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
			".png",
			sep=""
		);

		ggsave(
			paste(graph_file_path, folder_name, "/", "IntertransitionDist", file_name_ending, sep=""),
			plot = pl,
			width=22, height=5, limitsize=FALSE, dpi=50
		); # height=4
	}
}

graphics.off();
print(Sys.time()-start_time);
