# # Summary Plot for the Topology
# rm(list=ls())
#
# user_name <- Sys.info()[8][[1]];
# source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));
#
# folder_name <- paste("Change_point_application", sep="");
# dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);
#
# # source("Processing_Parameter_Table.R");
# Params <- fread(parameter_file, colClasses="character");
#
# start_time <- Sys.time();
#
# metrics_to_plot <- c();
#
# for(the_size in c(10000, 40000)) #
# {
#   Parameters <- Params[num(beta)==1 & num(init_prop)==0.05 & num(size)==num(the_size)];
#   Parameters[, c("V1", "risk"):=NULL];
#
# 	if(the_size == 40000)
# 	{
# 		Parameters <- Parameters[num(norm)==0 & init_prop==0.05 & infec==0.2];
# 		metrics_to_plot <-c("NV_join_count_mean", "mutual_info_mean", "VV_join_count_mean", "NN_join_count_mean", "moran_i_mean", "geary_c_mean")
# 		# "phys_S_mean", "phys_I_mean", "phys_R_mean", "phys_V_mean", "soc_N_mean", "soc_V_mean");
# 	}
# 	if(the_size == 10000) {
# 	  Parameters <- Parameters[num(norm)==0];
# 	   metrics_to_plot <- c("NV_join_count_mean", "mutual_info_mean", "VV_join_count_mean", "NN_join_count_mean", "moran_i_mean", "geary_c_mean");
# 	}
# 	Parameters <- unique(Parameters);
#
# 	for(line_index in 1:nrow(Parameters))
# 	{
# 		this_tuple <- Parameters[line_index];
#
# 		beta <- this_tuple[, beta];
# 		initial_vacc_prop <- this_tuple[, init_prop];
# 		size <- this_tuple[, size];
# 		duration <- this_tuple[, duration];
# 		infection_prob <- this_tuple[, infec];
# 		importation_rate <- this_tuple[, import];
# 		replenishment_rate <-  this_tuple[, birth_death];
# 		physical_topology <- this_tuple[, phys_top];
# 		physical_degree <- this_tuple[, phys_deg];
# 		social_topology <- this_tuple[, soc_top];
# 		social_degree <- this_tuple[, soc_deg];
# 		random_opinion_switch <- this_tuple[, random_switch];
# 		social_norm <- this_tuple[, norm];
# 		num_size <- as.numeric(size);
#
# 	    Prop_vs_Risk <- data.table();
#
# 		filenames <- Sys.glob(paste(
# 			CSV_file_path,
# 			"Summary",
# 			"_N_", size,
# 			"_dur_", duration,
# 			"_beta_", beta,
# 			"_vaccprop_", initial_vacc_prop,
# 			"_risk_", "*",
# 			"_inf_", infection_prob,
# 			"_imp_", importation_rate,
# 			"_rep_", replenishment_rate,
# 			"_ptop_", physical_topology,
# 			"_pdeg_", physical_degree,
# 			"_stop_", social_topology,
# 			"_sdeg_", social_degree,
# 			"_norm_", social_norm,
# 			"_switch_", random_opinion_switch,
# 			".csv",
# 			sep=""
# 		));
#
# 		if(length(filenames) == 0){ next; }
#
# 		for(file in filenames)
# 		{
# 			temp <- data.table(read.csv(file));
# 			Prop_vs_Risk <- rbind(Prop_vs_Risk, temp[is.na(instance), ], fill=TRUE);
# 		}
#
# 		Trans_Table <- get_transitions(Prop_vs_Risk, initial_vacc_prop);
#
# 		file_name_ending <- paste(
# 			"_N_", size,
# 			"_dur_", duration,
# 			"_beta_", beta,
# 			"_vaccprop_", initial_vacc_prop,
# 			"_inf_", infection_prob,
# 			"_imp_", importation_rate,
# 			"_rep_", replenishment_rate,
# 			"_ptop_", physical_topology,
# 			"_pdeg_", physical_degree,
# 			"_stop_", social_topology,
# 			"_sdeg_", social_degree,
# 			"_norm_", social_norm,
# 			"_switch_", random_opinion_switch,
# 			"_short.png",
# 			sep=""
# 		);
#
# 		Proportion 	<- Prop_vs_Risk[initial_vacc_proportion==initial_vacc_prop, ];
#
# 		# if there are less than 5 points to be plotted, then cancel this plot and skip to the next parameter set
# 		if(nrow(Proportion)<5){ next; }
#
# 		for(metric in metrics_to_plot)
# 		{
# 			Time_Series_Table <- data.table(matrix(ncol=3, nrow=0));
# 			names(Time_Series_Table) <- c("risk_end", "norm_change", "norm_p");
#
# 			Series <- Proportion[, .SD, .SDcols=c("perceived_vaccine_risk", metric)];
#
# 			the_sets <- vector(mode="list", length=4);
#
# 			label_local <- function(x)
# 			{
# 				if(x == "NV_join_count_mean"){ temp <- expression('B{<N,V'['s']*'>}'); }
# 				if(x == "VV_join_count_mean"){ temp <- expression('B{<V'['s']*',V'['s']*'>}'); }
# 				if(x == "NN_join_count_mean"){ temp <- expression('B{<N,N>}'); }
# 				if(x == "mutual_info_mean"){ temp <- expression('B{<I>}'); }
# 				if(x == "moran_i_mean"){ temp <- expression('B{<M>}'); }
# 				if(x == "geary_c_mean"){ temp <- expression('B{<C>}'); }
# 				# return( bquote(.(axis_label(Test))*'{'*.(axis_label(x))*'}') );
# 				return(temp);
# 			}
#
# 			# very fucking important line
# 			Series <- Series[order(perceived_vaccine_risk)];
#
# 			Test_Vector <- c("buishand_r", "lanzante", "pettitt", "snh");
# 			# colours <- brewer.pal(length(Tesst_Vector), name="Dark2")
# 			colours <- c('red', 'blue', 'green', 'purple');
# 			names(colours) <- Test_Vector;
#
# 			text_size <- 45;
#
# 			pl <- ggplot() +
# 				theme(
# 					panel.background = element_blank(),
# 					axis.title.y = element_text(size=text_size),
# 					axis.title.x = element_text(size=text_size),
# 					axis.text = element_text(size=text_size),
# 					axis.line = element_line(color="black", size = 1)
# 				) +
# 				geom_hline(yintercept=Trans_Table[, soc_intersec_risk], size=1.5, linetype='dashed', col="blue") +
# 				labs(x="number of terms", y=label_local(metric)) +
# 				ylim(c(-1,0.03));
#
# 			for(jk in 1:length(Test_Vector))
# 			{
# 				the_test <- Test_Vector[jk];
#
# 				norm_significant   <- data.table(risk=numeric(), number=numeric());
# 				norm_insignificant <- data.table(risk=numeric(), number=numeric());
#
# 				for(i in 3:nrow(Series))
# 				{
# 					# if we've already reached the transition or are beyong it, the time for predictions has passed
# 					if(Series[i, perceived_vaccine_risk] > Trans_Table[, soc_intersec_risk]) break;
#
# 					norm_set <- data.table();
# 					Tab_Here <- Series[1:i];
#
# 					# writen explicitly to test increasing lengths of the kappa series
# 				   if(the_test == "buishand_r"){ norm_set <- br.test(Tab_Here[[metric]], m=50000); }
# 				   if(the_test == "lanzante")  { norm_set <- lanzante.test(Tab_Here[[metric]]); }
# 				   if(the_test == "pettitt")   { norm_set <- pettitt.test(Tab_Here[[metric]]); }
# 				   if(the_test == "snh")       { norm_set <- snh.test(Tab_Here[[metric]]); }
#
# 					if(isTRUE( norm_set$p.value < 0.05 )){
# 						norm_significant <- rbind(norm_significant, list(Tab_Here[norm_set$estimate[[1]], perceived_vaccine_risk], i))
# 					} else {
# 						norm_insignificant <- rbind(norm_insignificant, list(Tab_Here[norm_set$estimate[[1]], perceived_vaccine_risk], i));
# 					}
# 				}
#
# 				pl <- pl +
# 					geom_point(data=norm_insignificant, aes(x=number, y=risk), color=colours[jk], size=5, shape=19) +
# 					geom_line( data=norm_significant,   aes(x=number, y=risk, colour=the_test), size=2, color=colours[jk])
# 			}
#
# 			pl <- pl +
# 				scale_colour_manual(name="Test used", values=c(myline1=colours[1], myline2=colours[2], myline3=colours[3], myline4=colours[4]), labels=as.vector(sapply(Test_Vector, function(x) test_name(x, "long")))) +
# 				theme(legend.position="bottom", legend.justification="center", legend.direction="vertical");
#
# 			ggsave(
# 				paste(graph_file_path, folder_name, "/", "Changes_", metric, file_name_ending, sep=""),
# 				plot=pl, width=10, height=5, dpi=50, limitsize=FALSE
# 			);
# 			dev.off()
# 		}
# 	}
# }

legend_plot_name <- paste(graph_file_path, folder_name, "/", "Changes_legend", file_name_ending, sep="");
png(
	legend_plot_name,
	width=100, height=100, unit="cm", res=50
);
par(mar=c(1,1,1,1));
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1);
legend("bottom", horiz=TRUE, xpd=TRUE, col=unname(colours), pch=15, pt.cex=15, cex=8, bty='n', legend=as.vector(sapply(names(colours), function(x) test_name(x, "old"))));
dev.off();
system(sprintf("convert %s -trim %s", legend_plot_name, legend_plot_name));


print(Sys.time()-start_time);
