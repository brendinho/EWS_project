# Summary Plot for the Topology
rm(list=ls())

source("/home/b2philli/Dropbox/Processing/Parameter_Values.R");

folder_name <- paste("Heatmaps", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

metrics_to_plot <- c( "watts_strogatz_vacc_mean") #"vacc_avg_conn_comp_size_mean", "vacc_echo_chambers"); #
# metrics_to_plot <- c("soc_V_mean")
# "phys_V_mean", "soc_V_mean", "NV_join_count_mean", "mutual_info_mean") # "geary_c_mean") #, , "NN_join_count_mean", "VV_join_count_mean", "moran_i_mean", "ratio_sec_to_all_infs_at_the_end_mean") # ) # ) # , ); #, "join_counts", "nonvacc_echo_chambers", "mutual_info", "conn_chamb_numbers", "nonvacc_conn_comps", "vacc_conn_comps"); #, "modularity", "watts_strogatz",  "geary_c", "getis_ord", "infections", "opinion");

time_start <- Sys.time();
text_size = 40;

Params <- fread(parameter_file, colClasses="character");

for(the_size in c("10000")) # "10000", "562500"))
{
	Parameters <- Params[num(size)==the_size & num(beta)==1];

	if(the_size == "10000") Parameters <- Parameters[num(init_prop)==0.05 & num(infec)==0.2 & num(import)==0.00025];
	if(the_size == "40000") Parameters <- Parameters[num(init_prop)==0.05 & num(infec)==0.2];
	if(the_size == "562500") Parameters <- Parameters[num(init_prop)==0.05 & soc_top=="random" & num(random_switch)==0.1];

	Parameters[, c("V1", "risk", "norm", "time"):=NULL];
	Parameters <- unique(Parameters);

	# stop()

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
		random_opinion_switch <- this_tuple[ , random_switch];
		social_topology <- physical_topology;
		social_degree <- physical_degree;

		DT <- data.table();

		file_list <- Sys.glob(paste(
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

		for(filename in file_list)
		{
			temp <- data.table(read.csv(filename));
			DT <- rbind(DT, temp[is.na(instance)], fill=TRUE);
		}
		DT[, X:=NULL]
		DT <- unique(DT);
		if(nrow(DT) < 20){ next; }

		if(num(the_size) == 40000)
		{
			if(num(infection_prob) == 0.2)
			{
				DT <- unique(DT[
					num(perceived_vaccine_risk)<=0.4 &
					num(perceived_vaccine_risk)>=-1 &
					num(social_norm)>=0 &
					num(social_norm)<=2.4
				]);
			}
			if(num(infection_prob) == 0.8)
			{
				DT <- unique(DT[
					num(perceived_vaccine_risk)<=0.2 &
					num(perceived_vaccine_risk)>=-1 &
					num(social_norm)>=0 &
					num(social_norm)<=2.4
				]);
			}
		}
		if(num(the_size) == 10000)
		{
			# DT <- unique(DT[
			# 	num(perceived_vaccine_risk)<=0.2 &
			# 	num(perceived_vaccine_risk)>=-1 &
			# 	num(social_norm)>=0 &
			# 	num(social_norm)<=2.4
			# ]);
		}
		if(num(the_size) == 562500)
		{
		  DT <- unique(DT[
		    num(perceived_vaccine_risk)<=4 &
		      num(perceived_vaccine_risk)>=-2 &
		      num(social_norm)>=0 &
		      num(social_norm)<=2.5
		    ]);

			# DT <- unique(DT[
			# 	num(perceived_vaccine_risk)<=4 &
			# 	num(perceived_vaccine_risk)>=-2 &
			# 	num(social_norm)>=0 &
			# 	num(social_norm)<=4
			# ]);
		}

		for(metric in metrics_to_plot)
		{
			if(metric == "ratio_sec_to_all_infs_at_the_end_mean")
			{ # the smoothing process interp() doesn't take NAs, so set them all to zero
				set( DT, which(is.na( DT[[metric]] )), which(names(DT) == metric), 0);
			}
			if(metric == "NV_join_count_mean")
			{ # removing the highest values of the NV so that the contour map shows actual detail
				# remove the outliers
				outliers <- boxplot.stats(DT[[metric]])$out;
				set(DT,	which(DT[[metric]] %in% outliers), which(names(DT) == metric), -1);
			}
			if(metric == "geary_c_mean")
			{
				outliers <- boxplot.stats(DT[[metric]])$out;
				DT <- DT[-(which(DT[[metric]] %in% outliers))]
			}

			mains <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_V", "soc_N");
			mains_mean <- sapply(mains, function(x) sprintf("%s_mean",x));

			joins <- c("NN_join_count", "NV_join_count", "VV_join_count");
			joins_mean <- sapply(joins, function(x) sprintf("%s_mean",x));

			# not fixing up the sd'd here since not relevant
			DT[, (mains_mean):=lapply(.SD, "/", N), .SDcols=mains_mean];
			DT[, (joins_mean):=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean), .SDcols=joins_mean];

			DT.interp <- with(DT,
				interp(
					x = social_norm, y = perceived_vaccine_risk, z = get(metric),
					duplicate="median",
					xo=seq(min(DT$social_norm), max(DT$social_norm), length = 100),
					yo=seq(min(DT$perceived_vaccine_risk), max(DT$perceived_vaccine_risk), length = 100),
					extrap=FALSE, linear=FALSE
				)
			);
			# Regrouping this list to a 3-columns data.frame
			melt_x <- rep(DT.interp$x, times=length(DT.interp$y));
			melt_y <- rep(DT.interp$y, each=length(DT.interp$x));
			melt_z <- as.vector(DT.interp$z);
			DT.smooth <- data.table(na.omit(data.frame(social_norm=melt_x, perceived_vaccine_risk=melt_y, Z=melt_z)));

			# Table <- DT
			Table <- DT.smooth;
			names(Table)[3] <- "Z";

			DT_mat <- matrix(nrow=length(unique(Table$social_norm)), ncol=length(unique(Table$perceived_vaccine_risk)), NaN);
			rownames(DT_mat) <- sapply(sort(unique(Table$social_norm)), toString)
			colnames(DT_mat) <- sapply(sort(unique(Table$perceived_vaccine_risk)), toString)

			for(norm in sort(unique(Table$social_norm))){
			for(risk in sort(unique(Table$perceived_vaccine_risk)))
			{
				temp <- Table[num(social_norm)==num(norm) & num(perceived_vaccine_risk)==num(risk), Z]
				if(length(temp) == 0) next;
				if(metric %in% c("phys_V_mean", "soc_V_mean"))
				{
					if(temp < 0) temp <- 0;
					if(temp > num(size)) temp <- num(size);
				}
				else if(metric %in% c("ratio_sec_to_all_infs_at_the_end_mean"))
				{
					if(is.na(temp)) temp <- 0;
					if(temp < 0) temp <- 0;
					if(temp > 1) temp <- 1;
				}
				DT_mat[toString(norm), toString(risk)] <- temp
			}}

			file_name_ending <- paste(
				"_N_", size,
				"_vacc_", initial_vacc_prop,
				"_dur_", duration,
				"_inf_", infection_prob,
				"_imp_", importation_rate,
				"_rep_", replenishment_rate,
				"_ptop_", physical_topology,
				"_pdeg_", physical_degree,
				"_stop_", social_topology,
				"_sdeg_", social_degree,
				"_switch_", random_opinion_switch,
				sep=""
			);

			the_title <- list(family = "sans-serif", size=text_size, color = "black");
			the_tick <- list(family = "sans-serif", size=text_size, color = "black");

			x_ <- list(showticklabels = TRUE, tickfont = the_tick, titlefont = the_title,
				title='σ'
			);
			y_ <- list(showticklabels = TRUE, tickfont = the_tick, titlefont = the_title,
				title="κ"
			);
			l <- list(font = list(family="sans-serif", size=text_size, color="#000"), bgcolor="#E2E2E2", bordercolor="#FFFFFF", borderwidth=2);

			p <- plot_ly(
				x=rownames(DT_mat),
				y=colnames(DT_mat),
				z=DT_mat,
				autocontour=TRUE,
				line=list(smoothing=0),
				type='contour',
				colorscale='Jet',
				contours=list(start=0, end=size, size=10),
				contours=list(showlabels=TRUE),
				colorbar=list(tickfont=list(size=text_size)),
			) %>% layout(xaxis=x_, yaxis=y_) %>% colorbar(size=text_size, len=1) # legend = list(font = list(size(text_size)))

			print(sprintf("Contour_%s%s", metric, file_name_ending));

			# orca(p, paste(graph_file_path, folder_name, "/", "Contour_", metric, file_name_ending, ".jpg", sep=""))

			print(p)

			stop(Sys.time() - time_start);
		}
	}
}

# a <- list(
#   autotick = FALSE,
#   ticks = "outside",
#   tick0 = 0,
#   dtick = 0.25,
#   ticklen = 5,
#   tickwidth = 2,
#   tickcolor = toRGB("blue")
# )
# norm_range <- seq(0, 2, by=0.25)
# risk_range <- seq(-1, 1, by=0.25)
# # p <- plot_ly(z=DT_mat, type="contour") %>% layout(xaxis = a, yaxis = a)
# p <- plot_ly(x=~norm_range, z=~risk_range, type="contour") %>% add_markers() %>%
#   add_markers(y = ~rev(s)) %>%
#   layout(xaxis = a, yaxis = a)
#
# contour.mat <- ifelse(DT_mat < 20000, 0, DT_mat)
# contour.mat2 <- ifelse(DT_mat < 30000, 0, DT_mat)
#
# filled.contour(DT_mat, color = terrain.colors,
#            plot.axes = contour(contour.mat, levels = 1,
#                                drawlabels = FALSE, axes = FALSE,
#                                frame.plot = FFALSE, add = TRUE) +
# 							   contour(contour.mat2, levels = 1,
# 				                                   drawlabels = FALSE, axes = FALSE,
# 				                                   frame.plot = FFALSE, add = TRUE))
