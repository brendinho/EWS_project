# Brendon Phillips
# PhD candidate
# Department of Applied Mathematics
# University of Waterloo

rm(list=ls());

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Heatmaps", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character")
Params[, c("V1", "perceived_vaccine_risk", "social_norm"):=NULL];
Params <- Params[network_structure=="smallworld" & num(initial_vacc_proportion)==0.25];
Params <- unique(Params);

start_time <- Sys.time();

trend_width <- 25;
trend_height <- 6;
text_size <- 40;
title_size <- 40;

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

	the_means <- c("soc_H_mean"); # "phys_S_mean", "phys_I_mean", "phys_R_mean", "phys_V_mean", "soc_N_mean", "soc_V_mean", 
	the_take_aways <- c(hes_fixed_values, the_means);

	DT <- data.table(matrix(ncol=length(the_take_aways)));
	names(DT) <- the_take_aways;
	for(file in filenames)
	{
		temp <- data.table(read.csv(file));
		DT <- rbind(DT, temp[is.na(instance), .SD, .SDcols=the_take_aways], fill=TRUE);
	}
	DT <- na.omit(DT);

	# means <- "phys_R_mean";

	for(means in the_means)
	{
		# Table_Here <- DT[, .SD, .SDcols=c('social_norm', 'perceived_vaccine_risk', means)];
		# DM <- make_matrix(DT, 'perceived_vaccine_risk', 'social_norm', means);

		DT.loess <- loess(get(means) ~ social_norm * perceived_vaccine_risk, data=DT);
		normgrid <- seq(min(DT$social_norm), max(DT$social_norm), 0.03125);
		riskgrid <- seq(min(DT$perceived_vaccine_risk), max(DT$perceived_vaccine_risk), 0.03125);
		data.fit <- expand.grid(social_norm=normgrid, perceived_vaccine_risk=riskgrid);
		mtrx3d <-  predict(DT.loess, newdata = data.fit);

		mtrx.melt <- melt(mtrx3d, id.vars = c('social_norm', 'perceived_vaccine_risk'), measure.vars = 'dep')
		names(mtrx.melt) <- c('norm', 'risk', 'dep')

		mtrx.melt$norm <- as.numeric(str_sub(mtrx.melt$norm, str_locate(mtrx.melt$norm, '=')[1,1] + 1))
		mtrx.melt$risk <- as.numeric(str_sub(mtrx.melt$risk, str_locate(mtrx.melt$risk, '=')[1,1] + 1))

		the_title <- list(family = "sans-serif", size=text_size, color = "black");
		the_tick <- list(family = "sans-serif", size=text_size, color = "black");

		x_ <- list(showticklabels = TRUE, tickfont = the_tick, titlefont = the_title, title='σ');
		y_ <- list(showticklabels = TRUE, tickfont = the_tick, titlefont = the_title, title='κ');
		l <- list(font = list(family="sans-serif", size=text_size, color="#000"), bgcolor="#E2E2E2", bordercolor="#FFFFFF", borderwidth=2);

		DM <- make_matrix(data.table(mtrx.melt), "norm", "risk", "dep")

		vals <- unique(scales::rescale(c(DM)))
		o <- order(vals, decreasing = FALSE)
		# cols <- scales::col_numeric("Blues", domain = NULL)(vals)
		cols <- scales::col_numeric('Reds', domain = NULL)(vals)
		colz <- setNames(data.frame(vals[o], cols[o]), NULL)

		fig <- plot_ly(
			x = rownames(DM),
			y = colnames(DM),
			z = DM,
			type = 'heatmap',
			# colorscale = colz,
			width = 700,
			height = 500,
			colorbar=list(tickfont=list(size=text_size))
		) %>% layout(xaxis=x_, yaxis=y_) %>% colorbar(size=text_size, len=1);

		graph_file_name <-  paste("Hes_Contour_", means, "_N_", the_size, "_struct_", structure, "_deg_", degree, "_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_", proportion, ".png", sep="");
		# paste(graph_file_path, folder_name, "/Hes_Contour_", means, "_N_", the_size, "_struct_", structure, "_deg_", degree, "_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_", proportion, ".png", sep="");

		print(fig)
		# stop()
		readline(prompt=sprintf("\n\n%s", graph_file_name));
		# print(graph_file_name);


		# fig <- plot_ly( mtrx.melt,
		# 	x = ~norm,
		# 	y = ~risk,
		# 	z = ~dep,
		# 	type = 'heatmap',
		# 	# colorscale = colz,
		# 	width = 700,
		# 	height = 500,
		# 	colorbar=list(tickfont=list(size=text_size))
		# ) %>% layout(xaxis=x_, yaxis=y_) %>% colorbar(size=text_size, len=1);


		# fig <- plot_ly(DT,
		# 	x = ~social_norm,
		# 	y = ~perceived_vaccine_risk,
		# 	z = ~get(means),
		# fig <- plot_ly(
		# 	x = rownames(DM),
		# 	y = colnames(DM),
		# 	z = DM,
		# 	autocontour=TRUE,
		# 	line=list(smoothing=0),
		# 	type = 'contour',
		# 	colorscale='Jet',
		# 	width = 600,
		# 	height = 500,
		# 	contours=list(start=0, end=the_size, size=2000),
		# 	contours=list(showlabels=TRUE),
		# 	colorbar=list(tickfont=list(size=text_size))
		# ) %>% layout(xaxis=x_, yaxis=y_) %>% colorbar(size=text_size, len=1);

		# fig <- filled.contour(
		# 	x = 10*1:nrow(DM),
		# 	y = 10*1:ncol(DM),
		# 	z = DM,
		# 	color.palette=terrain.colors,
		# 	plot.title = { par(cex.main=50); title(xlab=expression(kappa), ylab=expression(sigma), cex=10) },
		# 	# plot.axes = {axis(1, seq(-1, 1, by=0.25)); axis(2, seq(0, 2.5, 0.5)) },
		# 	# key.title = title(main=axis_label(means)),
		# 	key.axes=axis(4, seq(0, the_size, by=2000))
		# );
	}
}

print(Sys.time() - start_time);

graphics.off()
