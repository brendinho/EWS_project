# Brendon  Phillips
# PhD candidate
# Department of Applied Mathematics
# Faculty of Mathematics
# University of Waterloo
# rm(list=ls())

# options(width=Sys.getenv("COLUMNS"));
options(bitmapType="cairo");

# source("find_colours.R")

if(!require(pacman)) install.packages('pacman')
pacman::p_load(
				akima,
				# Cairo,
				colorRamps,
				comprehenr,
				# ContourFunctions,
				# colorspace,
				cowplot,
				data.table,
				# doParallel,
				dplyr,
				e1071,
				# foreach,
				# graphics,
				# grDevices,
				gridExtra,
				ggformula,
				ggplot2,
				ggpol,
				# ggpubr,
				# ggthemes,
				gtools,
				# janitor,
				# laGP,
				latex2exp,
				magrittr,
				# Matching,
				# magrittr,
				# moments,
				# multipanelfigure,
				orca,
				plotly,
				plyr,
				png,
				# pracma,
				# processx,
				quantreg,
				reshape2,
				RColorBrewer,
				scales,
				# showtext,
				splines,
				stringr,
				svMisc,
				# threadpool,
				tidyverse,
				tools,
				trend #,
				# wrapr,
				# viridis,
				# xtable
);

x11 = function (...) grDevices::x11(...,type='cairo');

Sys.setenv('MAPBOX_TOKEN'='pk.eyJ1IjoiYjJwaGlsbGkiLCJhIjoiY2s1NjZ5bzJlMDI4dzNucGtkdjA2NjNydiJ9.CNeuv4-hVP3cRnfXtvP3Ew');
Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/anaconda3/bin", sep = .Platform$path.sep))

# setHook(packageEvent("grDevices", "onLoad"),
# function(...) grDevices::X11.options(type='cairo'))
# options(device='x11')

# font_add_google("Lobster")

user_name <- Sys.info()[8][[1]];

Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/anaconda3/bin", sep = .Platform$path.sep))

# home_folder <- sprintf("/media/%s/Simulations", user_name);
home_folder <- sprintf("/home/%s/Dropbox", user_name);
parameter_file <- sprintf("%s/CSV/params.csv", home_folder, user_name);
stored_data_path <- sprintf("/media/%s/Simulations", user_name);

threshold_for_declaring_high_proportion <- 0.9;
threshold_for_declaring_low_proportion <- 0.1;

number_of_time_steps_averaged <- 500;

vaxx_allowed_vector <- c(TRUE);
use_this_many_entries_for_ratio_sec_new <- 500;

# values that do not change between the time steps in each individual trial or not to be averaged
fixed_values <- c(
					"N", "instance", "social_topology", "social_degree", "physical_topology",
					"physical_degree", "infec_prob", "social_norm", "importation",
					"initial_vacc_proportion", "random_opinion_switch", "perceived_vaccine_risk",
					"beta_parameter"
				);

# the parameter values tha we want to reisnert to the data that were excluded above
reinserted_values <- setdiff(fixed_values, c("instance"));

user_name <- Sys.info()[8][[1]];

CSV_file_path <-   sprintf("%s/CSV/", home_folder);
graph_file_path <- sprintf("/home/%s/Dropbox/Apps/Overleaf/phd_thesis/Thesis_Graphs/", user_name);
# graph_file_path <- sprintf("/home/%s/Dropbox/Apps/Overleaf/warning_signals_paper/Graphs/", user_name);

maxx <- function(x) max(x, na.rm=T);
minn <- function(x) min(x, na.rm=T);

mode <- "long";

axissize <- 2;

high <- 10000;
wide <- 8000;

labelsize <- 3;
legend_linewidth <- 3;
linewidth <- 4;

paragraph <- c(4.25,5.5,4,1);
physical_transition_line_colour <- "black";
plot_border_width <- 3;

resolut <- 300;

social_transition_line_colour <- "black";

transition_line_type <- 3; # dashed
transition_linewidth <- 3;

num <- function(x) { return( as.numeric(x) ); }

composite <- function(f,g) function(...) f(g(...));

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot)
{
	tmp <- ggplot_gtable(ggplot_build(a.gplot));
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box");
	legend <- tmp$grobs[[leg]];
	return(legend)
}

test_name <- function(x, y="old")
{
    if(x == "buishand_r") return(if(y=="old") "Buishand R" else "Buishand Range Test");
    if(x == "pettitt") return(if(y=="old") "Pettitt" else "Pettitt Test");
    if(x == "lanzante") return(if(y=="old") "Lanzante" else "Lanzante Test");
    if(x == "snh") return(if(y=="old") "SNHT" else "Standard Normal Homogeneity Test");
}

axis_label <- function(x)
{
	if(x == "plain"){ return("Output"); }
	if(x == "number"){ return("Groups"); }
	if(x == "Dyn"){ return("S/P");}

	if(x == "buishand_r") return("BR");
	if(x == "pettitt") return("Pet");
	if(x == "lanzante") return("Lan");
	if(x == "snh") return("SNHT");

	if(x == "join_counts"){ return("Join counts"); }
	if(x %in% c("NN_join_count", "NN")){ return("[N,N]"); }
	# if(x %in% c("NV_join_count", "NV")){ return( expression('[N,V'*['s']*']') ); }
	if(x == "NN_join_count_mean"){ return("[N,N]"); }
	if(x == "NV_join_count_mean"){ return("[N,V]"); }
	if(x %in% c("VV_join_count_mean", "VV", "VV_join_count")){ return("[V,V]"); }

	if(x == "conn_comp_size"){ return(expression('<|Z'['*']*'|>')); }
	if(x == "conn_comp_number"){ return(expression('<#Z'['*']*'>')); }
	if(x == "echo_chamber_size"){ return(expression('<|J'['*']*'|>')); }
	if(x == "echo_chamber_number"){ return(expression('<#J'['*']*'>')); }
	if(x == "modularity"){ return(expression('<Q'['*']*'>')); }
	if(x == "opinion_change"){ return(expression('<'*Theta ['*']*'>')); }
	if(x == "prob_sick"){ return(expression('<'*Gamma ['*']*'>')); }
	if(x == "watts_strogatz"){ return(expression('<C'['*']*'>')); }

    if(x %in% c("vacc.modularity.change", "vacc_modularity")){ return(expression('<Q'['V']*'>')); }
    if(x %in% c("nonvacc.modularity.change", "nonvacc_modularity")){ return(expression('<Q'['N']*'>')); }
    if(x %in% c("total.modularity.change", "total_modularity")){ return(expression('<Q'[Sigma]*'>')); }

	if(x %in% c("nonvacc.avg.chamber.size.change", "nonvacc_avg_chamber_size")){ return(expression('<|J'['N']*'|>')); }
	if(x == "nonvacc.min.chamber.size.change"){ return(expression('min(|J'['N']*'|)')); }
	if(x == "nonvacc.max.chamber.size.change"){ return(expression('max(|J'['N']*'|)')); }
    if(x == "nonvacc_echo_chamber_size"){ return(expression('|J'['N']*'|')); }

	if(x %in% c("vacc.avg.chamber.size.change", "vacc_avg_chamber_size")){ return(expression('<|J'['V']*'|>')); }
	if(x == "vacc.min.chamber.size.change"){ return(expression('min(|J'['V']*'|)')); }
	if(x == "vacc.max.chamber.size.change"){ return(expression('max(|J'['V']*'|)')); }
    if(x == "vacc_echo_chamber_size"){ return(expression('|J'['V']*'|')); }

    if(x %in% c("nonvacc.number.chambers.change", "nonvacc_number_chambers")){ return(expression('<#J'['N']*'>')); }
    if(x %in% c("vacc.number.chambers.change", "vacc_number_chambers")){ return(expression('<#J'['V']*'>')); }

    if(x %in% c("vacc.avg.conn.comp.size.change", "vacc_avg_conn_comp_size")){ return(expression('<|Z'['V']*'|>')); }
    if(x %in% c("vacc.min.conn.comp.size.change")){ return(expression('min(|Z'['V']*'|)')); }
    if(x == "vacc.max.conn.comp.size.change"){ return(expression('max(|Z'['V']*'|)')); }

    if(x %in% c("nonvacc.avg.conn.comp.size.change", "nonvacc_avg_conn_comp_size")){ return(expression('<|Z'['N']*'|>')); }
    if(x == "nonvacc.min.conn.comp.size.change"){ return(expression('min(|Z'['N']*'|)')); }
    if(x == "nonvacc.max.conn.comp.size.change"){ return(expression('max(|Z'['N']*'|)')); }

    if(x %in% c("nonvacc.number.conn.comps.change", "nonvacc_number_conn_comps")){ return(expression('<#Z'['N']*'>')); }
    if(x %in% c("vacc.number.conn.comps.change", "vacc_number_conn_comps")){ return(expression('<#Z'['V']*'>')); }

    if(x == "VV.join.count.change"){ return( parse(text='"<N,N>"') ); }
    if(x == "NV.join.count.change"){ return( parse(text='"<N,V"[s]*">"') ); }
    if(x == "NN.join.count.change"){ return( parse(text='"<V"[s]*",V"[s]*">"') ); }

	if(x %in% c("phys_S", "phys_S_mean", "phys.S.change")){ return("<S>"); }
	if(x %in% c("phys_I", "phys_I_mean", "phys.I.change")){ return("<I>"); }
	if(x %in% c("phys_R", "phys_R_mean", "phys.R.change")){ return("<R>"); }
	if(x %in% c("phys_V", "phys_V_mean", "phys.V.change")){ return(parse(text='"<V"[p]*">"')); }
	if(x %in% c("soc_N",  "soc_N_mean" , "soc.N.change")){ return("<N>"); }
	if(x %in% c("soc_V",  "soc_V_mean", "soc.V.change")){ return(parse(text='"<V"[p]*">"')); }

	if(x %in% c("moran")){ return("Moran's I"); }
	if(x == "moran_i"){ return("Moran's I"); }
	if(x == "geary"){ return("Geary's C"); }
	if(x == "geary_c"){ return("Geary's C"); }
	if(x == "mutual"){ return("Mutual Inf."); }
	if(x %in% c("mutual_info", "mutual_info_mean")){ return("<M>"); }

	if(x %in% c("moran.i.change")){ return(expression('<I>')); }

	if(x == "geary.c.change"){ return(expression('<C>')); }
	if(x == "mutual.info.change"){ return(expression('<M>')); }

	if(x == "V_Conn_Size"){ return(quote('|Z'['V']*'|')); }
	if(x == "NV_Conn_Size"){ return(quote('|Z'['N']*'|')); }
    if(x == "Conn_Size"){ return(quote('|Z'['*']*'|')); }
	if(x == "Conn_Echo_Count"){ return(quote('<# '['*']*'>')); }

	if(x == "V_Echo_Size"){ return(quote('|J'['V']*'|')); }
	if(x == "NV_Echo_Size"){ return(quote('|J'['N']*'|')); }
    if(x == "Echo_Size"){ return(quote('|J'['*']*'|')) }

	if(x == "Prob_Sick"){ return(quote('<'*Gamma['*']*'>')); }
	if(x == "WattsGCC"){ return(quote('<C'['*']*'>')); }
	if(x == "Modu"){ return(quote('<Q'['*']*'>'));  }
    if(x == "Senti"){ return(quote('<'*Theta['*']*'>')); }

    if(x %in% c("prob.sick.nonvacc.change", "prob_sick_nonvacc")){ return(expression('<'*Gamma['N']*'>')); }
    if(x %in% c("prob.sick.vacc.change", "prob_sick_vacc")){ return(expression('<'*Gamma['V']*'>')); }

    if(x %in% c("watts.strogatz.all.change", "watts_strogatz_all_mean", "watts_strogatz_all")){ return(expression('<C'[Sigma]*'>')); }
    if(x %in% c("watts.strogatz.nonvacc.change", "watts_strogatz_nonvacc_mean", "watts_strogatz_nonvacc")){ return(expression('<C'['N']*'>')); }
    if(x %in% c("watts.strogatz.vacc.change", "watts_strogatz_vacc_mean", "watts_strogatz_vacc")){ return(expression('<C'['V']*'>')); }

    if(x %in% c("nonvacc.opinion.change.change", "nonvacc_opinion_change")){ return(expression('<'*Theta['N']*'>')); }
    if(x %in% c("vacc.opinion.change.change", "vacc_opinion_change")){ return(expression('<'*Theta['V']*'>')); }
    if(x %in% c("total.opinion.change.change", "total_opinion_change")){ return(expression('<'*Theta[Sigma]*'>')); }

	else return(x);
}

plot_these <- function(metric)
{
	if(metric == "mutual_info") 					{ plot_these <- c("mutual_info"); }
	else if(metric == "plain")						{ plot_these <- c("phys_R", "phys_V", "soc_N", "soc_V"); }
	else if(metric == "join_counts") 	 			{ plot_these <- c("NV_join_count", "NN_join_count", "VV_join_count"); }
	else if(metric == "NV_join_count")				{ plot_these <- c("NV_join_count"); }
	else if(metric == "moran_i") 					{ plot_these <- c("moran_i"); } #, "moran_i_calc"); }
	else if(metric == "geary_c")					{ plot_these <- c("geary_c"); } #, "geary_c_calc"); }
	else if(metric == "getis_ord")					{ plot_these <- c("getis_ord", "getis_ord_calc"); }
	else if(metric == "vacc_conn_comp_size") 		{ plot_these <- c("vacc_min_conn_comp_size", "vacc_max_conn_comp_size", "vacc_avg_conn_comp_size"); }
    else if(metric == "number")                     { plot_these <- c("vacc_number_conn_comps", "nonvacc_number_conn_comps", "vacc_number_chambers", "nonvacc_number_chambers"); }
	else if(metric == "nonvacc_conn_comp_size") 	{ plot_these <- c("nonvacc_min_conn_comp_size", "nonvacc_max_conn_comp_size", "nonvacc_avg_conn_comp_size"); }
	else if(metric == "conn_comp_number") 			{ plot_these <- c("vacc_number_conn_comps", "nonvacc_number_conn_comps"); }
	else if(metric == "vacc_echo_chamber_size") 	{ plot_these <- c("vacc_min_chamber_size", "vacc_max_chamber_size", "vacc_avg_chamber_size");  }
	else if(metric == "nonvacc_echo_chamber_size") 	{ plot_these <- c("nonvacc_min_chamber_size", "nonvacc_max_chamber_size", "nonvacc_avg_chamber_size");  }
	else if(metric == "echo_chamber_number")		{ plot_these <- c("vacc_number_chambers", "nonvacc_number_chambers"); }
	else if(metric == "echo_chamber_size")			{ plot_these <- c("vacc_avg_chamber_size", "nonvacc_avg_chamber_size"); }
	else if(metric == "conn_comp_size")				{ plot_these <- c("vacc_avg_conn_comp_size", "nonvacc_avg_conn_comp_size"); }
	else if(metric == "watts_strogatz")				{ plot_these <- c("watts_strogatz_vacc", "watts_strogatz_nonvacc", "watts_strogatz_all"); }
	else if(metric == "modularity")					{ plot_these <- c("vacc_modularity", "nonvacc_modularity", "total_modularity"); }
	else if(metric == "opinion_change")				{ plot_these <- c("nonvacc_opinion_change", "vacc_opinion_change", "total_opinion_change") }
	else if(metric == "prob_sick")					{ plot_these <- c("prob_sick_nonvacc", "prob_sick_vacc"); }
	# else if(metric == "vacc_echo_chamber_size")	{ plot_these <- c("vacc_min_chamber_size", "vacc_max_chamber_size", "vacc_avg_chamber_size"); }
	# else if(metric == "nonvacc_echo_chamber_size"){ plot_these <- c("nonvacc_min_chamber_size", "nonvacc_max_chamber_size", "nonvacc_avg_chamber_size"); }
	return(plot_these);
}

names.change <- function(name_vector)
{  # takes the names of variable for which we need lead times, and returns the names of the lead time table
	unlist(lapply(
		name_vector,
		# function(z) str_c(c(head(strsplit(z, '_')[[1]], -1), "change"), collapse=".")
		function(z) str_c(c(strsplit(strsplit(z, "_mean")[[1]], '_')[[1]], "change"), collapse=".")
	));
}

filter_these <- function(the_list, yes=c(), no=c())
{	# function to filter out any string from the_list containing any substring from not_these
	foo <- the_list;
	for(bad_thing in no){ if(length(bad_thing) != 0){ foo <- foo[!grepl(bad_thing, foo)]; }}
	for(good_thing in yes)
	{
		if(length(good_thing) == 0) next;
		foo <- foo[grepl(good_thing, foo)];
		if(length(foo) == 0) break
	}
	return(foo);
}

get_transitions <- function(Prop_vs_Risk_table, initial_vacc_prop)
{
	# since we're plotting the sequence with lines, we need to make sure that the perceived_vaccine_risks are in numerical order
	Prop_Table <- unique(Prop_vs_Risk_table[
			initial_vacc_proportion==initial_vacc_prop,
			lapply(.SD, num),
			.SDcols=c("perceived_vaccine_risk", "phys_V_mean", "phys_R_mean", "initial_vacc_proportion", "soc_V_mean", "soc_N_mean")
		][order(perceived_vaccine_risk)]);

	Trans_Table <- data.table(
		is_phys_trans = c(FALSE),
		is_soc_trans = c(FALSE),
		soc_intersec_risk = c(NaN),
		phys_intersec_risk = c(NaN),
		soc_intersec_prop = c(NaN),
		phys_intersec_prop = c(NaN),
		soc_trans_risks = vector("list", 1L),
		phys_trans_risks = vector("list", 1L)
	);

	if( length(unique(Prop_Table[, soc_V_mean-soc_N_mean] < 0)) != 1 )
	{
		Trans_Table[, is_soc_trans:= TRUE];

		# some of the curves have more than one transition, so we're taking the min value of the twp or three that come out
		all_the_social_transition_rows <- Prop_Table[, which(diff(soc_V_mean>soc_N_mean)!=0)];
		# using the location of the first transition
		soc_trans_row <- min(all_the_social_transition_rows);
		Trans_Table[, soc_trans_risks:=list(list(  Prop_Table[all_the_social_transition_rows, perceived_vaccine_risk] ))]

		soc_first_row <- soc_trans_row;
		soc_second_row <- soc_trans_row+1;

		x1 <- Prop_Table[soc_first_row, perceived_vaccine_risk];
		y1 <- Prop_Table[soc_first_row, soc_V_mean];

		x2 <- Prop_Table[soc_second_row, perceived_vaccine_risk];
		y2 <- Prop_Table[soc_second_row, soc_V_mean];

		x3 <- Prop_Table[soc_first_row, perceived_vaccine_risk];
		y3 <- Prop_Table[soc_first_row, soc_N_mean];

		x4 <- Prop_Table[soc_second_row, perceived_vaccine_risk];
		y4 <- Prop_Table[soc_second_row, soc_N_mean];

		soc_denom <- (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);

		soc_intersec_risk_numerator <- (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4);

		soc_intersec_prop_numerator <- (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4);

		Trans_Table[, ("soc_intersec_risk"):= soc_intersec_risk_numerator/soc_denom ];
		Trans_Table[, ("soc_intersec_prop"):= soc_intersec_prop_numerator/soc_denom ];
	}
	if( length(unique(Prop_Table[, phys_V_mean-phys_R_mean] < 0)) != 1 )
	{
		Trans_Table[, ("is_phys_trans"):= TRUE];

		# some of the curves have multiple transition points, so we'll just report all them
		all_the_physical_transition_rows <- Prop_Table[, which(diff(phys_R_mean>phys_V_mean)==1)];
		# using the first intersection of the curves for the physical transition
		phys_trans_row <- min(all_the_physical_transition_rows);
		Trans_Table[, phys_trans_risks:=list(list(  Prop_Table[all_the_physical_transition_rows, perceived_vaccine_risk]))]

		phys_first_row <- phys_trans_row;
		phys_second_row <- phys_trans_row+1;

		x1 <- Prop_Table[phys_first_row, perceived_vaccine_risk];
		y1 <- Prop_Table[phys_first_row, phys_R_mean];

		x2 <- Prop_Table[phys_second_row, perceived_vaccine_risk];
		y2 <- Prop_Table[phys_second_row, phys_R_mean];

		x3 <- Prop_Table[phys_first_row, perceived_vaccine_risk];
		y3 <- Prop_Table[phys_first_row, phys_V_mean];

		x4 <- Prop_Table[phys_second_row, perceived_vaccine_risk];
		y4 <- Prop_Table[phys_second_row, phys_V_mean];

		phys_denom <- (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);

		phys_intersec_risk_numerator <- (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4);

		phys_intersec_prop_numerator <- (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4);

		Trans_Table[, ("phys_intersec_risk"):= phys_intersec_risk_numerator/phys_denom ];
		Trans_Table[, ("phys_intersec_prop"):= phys_intersec_prop_numerator/phys_denom ];
	}

	return(Trans_Table);
}

change_point_test <- function(Prop_vs_Risk_table, metric, test_to_use="snh")
{
	mus <- data.table(kappa=c(NaN), p.value=c(NaN));

	# if(isTRUE( length(unique(Prop_vs_Risk_table$social_norm)) != 1 )){return( mus )}
	if(!isTRUE( metric %in% names(Prop_vs_Risk_table) )){return( mus )}
	Table_Here <- Prop_vs_Risk_table[order(perceived_vaccine_risk)];

	the_whole_series <- Table_Here[[metric]];
	if(length(the_whole_series)<3){ return(mus); }

	for(i in 3:length(the_whole_series))
	{
		test_data <- unique(the_whole_series[1:i]);
		if(sum(test_data) == 0) next;
		if(length(test_data) < 3) next;
			 if(test_to_use == "buishand_r"){ Test <- br.test(the_whole_series[1:i], m=5000); }
		else if(test_to_use == "lanzante"){ Test <- lanzante.test(the_whole_series[1:i]); }
		else if(test_to_use == "pettitt"){ Test <- pettitt.test(the_whole_series[1:i]); }
		else if(test_to_use == "snh"){ Test <- snh.test(the_whole_series[1:i]); }

		if(Test$p.value < 0.05)
		{
			mus <- data.table(kappa=c(Table_Here[Test$estimate[[1]], perceived_vaccine_risk]), p.value=c(Test$p.value))
			return(mus)
		};
	}
	return(mus)
}

as_sci <- function(vec){ formatC(vec, format = "f", digits = 1, width=2); }
