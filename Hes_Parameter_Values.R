# Brendon  Phillips
# PhD candidate
# Department of Applied Mathematics
# Faculty of Mathematics
# University of Waterloo
# rm(list=ls())

options(width=Sys.getenv("COLUMNS"));
options(bitmapType="cairo");

# source("find_colours.R")

if(!require(pacman)) install.packages('pacman')
pacman::p_load(
				akima,
				# Cairo,
				changepoint,
				colorRamps,
				comprehenr,
				# ContourFunctions,
				# colorspace,
				cowplot,
				data.table,
				doParallel,
				digest,
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
				pbmcapply,
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

fragments <- c("N", "struct", "deg", "risk", "infec", "norm", "imp", "init", "switch", "prop");
options(device='x11')

# Sys.setenv('MAPBOX_TOKEN'='pk.eyJ1IjoiYjJwaGlsbGkiLCJhIjoiY2s1NjZ5bzJlMDI4dzNucGtkdjA2NjNydiJ9.CNeuv4-hVP3cRnfXtvP3Ew');
# Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/anaconda3/bin", sep = .Platform$path.sep))

# setHook(packageEvent("grDevices", "onLoad"),
# function(...) grDevices::X11.options(type='cairo'))

# font_add_google("Lobster")

user_name <- "b2philli"; # Sys.info()[8][[1]];

Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/anaconda3/bin", sep = .Platform$path.sep))

# home_folder <- sprintf("/media/%s/Simulations", user_name);
home_folder <- sprintf("/home/%s/Dropbox", user_name);
parameter_file <- "/home/b2philli/Dropbox/CSV/Hes_Parameters.csv";
stored_data_path <- sprintf("/media/%s/Simulations/Hesitance", user_name);

threshold_for_declaring_high_proportion <- 0.9;
threshold_for_declaring_low_proportion <- 0.1;

number_of_time_steps_averaged <- 500;

vaxx_allowed_vector <- c(TRUE);
use_this_many_entries_for_ratio_sec_new <- 500;

# values that do not change between the time steps in each individual trial or not to be averaged
hes_fixed_values <- c("N", "network_structure", "mean_degree", "perceived_vaccine_risk", "infec_prob", "social_norm", "importation", "initial_vacc_proportion", "random_opinion_switch", "proportion_of_nodes");

# the parameter values tha we want to reisnert to the data that were excluded above
# reinserted_values <- setdiff(hes_fixed_values, c("instance"));

user_name <- Sys.info()[8][[1]];

CSV_file_path <-   sprintf("%s/CSV/", home_folder);
graph_file_path <- sprintf("/home/%s/Dropbox/Apps/Overleaf/phd_thesis/Thesis_Graphs/", user_name);
# graph_file_path <- sprintf("/home/%s/Dropbox/Apps/Overleaf/warning_signals_paper/Graphs/", user_name);

maxx <- function(x) max(x, na.rm=T);
minn <- function(x) min(x, na.rm=T);

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
	if(x == "number"){ return(quote('<# '['*']*'>')); }
	if(x == "Dyn"){ return("S/P");}

	#############################################################################################################

	if(x == "lanzante"){ return("Lan"); }
	if(x == "pettitt"){ return("Pet"); }
	if(x == "snh"){ return("SNHT"); }
	if(x == "buishand_r"){ return("BR"); }

	#############################################################################################################

	# time series
	if(x == "outputs"){ return("Outputs"); }
	if(x %in% c("phys_S")){ return("[S]"); }
	if(x %in% c("phys_I")){ return("[I]"); }
	if(x %in% c("phys_R")){ return("[R]"); }
	if(x %in% c("phys_V")){ return(quote('[V'['p']*']')); }
	if(x %in% c("soc_N")){ return("[N]"); }
	if(x %in% c("soc_H")){ return("[H]"); }
	if(x %in% c("soc_V")){ return(quote('[V'['s']*']')); }

	# if(x == "NN_join_count"){ return("[N,N]"); }
	# if(x == "NV_join_count"){ return( parse(text='"["*scriptstyle(N,)*","*scriptstyle(V[s])*"]"') ); }
	# if(x == "VV_join_count"){ return("[V,V]"); }

	#############################################################################################################

	if(x %in% c("phys.S.change", "phys_S_mean")){ return("<S>"); }
	if(x %in% c("phys.I.change", "phys_I_mean")){ return("<I>"); }
	if(x %in% c("phys.R.change", "phys_R_mean")){ return("<R>"); }
	if(x %in% c("phys.V.change", "phys_V_mean")){ return(quote('<V'['p']*'>')); }
	if(x %in% c("soc.V.change", "soc_V_mean")){ return(quote('<V'['s']*'>')); }
	if(x %in% c("soc.N.change", "soc_N_mean")){ return("<N>"); }
	if(x %in% c("soc.H.change", "soc_H_mean")){ return("<H>"); }

	if(x == "watts_strogatz"){ return(quote('<C'['*']*'>')); }
	if(x %in% c("watts.strogatz.all.change", "watts_strogatz_all")){ return(quote('<C'[Sigma]*'>')); }
	if(x %in% c("watts.strogatz.nonvacc.change", "watts_strogatz_nonvacc")){ return(quote('<C'['N']*'>')); }
	if(x %in% c("watts.strogatz.vacc.change", "watts_strogatz_vacc")){ return(quote('<C'['V']*'>')); }
	if(x %in% c("watts.strogatz.hesitant.change", "watts_strogatz_hesitant")){ return(quote('<C'['H']*'>')); }

	if(x == "join_counts"){ return(quote('<*,*>')); }
	if(x == "hesitant_joins"){ return(quote('<H,*>')); }
	if(x %in% c("NN.join.count.change", "NN_join_count_mean", "NN_join_count")){ return(quote('<'*'N,N'*'>') ); }
	if(x %in% c("NV.join.count.change", "NV_join_count_mean", "NV_join_count")){ return(quote('<'*'N,V'['s']*'>')); }
	if(x %in% c("VV.join.count.change", "VV_join_count_mean", "VV_join_count")){ return(quote('<'*'V'['s']*',V'['s']*'>') ); }
	if(x %in% c("HH.join.count.change", "HH_join_count_mean", "HH_join_count")){ return(quote('<'*'H,H>')); }
	if(x %in% c("HV.join.count.change", "HV_join_count_mean", "HV_join_count")){ return(quote('<'*'H,V'['s']*'>') ); }
	if(x %in% c("HN.join.count.change", "HN_join_count_mean", "HN_join_count")){ return(quote('<'*'H,N>')); }

	if(x == "opinion_change"){ return(quote('<'*Theta ['*']*'>')); }
	if(x %in% c("nonvacc.opinion.change.change", "nonvacc_opinion_change")){ return(quote('<'*Theta['N']*'>')); }
	if(x %in% c("vacc.opinion.change.change", "vacc_opinion_change")){ return(quote('<'*Theta['V']*'>')); }
	if(x %in% c("hesitant.opinion.change.change", "hesitant_opinion_change")){ return(quote('<'*Theta['H']*'>')); }
	if(x %in% c("total.opinion.change.change", "total_opinion_change")){ return(quote('<'*Theta[Sigma]*'>')); }

	if(x == "conn_comp_size"){ return(quote('|Z'['*']*'|')); }

	if(x == "vacc_conn_comp_size"){ return(quote('|Z'['V']*'|')); }
	if(x %in% c("vacc.avg.conn.comp.size.change", "vacc_avg_conn_comp_size")){ return(quote('<|Z'['V']*'|>')); }
	if(x %in% c("vacc.min.conn.comp.size.change", "vacc_min_conn_comp_size")){ return(quote('min(|Z'['V']*'|)')); }
	if(x %in% c("vacc.max.conn.comp.size.change", "vacc_max_conn_comp_size")){ return(quote('max(|Z'['V']*'|)')); }

	if(x == "nonvacc_conn_comp_size"){ return(quote('|Z'['N']*'|')); }
	if(x %in% c("nonvacc.avg.conn.comp.size.change", "nonvacc_avg_conn_comp_size")){ return(quote('<|Z'['N']*'|>')); }
	if(x %in% c("nonvacc.min.conn.comp.size.change", "nonvacc_min_conn_comp_size")){ return(quote('min(|Z'['N']*'|)')); }
	if(x %in% c("nonvacc.max.conn.comp.size.change", "nonvacc_max_conn_comp_size")){ return(quote('max(|Z'['N']*'|)')); }

	if(x == "hesitant_conn_comp_size"){ return(quote('|Z'['H']*'|')); }
	if(x %in% c("hesitant.avg.conn.comp.size.change", "hesitant_avg_conn_comp_size")){ return(quote('<|Z'['H']*'|>')); }
	if(x %in% c("hesitant.min.conn.comp.size.change", "hesitant_min_conn_comp_size")){ return(quote('min(|Z'['H']*'|)')); }
	if(x %in% c("hesitant.max.conn.comp.size.change", "hesitant_max_conn_comp_size")){ return(quote('max(|Z'['H']*'|)')); }

	if(x == "echo_chamber_size"){ return(quote('|J'['*']*'|')); }

	if(x == "nonvacc_echo_chamber_size"){ return(quote('|J'['N']*'|')); }
	if(x %in% c("nonvacc.avg.chamber.size.change", "nonvacc_avg_chamber_size")){ return(quote('<|J'['N']*'|>')); }
	if(x %in% c("nonvacc.min.chamber.size.change", "nonvacc_min_chamber_size")){ return(quote('min(|J'['N']*'|)')); }
	if(x %in% c("nonvacc.max.chamber.size.change", "nonvacc_max_chamber_size")){ return(quote('max(|J'['N']*'|)')); }

	if(x == "vacc_echo_chamber_size"){ return(quote('|J'['V']*'|')); }
	if(x %in% c("vacc.avg.chamber.size.change", "vacc_avg_chamber_size")){ return(quote('<|J'['V']*'|>')); }
	if(x %in% c("vacc.min.chamber.size.change", "vacc_min_chamber_size")){ return(quote('min(|J'['V']*'|)')); }
	if(x %in% c("vacc.max.chamber.size.change", "vacc_max_chamber_size")){ return(quote('max(|J'['V']*'|)')); }

	if(x == "hesitant_echo_chamber_size"){ return(quote('|J'['H']*'|')); }
	if(x %in% c("hesitant.avg.chamber.size.change", "hesitant_avg_chamber_size")){ return(quote('<|J'['H']*'|>')); }
	if(x %in% c("hesitant.min.chamber.size.change", "hesitant_min_chamber_size")){ return(quote('min(|J'['H']*'|)')); }
	if(x %in% c("hesitant.max.chamber.size.change", "hesitant_max_chamber_size")){ return(quote('max(|J'['H']*'|)')); }

	if(x == "conn_comp_number"){ return(quote('<#Z'['*']*'>')); }
	if(x %in% c("nonvacc.number.conn.comps.change", "nonvacc_number_conn_comps")){ return(quote('<#Z'['N']*'>')); }
	if(x %in% c("vacc.number.conn.comps.change", "vacc_number_conn_comps")){ return(quote('<#Z'['V']*'>')); }
	if(x %in% c("hesitant.number.conn.comps.change", "hesitant_number_conn_comps")){ return(quote('<#Z'['H']*'>')); }

	if(x == "echo_chamber_number"){ return(quote('<#J'['*']*'>')); }
	if(x %in% c("nonvacc.number.chambers.change", "nonvacc_number_chambers")){ return(quote('<#J'['N']*'>')); }
	if(x %in% c("vacc.number.chambers.change", "vacc_number_chambers")){ return(quote('<#J'['V']*'>')); }
	if(x %in% c("hesitant.number.chambers.change", "hesitant_number_chambers")){ return(quote('<#J'['H']*'>')); }

	if(x == "net_diameter"){ return(quote('<'*Omega['*']*'>')); }
	if(x %in% c("hesitant.diameter.change", "hesitant_diameter")){ return(quote('<'*Omega['H']*'>')); }
	if(x %in% c("nonvacc.diameter.change", "nonvacc_diameter")){ return(quote('<'*Omega['N']*'>')); }
	if(x %in% c("vacc.diameter.change", "vacc_diameter")){ return(quote('<'*Omega['V']*'>')); }

	if(x == "num_triads"){ return(quote('<'*Delta['*']*'>')); }
	if(x %in% c("num.hesitant.triads.change", "num_hesitant_triads")){ return(quote('<'*Delta['H']*'>')); }
	if(x %in% c("num.nonvacc.triads.change", "num_nonvacc_triads")){ return(quote('<'*Delta['N']*'>')); }
	if(x %in% c("num.vacc.triads.change", "num_vacc_triads")){ return(quote('<'*Delta['V']*'>')); }
	if(x %in% c("total.num.triads.change", "total_num_triads")){ return(quote('<'*Delta[Sigma]*'>')); }

	if(x == "prob_sick"){ return(quote('<'*Gamma['*']*'>')); }
	if(x %in% c("prob.sick.nonvacc.change", "prob_sick_nonvacc")){ return(quote('<'*Gamma['N']*'>')); }
	if(x %in% c("prob.sick.vacc.change", "prob_sick_vacc")){ return(quote('<'*Gamma['V']*'>')); }
	if(x %in% c("prob.sick.hesitant.change", "prob_sick_hesitant")){ return(quote('<'*Gamma['H']*'>')); }

	if(x %in% c("mutual.info.change", "mutual", "mutual_info")){ return(quote('<M>')); }

	else return(x);

	# if(x == "moran_i"){ return('<I>'); }
	# if(x == "geary_c"){ return('<C>'); }
	# if(x == "mutual_info"){ return('<M>'); }
	# if(x == "getis_ord"){ return("Getis-Ord"); }
	# if(x == "outputs"){ return("Outputs"); }
	#
	# if(x == "modularity"){ return(quote('<Q'['*']*'>')); }
	# if(x == "prob_sick"){ return(quote('<'*Gamma ['*']*'>')); }
	# if(x == "watts_strogatz"){ return(quote('<C'['*']*'>')); }
	#
	# # if(x == "ratio_sick_nbr"){ return(expression('<'*Xi*'>')); }
	#
    # if(x == "vacc.modularity.change"){ return(expression('<Q'['V']*'>')); }
    # if(x == "nonvacc.modularity.change"){ return(expression('<Q'['N']*'>')); }
    # if(x == "total.modularity.change"){ return(expression('<Q'[Sigma]*'>')); }
	#
	# if(x == "moran.i.change"){ return(quote('<I>')); }
	# if(x == "geary.c.change"){ return(quote('<C>')); }
	#
	# if(x == "V_Conn_Size"){ return(quote('|Z'['V']*'|')); }
	# if(x == "NV_Conn_Size"){ return(quote('|Z'['N']*'|')); }
    # if(x == "Conn_Size"){ return(quote('|Z'['*']*'|')); }
	# if(x == "Conn_Echo_Count"){ return(quote('<# '['*']*'>')); }
	#
	# if(x == "V_Echo_Size"){ return(quote('|J'['V']*'|')); }
	# if(x == "NV_Echo_Size"){ return(quote('|J'['N']*'|')); }
    # if(x == "Echo_Size"){ return(quote('|J'['*']*'|')) }
	#
	# if(x == "Prob_Sick"){ return(quote('<'*Gamma['*']*'>')); }
	# if(x == "WattsGCC"){ return(quote('<C'['*']*'>')); }
	# if(x == "Modu"){ return(quote('<Q'['*']*'>'));  }
    # if(x == "Senti"){ return(quote('<'*Theta['*']*'>')); }
	#
	# if(x %in% c("net_diameter", "net_diameter_mean", "net.diameter.change")){ return(quote('<'*Omega['*']*'>')); }
	# if(x %in% c("num_triads", "num_triads_mean", "num.triads.change")){ return(quote('<'*Delta['*']*'>')); }
	#
	# if(x %in% c("hesitant_joins", "hesitant_joins_mean", "hesitant.joins.change")){ return(quote('[H,*]')); }
}

plot_these <- function(metric)
{
	if(metric %in% c("mutual_info", "mutual")) 		{ cest_voila <- c("mutual_info"); }
	else if(metric == "outputs")					{ cest_voila <- c("phys_R", "phys_V", "phys_S", "phys_I", "soc_N", "soc_V", "soc_H"); }
	else if(metric == "plain")						{ cest_voila <- c("phys_R", "phys_V", "soc_N", "soc_V", "soc_H"); }
	else if(metric == "join_counts") 	 			{ cest_voila <- c("NV_join_count", "NN_join_count", "VV_join_count", "HH_join_count", "HN_join_count", "HV_join_count"); }
	else if(metric == "NV_join_count")				{ cest_voila <- c("NV_join_count"); }
	else if(metric == "hesitant_joins") 	 		{ cest_voila <- c("HH_join_count", "HN_join_count", "HV_join_count"); }
	else if(metric == "vacc_conn_comp_size") 		{ cest_voila <- c("vacc_min_conn_comp_size", "vacc_max_conn_comp_size", "vacc_avg_conn_comp_size"); }
	else if(metric == "nonvacc_conn_comp_size") 	{ cest_voila <- c("nonvacc_min_conn_comp_size", "nonvacc_max_conn_comp_size", "nonvacc_avg_conn_comp_size"); }
	else if(metric == "hesitant_conn_comp_size")	{ cest_voila <- c("hesitant_min_conn_comp_size", "hesitant_max_conn_comp_size", "hesitant_avg_conn_comp_size"); }
	else if(metric == "conn_comp_size")				{ cest_voila <- c("vacc_avg_conn_comp_size", "nonvacc_avg_conn_comp_size", "hesitant_avg_conn_comp_size");}  # "vacc_min_conn_comp_size", "vacc_max_conn_comp_size", "nonvacc_min_conn_comp_size", "nonvacc_max_conn_comp_size", "hesitant_min_conn_comp_size", "hesitant_max_conn_comp_size",
	else if(metric == "number")                     { cest_voila <- c("vacc_number_conn_comps", "nonvacc_number_conn_comps", "hesitant_number_conn_comps", "vacc_number_chambers", "nonvacc_number_chambers", "hesitant_number_chambers"); }
	else if(metric == "conn_comp_number") 			{ cest_voila <- c("vacc_number_conn_comps", "nonvacc_number_conn_comps", "hesitant_number_conn_comps"); }
	else if(metric == "vacc_echo_chamber_size") 	{ cest_voila <- c("vacc_min_chamber_size", "vacc_max_chamber_size", "vacc_avg_chamber_size");  }
	else if(metric == "nonvacc_echo_chamber_size") 	{ cest_voila <- c("nonvacc_min_chamber_size", "nonvacc_max_chamber_size", "nonvacc_avg_chamber_size");  }
	else if(metric == "hesitant_echo_chamber_size")	{ cest_voila <- c("hesitant_min_chamber_size", "hesitant_max_chamber_size", "hesitant_avg_chamber_size"); }
	else if(metric == "echo_chamber_number")		{ cest_voila <- c("vacc_number_chambers", "nonvacc_number_chambers", "hesitant_number_chambers"); }
	else if(metric == "echo_chamber_size")			{ cest_voila <- c("vacc_avg_chamber_size", "nonvacc_avg_chamber_size", "hesitant_avg_chamber_size"); }
	else if(metric == "conn_comp_size")				{ cest_voila <- c("vacc_avg_conn_comp_size", "nonvacc_avg_conn_comp_size", "hesitant_avg_conn_comp_size"); }
	else if(metric == "conn_echo_count")			{ cest_voila <- c("vacc_number_conn_comps_mean", "nonvacc_number_conn_comps_mean", "hesitant_number_conn_comps_mean", "vacc_number_chambers_mean", "nonvacc_number_chambers_mean", "hesitant_number_chambers_mean"); }
	else if(metric == "watts_strogatz")				{ cest_voila <- c("watts_strogatz_vacc", "watts_strogatz_nonvacc", "watts_strogatz_all", "watts_strogatz_hesitant"); }
	else if(metric == "opinion_change")				{ cest_voila <- c("nonvacc_opinion_change", "vacc_opinion_change", "hesitant_opinion_change"); } # , "total_opinion_change"
	else if(metric == "ratio_sick_nbr")				{ cest_voila <- c("ratio_hesitant_infected_neighbours", "ratio_vacc_infected_neighbours", "ratio_nonvacc_infected_neighbours"); }
	else if(metric == "prob_sick")					{ cest_voila <- c("prob_sick_vacc", "prob_sick_nonvacc", "prob_sick_hesitant"); }
	else if(metric == "num_triads")					{ cest_voila <- c("num_hesitant_triads", "num_vacc_triads", "num_nonvacc_triads", "total_num_triads"); } #
	else if(metric == "net_diameter")				{ cest_voila <- c("hesitant_diameter", "nonvacc_diameter", "vacc_diameter"); }
	else if(metric == "modularity")					{ cest_voila <- c("vacc_modularity", "nonvacc_modularity", "total_modularity"); }
	if( !exists("cest_voila") ){ print("no plotting variables specified for the requested plot");}
	return(cest_voila);
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
			.SDcols=c("perceived_vaccine_risk", "phys_V_mean", "phys_R_mean", "initial_vacc_proportion", "soc_V_mean", "soc_N_mean", "soc_H_mean")
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

	if( length(unique(Prop_Table[, soc_V_mean-(soc_N_mean + soc_H_mean)] < 0)) != 1 )
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
		y3 <- Prop_Table[soc_first_row, soc_N_mean+soc_H_mean];

		x4 <- Prop_Table[soc_second_row, perceived_vaccine_risk];
		y4 <- Prop_Table[soc_second_row, soc_N_mean+soc_H_mean];

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

the_series <- c(0.0048653955, 0.0048992766, 0.0049031573,0.0049565017, 0.0050303633, 0.0050837171, 0.0050757092, 0.0051124101, 0.0051745790, 0.0052399448, 0.0052638722, 0.0053621846, 0.0053520912, 0.0053784857, 0.0054950007, 0.0055586844, 0.0055813591, 0.0055697585, 0.0057371584, 0.0058000484, 0.0057900996, 0.0059269815, 0.0058087678, 0.0059436082, 0.0059905875, 0.0059627250, 0.0060330715, 0.0060017463, 0.0061239627, 0.0061141295, 0.0061637498, 0.0021726667, 0.0062169640, 0.0062954884, 0.0063006051, 0.0062997526, 0.0024059693, 0.0064589441, 0.0063750514, 0.0064429451, 0.0064957918, 0.0022827610, 0.0065731248, 0.0065746145, 0.0066977065, 0.0066663136, 0.0240787354, 0.3837142934, 0.3040702428, 0.0016807157, 0.0016682305, 0.0007126583, 0.0016262820, 0.0016543582, 0.0016555937, 0.0016355423, 0.0007749936, 0.0016178569, 0.0016303006, 0.0015658200, 0.0015516826, 0.0006643949, 0.0015889219, 0.0015696159, 0.0015092203, 0.0015086613, 0.0014730193, 0.0014673662, 0.0014594937, 0.0014357535, 0.0012101372);

change_point_test <- function(Prop_vs_Risk_table, metric, test_to_use="snh")
{
	mus <- data.table(kappa=c(NaN), p.value=c(NaN));

	if(!isTRUE( metric %in% names(Prop_vs_Risk_table) )){return( mus )}
	Table_Here <- Prop_vs_Risk_table[order(perceived_vaccine_risk)];

	the_whole_series <- Table_Here[[metric]];
	# print(the_whole_series); stop();
	if(length(the_whole_series)<3){ return(mus); }

	if(test_to_use == "binseg")
	{
		temp <- c();
		for(i in 3:length(the_whole_series))
		{
			temp <- cpt.mean(the_whole_series[1:i], method="BinSeg", Q=1);
			if(length(cpts(temp)) != 0) break;
		}

	} else
	{
		for(i in 3:length(the_whole_series))
		{
			go_to_next <-FALSE;
			tryCatch(
				{
					if(test_to_use == "buishand_r"){ Test <- br.test(the_whole_series[1:i], m=5000); }
					else if(test_to_use == "lanzante"){ Test <- lanzante.test(the_whole_series[1:i]); }
					else if(test_to_use == "pettitt"){ Test <- pettitt.test(the_whole_series[1:i]); }
					else if(test_to_use == "snh"){ Test <- snh.test(the_whole_series[1:i]); }
				},
				error=function(e){ go_to_next<-TRUE; }
			);
			if(go_to_next){ next; }
			if(exists("Test"))
			{
				# print(Test)
				if(!is.na(Test$p.value)){ if(Test$p.value < 0.05)
				{
					mus <- data.table(kappa=c(Table_Here[Test$estimate[[1]], perceived_vaccine_risk]), p.value=c(Test$p.value))
					return(mus);
				}}
			}
			rm(Test);
		}
	}
	return(mus)
}

as_sci <- function(vec){ formatC(vec, format = "f", digits = 1, width=2); }

convert_if_number <- function(x)
{
	if(!is.na(suppressWarnings(as.numeric(x)))) return(as.numeric(x))
	return(x)
}

make_matrix <- function(the_table, row_names, col_names, the_vals)
{
	rowvals <- sort(the_table[, unique(get(row_names))]);
	colvals <- sort(the_table[, unique(get(col_names))]);
	DM <- matrix( 0, nrow=length(rowvals), ncol=length(colvals) );
	for(row_i in seq_along(rowvals)){
		for(col_i in seq_along(colvals)){
			prospect <- the_table[get(row_names)==rowvals[row_i] & get(col_names)==colvals[col_i], get(the_vals)]
			DM[row_i, col_i] <- if(length(prospect) == 0) NaN else prospect;
		}
	}
	rownames(DM) <- sapply(rowvals, toString);
	colnames(DM) <- sapply(colvals, toString);
	return( DM );
}

matrix.axes <- function(data) {
	x <- (1:dim(data)[1] - 1) / (dim(data)[1] - 1);
	axis(side=1, at=x, labels=rownames(data), las=2);
	x <- (1:dim(data)[2] - 1) / (dim(data)[2] - 1);
	axis(side=2, at=x, labels=colnames(data), las=2);
}

# change_point_test <- function(Prop_vs_Risk_table, metric, test_to_use="snh")
# {
# 	mus <- data.table(kappa=c(NaN), p.value=c(NaN));
#
# 	# if(isTRUE( length(unique(Prop_vs_Risk_table$social_norm)) != 1 )){return( mus )}
# 	if(!isTRUE( metric %in% names(Prop_vs_Risk_table) )){return( mus )}
# 	Table_Here <- Prop_vs_Risk_table[order(perceived_vaccine_risk)];
#
# 	the_whole_series <- Table_Here[[metric]];
# 	if(length(the_whole_series)<3){ return(mus); }
#
# 	# parallel_kappa_series <- Table_Here$perceived_vaccine_risk;
# 	# the_non_na_entries <- which(!is.na(the_whole_series));
# 	# the_whole_series <- the_whole_series[the_non_na_entries];
# 	# parallel_kappa_series <- parallel_kappa_series[the_non_na_entries];
#
# 	for(i in 3:length(the_whole_series))
# 	{
# 		go_to_next <-FALSE;
# 		tryCatch(
# 			{
# 				if(test_to_use == "buishand_r"){ Test <- br.test(the_whole_series[1:i], m=5000); }
# 				else if(test_to_use == "lanzante"){ Test <- lanzante.test(the_whole_series[1:i]); }
# 				else if(test_to_use == "pettitt"){ Test <- pettitt.test(the_whole_series[1:i]); }
# 				else if(test_to_use == "snh"){ Test <- snh.test(the_whole_series[1:i]); }
# 			},
# 			warning=function(w){ go_to_next <- TRUE; },
# 			error=function(e){ go_to_next<-TRUE; }
# 		);
# 		if(go_to_next){ next; }
# 		if(exists("Test")){ if(Test$p.value < 0.05)
# 		{
# 			mus <- data.table(kappa=c(Table_Here[Test$estimate[[1]], perceived_vaccine_risk]), p.value=c(Test$p.value))
# 			return(mus)
# 		}}
# 	}
# 	return(mus)
# }
