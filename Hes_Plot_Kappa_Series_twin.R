# This is the file for the first diagram tableau in the figure - the oen for multiple transitions is in another file
rm(list=ls())

options(scipen=10000)

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Kappa_Series", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

# source("Processing_Parameter_Table.R");
Parameters <- fread(parameter_file, colClasses="character");
Parameters <- Parameters[num(initial_vacc_proportion)==0.25 & network_structure=="smallworld" & num(proportion_of_nodes)==1];
Parameters[, c("V1","perceived_vaccine_risk", "social_norm"):=NULL];
Parameters <- unique(Parameters);

image_height <- 4;
image_width <- 25;

norms_to_use <- c(0, 0.25);

mains <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_V", "soc_N", "soc_H");
mains_mean <- sapply(mains, function(x) sprintf("%s_mean",x));
mains_sd <- sapply(mains, function(x) sprintf("%s_sd",x));

joins <- c("NN_join_count", "NV_join_count", "VV_join_count", "HH_join_count", "HV_join_count", "HN_join_count");
joins_mean <- sapply(joins, function(x) sprintf("%s_mean",x));
joins_sd <- sapply(joins, function(x) sprintf("%s_sd",x));

metrics_to_plot <- c("plain") # "join_counts", "mutual_info", "hesitant_conn_comp_size", "hesitant_echo_chamber_size", "prob_sick", "plain", "num_triads", "hesitant_joins", "NV_join_count", "conn_comp_size", "echo_chamber_size", "number", "watts_strogatz", "opinion_change", "echo_chamber_number", "conn_comp_number", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "vacc_conn_comp_size", "nonvacc_conn_comp_size", "net_diameter");

time_start <- Sys.time();

for(line_index in 1:nrow(Parameters))
{
    this_tuple <- Parameters[line_index];

    initial_vacc_prop <- convert_if_number(this_tuple$initial_vacc_prop);
    the_size <- convert_if_number(this_tuple$N);
    infection_prob <- convert_if_number(this_tuple$infec_prob);
    importation_rate <- "2.5e-05"; # convert_if_number(this_tuple$importation);
    structure <- convert_if_number(this_tuple$network_structure);
    degree <- convert_if_number(this_tuple$mean_degree);
    random_opinion_switch <- "1e-04"; # convert_if_number(this_tuple$random_opinion_switch);
    proportion <- convert_if_number(this_tuple$proportion_of_nodes);

    Proportion <- data.table();
    Transitions <-  data.table(norm=numeric(), soc=numeric(), phys=numeric());

    for(the_norm in norms_to_use)
    {
        Proportion_Here <- data.table();
        filenames <- Sys.glob(paste(CSV_file_path, "Hes_N_", the_size, "_struct_", structure, "_deg_", degree,  "_risk_*", "_infec_", infection_prob,
                            "_norm_", the_norm,  "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_", proportion,
                            "_summary.csv", sep=""));

        if(length(filenames) == 0){ next; }

        for(file in filenames)
        {
            temp <- data.table(read.csv(file));
            Proportion_Here <- rbind(Proportion_Here, temp[is.na(instance)], fill=TRUE);
        }
        Proportion_Here[, c("X", "instance"):=NULL];

        Trans_Table <- get_transitions(Proportion_Here, initial_vacc_prop);
        Transitions <- rbind(Transitions, list(the_norm, Trans_Table$soc_intersec_risk, Trans_Table$phys_intersec_risk));

        # Proportion_Here[, ("phys_trans"):=Trans_Table$phys_intersec_risk];
        # Proportion_Here[, ("soc_trans"):=Trans_Table$soc_intersec_risk];

        Proportion_Here[, (mains_sd)  :=lapply(.SD, "/", the_size), .SDcols=mains_sd  ];
        Proportion_Here[, (mains_mean):=lapply(.SD, "/", the_size), .SDcols=mains_mean];

        # do the sd first before normalising the means, since normalising the means changes to total 1
        # normalising the sd's by one doesn't do fuck all, then
        Proportion_Here[, (joins_sd)  :=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean+HN_join_count_mean+HH_join_count_mean+HV_join_count_mean), .SDcols=joins_sd ];
        Proportion_Here[, (joins_mean):=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean+HN_join_count_mean+HH_join_count_mean+HV_join_count_mean), .SDcols=joins_mean];

        Proportion <- rbind(Proportion, Proportion_Here[abs(num(perceived_vaccine_risk))<=0.1][order(perceived_vaccine_risk)]);
    }

    # if there are less than 5 points to be plotted, then cancel this plot and skip to the next parameter set
    if(nrow(Proportion)<5){ next; }

    for(metric in metrics_to_plot)
    {
        if(metric != "plain"){
        for(the_thing in plot_these(metric))
        {
            if(!( the_thing %in% names(metric) )) next;
        }}

        print(metric);
        plots_to_make <- plot_these(metric);
        DT <- data.table(grp=character(), norm=numeric(), risk=numeric(), val=numeric(), lwr=numeric(), upr=numeric());

        for(the_norm in norms_to_use){ for(index in seq_along(plots_to_make))
        {
            Proportion_Here <- Proportion[social_norm==the_norm];

            prefix <- plot_these(metric)[index];
            tmean <- sprintf("%s_mean", prefix);
            tsd <- sprintf("%s_sd", prefix);

            temp_table <- data.table(
                grp = rep(prefix, nrow(Proportion_Here)),
                norm = rep(the_norm, nrow(Proportion_Here)),
                risk = Proportion_Here[, perceived_vaccine_risk],
                val = Proportion_Here[, get(tmean)],
                lwr = Proportion_Here[, get(tmean)-get(tsd)],
                upr = Proportion_Here[, get(tmean)+get(tsd)]
            );

            DT <- data.table(do.call(smartbind, list(DT, temp_table)));
        }}

        text_size <- 30;

        scientific_10 <- function(x) {
            parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
        }

        pl <- ggplot(data=DT, aes(x=risk, y=val, ymin=lwr, ymax=upr, fill=grp, colour=grp)) + #
            geom_line(size=2) +
            geom_ribbon(alpha=.3, linetype=0) +
            xlab(expression(kappa)) + ylab(axis_label(metric)) +
            theme(
                    # legend.position="bottom",
                    # legend.direction="horizontal",
                    # legend.key.width=unit(4, 'cm'),
                    panel.background=element_blank(),
                    axis.text=element_text(size=text_size),
                    legend.text=element_text(size=text_size+5),
                    strip.text.x = element_text(size = text_size),
                    axis.title=element_text(size=text_size+10),
                    panel.spacing = unit(3, "lines"),
                    legend.title=element_blank(),
                    legend.key.width=unit(1, 'cm'),
                    # legend.spacing.x=unit(1, 'cm'),
                    legend.direction="vertical",
                    legend.justification='center',
                    legend.position="right",
            ) +
            scale_colour_discrete(labels=as.vector(sapply(plots_to_make, axis_label))) +
            scale_fill_discrete(labels=as.vector(sapply(plots_to_make, axis_label))) +
            geom_vline(aes(xintercept=phys), data=Transitions, size=0.5, linetype="dashed") +
            geom_vline(aes(xintercept=soc),  data=Transitions, size=0.5, linetype="dotted") +
            facet_wrap(~norm, nrow=1, scales='free', labeller=label_bquote(sigma*' = '*.(norm)));

        if(metric %in% c("vacc_conn_comp_size", "nonvacc_conn_comp_size", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "num_triads")){ pl <- pl + scale_y_continuous(labels=scientific_10); }

        plot_file_name <- paste(graph_file_path, folder_name, "/Hes_twin_Kappa_", metric, "_N_", the_size, "_struct_", structure, "_deg_", degree,
            "_infec_", infection_prob, "_imp_", importation_rate, "_init_", initial_vacc_prop, "_switch_", random_opinion_switch, "_prop_", proportion, ".png", sep=""
        );

        ggsave(
            plot_file_name,
            plot=pl, width=image_width, height=image_height, limitsize=FALSE, dpi=50 #  height=5,
        )
        dev.off();

        system(sprintf("convert %s -trim %s", plot_file_name, plot_file_name));
        # stop()

    }

}

print(Sys.time()-time_start);

# make_mean_and_sd <- function(data, average_by, dependent_var)
# {
# 	mains <- setdiff(names(Here), c(average_by, dependent_var))
#
# 	mains_mean 	<- sprintf("%s_mean", mains);
# 	mains_sd 	<- sprintf("%s_sd", mains);
#
# 	Table_Here <- data.table( dependent = unique(Here[[dependent_var]]) )
# 	Table_Here[, (mains_mean):=data[, lapply(.SD, mean), by=time_step][, .SD, .SDcols=setdiff(names(data), c(average_by, dependent_var))]];
# 	Table_Here[, (mains_sd  ):=data[, lapply(.SD, sd  ), by=time_step][, .SD, .SDcols=setdiff(names(data), c(average_by, dependent_var))]];
#
# 	DT <- data.table(dep=numeric(), name=character(), value=numeric(), lower=numeric(), upper=numeric());
#
# 	for( dep_var in unique(data[[dependent_var]]) ){ for( prefix in mains )
# 	{
# 		DT <- rbind(DT, list(
# 			dep_var,
# 			prefix,
# 			Table_Here[dependent==dep_var, get(sprintf("%s_mean", prefix))],
# 			Table_Here[dependent==dep_var, get(sprintf("%s_mean", prefix)) - get(sprintf("%s_sd", prefix))],
# 			Table_Here[dependent==dep_var, get(sprintf("%s_mean", prefix)) + get(sprintf("%s_sd", prefix))]
# 		));
# 	}}
#
# 	return(DT);
# }
#
# RAW_DATA <- data.table(fread("Model_Output.csv"));
#
# parameters <- c("B_H", "B_C", "B_0", "not_same_school", "childless_teachers", "background_inf", "R_init", "prob_symptomatic", "max_capacity", "E_to_P", "P_to_Infected", "I_to_R", "A_to_R")
#
# DT <- RAW_DATA[, .SD, .SDcols=setdiff(names(RAW_DATA), parameters)];
#
# # Trends_Infected <- DT[time_step<=20, .SD, .SDcols=c("child_teacher_ratio", "class_grouping", "instance", "time_step" , "prop_I", "prop_I_in_school", "inf_background", "inf_home", "inf_class", "inf_commons")];
#
# Trends_Infected <- DT[, .SD, .SDcols=c("child_teacher_ratio", "class_grouping", "teacher_replacement", "instance", "time_step" , "prop_I", "prop_I_in_school", "inf_background", "inf_home", "inf_class", "inf_commons")];
#
# for(replace in unique(Trends_Infected$teacher_replacement)){
# for(ratio in unique(Trends_Infected$child_teacher_ratio)){
# for(groups in unique(Trends_Infected$class_grouping)){
#
# 	if(replace == 0){ replacement_string <- "not replacing teachers"; }
# 	else { replacement_string <- "replacing teachers"; }
#
# 	Here <- Trends_Infected[(child_teacher_ratio==ratio) & (class_grouping==groups) & (teacher_replacement==replace)];
# 	Here[, c("child_teacher_ratio", "class_grouping", "prop_I", "prop_I_in_school", "teacher_replacement"):=NULL];
#
# 	pl_here <- ggplot(
# 			make_mean_and_sd(data=Here, average_by="instance", dependent="time_step"),
# 			aes(x=dep, y=value, colour=name, fill=name, ymin=lower, ymax=upper)
# 		) +
# 		geom_line(size=0.75) +
# 		geom_ribbon(alpha=.3, linetype=0) +
# 		scale_x_continuous(expand = c(0,0)) +
# 		scale_y_continuous(expand = c(0,0)) +
# 		labs(x="Day", y="Number of Infections", colour=sprintf("Locale, %s, %s,\n%s", ratio, groups, replacement_string)) +
# 		theme(
# 			axis.text=element_text(size=15),
# 			axis.title=element_text(size=15),
# 			legend.text=element_text(size=15),
# 			legend.title=element_text(size=15)
# 		) +
# 		scale_colour_discrete(breaks=sprintf("inf_%s", c("background", "home", "class", "commons")), labels=c('Background','Household','Classroom','Common areas')) +
# 		scale_fill_discrete(name=sprintf("Locale, %s, %s", ratio, groups), breaks=sprintf("inf_%s", c("background", "home", "class", "commons")), labels=c('Background','Household','Classroom','Common areas')) +
# 		theme_bw();
#
# 	ggsave(plot=pl_here, file=sprintf("covid_locales_ratio_%s_arrangement_%s_teacher_replacement_%s.png", ratio, groups, replace), width=12, height=6, units='in');
#
# 	stop()
#
# }}}
#
# scale_colour_discrete(name=sprintf("Locale, %s, %s", ratio, groups), breaks=sprintf("inf_%s", c("background", "home", "class", "commons")), labels=c('Background','Household','Classroom','Common areas')) +
# scale_fill_discrete(name=sprintf("Locale, %s, %s", ratio, groups), breaks=sprintf("inf_%s", c("background", "home", "class", "commons")), labels=c('Background','Household','Classroom','Common areas')) +
#
