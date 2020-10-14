# This is the file for the first diagram tableau in the figure - the oen for multiple transitions is in another file
rm(list=ls())

options(scipen=10000)

source("/home/b2philli/Dropbox/Processing/Hes_Parameter_Values.R");

folder_name <- paste("Kappa_Series", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Parameters <- fread(parameter_file, colClasses="character");
Parameters <- Parameters[num(initial_vacc_proportion)==0.25 & network_structure=="smallworld" & num(social_norm)%in%c(0, 0.25) &  num(proportion_of_nodes)==1];
Parameters[, c("V1","perceived_vaccine_risk"):=NULL];
Parameters <- unique(Parameters);

image_height <- 6.5;
image_width <- 22;

metrics_to_plot <- c("join_counts", "opinion_change", "mutual_info", "hesitant_conn_comp_size", "hesitant_echo_chamber_size", "prob_sick", "plain", "num_triads", "hesitant_joins", "NV_join_count", "conn_comp_size", "echo_chamber_size", "number", "watts_strogatz", "echo_chamber_number", "conn_comp_number", "vacc_echo_chamber_size", "nonvacc_echo_chamber_size", "vacc_conn_comp_size", "nonvacc_conn_comp_size", "net_diameter");

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
    the_norm <- convert_if_number(this_tuple$social_norm);
    proportion <- convert_if_number(this_tuple$proportion_of_nodes);

    filenames <- Sys.glob(paste(CSV_file_path, "Hes_N_", the_size, "_struct_", structure, "_deg_", degree,  "_risk_*", "_infec_", infection_prob,
                                "_norm_", the_norm,  "_imp_", importation_rate, "_init_", initial_vacc_prop,  "_switch_", random_opinion_switch, "_prop_", proportion,
                                "_summary.csv", sep=""));

    if(length(filenames) == 0){ next; }

    Proportion <- data.table();

    for(file in filenames)
    {
        temp <- data.table(read.csv(file));
        Proportion <- rbind(Proportion, temp[is.na(instance), ], fill=TRUE);
    }
    Proportion[, c("X", "instance"):=NULL];

    Trans_Table <- get_transitions(Proportion, initial_vacc_prop);

    mains <- c("phys_S", "phys_I", "phys_R", "phys_V", "soc_V", "soc_N", "soc_H");
    mains_mean <- sapply(mains, function(x) sprintf("%s_mean",x));
    mains_sd <- sapply(mains, function(x) sprintf("%s_sd",x));

    Proportion[, (mains_sd)  :=lapply(.SD, "/", the_size), .SDcols=mains_sd  ];
    Proportion[, (mains_mean):=lapply(.SD, "/", the_size), .SDcols=mains_mean];

    joins <- c("NN_join_count", "NV_join_count", "VV_join_count", "HH_join_count", "HV_join_count", "HN_join_count");
    joins_mean <- sapply(joins, function(x) sprintf("%s_mean",x));
    joins_sd <- sapply(joins, function(x) sprintf("%s_sd",x));

    # do the sd first before normalising the means, since normalising the means changes to total 1
    # normalising the sd's by one doesn't do fuck all, then
    Proportion[, (joins_sd)  :=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean+HN_join_count_mean+HH_join_count_mean+HV_join_count_mean), .SDcols=joins_sd ];
    Proportion[, (joins_mean):=lapply(.SD, "/", NN_join_count_mean+NV_join_count_mean+VV_join_count_mean+HN_join_count_mean+HH_join_count_mean+HV_join_count_mean), .SDcols=joins_mean];

    # Proportion <- Proportion[abs(num(perceived_vaccine_risk))<=0.2][order(perceived_vaccine_risk)];
    Proportion <- Proportion[abs(num(perceived_vaccine_risk))<=0.1][order(perceived_vaccine_risk)];

    # if there are less than 5 points to be plotted, then cancel this plot and skip to the next parameter set
    if(nrow(Proportion)<5){ next; }

    for(metric in metrics_to_plot)
    {
        if(metric != "plain"){
            for(the_thing in plot_these(metric))
            {
                if(!( the_thing %in% names(metric) )) next;
            }}

        print(metric)

        plots_to_make <- plot_these(metric);

        DT <- data.table(grp=character(), dep=numeric(), val=numeric(), lwr=numeric(), upr=numeric()); # dep means dependent

        for(index in seq_along(plots_to_make))
        {
            prefix <- plot_these(metric)[index];
            tmean <- sprintf("%s_mean", prefix);
            tsd <- sprintf("%s_sd", prefix);

            temp_table <- data.table(
                grp = rep(prefix, nrow(Proportion)),
                dep = Proportion[, perceived_vaccine_risk],
                val = Proportion[, get(tmean)],
                lwr = Proportion[, get(tmean)-get(tsd)],
                upr = Proportion[, get(tmean)+get(tsd)]
            );

            DT <- do.call(smartbind, list(DT, temp_table));
        }

        DT <- data.table(DT);

        getPalette <- colorRampPalette(brewer.pal(8, "Dark2"));
        colours_here <- getPalette(length(plots_to_make));
        names(colours_here) <- plots_to_make;

        pl <- ggplot(data=DT, aes(x=dep, y=val, ymin=lwr, ymax=upr, fill=grp, colour=grp)) + #
            geom_line(size=2) +
            geom_ribbon(alpha=.3, linetype=0) +
            xlab(expression(kappa)) + ylab(axis_label(metric)) +
            theme(
                    legend.position="bottom",
                    panel.background=element_blank(),
                    axis.text=element_text(size=70),
                    axis.title=element_text(size=80),
                    axis.title.y=element_text(vjust=4) # hjust=-10,
            ) +
            scale_colour_discrete(labels=as.vector(sapply(plots_to_make, axis_label))) +
            scale_fill_discrete(labels=as.vector(sapply(plots_to_make, axis_label))) +
            theme(
                legend.title=element_blank(),
                legend.text=element_text(size=60),
                legend.key.width=unit(4, 'cm'),
                legend.spacing.x=unit(1, 'cm'),
                legend.direction="horizontal"
            )
            guides(colour = guide_legend(nrow=1, parse=TRUE));
            # guides(fill = guide_legend(nrow=1, parse=TRUE));


        if(isTRUE( Trans_Table$is_soc_trans ))
        {
            pl <-  pl + geom_vline(xintercept=Trans_Table[, soc_intersec_risk], size=1, linetype="dotted");
        }
        if(isTRUE( Trans_Table$is_phys_trans ))
        {
            pl <-  pl + geom_vline(xintercept=Trans_Table[, phys_intersec_risk], size=1, linetype="dashed");
        }

        plot_file_name <- paste(graph_file_path, folder_name, "/Hes_Kappa_", metric, "_N_", the_size, "_struct_", structure, "_deg_", degree,
                                "_infec_", infection_prob, "_norm_", the_norm,  "_imp_", importation_rate, "_init_", initial_vacc_prop,
                                "_switch_", random_opinion_switch, "_prop_", proportion, ".png", sep="");

        ggsave(
            plot_file_name,
            plot=pl, width=image_width, height=image_height, units="in", limitsize=FALSE #  height=5,
        )
        dev.off();

        # graphics.off(); stop();

    }

}

print(Sys.time()-time_start);
