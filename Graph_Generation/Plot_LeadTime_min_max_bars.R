# Dear God, why have I retooled this code so many fucking times???

# Summary Plot for the Topology
rm(list=ls());

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Lead_Time", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

Params <- fread(parameter_file, colClasses="character")
Params[, c("V1", "risk", "norm"):=NULL];

start_time <- Sys.time();
#
for(the_pair in c( list(c(10000, "new")) )) # list(c(10000, "new")), list(c(10000, "old")),list(c(40000, "old")),
{
    the_size <- num(the_pair[[1]])
    old_or_new <- the_pair[[2]]

    Parameters <- unique(Params[phys_top=="random" & num(init_prop)==0.05 & num(beta)==1 & num(size)==the_size]);

    if(num(the_size)==10000){ Parameters <- Parameters[num(import)==0.00025]; }
    else if(num(the_size)==40000){ Parameters <- Parameters[num(infec)==0.2]; }
    else if(num(the_size)==562500){ Parameters <- Parameters[]; }

    Parameters <- unique(Parameters);

    for(line_index in 1:nrow(Parameters))
    {
        this_tuple <- Parameters[line_index];

        initial_vacc_prop <- this_tuple[, init_prop]; size <- this_tuple[, size]; beta <- this_tuple[, beta]; duration <- this_tuple[, duration];
        infection_prob <- this_tuple[, infec]; importation_rate <- this_tuple[, import]; replenishment_rate <-  this_tuple[, birth_death];
        random_opinion_switch <- this_tuple[, random_switch]; physical_topology <- this_tuple[, phys_top]; physical_degree <- this_tuple[, phys_deg];
        social_topology <- this_tuple[, soc_top]; social_degree <- this_tuple[, soc_deg];

        filenames <- Sys.glob(paste(CSV_file_path, "Summary", "_N_", size, "_dur_", duration, "_beta_", beta, "_vaccprop_", initial_vacc_prop, "_risk_", "*", "_inf_",
                                    infection_prob, "_imp_", importation_rate, "_rep_", replenishment_rate, "_ptop_", physical_topology, "_pdeg_", physical_degree, "_stop_",
                                    social_topology, "_sdeg_", social_degree, "_norm_", "*", "_switch_", random_opinion_switch, ".csv", sep=""));

        if(length(filenames) < 10){ next;   }

        writeLines("\n"); writeLines(paste("\treplenishment ratio: ", replenishment_rate, sep="")); writeLines(paste("\timportation rate: ", importation_rate, sep=""));
        writeLines(paste("\tnetwork size: ", size, sep="")); writeLines(paste("\tinfection probability: ", infection_prob, sep=""));
        writeLines(paste("\tphysical topology: ", physical_topology, sep=""));  writeLines(paste("\tphysical degree: ", physical_degree, sep=""));
        writeLines(paste("\tsocial topology: ", social_topology, sep="")); writeLines(paste("\tsocial degree: ", social_degree, sep=""));
        writeLines(paste("\topinion switching rate: ", random_opinion_switch, sep=""));

        Prop_vs_Risk <- data.table();
        for(file in filenames){ temp <- data.table(read.csv(file)); Prop_vs_Risk <- rbind(Prop_vs_Risk, temp[is.na(instance)], fill=TRUE); }
        Prop_vs_Risk[, ("X"):=NULL];
        write.csv(Prop_vs_Risk, "temp_risk.csv");

        # Prop_vs_Risk <- data.table(read.csv("temp_risk.csv")); Prop_vs_Risk[, X:=NULL];

        metrics <- c();
        if(the_size == 40000)
        {
            metrics <- filter_these(names(Prop_vs_Risk), no=c("sd", fixed_values[fixed_values!="N"], "dist", "chamber", "conn", "mod", "watts", "ratio", "prob",
                                                              "opinion", "phys", "soc", "getis"));
        }
        else if(the_size == 10000)
        {
            if(old_or_new == "old"){
                metrics <- filter_these(names(Prop_vs_Risk), no=c("sd", fixed_values[fixed_values!="N"], "dist", "chamber", "conn", "mod", "watts", "ratio", "prob",
                                                                  "opinion", "phys", "soc", "getis"));
            } else {
                metrics <- filter_these(names(Prop_vs_Risk), no=c(fixed_values, "sd", "ratio", "moran", "getis", "join", "geary", "mutual", "phys", "soc"));
            }
        }
        metrics <- metrics[metrics!="N"];

        ordered_rows <- numeric(); # order the table so that the bar chart places related metrics together instead of spreading them around the graph
        the_names_to_group_together <- c();
        if(size == "10000"){
            if(old_or_new == "old"){
                the_names_to_group_together <- c("join.count", "phys", "soc", "moran", "geary", "mutual");
            }
            else{
                the_names_to_group_together <- c("phys", "soc", "watts", "avg.chamber.size", "min.chamber.size", "max.chamber.size", "number.chambers",
                                                 "avg.conn.comp", "min.conn.comp", "max.conn.comp", "number.conn", "prob", "modularity", "opinion");
            }
        } else if(size == "40000"){
            the_names_to_group_together <- c("join.count", "phys", "soc", "moran", "geary", "mutual");
        }

        col_names <- c("norm", "phys_trans", "soc_trans", names.change(metrics));

        for(CHANGE_TEST in c("snh", "lanzante", "pettitt", "buishand_r")) #
        {
            writeLines(sprintf("\t\t %s test", CHANGE_TEST));
            file_name_ending <- paste("_N_", size, "_dur_", duration, "_beta_", beta, "_vaccprop_", initial_vacc_prop, "_inf_", infection_prob,
                                      "_imp_", importation_rate, "_rep_", replenishment_rate, "_top_", physical_topology, "_deg_", physical_degree,
                                      "_switch_", random_opinion_switch, ".png", sep=""
            );

            Change_Predictions <- data.table(matrix(ncol=length(col_names)));
            names(Change_Predictions) <- col_names;
            for(norm_index in seq_along(Prop_vs_Risk[, unique(social_norm)]))
            {
                soc_norm <- Prop_vs_Risk[, sort(unique(social_norm))][norm_index];
                Prop_Table_Here <- Prop_vs_Risk[social_norm==soc_norm][order(perceived_vaccine_risk)];
                transitions <- get_transitions(Prop_Table_Here, initial_vacc_prop);
                change_points <- c();
                for(metric in metrics) change_points[[metric]] <- change_point_test(Prop_Table_Here, metric, CHANGE_TEST)$kappa;
                Change_Predictions <- rbind(
                    Change_Predictions,
                    as.list(c(
                        soc_norm, transitions$phys_intersec_risk, transitions$soc_intersec_risk,
                        data.table(change_points)[, change_points]
                    ))
                );
            }
            Change_Predictions <- Change_Predictions[rowSums(is.na(Change_Predictions)) != ncol(Change_Predictions)][order(norm)];
            write.csv(Change_Predictions, "temp_change_predictions.csv");

            # Change_Predictions <- data.table(read.csv("temp_change_predictions.csv")); Change_Predictions[, X:=NULL];

            if(the_size == 40000){ Change_Predictions <- Change_Predictions[norm<=2.625]; }

            All_Lead_Times <- Change_Predictions$soc_trans - Change_Predictions[, .SD, .SDcols=names.change(metrics)];
            All_Lead_Times <- cbind(Change_Predictions$norm, All_Lead_Times); names(All_Lead_Times)[1] <- "norm";
            # All_Lead_Times <- All_Lead_Times %>% select_if(function(x) any(!is.na(x)))

            Max_or_Min <- data.table(matrix(nrow=nrow(All_Lead_Times), ncol=ncol(All_Lead_Times)-1, "neither"));
            Max_or_Min <- cbind(Change_Predictions$norm, Max_or_Min);
            names(Max_or_Min) <- names(All_Lead_Times);

            for(num_norm in unique(All_Lead_Times$norm)) { for(name in setdiff(names(All_Lead_Times), "norm")) {
                if(is.na( All_Lead_Times[norm==num_norm, get(name)] )){ Max_or_Min[which(norm==num_norm), (name):="none"]; }
                else if(All_Lead_Times[norm==num_norm, get(name)]<0)
                {
                    Max_or_Min[which(norm==num_norm), (name):="fail"];
                    All_Lead_Times[which(norm==num_norm), (name):=NA];
                }
            }}

            for(i in 1:nrow(All_Lead_Times))
            {
                rows_with_max <- which(All_Lead_Times[i, ] == max(All_Lead_Times[i, -c("norm")], na.rm=TRUE));
                for(here in rows_with_max)
                {
                    if(Max_or_Min[i, ..here][[1]] == "neither"){ Max_or_Min[i, here] <- "max"; }
                }
                rows_with_min <- which(All_Lead_Times[i, ] == min(All_Lead_Times[i, -c("norm")], na.rm=TRUE));
                for(here in rows_with_min)
                {
                    if(Max_or_Min[i, ..here][[1]] == "max"){ Max_or_Min[i, here] <- "both"; }
                    else if(Max_or_Min[i, ..here][[1]] == "neither"){ Max_or_Min[i, here] <- "min"; }
                }
            }

            grid_image_width <- 30;
            grid_text_size <- 40;
            grid_title_size <-50;

            Max_Min_Grid <- Max_or_Min %>% gather(col_name, value, -norm) %>% ggplot(aes(factor(norm), col_name, fill=value)) +
                geom_tile(colour = 'black') +
                labs(x=expression(sigma), y=test_name(CHANGE_TEST, old_or_new)) +
                scale_x_discrete(breaks=seq(min(Max_or_Min$norm), max(Max_or_Min$norm), by=.25)) +
                scale_y_discrete(labels=as.vector(sapply(setdiff(names(Max_or_Min), "norm"), axis_label))) +
                scale_fill_manual(values = c('max'='green', 'neither'='grey', 'min'='red', 'both'='yellow', 'none'='white', 'fail'='black')) +
                theme(
                    legend.position="none",
                    axis.text.y=element_text(size=grid_title_size)
                );

            if(the_size == 40000){
                Max_Min_Grid <- Max_Min_Grid +
                    theme(
                        axis.text.x=element_text(size=grid_text_size),
                        axis.title=element_text(size=grid_title_size+5)
                );
                grid_height <- ncol(Max_or_Min)/2;
                grid_width <- grid_image_width;
            } else if(the_size == 10000){
                if(old_or_new == "new"){
                    Max_Min_Grid <- Max_Min_Grid +
                        theme(
                            axis.text.x=element_text(size=grid_text_size),
                            axis.text.y=element_text(size=grid_text_size),
                            axis.title=element_text(size=grid_title_size)
                        );
                    grid_height <- ncol(Max_or_Min)/1.8;
                    grid_width <- ncol(Max_or_Min);
                }
                else if(old_or_new == "old"){
                    Max_Min_Grid <- Max_Min_Grid +
                        theme(
                            axis.text.x=element_text(size=grid_text_size),
                            axis.text.y=element_text(size=grid_text_size),
                            axis.title=element_text(size=grid_title_size)
                        );
                    grid_height <- ncol(Max_or_Min)/3;
                    grid_width <- grid_image_width;
                }
            }

           grid_name <- paste(graph_file_path, folder_name, "/", "LeadTime_maxmin_grid", sep='');
           if(the_size == 10000){ grid_name <- paste(grid_name, "_", old_or_new, sep=''); }
           grid_name <- paste(grid_name, "_", CHANGE_TEST, file_name_ending, sep='');

            ggsave(
                grid_name, plot = Max_Min_Grid, width=grid_width, height=grid_height, limitsize=FALSE, dpi=50
            );
            dev.off();
            system(paste("convert ", grid_name, " -trim ", grid_name, sep=""));

            Joint <- data.table(EWS=character(), Maxima=numeric(), Minima=numeric());
            for(row_name in setdiff(names(Max_or_Min), "norm"))
            {
                Joint <- rbind(Joint,
                   list(
                        row_name,
                        length(which(Max_or_Min[, get(row_name)] == "max"))/nrow(Max_or_Min),
                        length(which(!(Max_or_Min[, get(row_name)] %in% c("max", "neither"))))/nrow(Max_or_Min)
                   )
                );
            }

            for(i in the_names_to_group_together)
            {
               ordered_rows <- c(ordered_rows, which(grepl(i, Joint$EWS))); # group all the like metrics together
            }

            Joint <- unique(Joint[ordered_rows]); # put the rows in order
            Joint[, Position:=1:nrow(Joint)]; # put in a row that preserves this order, so we can use the reorder() command

            if(the_size == 40000){
                bars_label_size <- 5;
                bars_text_size <- 40;
            } else if(the_size == 10000){
                if(old_or_new == "old"){
                   bars_label_size <- 7;
                   bars_text_size <- 40;
                } else if(old_or_new == "new") {
                   bars_label_size <- 15;
                   bars_text_size <- 45;
                }
            }

            bars_title_size <- 35;

            pl_max <- ggplot(Joint, aes(x=reorder(EWS, Position), y=Maxima)) +
                geom_bar(fill='green', stat='identity') +
                labs(title="maximal lead times")  +
                scale_x_discrete(labels=as.vector(sapply(Joint$EWS, axis_label)))  +
                geom_text(aes(label=sprintf("%i%%", round(100*abs(Maxima), 0))), size=bars_label_size, colour="black", hjust=-0.01) +
                coord_flip(clip="off") +
                theme(
                    legend.position='none',
                    plot.title=element_text(size=bars_title_size),
                    axis.text.y=element_text(size=bars_text_size, hjust=0.5),
                    axis.text.x=element_text(size=bars_text_size),
                    axis.ticks.y=element_blank(),
                    axis.line.y=element_blank(),
                    axis.title=element_blank(),
                    panel.grid.major=element_line(size = 0.50, linetype = 'solid', colour = 'grey'),
                    panel.grid.minor=element_line(size = 0.25, linetype = 'solid', colour = 'grey'),
                    plot.margin = margin(0, 2.1, 0, 0, "cm"),
                );

            pl_min <- ggplot(Joint, aes(x=reorder(EWS, Position), y=Minima)) +
                geom_col(fill='red') +
                labs(title=sprintf("minimal, undefined, negative lead times", axis_label(CHANGE_TEST))) +
                scale_y_continuous('', trans='reverse') +
                geom_text(aes(label=sprintf("%i%%", round(100*abs(Minima), 0))), size=bars_label_size, colour="black", hjust=1) +
                coord_flip(clip="off") +
                # expand_limits(x = 0.3) + # or some other arbitrarily large number
                theme(
                    legend.position='none',
                    plot.title=element_text(size=bars_title_size),
                    axis.text.x=element_text(size=bars_text_size),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.line.y=element_blank(),
                    axis.title=element_blank(),
                    panel.grid.major=element_line(size = 0.50, linetype = 'solid', colour = 'grey'),
                    panel.grid.minor=element_line(size = 0.25, linetype = 'solid', colour = 'grey'),
                    plot.margin = margin(0, 0, 0, 2.8, "cm"),
                );

            pl_performance <- grid.arrange(
                pl_min,
                pl_max,
                widths=c(0.46, 0.54),
                ncol=2,
                bottom=expression("Proportion of social("*sigma*") norm values")
            );

            # bars_image_height <- 30;
            bars_image_width <- 20;

            graph_name <- paste(graph_file_path, folder_name, "/", "LeadTime_performance_bars", sep='');
            if(the_size == 10000){ graph_name <- paste(graph_name, "_", old_or_new, sep=''); }
            graph_name <- paste(graph_name, "_", CHANGE_TEST, file_name_ending, sep='');

            number_of_bars <- length(unique(Joint$EWS));
            if(old_or_new == "new"){ bars_image_height <- length(unique(Joint$EWS))/1; }
            if(old_or_new == "old"){ bars_image_height <- length(unique(Joint$EWS))/2; }

            ggsave(
                graph_name, plot = pl_performance, width=bars_image_width, height=bars_image_height, limitsize=FALSE, dpi=50
            );
            dev.off();

            system(paste("convert ", graph_name, " -trim ", graph_name, sep=""));

#          # stop(Sys.time() - start_time);
        }
    }
}

graphics.off()

print(Sys.time()-start_time);
#
# good_vector <- c()
# bad_vector <- c()
#
# for(i in 1:nrow(Max_or_Min))
# {
#    good_vector <- c(good_vector,  names(Max_or_Min)[which(Max_or_Min[i] == "max")]);
#    bad_vector <- c(bad_vector, names(Max_or_Min)[which(!(Max_or_Min %in% c("neither", "max")))]);
# }
#
# pl <- ggplot(good_vector) + geom_bar()
