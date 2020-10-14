
for(line_index in 1:nrow(Params))
{


    for(CHANGE_TEST in c("snh", "pettitt", "buishand_r", "lanzante")) #
    {

       for(i in the_names_to_group_together)
       {
           ordered_rows <- c(ordered_rows, which(grepl(i, Joint$EWS))); # group all the like metrics together
       }

        Joint <- Joint[ordered_rows]; # put the rows in order
        Joint[, Position:=1:nrow(Joint)]; # put in a row that preserves this order, so we can use the reorder() command
        Joint[, Minima:=-(Joint$Minima)];

        min_name <- sprintf("min(%s)", axis_label(CHANGE_TEST));
        max_name <- sprintf("max(%s)", axis_label(CHANGE_TEST));

        names(Joint)[which(names(Joint)=="Maxima")] <- max_name;
        names(Joint)[which(names(Joint)=="Minima")] <- min_name;

        Joint <- melt(Joint, id.vars=c("EWS", "Position"), measure.vars=c(min_name, max_name));

        number_of_bars <- length(unique(Joint$EWS));
        if(old_or_new == "new"){ image_height <- length(unique(Joint$EWS))/3.5; }
        if(old_or_new == "old"){ image_height <- length(unique(Joint$EWS))/2; }

        pl_performance <- ggplot(Joint, aes(x=reorder(EWS, Position), y=value, fill=variable)) + #, color=rep(rainbow(length(metrics)-1), 2))
            geom_bar(stat = 'identity') +
            # scale_discrete_manual("point_colour", values=getPalette(number_of_bars), labels=as.vector(sapply(Joint$EWS, axis_label))) +
            facet_share(~variable, dir='h', scales='free', reverse_num=TRUE) +
            scale_x_discrete(labels=as.vector(sapply(Joint$EWS, axis_label))) +
            coord_flip() +
            xlab("Early Warning Signal") + ylab(bquote("Ratio of social norm ("*sigma*") values")) +
            geom_text(aes(label=sprintf("%i%%", round(100*abs(value), 0))), size=5, colour="black") + #, hjust=1
            # geom_hline(yintercept=c(0.25), linetype="dashed", colour="black") +
            theme(
                axis.text=element_text(size=text_size-2),
                axis.title.y=element_blank(),
                axis.title.x=element_text(size=text_size),
                panel.grid.major=element_line(size = 0.50, linetype = 'solid', colour = 'grey'),
                panel.grid.minor=element_line(size = 0.25, linetype = 'solid', colour = 'grey'),
                legend.position="none"
            );

        graph_name <- paste(graph_file_path, folder_name, "/", "LeadTime_performance_bars", sep='');
        if(the_size == 10000){ graph_name <- paste(graph_name, "_", old_or_new, sep=''); }
        graph_name <- paste(graph_name, "_", CHANGE_TEST, file_name_ending, sep='');

        ggsave(
            graph_name, plot = pl_performance, width=image_width, height=image_height, units="in", limitsize=FALSE
        );
        dev.off();

        system(paste("convert ", graph_name, " -trim ", graph_name, sep=""));

        # stop();
    }
}

graphics.off()

print(Sys.time()-start_time);
