# Summary Plot for the Topology
rm(list=ls())

user_name <- Sys.info()[8][[1]];
source(sprintf("/home/%s/Dropbox/Processing/Parameter_Values.R", user_name));

folder_name <- paste("Heatmaps", sep="");
dir.create(file.path(graph_file_path, folder_name), showWarnings=FALSE);

time_start <- Sys.time();
Params <- fread(parameter_file, colClasses="character");
text_size = 40;

for(the_size in c("40000")) # "40000", "562500"))
{
  Parameters <- Params[num(size)==the_size & num(beta)==1];

  if(the_size == "10000") Parameters <- Parameters[num(init_prop)==0.05 & num(beta)==1 & num(infec)==0.2 & num(import)==0.00025];
  if(the_size == "40000") Parameters <- Parameters[num(init_prop)==0.05 & num(beta)==1 & num(infec)==0.2];
  if(the_size == "562500") Parameters <- Parameters[num(init_prop)==0.05 & num(beta)==1 & soc_top=="random" & num(random_switch)==0.1];

  Parameters[, c("V1"):=NULL];
  Parameters <- unique(Parameters);

  data_file_path <- sprintf("%s/N_%s_cleaned/", stored_data_path, the_size);

  # Summary_Data <- data.table(matrix(ncol=4, nrow=length(list.files(data_file_path)), 0));
  # names(Summary_Data) <- c("instance", "risk", "norm", "time");

  Summary_Data <- data.table(matrix(ncol=3, nrow=nrow(Parameters), 0));
  names(Summary_Data) <- c("risk", "norm", "time");

  DT_index <- 1;

  for(line_index in 1:nrow(Parameters))
  {
    progress(line_index, max=nrow(Parameters), progress.bar=TRUE)
    
    this_tuple <- Parameters[line_index];

    beta <- this_tuple[, beta];
    initial_vacc_prop <- this_tuple[, init_prop];
    size <- this_tuple[, size];
    duration <- this_tuple[, duration];
    risk <- this_tuple[, risk];
    infection_prob <- this_tuple[, infec];
    importation_rate <- this_tuple[, import];
    replenishment_rate <-  this_tuple[, birth_death];
    physical_topology <- this_tuple[, phys_top];
    physical_degree <- this_tuple[, phys_deg];
    social_topology <- this_tuple[, soc_top];
    social_degree <- this_tuple[, soc_deg];
    random_opinion_switch <- this_tuple[, random_switch];
    social_norm <- this_tuple[, norm];

    data_file_list <- Sys.glob(
      paste(
        data_file_path,
        "stats*",
        "_N_", size,
        "_dur_", duration,
        "_beta_", beta,
        "_vaccprop_", initial_vacc_prop,
        "_pay_", risk,
        "_inf_", infection_prob,
        "_imp_", importation_rate,
        "_rep_", replenishment_rate,
        "_ptop_", physical_topology,
        "_pdeg_", physical_degree,
        "_stop_", social_topology,
        "_sdeg_", social_degree,
        "_norm_", social_norm,
        "_switch_", random_opinion_switch,
        "_inst*.bin",
        sep=""
      )
    );

    evolution_file_stub = paste(
      "_EVO",
      "_N_", size,
      "_dur_", duration,
      "_beta_", beta,
      "_vaccprop_", initial_vacc_prop,
      "_risk_", risk,
      "_in_", infection_prob,
      "_im_", importation_rate,
      "_rep_", replenishment_rate,
      "_pt_", physical_topology,
      "_pd_", physical_degree,
      "_st_", social_topology,
      "_sd_", social_degree,
      "_sn_", social_norm,
      "_sw_", random_opinion_switch,
      ".png",
      sep=""
    );

    if(length(data_file_list) == 0){ next; }

    list_of_times <- numeric();
    for(file_index in 1:length(data_file_list))
    {
      if(file.size(data_file_list[file_index]) == 0) { next; }
      temp <- fread(data_file_list[file_index]);
      if(nrow(temp) == 0){ next; }
      list_of_times <- c(list_of_times, max(temp$time));
      # Summary_Data[DT_index,] <- list(unique(temp$instance), unique(temp$perceived_vaccine_risk), unique(temp$social_norm), max(temp$time));
      # DT_index <- DT_index+1;
    }
    Summary_Data[DT_index,] <- list(unique(temp$perceived_vaccine_risk), unique(temp$social_norm), mean(list_of_times));
    DT_index <- DT_index + 1;
  }
  
  stop();

  Summary_Data <- Summary_Data[time!=0 & risk>=-1 & risk <=1 & norm<=3];
  write.csv(Summary_Data, sprintf("convergenge_data_size_%s_averaged.csv", the_size));

  all_the_norms <- sort(unique(Summary_Data$norm));
  all_the_risks <- sort(unique(Summary_Data$risk));

  DM <- matrix(nrow=length(all_the_norms), ncol=length(all_the_risks));
  rownames(DM) <- all_the_norms;
  colnames(DM) <- all_the_risks;

  for(temp_norm in all_the_norms){
    for(temp_risk in all_the_risks){
      the_time_taken <- Summary_Data[risk==temp_risk & norm==temp_norm, time];
      if(length(the_time_taken)) DM[which(rownames(DM)==temp_norm), which(colnames(DM)==temp_risk)] <- log(the_time_taken);
    }}

  time_end <- Sys.time();
  print(time_end - time_start);

  the_title <- list(family = "sans-serif", size=text_size, color = "black");
  the_tick  <- list(family = "sans-serif", size=text_size, color = "black");

  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  x <- list(
    title = "x Axis",
    titlefont = f
  )
  y <- list(
    title = "y Axis",
    titlefont = f
  )
  
  fig <- plot_ly(z = ~DM) %>% add_surface(
    contours = list(
      z = list(
        show=TRUE,
        usecolormap=TRUE,
        highlightcolor="#ff0000",
        project=list(z=TRUE)
      )
    )
  );
  fig <- fig %>% layout(
    # xaxis = list(showticklabels = TRUE, tickfont = the_tick, titlefont = the_title, title = "σ"),
    # yaxis = list(showticklabels = TRUE, tickfont = the_tick, titlefont = the_title, title = "κ")
    # xaxis = x,
    # yaxis = y
    scene = list(
      camera=list(
        eye = list(x=1.87, y=0.88, z=-0.64)
      )
    )
  ); 
  fig <- fig %>% colorbar(size=60, len=0.9);
  
  stop()
  
}
