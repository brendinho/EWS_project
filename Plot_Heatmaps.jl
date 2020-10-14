using Plots, DataFrames, CSV, Query, Glob, Printf, LaTeXStrings

folder_name = "Heatmaps"
full_folder_name = "/home/b2philli/Dropbox/Apps/Overleaf/phd_thesis/Thesis_Graphs/"*folder_name*"/"
CSV_file_path = "/home/b2philli/Dropbox/CSV/"
parameter_file = "/home/b2philli/Dropbox/CSV/params.csv"

isdir(full_folder_name) ||  mkdir(full_folder_name)

metrics_to_plot = ["soc_V_mean", "phys_V_mean"] # "watts_strogatz_vacc_mean"] #"vacc_avg_conn_comp_size_mean", "vacc_echo_chambers"); #
# metrics_to_plot = c("soc_V_mean")
# "phys_V_mean", "soc_V_mean", "NV_join_count_mean", "mutual_info_mean") # "geary_c_mean") #, , "NN_join_count_mean", "VV_join_count_mean", "moran_i_mean", "ratio_sec_to_all_infs_at_the_end_mean") # ) # ) # , ); #, "join_counts", "nonvacc_echo_chambers", "mutual_info", "conn_chamb_numbers", "nonvacc_conn_comps", "vacc_conn_comps"); #, "modularity", "watts_strogatz",  "geary_c", "getis_ord", "infections", "opinion");

Params = CSV.read(parameter_file);

for the_size in [40000]

	global Parameters = Params[(Params.size .== the_size) .& (Params.beta .== 1) .& (Params.init_prop .== 0.05), :]; #

	if the_size == 10000
		Parameters = Parameters[(Parameters.infec .== 0.2) .& (Parameters.import .== 0.00025), :]
	elseif the_size == 40000
		Parameters = Parameters[(Parameters.infec .== 0.2), :]
	elseif the_size == 562500
		Parameters = Parameters[(Parameters.soc_top .== "random") & (Parameters.random_switch .== 0.1), :]
	end

	Parameters = unique(select(Parameters, Not([:risk, :norm, :Column1, :time])))

	for this_tuple in eachrow(Parameters)

		global DT = DataFrame()

		global file_name =join([
			"Summary_N_", this_tuple.size,
			"_dur_", this_tuple.duration,
			"_beta_", @sprintf("%.6f", this_tuple.beta),
			"_vaccprop_", this_tuple.init_prop,
			"_risk_", "*",
			"_inf_", this_tuple.infec,
			# "_imp_", this_tuple.import,
			"_imp_", @sprintf("%.1e", this_tuple.import),
			"_rep_", this_tuple.birth_death,
			"_ptop_", this_tuple.phys_top,
			"_pdeg_", this_tuple. phys_deg,
			"_stop_", this_tuple.soc_top,
			"_sdeg_", this_tuple.soc_deg,
			"_norm_", "*",
			"_switch_", @sprintf("%.0e", this_tuple.random_switch),
			".csv",
		])

		global filelist = glob(file_name, CSV_file_path);

		for the_file in filelist
			global temp = CSV.read(the_file)
			# DT = join(DT, temp[(temp.instance .== "NA") .| isna(temp.instance), :], kind=:outer, on=intersect(names(DT), names(temp)))
			# DT = [ DT ; temp[(temp.instance .== "NA") .| isna(temp.instance), :] ]
			# DT = vcat(DT, temp[(temp.instance .== "NA") .| isna(temp.instance), :])
			append!(DT, temp[(temp.instance .== "NA") .| isna(temp.instance), :])
		end

		global norms = sort(unique(DT.social_norm))
		global risks = sort(unique(DT.perceived_vaccine_risk))

		break

		for metric in metrics_to_plot

			println(metric)

			global Mat = zeros(length(norms), length(risks))
			for i=1:length(norms), j=1:length(risks)
				global element = DT[(DT.social_norm .== norms[i]) .& (DT.perceived_vaccine_risk .== risks[j]), Symbol(metric)][1]
				Mat[i,j] = (typeof(element) == String) ? parse(Float64, element) : element
			end

			global file_name_ending = join([
				"_N_", this_tuple.size,
				"_vacc_", this_tuple.init_prop,
				"_dur_", this_tuple.duration,
				"_inf_", this_tuple.infec,
				"_imp_", this_tuple.import,
				"_rep_", this_tuple.birth_death,
				"_ptop_", this_tuple.phys_top,
				"_pdeg_", this_tuple.phys_deg,
				"_stop_", this_tuple.soc_top,
				"_sdeg_", this_tuple.soc_deg,
				"_switch_", this_tuple.random_switch,
				"_julia.png"
			])

			println(join([full_folder_name, "Contour_", metric, file_name_ending]))

			# png(heatmap(risks, norms, Mat, xlabel=L"\sigma", ylabel=L"\kappa", xtickfont=font(36, "Courier")), join([full_folder_name, "Contour_", metric, file_name_ending]))

			# "/home/b2philli/Dropbox/trial.png")
			break
		end
		break
	end
	break
end
