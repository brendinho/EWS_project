# Brendon Phillips
# PhD candidate
# Department of Applied Mathematics
# University of Waterloo

using DataFrames, CSV, ProgressMeter

Structures = ["smallworld"]
Sizes = [10000];
Degrees = [50];
Ensemble = [1:5;];
Infectivity = [0.2];
Risk = [-1:0.03125:1;];
Birth_Death = [2.4e-4];
Importation = [2.5e-5];
Norm = [0:0.125:2.5;];
Init_Vacc_Prop = [0.25];
Random_Switch = [0.0001];
Tested_Proportion = [0.8];

done = CSV.read("/home/b2philli/Dropbox/CSV/Hes_Parameters_instances.csv");
select!(done, Not(:Column1));

Unfinished = DataFrame(N=Int[], structure=String[], degree=Int[], risk=Float64[], infectivity=Float64[], norm=Float64[], importation=Float64[], vaccprop=Float64[], switch=Float64[], proportion=Float64[], instance=Int[]);

possible_combinations = Base.Iterators.product(Sizes, Structures, Degrees, Risk, Infectivity, Norm, Importation, Init_Vacc_Prop, Random_Switch, Tested_Proportion, Ensemble) |> collect; # Birth_Death,

@showprogress 1 "Checking..." for index in 1:length(possible_combinations)
	combo = possible_combinations[index];
	trials = filter(i->(collect(combo) == Array(done[i,:])), 1:nrow(done));
	if length(trials)==0
		push!(Unfinished, combo);
	end
end
sort!(Unfinished, order(:risk, by=abs));

CSV.write("/home/b2philli/Dropbox/Processing/jobs_to_be_done.csv", Unfinished);

outfile = open("/home/b2philli/Dropbox/Processing/job_table_unfinished.txt", "w")

for index in 1:nrow(Unfinished)
	row = Unfinished[index, :];
	write(outfile, "$(row.structure) $(row.N) $(row.degree) $(row.instance) $(row.infectivity) $(row.risk) 2.4e-4 $(row.importation) $(row.norm) $(row.vaccprop) $(row.switch) $(row.proportion)\n");
end

# write(outfile, param_string)
close(outfile)
