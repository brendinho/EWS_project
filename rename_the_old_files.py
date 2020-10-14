# Brendon Phillips
# Ph.D. candidate
# Department of Applied Mathematics
# Faculty of Mathematics
# University of Waterloo

import os
import glob

# all_the_files = os.listdir("/media/bren/Research_Files/All_Data");
all_the_files = glob.glob("/media/bren/Research_Files/All_Data/*vacc_true*duration*bin");

for filename in [all_the_files[1]]:

	elements = filename.split("_");

	size = elements[elements.index('N')+1];
	duration = elements[elements.index("duration")+1];
	latency = elements[elements.index("latency")+1];




# print(len(all_the_files))
