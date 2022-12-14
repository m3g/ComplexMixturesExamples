#
# This script will run all examples with the complete trajectories, # to test the outputs.
#
# Normally, run the script and save the results to a new branch
# of the repository. Checking the differences to the main branch allows
# visualization of the differences in the resulting images.
#
using ComplexMixtures
script_dir = @__DIR__

#
# Download trajectory files
#
found_all_downloads = isfile("$script_dir/../Test/trajectories/traj_Glyc.dcd")
found_all_downloads = isfile("$script_dir/../Test/trajectories/traj_POPC.dcd")
found_all_downloads = isfile("$script_dir/../Test/trajectories/traj_Polyacry.dcd")
found_all_downloads = isfile("$script_dir/../Test/trajectories/glyc50.dcd")

if !found_all_downloads
    error("""

    Missing full trajectory files at $script_dir/trajectories
    
    Download file (200Mb): https://drive.google.com/file/d/1BuXJ8AjBeduMSD2CkDJLDNxAAD2QNNg6/view?usp=sharing
    Download file (280Mb): https://drive.google.com/file/d/1ug43ncCLsBATaJrT9zlbaqK6AORVvhhx/view?usp=sharing
    Download file (365Mb): https://drive.google.com/file/d/12TT5tblkFp1NtFOAQgjjGhmnYaXA8vQi/view?usp=sharing
    Download file (3Gb): https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing

    """)
end

