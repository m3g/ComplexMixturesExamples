#
# This script will run all examples with the complete trajectories, # to test the outputs.
#
# Normally, run the script and save the results to a new branch
# of the repository. Checking the differences to the main branch allows
# visualization of the differences in the resulting images.
#
test_dir = @__DIR__

#
# Download trajectory files
#
found_all_downloads = isfile("$test_dir/../Test/trajectories/traj_Glyc.dcd")
found_all_downloads = isfile("$test_dir/../Test/trajectories/traj_POPC.dcd")
found_all_downloads = isfile("$test_dir/../Test/trajectories/traj_Polyacry.dcd")
found_all_downloads = isfile("$test_dir/../Test/trajectories/glyc50.dcd")

if !found_all_downloads
    error("""

    Missing full trajectory files at $test_dir/trajectories
    
    Download file (200Mb): https://drive.google.com/file/d/1BuXJ8AjBeduMSD2CkDJLDNxAAD2QNNg6/view?usp=sharing
    Download file (280Mb): https://drive.google.com/file/d/1ug43ncCLsBATaJrT9zlbaqK6AORVvhhx/view?usp=sharing
    Download file (365Mb): https://drive.google.com/file/d/12TT5tblkFp1NtFOAQgjjGhmnYaXA8vQi/view?usp=sharing
    Download file (3Gb): https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing

    """)
end

# Glycerol in water
dir = "$test_dir/../Glycerol_in_Water"
include("$dir/compute_mddf.jl")
include("$dir/plots.jl")

# Polyacrylamide in DMF
dir = "$test_dir/../Polyacrylamide_in_DMF"
include("$dir/compute_mddf.jl")
include("$dir/plots.jl")

# POPC in Water/Ethanol
dir = "$test_dir/../POPC_in_Water-Ethanol"
include("$dir/compute_mddf.jl")
include("$dir/plots.jl")

# Protein in Glycerol 
dir = "$test_dir/../Protein_in_Glycerol"
include("$dir/MDDF/mddf.jl")
include("$dir/Density2D/density2D.jl")
include("$dir/Density3D/density3D.jl")
cd("$dir/Density3D")
run(`vmd -dispdev text -e grid_render.vmd`)
run(`convert grid.tga grid.png`)

cd(test_dir)





