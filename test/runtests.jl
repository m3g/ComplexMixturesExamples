using TestItemRunner
@run_package_tests

@testitem  "Glycerol in water" begin
    test_dir = @__DIR__
    if !isfile("$test_dir/trajectories/traj_Glyc.dcd")
        println("""
    
        WARNING: Will run tests without recomputing the data: only test trajectory files were found for all or some examples.
    
        Missing full trajectory files at $test_dir/trajectories
        
        Download file (200Mb): https://drive.google.com/file/d/1BuXJ8AjBeduMSD2CkDJLDNxAAD2QNNg6/view?usp=sharing
        """)
    end
    dir = "$test_dir/../Glycerol_in_Water"
    include("$dir/compute_mddf.jl")
    include("$dir/plots.jl")
end

@testitem "Polyacrylamide in DMF" begin
    test_dir = @__DIR__
    if !isfile("$test_dir/../Test/trajectories/traj_Polyacry.dcd")
        println("""
    
        WARNING: Will run tests without recomputing the data: only test trajectory files were found for all or some examples.
    
        Missing full trajectory files at $test_dir/trajectories
        
        Download file (280Mb): https://drive.google.com/file/d/1ug43ncCLsBATaJrT9zlbaqK6AORVvhhx/view?usp=sharing
        """)
    end
    dir = "$test_dir/../Polyacrylamide_in_DMF"
    include("$dir/compute_mddf.jl")
    include("$dir/plots.jl")
end

@testitem "POPC in Water/Ethanol" begin
    test_dir = @__DIR__
    if !isfile("$test_dir/../Test/trajectories/traj_POPC.dcd")
        println("""
    
        WARNING: Will run tests without recomputing the data: only test trajectory files were found for all or some examples.
    
        Missing full trajectory files at $test_dir/trajectories
        
        Download file (365Mb): https://drive.google.com/file/d/12TT5tblkFp1NtFOAQgjjGhmnYaXA8vQi/view?usp=sharing
        """)
    end
    dir = "$test_dir/../POPC_in_Water-Ethanol"
    include("$dir/compute_mddf.jl")
    include("$dir/plots.jl")
end

@testitem "Protein in Glycerol" begin
    test_dir = @__DIR__
    if !isfile("$test_dir/../Test/trajectories/glyc50_complete.dcd")
        println("""
    
        WARNING: Will run tests without recomputing the data: only test trajectory files were found for all or some examples.
    
        Missing full trajectory files at $test_dir/trajectories
        
        Download file (3Gb): https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing
        """)
    end
    dir = "$test_dir/../Protein_in_Glycerol"
    include("$dir/MDDF/mddf.jl")
    include("$dir/Density2D/density2D.jl")
    include("$dir/Density3D/density3D.jl")
    cd("$dir/Density3D")
    run(`vmd -dispdev text -e grid_render.vmd`)
    run(`convert grid.tga grid.png`)
end
