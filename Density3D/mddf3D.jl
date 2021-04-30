"""

This function builds the grid of the 3D density function and fills an array of
mutable structures of type Atom, containing the position of the atoms of 
grid, the closest atom to that position, and distance. 

"""
function mddf_setgrid3D(solute,solute_atoms,mddf_result;
                        dmax=7.0,step=0.5,dmin=1.0)
  
  function interpolate(x₁,x₂,y₁,y₂,xₙ)
    dydx = (y₂-y₁)/(x₂-x₁)
    yₙ = y₁ + dydx*(xₙ-x₁) 
    return yₙ
  end

  # Maximum and minimum coordinates of the solute
  lims = maxmin(solute_atoms)
  n = @. ceil(Int,(lims.xlength + dmax)/step + 1)

  # Building the grid with the nearest solute atom information
  igrid = 0 
  grid = PDBTools.Atom[]
  for ix in 1:n[1]; x = lims.xmin[1] + step*(ix-1)
  for iy in 1:n[2]; y = lims.xmin[2] + step*(iy-1)
  for iz in 1:n[3]; z = lims.xmin[3] + step*(iz-1)
    rgrid = -1
    iat, r = PDBTools.closest(x,y,z,solute_atoms)
    if (dmin < r < dmax)  
      if rgrid < 0 || r < rgrid
        at = solute_atoms[iat]
        # Get contribution of this atom to the MDDF
        c = ComplexMixtures.contrib(solute,mddf_result.solute_atom,[at.index])
        # Interpolate c at the current distance
        iright = findfirst(d -> d > r, mddf_result.d)
        ileft = iright - 1
        cᵣ = interpolate(mddf_result.d[ileft],mddf_result.d[iright],
                         c[ileft],c[iright],r)
        gridpoint = Atom(index=at.index,index_pdb=at.index_pdb,
                         name=at.name,chain=at.chain,
                         resname=at.resname,resnum=at.resnum,
                         x=x,y=y,z=z,
                         occup=r, b=cᵣ,
                         model=at.model,segname=at.segname,)
        if rgrid < 0
          igrid += 1
          push!(grid,gridpoint)
        elseif r < rgrid
          grid[igrid] = gridpoint
        end
        rgrid = r
      end
    end
  end #iz
  end #iy
  end #ix

  # Now will scale the density to be between 0 and 99.9 in the temperature
  # factor column, such that visualization is good enough
  bmin, bmax = +Inf, -Inf
  for gridpoint in grid
    bmin = min(bmin,gridpoint.b)
    bmax = max(bmax,gridpoint.b)
  end
  for gridpoint in grid
    gridpoint.b = (gridpoint.b - bmin)/(bmax-bmin)
  end

  return grid
end

"""
#
# Function that writes the gmd3D to a pdb file
#

# If called already with the grid of the gmd3D computed

  file = open(output,"w")
  println(file,"REMARK  3D representation of MDDF density ")
  println(file,"REMARK  ")
  println(file,"REMARK  Occupancy column contains distance to solute. ")
  println(file,"REMARK  B-factor column contains the scaled density, such that 99.9 is the maximum ")
  println(file,"REMARK  density obsered and 0. is the minimum density observed. The actual limits are: ")
  println(file,"REMARK  Minimum density: \$minrho")
  println(file,"REMARK  Maximum density: \$maxrho")
  for i in 1:n
    iat = grid[i].atom
    name = mysim.atoms[iat].name
    resname = mysim.atoms[iat].resname
    chain = mysim.atoms[iat].chain
    resnum = mysim.atoms[iat].resid
    x = grid[i].x[1]
    y = grid[i].x[2]
    z = grid[i].x[3]
    b = 99.9*(grid[i].rho - minrho)/(maxrho-minrho)
    occup = grid[i].dmin
    model = 0
    atom = PDBTools.Atom(i,name,resname,chain,resnum,x,y,z,b,occup,model)
    println(file,PDBTools.write_atom(atom))
  end
  close(file)

end

# If called directly with the solute selection and gmd_solute file name

function gmd3D_write(mysim :: Simulation, 
                     solute :: Vector{Int64}, 
                     gmd_solute_file :: String, 
                     output :: String;
                     dmax = 7.0, dmin = 1.0, step = 0.5, scale=nothing)

  grid = gmd3D(mysim, solute, gmd_solute_file, 
               dmax = 7.0, dmin = 1.0, step = 0.5)

  gmd3D_write(mysim, grid, output, scale=scale)

end

"""
