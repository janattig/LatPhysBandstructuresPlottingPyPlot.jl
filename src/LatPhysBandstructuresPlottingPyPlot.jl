################################################################################
#
#   module LatPhysBandstructuresPlottingPyPlot
#   -> PyPlot
#   -> LatPhysBase
#   -> LatPhysReciprocal
#   -> LatPhysBandstructures
#   -> LinearAlgebra
#
#   --> PLOTTING OF BANDSTRUCTURES AND FERMI SURFACES
#           - 2D lattices
#           - 3D lattices
#
################################################################################


# start of module
module LatPhysBandstructuresPlottingPyPlot

# used libraries
using LatPhysBase
using LatPhysReciprocal
using LatPhysBandstructures
using LinearAlgebra
using LatPhysPlottingPyPlot
using LatPhysReciprocalPlottingPyPlot
using LaTeXStrings
using PyPlot

# explicitly import PyPlot.plot to overwrite
import PyPlot.plot



# include plotting of bandstructures
include("plotting_bandstructures/bandstructure_plotting.jl")

# include plotting of energy manifoldss
include("plotting_energy_manifolds/energy_manifold_plotting.jl")


end # module
