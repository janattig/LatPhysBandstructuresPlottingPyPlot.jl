



#####PLOTTING THE FERMI SURFACE


function plotFermiSurface2D(
            k_values::Array{Float64, 2},
            brillouin_zone:: BZ,
            kappa::Float64;
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )   where {L,NB, BZ<:AbstractBrillouinZone{R} where {
                 R <: AbstractReciprocalUnitcell{P,B} where {
                 P <: AbstractReciprocalPoint{D} where {D},
                 B<:AbstractBond{L,NB},
}}}


    ###########################
    #   INITIAL SETTINGS
    ###########################

    # configure plot environment
    rc("font", family="serif")

    # create a new figure
    fig = figure(figsize=figsize)



    ###########################
    #   PLOT FERMI SURFACE
    ###########################

    # if brillouin zone not empty, plot it as well
     if length(corners(brillouin_zone)) > 0
        plotBrillouinZone(brillouin_zone, new_figure=false)
    end

    # scatter the points
    
    scatter(k_values[:,1], k_values[:,2], color=plot_color)



    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # set the title
    if plot_title == "AUTO"
        # set the title to an automatically generated title
        title("Fermi surface")
    elseif plot_title == ""
        title("κ=$(kappa)")
    else
        # set the title to the given title
        title(plot_title)
    end

    # equal axis
    gca()[:set_aspect]("equal")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # save the plot
    if save_filename != "NONE"
        # make sure the directory exists
        if contains(save_filename, "/")
    		# get the containing folder
    		folder = save_filename[1:findlast(save_filename, '/')]
    		# build the path to that folder
    		mkpath(folder)
    	end
        # save the plot
        savefig(save_filename)
    end

    # maybe show the plot
    if showPlot
        show()
    end

    # return the figure object
    return fig

end

function plotFermiSurface2D(
            unitcell:: U,
            bondInteractionMatrix:: HB,
            N_points::Int64,
            kappa::Float64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(4),
            bounds_upper::Array{Float64,1}=2*pi.*ones(4),
            refold_to_first_BZ::Bool=true,
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        )   where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B},HB<:AbstractBondHamiltonian{L,NS} where {NS}}

    # Calculate the Brillouin Zone
     
    
        brillouin_zone = getBrillouinZone(getReciprocalUnitcell(unitcell))
    

    # calculate the Fermi surface first
    fermi_surface = getFermiSurface2D(
            unitcell,
            bondInteractionMatrix,
            N_points,
            kappa,
            fermi_energy=fermi_energy,
            enforce_hermitian=enforce_hermitian,
            epsilon=epsilon,
            epsilon_k=epsilon_k,
            bounds_lower=bounds_lower,
            bounds_upper=bounds_upper,
            refold_to_first_BZ=refold_to_first_BZ
        )

    # plot the fermi surface
    plotFermiSurface2D(
            fermi_surface,
            brillouin_zone,
            kappa,
            plot_title=plot_title,
            plot_color=plot_color,
            figsize=figsize,
            showPlot=showPlot,
            save_filename=save_filename
        )
    
end

# plot the Fermi surface (3D)
function plotFermiSurface3D(
            k_values::Array{Float64, 2},
            brillouin_zone:: BZ,
            kappa::Float64;
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        ) where {L,NB, BZ<:AbstractBrillouinZone{R} where {
                 R <: AbstractReciprocalUnitcell{P,B} where {
                 P <: AbstractReciprocalPoint{D} where {D},
                 B<:AbstractBond{L,NB},
}}}


    ###########################
    #   INITIAL SETTINGS
    ###########################

    # configure plot environment
    rc("font", family="serif")

    # create a new figure
    fig = figure(figsize=figsize)



    ###########################
    #   PLOT FERMI SURFACE
    ###########################

    # if brillouin zone not empty, plot it as well
   if length(corners(brillouin_zone)) > 0
      plotBrillouinZone(brillouin_zone, new_figure=false)
    end

    # scatter the points
    scatter3D(k_values[:,1], k_values[:,2], k_values[:,3], color=plot_color)
    #view3D(k_values[:,1], k_values[:,2], k_values[:,3], color=plot_color)



    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # set the title
    if plot_title == "AUTO"
        # set the title to an automatically generated title
        title("Fermi surface")
    elseif plot_title == ""
         title("κ=$(kappa)")
    else
        # set the title to the given title
        title(plot_title)
    end

    # equal axis
    gca()[:set_aspect]("equal")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # save the plot
    if save_filename != "NONE"
        # make sure the directory exists
        if contains(save_filename, "/")
    		# get the containing folder
    		folder = save_filename[1:findlast(save_filename, '/')]
    		# build the path to that folder
    		mkpath(folder)
    	end
        # save the plot
        savefig(save_filename)
    end

    # maybe show the plot
    if showPlot
        show()
    end

    # return the figure object
    return fig

end
function plotFermiSurface3D(
            unitcell:: U,
            bondInteractionMatrix:: HB,
            N_points::Int64,
            kappa::Float64;
            fermi_energy::Float64=0.0,
            enforce_hermitian::Bool=false,
            epsilon::Float64=1e-10,
            epsilon_k::Float64=1e-10,
            bounds_lower::Array{Float64,1}=-2*pi.*ones(3),
            bounds_upper::Array{Float64,1}=2*pi.*ones(3),
            refold_to_first_BZ::Bool=true,
            plot_title::String="",
            plot_color="b",
            figsize::Tuple=(6,6),
            showPlot::Bool=true,
            save_filename::String="NONE"
        ) where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B},HB<:AbstractBondHamiltonian{L,NS} where {NS}}
               
    # Calculate the BZ
   
       brillouin_zone = getBrillouinZone(getReciprocalUnitcell(unitcell))
   

    # calculate the Fermi surface first
    fermi_surface = getFermiSurface3D(
            unitcell,
            bondInteractionMatrix,
            N_points,
            kappa,
            fermi_energy=fermi_energy,
            enforce_hermitian=enforce_hermitian,
            epsilon=epsilon,
            epsilon_k=epsilon_k,
            bounds_lower=bounds_lower,
            bounds_upper=bounds_upper,
            refold_to_first_BZ=refold_to_first_BZ
        )

    # plot the fermi surface
    plotFermiSurface3D(
            fermi_surface,
            brillouin_zone,
            kappa,
            plot_title=plot_title,
            plot_color=plot_color,
            figsize=figsize,
            showPlot=showPlot,
            save_filename=save_filename
        )

end




