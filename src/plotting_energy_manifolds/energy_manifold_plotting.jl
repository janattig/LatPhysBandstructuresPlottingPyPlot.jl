################################################################################
#
#   ENERGY MANIFOLD PLOTTING
#
################################################################################

# function to plot an energy manifold (2D)
function plotEnergyManifold(
            em :: AEM
            ;
            new_figure :: Bool = true,
            figsize :: Tuple = (6,6),
            color :: Vector{<:Integer} = [0,0,255],
            plot_bz :: Bool = true,
            kwargs...
        ) where {LS,S<:AbstractSite{LS,2},L,B,UC<:AbstractUnitcell{S,B},HB,H<:AbstractHamiltonian{L,UC,HB}, AEM <: AbstractEnergyManifold{H}}


    ###########################
    #   INITIAL SETTINGS
    ###########################

    # create a new figure
    if plot_bz && new_figure
        # plot the brillouin zone
        plotBrillouinZone(getBrillouinZone(unitcell(hamiltonian(em))), figsize=figsize)
        # get the figure
        fig = gcf()
    elseif new_figure
        # only new figure
        fig = figure(figsize=figsize)
    else
        fig = gcf()
    end



    ###########################
    #   PLOT POINTS
    ###########################

    # obtain the points
    k_points = kpoints(em)
    kx = [k[1] for k in k_points]
    ky = [k[2] for k in k_points]

    # scatter all points
    scatter(kx,ky, color=color./255)



    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # equal axis
    gca().set_aspect("equal")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()


    # return the figure object
    return fig
end

# function to plot an energy manifold (3D)
function plotEnergyManifold(
            em :: AEM
            ;
            new_figure :: Bool = true,
            figsize :: Tuple = (6,6),
            color :: Vector{<:Integer} = [0,0,255],
            plot_bz :: Bool = true,
            kwargs...
        ) where {LS,S<:AbstractSite{LS,3},L,B,UC<:AbstractUnitcell{S,B},HB,H<:AbstractHamiltonian{L,UC,HB}, AEM <: AbstractEnergyManifold{H}}


    ###########################
    #   INITIAL SETTINGS
    ###########################

    # create a new figure
    if plot_bz && new_figure
        # plot the brillouin zone
        plotBrillouinZone(getBrillouinZone(unitcell(hamiltonian(em))), figsize=figsize)
        # get the figure
        fig = gcf()
    elseif new_figure
        # only new figure
        fig = figure(figsize=figsize)
    else
        fig = gcf()
    end



    ###########################
    #   PLOT POINTS
    ###########################

    # obtain the points
    k_points = kpoints(em)
    kx = [k[1] for k in k_points]
    ky = [k[2] for k in k_points]
    kz = [k[3] for k in k_points]

    # scatter all points
    scatter3D(kx,ky,kz, color=color./255)



    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # equal axis
    gca().set_aspect("equal")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()


    # return the figure object
    return fig
end

# export plotting function
export plotEnergyManifold
