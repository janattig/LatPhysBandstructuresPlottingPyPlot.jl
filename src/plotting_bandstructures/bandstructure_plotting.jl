################################################################################
#
#   BAND STRUCTURE PLOTTING
#
################################################################################

# function to plot a bandstructure
function plotBandstructure(
            bs :: BS
            ;
            new_figure :: Bool = true,
            figsize :: Tuple = (6,4),
            color :: Vector{<:Integer} = [100,120,255],
            kwargs...
        ) where {RP, P<:AbstractReciprocalPath{RP}, L,UC,HB,H<:AbstractHamiltonian{L,UC,HB}, BS <: AbstractBandstructure{P,H}}

    ###########################
    #   INITIAL SETTINGS
    ###########################

    # configure plot environment
    #rc("font", family="serif")

    # create a new figure
    if new_figure
        fig = figure(figsize=figsize)
    else
        fig = gcf()
    end




    ###########################
    #   PLOT BANDS
    ###########################

    # collect the segment breaks
    k_point_indices = zeros(Int64, length(energies(bs))+1)
    # get the labels of these breaks
    k_point_labels = label.(path(bs))

    # plot the band structure
    for s in 1:length(energies(bs))
        # push the next k point as index into the array
        k_point_indices[s+1] = length(energies(bs)[s][1]) + k_point_indices[s] - 1
        # plotting segment s, collecting x values for bands
        xvals = range(k_point_indices[s], stop=k_point_indices[s+1], length=length(energies(bs)[s][1]))
        # plot all bands
        bands = zeros(Float64, length(energies(bs)[s][1]), length(energies(bs)[s]))
        for b in 1:length(energies(bs)[s])
            bands[:,b] .= energies(bs)[s][b]
        end
        plot(xvals,bands, color=color./255)
    end



    ###########################
    #   SET ALL TICKS (POINTS)
    ###########################

    # get the current axis
    ax = gca()
    axx = ax.get_xaxis()

    # configure tick labels on x axis
    xticks(k_point_indices, k_point_labels)
    # momentum limits (x axis)
    xlim(0, k_point_indices[end]+1)

    # configure ticks on x axis
    axx.set_tick_params(which="both", direction="out")
    axx.set_tick_params(which="top", color="none")

    # configure ticks on x axis
    axy = ax.get_yaxis()
    axy.set_tick_params(which="both", direction="out")

    # plot vertical lines for each point
    for p in k_point_indices
        axvline(p,color=[0.6, 0.6, 0.6], linestyle="--")
    end


    ###########################
    #   CONFIGURE AXIS
    ###########################

    # label the axis
    xlabel("momentum")
    ylabel("energy")



    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # return the figure object
    return fig
end

# pass unknown arguments directly to construction of band structure
function plotBandstructure(
        args...
        ;
        kwargs...
    )

    # create and plot a bandstructure
    plotBandstructure(getBandstructure(args...); kwargs...)
end

# export the function
export plotBandstructure
