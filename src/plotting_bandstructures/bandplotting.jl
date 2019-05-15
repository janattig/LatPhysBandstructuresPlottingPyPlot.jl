################################################################################
#
#   BAND STRUCTURE PLOTTING
#
################################################################################
using PyPlot

"""
    plotBandstructure(
            bandstructure::Bandstructure
          [;  limits_energy="AUTO",
            plot_title::String="",
            plot_color_valid="b",
            plot_color_invalid="r",
            figsize::Tuple=(6,4),
            showPlot::Bool=true,
            save_filename::String="NONE" ]
        )
    plotLTBandstructure(
            unitcell::Unitcell,
            path::Path,
         [  bondInteractionMatrix::Function
          ; resolution::Int64=-1,
            enforce_hermitian::Bool=false,
            ... ]
        )
Plots the band struture of a passed `Bandstructure` object along some its path
and returns the plot as a `PyPlot.Figure` object.

Additional options include plotting related options of `PyPlot` as well as determining if the plot is saved or shown.
# Examples
```julia-repl
julia> plotLTBandstructure(unitcell, path)
PyPlot.Figure(...)
julia> plotLTBandstructure(unitcell, path, c->diagm([ c[3] ]))
PyPlot.Figure(...)
julia> plotLTBandstructure(unitcell, path, showPlot=false)
PyPlot.Figure(...)
julia> plotLTBandstructure(unitcell, save_filename="myplot.pdf")
PyPlot.Figure(...)
julia> plotLTBandstructure(bandstructure)
PyPlot.Figure(...)
```
"""
function plotBandstructure(
            bandstructure:: Bandstructure{RPA},
            segment_resolution::Array{Int64, 1},
            kappa::Float64;
            limits_energy="AUTO",
            plot_title::String="",

            plot_color_valid="b",
            plot_color_invalid="r",
            figsize::Tuple=(6,4),
            showPlot::Bool=false,
            save_filename::String="NONE"
        ) where {RPA<: AbstractReciprocalPath{P} where {P <: AbstractReciprocalPoint{D} where {D}}}

    ###########################
    #   INITIAL SETTINGS
    ###########################



    # get the path from the bandstructure
    path = bandstructure.path

    # configure plot environment
    rc("font", family="serif")

    # create a new figure
    fig = figure(figsize=figsize)




    ###########################
    #   PLOT BANDS
    ###########################

    # plot the band structure

    for s in 1:length(bandstructure.bands)
        # plot the segment (only invalid stuff)
        for b in 1:length(bandstructure.bands[s])
            # xvalues
            xvals = collect(1:segment_resolution[s]) .+ sum(segment_resolution[1:s-1])
            yvals = bandstructure.bands[s][b]

            # plot everything
            plot(
                xvals,yvals,color = [color[cnum].r,color[cnum].g,color[cnum].b])

        end


    end



    ###########################
    #   SET ALL TICKS (POINTS)
    ###########################

    # get the current axis
    ax = gca()
    axx = ax[:get_xaxis]()
    # compile tick positions and labels
    point_pos = Int64[]
    push!(point_pos, 1)
    for l in 1:length(segment_resolution)
        push!(point_pos, sum(segment_resolution[1:l]))
    end

    point_labels = String[point(path,i).label_LaTeX for i in 1:numPoints(path)]
    # configure tick labels
    xticks(point_pos, point_labels)
    # configure ticks
    axx[:set_tick_params](which="both", direction="out")
    axx[:set_tick_params](which="top", color="none")
    axy = ax[:get_yaxis]()
    axy[:set_tick_params](which="both", direction="out")

    # plot vertical lines for each point
    for p in point_pos
        axvline(p,color=[0.6, 0.6, 0.6], linestyle="--")
    end


    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # label the axis
    xlabel("momentum")
    ylabel("energy")


    # energy limits
    # check if specific boundaries are desired
    if !(limits_energy == "AUTO")
        ylim(limits_energy[1], limits_energy[2])
    end

    # momentum limits (x axis)
    xlim(0, maximum(point_pos)+1)

    # set the title
    if plot_title == "AUTO"
        # set the title to an automatically generated title
        title("energy spectrum")
    elseif plot_title == ""
        # do nothing title related
    else
        # set the title to the given title
        title(plot_title)
    end





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


function plotBandstructure(
            unitcell :: U,
            path :: RPA,
            bondInteractionMatrix:: HB,
            segment_resolution::Array{Int64, 1};
            enforce_hermitian::Bool=false,
            limits_energy="AUTO",
            plot_title::String="",
            plot_color_valid="b",
            plot_color_invalid="r",
            figsize::Tuple=(6,4),
            showPlot::Bool=true,
            save_filename::String="NONE"
        ) where {L,NB,S,B<:AbstractBond{L,NB},U<:AbstractUnitcell{S,B},HB<:AbstractBondHamiltonian{L,NS} where {NS}, RPA<: AbstractReciprocalPath{P} where {P <: AbstractReciprocalPoint{D} where {D}}}

    # calculate the bandstructure
    bandstructure = getBands(unitcell, path, bondInteractionMatrix, segment_resolution, enforce_hermitian=enforce_hermitian)
    # call the respective function
    return plotBandstructure(
                bandstructure;
                limits_energy=limits_energy,
                plot_title=plot_title,
                plot_color_valid=plot_color_valid,
                plot_color_invalid=plot_color_invalid,
                figsize=figsize,
                showPlot=showPlot,
                save_filename=save_filename
            )
end
export plotBandstructure
