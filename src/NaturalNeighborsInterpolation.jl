module NaturalNeighborsInterpolation

export NaturalNeighborsInterpolator

using VoronoiDelaunay
using VoronoiCells
using Plots
using Gadfly
import Cairo, Fontconfig
using GeometricalPredicates
using Combinatorics
using Polyhedra, CDDLib
include("utils.jl")

"""    interpolate(pointlist, values, interpolation_point)

Interpolate at `interpolation_point` using the known `values` at coordinate `pointlist`.
`points` is a `Vector` of `Point2D`, and `interpolation_point` a single `Point2D`.
"""
function interpolate(points, values, interpolation_point, tess, rect, convex_hull, fallback; attempt=0)
    # If point is on grid of known values, return that kown value
    if interpolation_point in points
        idx = findall(x -> x == interpolation_point, points)[1]
        return values[idx]
    end

    if fallback in ["nearest", "nan"] && !inpolygon(convex_hull, interpolation_point)
        if fallback == "nan"
            return NaN64
        elseif fallback=="nearest"
            return values[findNearest(interpolation_point, points)]
        end
    end
    #= Plot
    scatter([getx(p) for p in points], [gety(p) for p in points])
    Plots.plot!(tess)
    savefig("OriginalVoronoi.png")
    # end Plot # =#

    allpoints = copy(points)
    push!(allpoints, interpolation_point)
    tess2 = voronoicells(allpoints, rect)

    #= Plot 
    scatter([getx(p) for p in allpoints], [gety(p) for p in allpoints])
    Plots.plot!(tess2)
    savefig("ExtendedVoronoi.png")
    # end Plot # =#

    # Find the newly made cell of the interpolation point
    interpolant_cell = nothing
    for cell in tess2.Cells
        if inpolygon(Polygon(cell...), interpolation_point)
            interpolant_cell = cell
            break
        end
    end
    if isnothing(interpolant_cell)
        if attempt>10
            return values[findNearest(interpolation_point, points)]
        elseif attempt <= 10
            return interpolate(points, values, interpolation_point, tess, rect, convex_hull, fallback, attempt=attempt+1)
        end
    end
    interpolant_cell = reduce(interpolant_cell)
    
    #= Plot
    scatter([getx(p) for p in points], [gety(p) for p in points])
    Plots.plot!(tess)
    Plots.plot!([getx(p) for p in interpolant_cell], [gety(p) for p in interpolant_cell])
    savefig("Overlayedoronoi.svg")
    # end Plot # =#

    interpolant_polygon = Polygon(interpolant_cell...)
    interpolant_edges = getEdges(interpolant_cell)
    neighboring_cells = []
    areas = []
    relevant_values = []
    relevant_points = []
    for cell in tess.Cells
        # List of points of that cell that are inside the interpolants cell
        inner_points = []
        for cellpoint in cell
            if inpolygon(interpolant_polygon, cellpoint)||
               any([isBetween(edge.a, edge.b, cellpoint) for edge in interpolant_edges])
                push!(inner_points, cellpoint)
            end
        end
        # No inner points means the cell is no natural neighbor; skip it
        if length(inner_points) ==0
            continue
        end
        # Now the other way around: Find the points of interpolant_cell living on
        # the edge of cell
        polygon = Polygon(cell...)
        for itp_edge in interpolant_edges
            if inpolygon(polygon, midPoint(itp_edge))
                push!(inner_points, itp_edge.a)
                push!(inner_points, itp_edge.b)
            end
        end

        A = getAreaSimple(sort_angular(inner_points))
        push!(areas, A)
        push!(neighboring_cells, inner_points)
        for p_i in 1:length(points)
            if inpolygon(polygon, points[p_i])
                push!(relevant_values, values[p_i])
                push!(relevant_points, points[p_i])
            end
        end
    end

    λ = areas ./ sum(areas)

    #= Plot
    colors = [:black, :red, :yellow, :green, :blue, :orange, :darkred]
    for (i,cell) in enumerate(neighboring_cells)
        Plots.plot(tess)
        scatter!([getx(p) for p in cell], [gety(p) for p in cell], color=colors[i])
        savefig("plots/Cell_"*string(interpolation_point)*"_$i.svg")
    end
    Plots.plot(tess)
    for (i,cell) in enumerate(neighboring_cells)
        scatter!([getx(p) for p in cell], [gety(p) for p in cell], color=:black)
    end
    scatter!([getx(p) for p in relevant_points], [gety(p) for p in relevant_points])
    Plots.plot!([getx(p) for p in interpolant_cell], [gety(p) for p in interpolant_cell])
    annotate!([getx(p) for p in relevant_points], [gety(p) for p in relevant_points],
     [text(string(λ[k]), 3,  :top) for k in eachindex(λ)])
    savefig("plots/InnerPoints_"*string(interpolation_point)*".svg")
    # end plot =#

    #== Test if λ is correct via the Local Coordinates Property ==#
    #supposed_x = sum([λ[k] * getx(relevant_points[k]) for k in 1:length(λ)])
    #supposed_y = sum([λ[k] * gety(relevant_points[k]) for k in 1:length(λ)])
    #println(supposed_x,",", getx(interpolation_point))
    #println(supposed_y,",", gety(interpolation_point))

    interpolated_value = sum([λ[k] * relevant_values[k] for k in 1:length(λ)])
    return interpolated_value
end

"""    NaturalNeighborsInterpolator(pointlist, values; fallback= \"nearest\")

Return the interpolatant object.

`pointlist` is an array of `Point2D` or an array of size `(2,n)`, representing the `n` coordinates
at which `values` of the function are known. Points outside the convex
hull of the `pointlist` are treated with the method specified by `fallback`:
Return the value of the nearest neighbor (`\"nearest\"`), use the natural neighbor
interpolation (`\"natural\"`) or return NaN64 (\"nan\").
"""
function NaturalNeighborsInterpolator(pointlist, values; fallback="nearest")
    if eltype(pointlist)<:Real
        pointlist = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:size(pointlist)[2]]
    end
    @assert fallback in ["nearest", "natural", "nan"] "fallback has to be one of [\"nearest\", \"natural\", \"nan\"] but was \"$fallback\""
    rect = Rectangle(Point(1, 1), Point(2, 2))
    tess = voronoicells(pointlist, rect)

    A = vrep([[getx(p), gety(p)] for p in pointlist])
    P = polyhedron(A, CDDLib.Library())
    Pch = convexhull(P, P)
    removevredundancy!(Pch)
    hullpoints = [Point2D(e[1], e[2]) for e in points(Pch)]
    hull = Polygon(sort_angular(hullpoints)...)

    interpolator = let pointlist=pointlist, values=values, tess=tess, rect=rect, hull=hull, fallback=fallback
        function interpolator(interpolation_point)
            return interpolate(pointlist, values, interpolation_point, tess, rect, hull, fallback)
        end
    end
    return interpolator
end

end # module