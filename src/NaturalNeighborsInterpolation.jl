module NaturalNeighborsInterpolation
using VoronoiDelaunay
using VoronoiCells
using Plots
using Gadfly
import Cairo, Fontconfig
using Random
using GeometricalPredicates
using Combinatorics

"""    Edge(a,b)

An undirected egde between `Point2D` `a` and `Point2D` `b`. 
"""
struct Edge
    a::Point2D
    b::Point2D
end

function equals(e1::Edge, e2::Edge)
    if ((e1.a == e2.a) && (e1.b == e2.b)) || ((e1.a == e2.b) && (e1.b == e2.a))
        return true
    end
    return false
end

"""    sharePoint(e1::Edge, e2::Edge)

Return the shared point if `e1` and `e2` share a common point, else `nothing`.
"""
function sharePoint(e1::Edge, e2::Edge)
    a,b = e1.a, e1.b
    c,d = e2.a, e2.b
    if a in [c,d]
        return a
    elseif b in [c,d] 
        return b
    end
    return nothing
end


function getplotxy(t::VoronoiDelaunay.DelaunayTriangle{Point2D})
    x = [getx(getpoint(t)) for getpoint in [geta, getb, getc, geta]]
    y = [gety(getpoint(t)) for getpoint in [geta, getb, getc, geta]]
    return x, y
end

function getBowyerWatsonEnvelope(tesselation, interpolation_point)
    loc = locate(tesselation, Point(getx(interpolation_point), gety(interpolation_point)))
    triangles = Array{VoronoiDelaunay.DelaunayTriangle{Point2D},1}(undef, 4)
    triangles[1] = loc
    x, y = getplotxy(loc)
    removed_edges = [Edge(Point(x[i], y[i]), Point(x[i+1], y[i+1])) for i in 1:length(x)-1]
    envelope_edges = []
    circumcenters = Array{Point2D,1}(undef, 4)
    circumcenters[1] = circumcenter(loc)
    for (i, move) in enumerate([movea, moveb, movec])
        triangle = move(tesselation, loc)
        triangles[i+1] = triangle
        circumcenters[i+1] = circumcenter(triangle)
        ps = Array{Point2D,1}(undef, 4)
        for (i, getpoint) in enumerate([geta, getb, getc, geta])
            p = getpoint(triangle)
            ps[i] = p
        end
        edges = [Edge(ps[i], ps[i+1]) for i in 1:length(ps)-1]
        push!(envelope_edges, edges...)
    end
    # Final edges are the envelope edges without the removed edges
    final_edges = []
    for e_edge in envelope_edges
        removed = false
        for r_edge in removed_edges
            if equals(e_edge, r_edge)
                removed = true
            end
        end
        if !removed
            push!(final_edges, e_edge)
        end
    end
    #= 
    Find the removed points. Sometimes a two edges are removed at a point that was only 
    connected to these two edges. If such a thing happens, the construction of the envelope
    was not successful. This following part selects points that are possibly removed. The final_edges
    check happens in `interpolate`.
    =#
    removed_points = []
    for (r1, r2) in combinations(removed_edges, 2)
        shared_point = sharePoint(r1,r2)
        if !isnothing(shared_point)
            push!(removed_points, shared_point)
        end
    end

    return final_edges, triangles, circumcenters, removed_points
end

function sort_angular(points)
    xs = [getx(p) for p in points]
    ys = [gety(p) for p in points]
    midpoint = centerOfWeight(points)
    mid_x = getx(midpoint)
    mid_y = gety(midpoint)
    angles = atan.(ys .- mid_y, xs .- mid_x)
    sortind = sortperm(angles)
    points = points[sortind]
    return points
end

function isConvex(points)
    n = length(points)
    crossproducts = Array{Float64,1}(undef, n)
    for i in 1:n
        x1, y1 = getx(points[i]), gety(points[i])
        if i < n - 1
            x2, y2 = getx(points[i+1]), gety(points[i+1])
            x3, y3 = getx(points[i+2]), gety(points[i+2])
        elseif i == n - 1
            x2, y2 = getx(points[i+1]), gety(points[i+1])
            x3, y3 = getx(points[1]), gety(points[1])
        elseif i == n
            x2, y2 = getx(points[1]), gety(points[1])
            x3, y3 = getx(points[2]), gety(points[2])
        end
        dx1 = x2 - x1
        dy1 = y2 - y1
        dx2 = x3 - x2
        dy2 = y3 - y2
        crossproducts[i] = dx1 * dy2 - dy1 * dx2
    end
    println(crossproducts)
    return all(crossproducts .<= 0) || all(crossproducts .>= 0)
end

"""    getEdges(points)

Return the edges that connect the `points`.
"""
function getEdges(points)
    n = length(points)
    lines = [Edge(points[i], points[i+1]) for i in 1:n-1]
    push!(lines, Edge(points[n], points[1]))
    return lines
end

"""    twoLinesIntersect(a,b,c,d)

Find the intersection of two lines, one going from pont `a` to point `b`, the other from `c` to `d`.


Returns two parameters `t2, t1` which describe how far along the second / first line the intersection is. 
Line segments intersect if `∀ t ∈ [t1, t2]: 0 < t < 1 `.
"""
function twoLinesIntersect(a, b, c, d)
    a_x, a_y = getx(a), gety(a)
    b_x, b_y = getx(b), gety(b)
    c_x, c_y = getx(c), gety(c)
    d_x, d_y = getx(d), gety(d)
    bxax = (b_x - a_x)
    byay = (b_y - a_y)
    cyay = (c_y - a_y)
    aycy = (a_y - c_y)
    cxax = (c_x - a_x)
    dycy = (d_y - c_y)
    dxcx = (d_x - c_x)
    num = (bxax * aycy + byay * cxax)
    denom = (bxax * dycy - byay * dxcx)
    t2 = num / denom
    if bxax==0
        t1 = (cyay + dycy * t2)/byay
    else
        t1 = (cxax + dxcx*t2)/bxax
    end
    return t2, t1
end

"""    isSelfIntersecting(points)

Return `true` if the polygon defined by `points` is self-intersecting.
"""
function isSelfIntersecting(points)
    lines = getEdges(points)
    for (line1, line2) in combinations(lines, 2)
        a, b = line1.a, line1.b
        c, d = line2.a, line2.b
        t2, t1 = twoLinesIntersect(a, b, c, d)
        if (0 + eps(Float64) < t2 < 1 - eps(Float64)) && (0 + eps(Float64) < t1 < 1 - eps(Float64))
            return true
        end
    end
    return false
end

"""    decomposeIntersection(points)

Decompose a self intersecting polygon into its non-intersecting parts.
Assumes `points` is already confirmed self-intersecting, otherwise this will 
do ineffective calculations.
"""
function decomposeIntersection(points)
    lines = getEdges(points)
    intersection_points = []
    for (l1, l2) in combinations(lines, 2)
        a, b = l1.a, l1.b
        c, d = l2.a, l2.b
        t, t1 = twoLinesIntersect(a, b, c, d)
        if !((0 + eps(Float64) < t < 1 - eps(Float64)) && (0 + eps(Float64) < t1 < 1 - eps(Float64)))
            continue
        end
        t,_ = twoLinesIntersect(a, b, c, d)
        a = [getx(a), gety(a)]
        c = [getx(c), gety(c)]
        d = [getx(d), gety(d)]
        push!(intersection_points, (c, c .+ (d .- c) .* t))
        push!(intersection_points, (a, c .+ (d .- c) .* t))
    end
    newpoints = []
    cross_indices = []
    for line in lines
        push!(newpoints, line.a)
        for intersection in intersection_points
            a = Point(intersection[1]...)
            if line.a == a
                push!(newpoints, Point(intersection[2]...))
                push!(cross_indices, length(newpoints))
            end
        end
    end
    n_polygons = length(cross_indices)
    newpoints = circshift(newpoints, 1 - cross_indices[1])
    cross_indices .+= 1 - cross_indices[1]
    polys = []
    for n in 1:n_polygons-1
        poly = newpoints[cross_indices[n]:cross_indices[n+1]-1]
        push!(polys, poly)
    end
    poly = newpoints[cross_indices[end]:end]
    push!(polys, poly)
    return polys
end


function polygonOrientation(points)
    n = length(points)
    xs = [getx(p) for p in points]
    ys = [gety(p) for p in points]
    min_idx = -1
    min_idx_x = findall(p->p==minimum(xs), xs)
    if length(min_idx_x)>1
        min_y = maximum(ys)
        for i in min_idx_x
            if ys[i]<min_y
                min_idx = i
                min_y = ys[i]
            end
        end
    else
        min_idx = min_idx_x[1]
    end
    # Point at min_idx is part of the convex hull
    b = points[min_idx]
    if min_idx==1
        a = points[n]
        c = points[2]
    elseif min_idx==n
        a = points[n-1]
        c = points[1]
    else
        a = points[min_idx-1]
        c = points[min_idx+1]
    end
    # Determinant of these locally convex three points
    det = getx(b)*gety(c) + getx(a)*gety(b) + gety(a)*getx(c) - (gety(a)*getx(b) + gety(b)*getx(c) + getx(a)*gety(c))
    return sign(det)
end

function centerOfWeight(points)
    n = length(points)
    x = [getx(p) for p in points]
    y = [gety(p) for p in points]
    return Point(sum(x)/n, sum(y)/n)
end

function getAreaSimple(points)
    @assert !isSelfIntersecting(points)
    A = zero(Float64)
    n = length(points)
    xs = [getx(p) for p in points]
    ys = [gety(p) for p in points]
    mid_x = sum(xs) / n
    mid_y = sum(ys) / n
    midpoint = Point(mid_x, mid_y)
    for q in 1:length(points)-1
        temp_triangle = Primitive(midpoint, points[q], points[q+1])
        A += area(temp_triangle)
    end
    A += area(Primitive(midpoint, points[1], points[end]))
end

"""    getArea(points)


Calculate the area of a polygon
"""
function getArea(points, envelope)
    A = zero(Float64)
    if !isSelfIntersecting(points)
        points = sort_angular(points)
        n = length(points)
        xs = [getx(p) for p in points]
        ys = [gety(p) for p in points]
        mid_x = sum(xs) / n
        mid_y = sum(ys) / n
        midpoint = Point(mid_x, mid_y)
        for q in 1:length(points)-1
            temp_triangle = Primitive(midpoint, points[q], points[q+1])
            A += area(temp_triangle)
        end
        A += area(Primitive(midpoint, points[1], points[end]))
    else
        #name = randstring(10)
        #plotPoly(points, name*".svg")
        #println(name)
        polygons = decomposeIntersection(points)
        for (i,poly) in enumerate(polygons)
            if inpolygon(envelope, centerOfWeight(poly))
                A += getArea(poly, envelope)
            end
        end
    end
    return abs(A)
end

"""    findNearest(p, points)

Find nearest point in `points` to point `p`.

Returns index of nearest point.
"""
function findNearest(p, points)
    function distSq(p, q)
        return (getx(p) - getx(q))^2 + (gety(p) - gety(q))^2
    end
    min_distance = 0.0
    min_index = -1
    for (i,q) in enumerate(points)
        d = distSq(p, q)
        if d<min_distance || min_index==-1
            min_index = i
            min_distance = d
        end
    end
    return min_index
end

"""    areCollinear(a,b,c)

If points `a`, `b` and `c` are collinear, return `true`
"""
function areCollinear(a,b,c)
    dx1 = getx(a) - getx(b)
    dx2 = getx(c) - getx(b)
    dy1 = gety(a) - gety(b)
    dy2 = gety(c) - gety(b)
    crossproduct = dx1 * dy2 - dy1 * dx2
    return abs(crossproduct) < eps(Float64)
end

"""    isBetween(a,b,c)

Check if `c` is on a line between `a` and `b`.
"""
function isBetween(a,b,c)
    if !areCollinear(a,b,c)
        return false
    end
    line = [getx(a) - getx(b), gety(a) - gety(b)]
    extremum1 = line[1] * getx(a) + line[2] * gety(a)
    extremum2 = line[1] * getx(b) + line[2] * gety(b)
    value = line[1] * getx(c) + line[2] * gety(c)
    if extremum1 > extremum2
        return extremum1 >= value >= extremum2
    else
        return extremum1 <= value <= extremum2
    end
end

function plotPoly(poly, name)
    x = [getx(p) for p in poly]
    y = [gety(p) for p in poly]
    tlayer = layer(x=x,y=y, Geom.path)
    myplot = plot(tlayer, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.0, maxvalue=2.0))
    draw(SVG(name, 8inch, 8inch), myplot)
end

"""    interpolate(pointlist, values, interpolation_point)

Interpolate at `interpolation_point` using the known `values` at coordinate `pointlist`.
`points` is a `Vector` of `Point2D`, and `interpolation_point` a single `Point2D`.
`tess` is the `DelaunayTessellation` of `points`.
"""
function interpolate(points, values, interpolation_point, tess)
    # If point is on grid of known values, return that kown value
    if interpolation_point in points
        idx = findall(x -> x == interpolation_point, points)[1]
        return values[idx]
    end

    final_edges, enclosed_triangles, circumcenters, removed_points = getBowyerWatsonEnvelope(tess, interpolation_point)
    
    n_edges = zeros(Int, length(removed_points))
    for e in delaunayedges(tess)
        a,b = geta(e), getb(e)
        for i in 1:length(removed_points)
            if n_edges[i] <3 && (a==removed_points[i] || b==removed_points[i])
                n_edges[i] += one(Int)
            end
        end
    end
    
    #== Plotting ==
    # Plot some triangles
    loc = locate(tess, interpolation_point)
    x,y = getplotxy(loc)
    trianglelayers = []
    for move in [movea, moveb, movec]
        triangle = move(tess, loc)
        x,y = getplotxy(triangle)
        push!(trianglelayers, layer(x=x, y=y, Geom.path, color=[colorant"hotpink"]))
    end

    x,y = VoronoiDelaunay.getplotxy(delaunayedges(tess))
    delaunay_layer = layer(x=x,y=y, Geom.path)
    point_layer = layer(x=[getx(p) for p in points], y=[gety(p) for p in points], color=values)
    circumcenter_layer = layer(x=[getx(c) for c in circumcenters], y=[gety(c) for c in circumcenters], color=[colorant"green"])
    plot_egdes_layer = []
    for f_edge in final_edges
        x = [getx(f_edge.a), getx(f_edge.b)]
        y = [gety(f_edge.a), gety(f_edge.b)]
        push!(plot_egdes_layer, layer(x=x, y=y, Geom.path, color=[colorant"black"]))
    end
    interp_layer = layer(x=[getx(interpolation_point)], y=[gety(interpolation_point)], shape=[Shape.xcross])
    myplot = plot(interp_layer, point_layer,circumcenter_layer,plot_egdes_layer...,delaunay_layer, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.0, maxvalue=2.0))
    draw(SVG("delaunay.svg", 8inch, 8inch), myplot)# =#


    # Compute areas T of original cells
    points_in_envelope = [f_e.a for f_e in final_edges]
    relevant_values = Array{eltype(values),1}(undef, length(points_in_envelope))
    points_in_envelope = sort_angular(points_in_envelope)
    for (i, p) in enumerate(points_in_envelope)
        idx = findall(x -> x == p, points)
        if isempty(idx)
            # Hopefully getting the intersection speeds up the nearest neighbor search
            valid_envelope = intersect(points, points_in_envelope)
            idx = [findNearest(p, valid_envelope)]
            idx = findfirst(valid_envelope[idx] .== points)
            # Check if construction of envelope was successful
            if length(unique(points_in_envelope)) < 6 || minimum(n_edges) <= 2
                idx = [findNearest(interpolation_point, valid_envelope)]
                idx = findfirst(valid_envelope[idx] .== points)
                # If not, fall back to nearest neighbor interpolation
                return values[idx[1]]
            end
        end
        relevant_values[i] = values[idx[1]]
    end
    advanced_points = Array{Point2D,1}(undef, length(points_in_envelope))
    retarded_points = Array{Point2D,1}(undef, length(points_in_envelope))
    circshift!(advanced_points, points_in_envelope, -1)
    circshift!(retarded_points, points_in_envelope, 1)
    # Midpoints between tesselation points
    m1s = [Point2D(0.5 * (getx(p) + getx(q)), 0.5 * (gety(p) + gety(q))) for (p, q) in zip(points_in_envelope, retarded_points)]
    m2s = [Point2D(0.5 * (getx(p) + getx(q)), 0.5 * (gety(p) + gety(q))) for (p, q) in zip(points_in_envelope, advanced_points)]
    # Circumcenters added by insertion of interpolation points
    g1s = Array{Point2D, 1}(undef, length(points_in_envelope))
    g2s = Array{Point2D, 1}(undef, length(points_in_envelope))
    for i in 1:length(points_in_envelope)
        p = points_in_envelope[i]
        r = retarded_points[i]
        a = advanced_points[i]
        # Sometimes the point are on a line, thus not forming a triangle. Then the 
        # next best thing to the circumcenter of all three points is the circumcenter
        # of the two extremal points. To find them, the points are projected on their common line.
        if areCollinear(p, r, interpolation_point)
            direction = [getx(p) - getx(r), gety(p)-gety(r)]
            temp_points= [p, r, interpolation_point]
            projections = [direction[1]*getx(q) + direction[2]*gety(q) for q in temp_points]
            min_idx = argmin(projections)
            max_idx = argmax(projections)
            min_point = temp_points[min_idx]
            max_point = temp_points[max_idx]
            g1s[i] = Point(0.5*(getx(min_point) + getx(max_point)),
                            0.5*(gety(min_point) + gety(max_point)))
        else
            g1s[i] = circumcenter(Primitive(p, r, interpolation_point))
        end

        if areCollinear(p, a, interpolation_point)
            direction = [getx(p) - getx(a), gety(p)-gety(a)]
            temp_points= [p, a, interpolation_point]
            projections = [direction[1]*getx(q) + direction[2]*gety(q) for q in temp_points]
            min_idx = argmin(projections)
            max_idx = argmax(projections)
            min_point = temp_points[min_idx]
            max_point = temp_points[max_idx]
            g2s[i] = Point(0.5*(getx(min_point) + getx(max_point)),
                            0.5*(gety(min_point) + gety(max_point)))
        else
            g2s[i] = circumcenter(Primitive(p, a, interpolation_point))
        end
    end
    # Area of unmodified tesselation
    Ts = zeros(Float64, length(points_in_envelope))
    # Area of modified tesselation
    As = zeros(Float64, length(points_in_envelope))
    envelope = Polygon(points_in_envelope...)
    for (j, p) in enumerate(points_in_envelope)
        relevant_cirumcenters = []
        for (i, t) in enumerate(enclosed_triangles)
            if p in [getpoint(t) for getpoint in [geta, getb, getc]]
                push!(relevant_cirumcenters, circumcenters[i])
            end
        end
        # These are not guarenteed in the right (orientied) order!
        # Since they are also never self-intersecting, they can be sorted beforehand
        relevant_points = sort_angular([p, m1s[j], relevant_cirumcenters..., m2s[j]])
        Ts[j] += getArea(relevant_points, envelope)
        relevant_points = [p, m1s[j], g1s[j], g2s[j], m2s[j]]
        As[j] += getArea(relevant_points, envelope)
    end
    w_k = abs.(Ts .- As)
    λ = w_k ./ sum(w_k)

    #== Test if λ is correct via the Local Coordinates Property ==#
    supposed_x = sum([λ[k] * getx(points_in_envelope[k]) for k in 1:length(λ)])
    supposed_y = sum([λ[k] * gety(points_in_envelope[k]) for k in 1:length(λ)])
    #println(supposed_x,",", getx(interpolation_point))
    #println(supposed_y,",", gety(interpolation_point))
    interpolated_value = sum([λ[k] * relevant_values[k] for k in 1:length(λ)])
    #println(interpolated_value)
    return interpolated_value
end

function interpolate2(points, values, interpolation_point)
    # If point is on grid of known values, return that kown value
    if interpolation_point in points
        idx = findall(x -> x == interpolation_point, points)[1]
        return values[idx]
    end
    rect = Rectangle(Point(1, 1), Point(2, 2))
    tess = voronoicells(points, rect)

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
        polygon = Polygon(cell...)
        if inpolygon(polygon, interpolation_point)
            interpolant_cell = cell
            break
        end
    end
    
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
            if inpolygon(interpolant_polygon, cellpoint) ||
               any([isBetween(edge.a, edge.b, cellpoint) for edge in interpolant_edges])
                push!(inner_points, cellpoint)
            end
        end
        # No inner points means the cell is no natural neighbor; skip it
        if length(inner_points) ==0
            continue
        end
        # Now the other way around: Find the points of interpolant_cell livin on
        # the edge of cell
        celledges = getEdges(cell)
        for point in interpolant_cell
            for edge in celledges
                if areCollinear(edge.a, edge.b, point)
                    push!(inner_points, point)
                end
            end
        end
        A = getAreaSimple(sort_angular(inner_points))
        push!(areas, A)
        push!(neighboring_cells, inner_points)
        polygon = Polygon(cell...)
        for p_i in 1:length(points)
            if inpolygon(polygon, points[p_i])
                push!(relevant_values, values[p_i])
                push!(relevant_points, points[p_i])
            end
        end
    end
    #= Plot
    Plots.plot(tess)
    for cell in neighboring_cells
        scatter!([getx(p) for p in cell], [gety(p) for p in cell])
    end
    savefig("InnerPoints.png")
    # end plot =#
    λ = areas ./ sum(areas)

    #== Test if λ is correct via the Local Coordinates Property ==#
    supposed_x = sum([λ[k] * getx(relevant_points[k]) for k in 1:length(λ)])
    supposed_y = sum([λ[k] * gety(relevant_points[k]) for k in 1:length(λ)])
    #println(supposed_x,",", getx(interpolation_point))
    #println(supposed_y,",", gety(interpolation_point))

    interpolated_value = sum([λ[k] * relevant_values[k] for k in 1:length(λ)])
    return interpolated_value
end

"""    NaturalNeighborsInterpolator(pointlist, values)

Return the interpolatant object.

`pointlist` is an array of `Point2D`, representing the Coordinates
at which `values` of the function are known.
"""
function NaturalNeighborsInterpolator(pointlist, values)
    n = length(pointlist)
    tess = DelaunayTessellation(n)
    # Pushing into tess modifies points
    push!(tess, copy(pointlist))
    interpolator = let pointlist=pointlist, values=values, tess=tess
        function interpolator(interpolation_point)
            return interpolate(pointlist, values, interpolation_point, tess)
        end
    end
    return interpolator
end

end # module