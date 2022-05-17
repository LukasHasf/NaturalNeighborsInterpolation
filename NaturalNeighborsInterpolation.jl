module NaturalNeighborsInterpolation
using VoronoiDelaunay
using Gadfly
import Cairo, Fontconfig
using Random
using GeometricalPredicates
using Combinatorics

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


function getplotxy(t::VoronoiDelaunay.DelaunayTriangle{Point2D})
    x = [getx(getpoint(t)) for getpoint in [geta, getb, getc, geta]]
    y = [gety(getpoint(t)) for getpoint in [geta, getb, getc, geta]]
    return x,y
end

function getBowyerWatsonEnvelope(tesselation, interpolation_point)
    loc = locate(tesselation, Point(getx(interpolation_point), gety(interpolation_point)))
    triangles = Array{VoronoiDelaunay.DelaunayTriangle{Point2D},1}(undef, 4)
    triangles[1] = loc
    x,y = getplotxy(loc)
    removed_edges = [Edge(Point(x[i], y[i]), Point(x[i+1], y[i+1])) for i in 1:length(x)-1]
    envelope_edges = []
    circumcenters = Array{Point2D, 1}(undef, 4)
    circumcenters[1] = circumcenter(loc)
    for (i,move) in enumerate([movea, moveb, movec])
        triangle = move(tesselation, loc)
        triangles[i+1] = triangle
        circumcenters[i+1] = circumcenter(triangle)
        ps = Array{Point2D, 1}(undef, 4)
        for (i,getpoint) in enumerate([geta, getb, getc, geta])
            p = getpoint(triangle)
            ps[i] = p
        end
        edges = [Edge(ps[i], ps[i+1]) for i in 1:length(ps)-1]
        push!(envelope_edges, edges...)
    end
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
    return final_edges, triangles, circumcenters
end

function sort_angular(points)
    n = length(points)
    xs = [getx(p) for p in points]
    ys = [gety(p) for p in points]
    mid_x = sum(xs)/n
    mid_y = sum(ys)/n
    angles = atan.(ys .- mid_y, xs .- mid_x)
    sortind = sortperm(angles)
    points = points[sortind]
    return points
end

function isConvex(points)
    n = length(points)
    crossproducts = Array{Float64, 1}(undef, n)
    for i in 1:n
        x1, y1 = getx(points[i]), gety(points[i])
        if i<n-1
            x2, y2 = getx(points[i+1]), gety(points[i+1])
            x3, y3 = getx(points[i+2]), gety(points[i+2])
        elseif i==n-1
            x2, y2 = getx(points[i+1]), gety(points[i+1])
            x3, y3 = getx(points[1]), gety(points[1])
        elseif i==n
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


Returns a parameter `t` which describes how far along the second line the intersection is. 
Line segments intersect if `0 < t < 1`.
"""
function twoLinesIntersect(a,b,c,d)
    a_x, a_y = getx(a), gety(a)
    b_x, b_y = getx(b), gety(b)
    c_x, c_y = getx(c), gety(c)
    d_x, d_y = getx(d), gety(d)
    bxax = (b_x-a_x)
    byay = (b_y-a_y)
    aycy = (a_y-c_y)
    cxax = (c_x-a_x)
    dycy = (d_y-c_y)
    dxcx = (d_x-c_x)
    num =   (bxax*aycy + byay*cxax)
    denom = (bxax*dycy - byay*dxcx)
    t =  num / denom 
    return t
end

function isSelfIntersecting(points)
    lines = getEdges(points)
    for line1 in lines
        for line2 in lines
            if equals(line1, line2)
                continue
            end
            a, b = line1.a, line1.b
            c, d = line2.a, line2.b
            t = twoLinesIntersect(a,b,c,d)
            if 0+eps(Float64) < t < 1-eps(Float64)
                return true
            end
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
    for (l1,l2) in combinations(lines, 2)
            if equals(l1, l2)
                continue
            end
            a,b = l1.a, l1.b
            c,d = l2.a, l2.b
            t = twoLinesIntersect(a,b,c,d)
            if !(0+eps(Float64) < t < 1-eps(Float64))
                continue
            end
            a = [getx(a), gety(a)]
            b = [getx(b), gety(b)]
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
    display(newpoints)
    display(cross_indices)
    display(intersection_points)


    return [points]
end

"""    getArea(points)


Calculate the area of a polygon.
"""
function getArea(points)
    A = zero(Float64)
    if !isSelfIntersecting(points)
        points = sort_angular(points)
        n = length(points)
        xs = [getx(p) for p in points]
        ys = [gety(p) for p in points]
        mid_x = sum(xs)/n
        mid_y = sum(ys)/n
        midpoint = Point(mid_x, mid_y)
        for q in 1:length(points)-1
            temp_triangle = Primitive(midpoint, points[q], points[q+1])
            A += area(temp_triangle)
        end
        A += area(Primitive(midpoint, points[1], points[end]))
    else
        polygons = decomposeIntersection(points)
        for poly in polygons
            A += getArea(poly)
        end
    end
    return A
end

function interpolate(pointlist, values, interpolation_point)
    n = size(pointlist)[2]
    points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:n]
    if interpolation_point in points
        idx = findall(x->x==interpolation_point, points)[1]
        return values[idx]
    end
    tess = DelaunayTessellation(n)
    push!(tess, copy(points))
    

    # Plot some triangles
    loc = locate(tess, interpolation_point)
    x,y = getplotxy(loc)
    trianglelayers = []
    for move in [movea, moveb, movec]
        triangle = move(tess, loc)
        x,y = getplotxy(triangle)
        push!(trianglelayers, layer(x=x, y=y, Geom.path, color=[colorant"hotpink"]))
    end
    final_edges, enclosed_triangles, circumcenters = getBowyerWatsonEnvelope(tess, interpolation_point)
    
    #== Plotting ==#

    x,y = VoronoiDelaunay.getplotxy(delaunayedges(tess))
    delaunay_layer = layer(x=x,y=y, Geom.path)
    point_layer = layer(x=pointlist[1,:], y=pointlist[2,:], color=values)
    circumcenter_layer = layer(x=[getx(c) for c in circumcenters], y=[gety(c) for c in circumcenters], color=[colorant"green"])
    plot_egdes_layer = []
    for f_edge in final_edges
        x = [getx(f_edge.a), getx(f_edge.b)]
        y = [gety(f_edge.a), gety(f_edge.b)]
        push!(plot_egdes_layer, layer(x=x, y=y, Geom.path, color=[colorant"black"]))
    end
    myplot = plot(point_layer,circumcenter_layer,plot_egdes_layer...,delaunay_layer, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.0, maxvalue=2.0))
    draw(SVG("delaunay.svg", 8inch, 8inch), myplot)
    
    
    # Compute areas T of original cells
    points_in_envelope = [f_e.a for f_e in final_edges]
    relevant_values = Array{eltype(values), 1}(undef, length(points_in_envelope))
    points_in_envelope = sort_angular(points_in_envelope) # stable
    for (i,p) in enumerate(points_in_envelope)
        idx = findall(x->x==p, points)
        relevant_values[i] = values[idx[1]]
    end
    advanced_points = circshift(points_in_envelope, -1)
    retarded_points = circshift(points_in_envelope, 1)
    # Midpoints between tesselation points
    m1s = [Point2D(0.5*(getx(p) + getx(q)), 0.5*(gety(p) + gety(q))) for (p,q) in zip(points_in_envelope, retarded_points)]
    m2s = [Point2D(0.5*(getx(p) + getx(q)), 0.5*(gety(p) + gety(q))) for (p,q) in zip(points_in_envelope, advanced_points)]
    # Circumcenters added by insertion of interpolation points
    g1s = [circumcenter(Primitive(points_in_envelope[i], retarded_points[i], interpolation_point)) for i in 1:length(points_in_envelope)]
    g2s = [circumcenter(Primitive(points_in_envelope[i], advanced_points[i], interpolation_point)) for i in 1:length(points_in_envelope)]
    # Area of unmodified tesselation
    Ts = zeros(Float64, length(points_in_envelope))
    # Area of modified tesselation
    As = zeros(Float64, length(points_in_envelope))
    for (j,p) in enumerate(points_in_envelope)
        relevant_cirumcenters = []
        for (i,t) in enumerate(enclosed_triangles)
            if p in [getpoint(t) for getpoint in [geta, getb, getc]]
                push!(relevant_cirumcenters, circumcenters[i])
            end
        end
        # These are not guarenteed in the right (orientied) order!
        relevant_points = [p, m1s[j], relevant_cirumcenters..., m2s[j]] # correct
        Ts[j] += getArea(relevant_points)
        relevant_points = [p,m1s[j], g1s[j], g2s[j], m2s[j]] # consistent
        display(relevant_points)
        display(isSelfIntersecting(relevant_points))
        As[j] += getArea(relevant_points)
    end
    w_k = Ts .- As
    λ = w_k ./ sum(w_k)
    display(As)

    #== Test if λ is correct via the Local Coordinates Property ==#
        supposed_x = sum([λ[k] * getx(points_in_envelope[k]) for k in 1:length(λ)])
        supposed_y = sum([λ[k] * gety(points_in_envelope[k]) for k in 1:length(λ)])
        println(supposed_x,",", getx(interpolation_point))
        println(supposed_y,",", gety(interpolation_point))
        interpolated_value = sum([λ[k] * relevant_values[k] for k in 1:length(λ)])
        println(interpolated_value)
end

function BowyerWatson(pointlist)
    if size(pointlist)[1]==2
        return BowyerWatson2D(pointlist)
    end
end

end # module