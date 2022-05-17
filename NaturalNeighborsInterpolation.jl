module NaturalNeighborsInterpolation
using VoronoiDelaunay
using Gadfly
import Cairo, Fontconfig
using Random
using GeometricalPredicates

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
            ps[i] = getpoint(triangle)
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

function getArea(points)
    # TODO: Add test
    n = length(points)
    xs = [getx(p) for p in points]
    ys = [gety(p) for p in points]
    mid_x = sum(xs)/n
    mid_y = sum(ys)/n
    midpoint = Point(mid_x, mid_y)
    angles = atan.(ys .- mid_y, xs .- mid_x)
    sortind = sortperm(angles)
    points = points[sortind]
    A = zero(Float64)
    for q in 1:length(points)-1
        temp_triangle = Primitive(midpoint, points[q], points[q+1])
        A += area(temp_triangle)
    end
    A += area(Primitive(midpoint, points[1], points[end]))
    return A
end

function interpolate(pointlist, values, interpolation_point)
    n = size(pointlist)[2]
    points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:n]
    tess = DelaunayTessellation(n)
    push!(tess, points)
    

    # Find the triangle where the interpolated point is located
    loc = locate(tess, interpolation_point)
    x,y = getplotxy(loc)
    trianglelayers = []
    for move in [movea, moveb, movec]
        triangle = move(tess, loc)
        x,y = getplotxy(triangle)
        push!(trianglelayers, layer(x=x, y=y, Geom.path, color=[colorant"hotpink"]))
        ps = Array{Point2D, 1}(undef, 4)
        for (i,getpoint) in enumerate([geta, getb, getc, geta])
            ps[i] = getpoint(triangle)
        end
    end
    final_edges, enclosed_triangles, circumcenters = getBowyerWatsonEnvelope(tess, interpolation_point)
    # Compute areas T of original cells
    points_in_envelope = [f_e.a for f_e in final_edges]
    relevant_values = Array{eltype(values), 1}(undef, length(points_in_envelope))
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
        relevant_points = [p, m1s[j], relevant_cirumcenters..., m2s[j]]
        Ts[j] += getArea(relevant_points)
        relevant_points = [p,m1s[j], g1s[j], g2s[j], m2s[j]]
        As[j] += getArea(relevant_points)
    end
    display(Ts)
    display(As)
    w_k = Ts .- As
    λ = w_k ./ sum(w_k)
    display(λ)

    #== Test if λ is correct via the Local Coordinates Property ==#
        supposed_x = sum([λ[k] * getx(points_in_envelope[k]) for k in 1:length(λ)])
        supposed_y = sum([λ[k] * gety(points_in_envelope[k]) for k in 1:length(λ)])
        println(supposed_x,",", getx(interpolation_point))
        println(supposed_y,",", gety(interpolation_point))
        interpolated_value = sum([λ[k] * relevant_values[k] for k in 1:length(λ)])
        println(interpolated_value)
    #== Plotting ==#

    x,y = VoronoiDelaunay.getplotxy(delaunayedges(tess))
    delaunay_layer = layer(x=x,y=y, Geom.path)
    point_layer = layer(x=pointlist[1,:], y=pointlist[2,:], color=values)
    plot_egdes_layer = []
    for f_edge in final_edges
        x = [getx(f_edge.a), getx(f_edge.b)]
        y = [gety(f_edge.a), gety(f_edge.b)]
        push!(plot_egdes_layer, layer(x=x, y=y, Geom.path, color=[colorant"black"]))
    end
    myplot = plot(point_layer,plot_egdes_layer...,delaunay_layer, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.0, maxvalue=2.0))
    draw(SVG("delaunay.svg", 8inch, 8inch), myplot)
end

function BowyerWatson(pointlist)
    if size(pointlist)[1]==2
        return BowyerWatson2D(pointlist)
    end
end

end # module