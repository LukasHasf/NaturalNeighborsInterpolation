module NaturalNeighborsInterpolation
using VoronoiDelaunay
using Gadfly
import Cairo, Fontconfig
using Random

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
    loc = locate(tesselation, Point(interpolation_point[1], interpolation_point[2]))
    for move in [movea, moveb, movec]
        triangle = move(tesselation, loc)
        for getpoint in [geta, getb, getc]
            p = getpoint(triangle)
            if !(p in points_of_envelope)
                push!(points_of_envelope, p)
            end
        end
    end
end

function interpolate(pointlist)
    n = size(pointlist)[2]
    points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:n]
    tess = DelaunayTessellation(n)
    push!(tess, points)
    x,y = VoronoiDelaunay.getplotxy(delaunayedges(tess))
    delaunay_layer = layer(x=x,y=y, Geom.path)
    point_layer = layer(x=pointlist[1,:], y=pointlist[2,:])
    

    # Find the triangle where the interpolated point is located
    loc = locate(tess, Point(1.5, 1.5))
    x,y = getplotxy(loc)
    removed_edges = [Edge(Point(x[i], y[i]), Point(x[i+1], y[i+1])) for i in 1:length(x)-1]
    envelope_prestep_layer = layer(x=x, y=y, Geom.path, color=[colorant"red"])
    trianglelayers = []
    envelope_edges = []
    for move in [movea, moveb, movec]
        triangle = move(tess, loc)
        x,y = getplotxy(triangle)
        push!(trianglelayers, layer(x=x, y=y, Geom.path, color=[colorant"hotpink"]))
        ps = Array{Point2D, 1}(undef, 4)
        for (i,getpoint) in enumerate([geta, getb, getc, geta])
            ps[i] = getpoint(triangle)
        end
        edges = [Edge(ps[i], ps[i+1]) for i in 1:length(ps)-1]
        push!(envelope_edges, edges...)
    end
    display(envelope_edges)
    display(removed_edges)
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
    display(final_edges)
    plot_egdes_layer = []
    for f_edge in final_edges
        x = [getx(f_edge.a), getx(f_edge.b)]
        y = [gety(f_edge.a), gety(f_edge.b)]
        push!(plot_egdes_layer, layer(x=x, y=y, Geom.path, color=[colorant"black"]))
    end
    envelope_points_layer = layer(x=x, y=y, shape=[Shape.square], color=[1,1.5,2,4,5,6])
    myplot = plot(envelope_points_layer ,point_layer,plot_egdes_layer...,delaunay_layer, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.0, maxvalue=2.0))
    draw(SVG("delaunay.svg", 8inch, 8inch), myplot)
end

function BowyerWatson(pointlist)
    if size(pointlist)[1]==2
        return BowyerWatson2D(pointlist)
    end
end

end # module