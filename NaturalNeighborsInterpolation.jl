module NaturalNeighborsInterpolation
using VoronoiDelaunay

function BowyerWatson2D(pointlist)
    n = size(pointlist)[2]
    points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:n]
    tess = DelaunayTessellation(n)
    push!(tess, points)
    for edge in delaunayedges(tess)
        display(edge)
    end
end

function BowyerWatson(pointlist)
    if size(pointlist)[1]==2
        return BowyerWatson2D(pointlist)
    end
end

end # module