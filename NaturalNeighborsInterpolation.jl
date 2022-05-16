module NaturalNeighborsInterpolation
using VoronoiDelaunay
using Gadfly
import Cairo, Fontconfig

function BowyerWatson2D(pointlist)
    n = size(pointlist)[2]
    points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:n]
    tess = DelaunayTessellation(n)
    push!(tess, points)
    x,y = getplotxy(delaunayedges(tess))
    myplot = plot(x=x,y=y, Geom.path)
    draw(PNG("delaunay.jpg", 3inch, 3inch), myplot)
end

function BowyerWatson(pointlist)
    if size(pointlist)[1]==2
        return BowyerWatson2D(pointlist)
    end
end

end # module