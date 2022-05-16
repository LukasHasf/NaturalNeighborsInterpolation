module NaturalNeighborsInterpolation

function orientation(p, q, r)
    val = ((q[2] - p[2]) * (r[1] - q[1]) - 
            (q[1] - p[1]) * (r[2] - q[2]))
    if val==0
        return 0
    elseif val > 0
        return 1
    else
        return 2
    end
end

function distSq(p, q)
    return sum(abs2(p .- q))
end

function compare(p0, p1, p2)
    o = orientation(p0, p1, p2)
    if o==0
        if distSq(p0, p2) >= distSq(p0, p1)
            return -1
        else
            return 1
        end
    else
        if o==2
            return -1
        else
            return 1
        end
    end
end

function ConvexHull2D(pointlist)
    # Find the bottommost point. If two are found, prefer the leftmost point
    ymin = pointlist[2, 1]
    indexmin = 1
    n = size(pointlist)[2]
    for i in 1:n
        y = pointlist[2,i]
        if (y<ymin) || (ymin==y && pointlist[1,i] < pointlist[1,indexmin])
            ymin = y
            indexmin = i
        end
    end
    pointlist[:, 1], pointlist[:, indexmin] = pointlist[:, indexmin], pointlist[:, 1] 
    p0 = pointlist[:, 1]
    
    return pointlist
end

function BowyerWatson2D(pointlist)
    triangulation = []

end

function BowyerWatson(pointlist)
    if size(pointlist)[1]==2
        return BowyerWatson2D(pointlist)
    end
end

end # module