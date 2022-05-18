using GeometricalPredicates
using Combinatorics

struct Edge
    a::Point2D
    b::Point2D
end

function getEdges(points)
    n = length(points)
    lines = [Edge(points[i], points[i+1]) for i in 1:n-1]
    push!(lines, Edge(points[n], points[1]))
    return lines
end

function twoLinesIntersect(a, b, c, d)
    a_x, a_y = getx(a), gety(a)
    b_x, b_y = getx(b), gety(b)
    c_x, c_y = getx(c), gety(c)
    d_x, d_y = getx(d), gety(d)
    bxax = (b_x - a_x)
    byay = (b_y - a_y)
    aycy = (a_y - c_y)
    cxax = (c_x - a_x)
    dycy = (d_y - c_y)
    dxcx = (d_x - c_x)
    num = (bxax * aycy + byay * cxax)
    denom = (bxax * dycy - byay * dxcx)
    t2 = num / denom
    t1 = (cxax + dxcx*t2)/bxax
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

function decomposeIntersection(points)
    lines = getEdges(points)
    intersection_points = []
    loopcounter = 0
    for (l1, l2) in combinations(lines, 2)
        loopcounter += 1
        println("Loop ", loopcounter)
        a, b = l1.a, l1.b
        c, d = l2.a, l2.b
        t, t1 = twoLinesIntersect(a, b, c, d)
        if !((0 + eps(Float64) < t < 1 - eps(Float64)) && (0 + eps(Float64) < t1 < 1 - eps(Float64)))
            continue
        end
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


points = [
 Point2D(1.7, 1.55)
 Point2D(1.65, 1.525)
 Point2D(1.5984782608695651, 1.6280434782608695)
 Point2D(1.5902631578947368, 1.455526315789474)
 Point2D(1.65, 1.5750000000000002)]

# Intersection is at (1.6375, 1.55)

decomposeIntersection(points)