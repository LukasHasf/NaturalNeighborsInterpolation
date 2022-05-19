export Edge, equals, sharePoint, midPoint, getEdges, isBetween
export sort_angular, distSq, reduce, getAreaSimple

"""    Edge(a,b)

An undirected egde between `Point2D` `a` and `Point2D` `b`. 
"""
struct Edge
    a::Point2D
    b::Point2D
end

"""    equals(e1::Edge, e2::Edge)

Return `true` if `e1` and `e2` connect the same points.
"""
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

"""    midPoint(e::Edge)

Return the `Point2D` in the middle of `e`.
"""
function midPoint(e::Edge)
    return Point(0.5*(getx(e.a) + getx(e.b)), 0.5*(gety(e.a) + gety(e.b)))
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

"""    centerOfWeight(points)

Compute the center of weight of `points`.
"""
function centerOfWeight(points)
    n = length(points)
    x = [getx(p) for p in points]
    y = [gety(p) for p in points]
    return Point(sum(x)/n, sum(y)/n)
end

"""    sort_angular(points)

Sort `points` according to their angle.
"""
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

"""    distSq(p, q)

Distance squared between `Point2D`s `p` and `q`.
"""
function distSq(p, q)
    return (getx(p) - getx(q))^2 + (gety(p) - gety(q))^2
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

"""    reduce(cell)

Remove (almost) duplicate points from `cell`.
Threshold for points to be considered same is a distance less than `1e-14`.
"""
function reduce(cell)
    cell = unique(cell)
    reduced_cell = []
    checked_points = []
    for p1 in cell
        isOkay = true
        for p2 in cell
            if p1==p2
                continue
            end
            if distSq(p1, p2) < 1e-28 && !(p2 in checked_points)
                isOkay = false
            end
        end
        if isOkay
            push!(reduced_cell, p1)
        end
        push!(checked_points, p1)
    end
    return reduced_cell
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

"""    getAreaSimple(points)

Calculate the area of a simple (i.e. non-self-intersecting) polygon with vertices `points`.    
"""
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

#== Unused functions ==#

"""    isConvex(points)

Return `true` if polygon defined by `points` is convex.
"""
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

"""    polygonOrientation(points)

Return the orientation of a polygon defined by `points`.
`-1` means clockwise, `1` means counter clockwise.
"""
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

"""    findNearest(p, points)

Find nearest point in `points` to point `p`.

Returns index of nearest point.
"""
function findNearest(p, points)
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