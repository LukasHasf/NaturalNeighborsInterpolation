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

function isSelfIntersecting(points)
    lines = getEdges(points)
    for (line1, line2) in combinations(lines, 2)
        a, b = line1.a, line1.b
        c, d = line2.a, line2.b
        t2, t1 = twoLinesIntersect(a, b, c, d)
        println(t2,", ", t1)
        if (0 + eps(Float64) < t2 < 1 - eps(Float64)) && (0 + eps(Float64) < t1 < 1 - eps(Float64))
            println(t1)
            println(a)
            println(b)
            println(c)
            println(d)
            return true
        end
    end
    return false
end


a= Point2D(1.4931818181818182, 1.4977272727272728)
b= Point2D(1.525, 1.45)
c= Point2D(1.6, 1.5)
d= Point2D(1.65, 1.525)
e= Point2D(1.6375, 1.55)
f= Point2D(1.50625, 1.55)
points = [a,b,c,d,e,f]


points = [Point(1.5, 1.75), Point(1.5, 1.25), Point(1.75, 1.5), Point(1.25, 1.5)]
#isSelfIntersecting(points)
twoLinesIntersect(Point(1.5, 1.75), Point(1.5, 1.25), Point(1.75, 1.5), Point(1.25, 1.5))