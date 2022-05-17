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

    n = length(points)
    isIntersecting = false
    lines = [Edge(points[i], points[i+1]) for i in 1:n-1]
    push!(lines, Edge(points[n], points[1]))
    display(lines)
    for line1 in lines
        for line2 in lines
            if equals(line1, line2)
                continue
            end
            a, b = line1.a, line1.b
            c,d = line2.a, line2.b
            #display(line1)
            #display(line2)
            t = twoLinesIntersect(a,b,c,d)
            if getx(a)==1.4916666666666667
                println(a,b,c,d)
                println(points)
                println(t)
            end
            if 0+eps(Float64) < t < 1-eps(Float64)
                return true
            end
        end
    end
    return false
end


a = Point(1.492, 1.608)
b = Point(1.425, 1.475)
c = Point(1.4, 1.55)
d = Point(1.5, 1.575)

points = [a,b,c,d]
n = length(points)
lines = [Edge(points[i], points[i+1]) for i in 1:n-1]
push!(lines, Edge(points[n], points[1]))
#display(lines)
isSelfIntersecting(points)