using GeometricalPredicates

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
    #println(bxax) # ok
    #println(byay) #ok
    #println(aycy) # ok
    #println(cxax) # ok
    #println(dycy) # ok
    #println(dxcx) # ok
    #println(num) # ok
    #println(denom) # ok
    t =  num / denom 
    println(t)
    return 0 < t <1
end

a = Point(1.492, 1.608)
b = Point(1.425, 1.475)
c = Point(1.4, 1.55)
d = Point(1.5, 1.575)

twoLinesIntersect(a,b,c,d)