using GeometricalPredicates
using Test
using Random

function dist(p, q)
    return sqrt((getx(p) - getx(q))^2 + (gety(p) - gety(q))^2)
end

function triArea(triangle)
    a,b,c = geta(triangle), getb(triangle), getc(triangle)
    ba = Point(getx(a) - getx(b), gety(a) - gety(b))
    bc = Point(getx(c) - getx(b), gety(c) - gety(b))
    area = 0.5 * (getx(ba) * gety(bc) - gety(ba) * getx(bc))
    return abs(area)
end

min_coord = 1.0 +eps(Float64)
max_coord = 2 - 2*eps(Float64)
width = max_coord - min_coord

square = [Point(min_coord, min_coord), Point(min_coord, max_coord),
            Point(max_coord, max_coord), Point(max_coord, min_coord)]
correctArea = (max_coord - min_coord)^2
calculatedArea = NaturalNeighborsInterpolation.getArea(square, nothing)
@test correctArea ≈ calculatedArea

#poly = [Point(1.5, 1.75), Point(1.5, 1.25), Point(1.75, 1.5), Point(1.25, 1.5)]
#correct_area = 1/16
#calculatedArea = NaturalNeighborsInterpolation.getArea(poly, nothing)
#@test correct_area ≈ calculatedArea
