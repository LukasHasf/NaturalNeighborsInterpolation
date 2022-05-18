using GeometricalPredicates

a = Point(1.492, 1.608)
b = Point(1.425, 1.475)
c = Point(1.4, 1.55)
d = Point(1.5, 1.575)

points = [a,b,c,d]
n = length(points)
lines = [NaturalNeighborsInterpolation.Edge(points[i], points[i+1]) for i in 1:n-1]
push!(lines, NaturalNeighborsInterpolation.Edge(points[n], points[1]))
@test NaturalNeighborsInterpolation.isSelfIntersecting(points) == true