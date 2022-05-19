# NaturalNeighborsInterpolation

A very work-in-progress Julia implementation of natural neighbor interpolation.

Currently not in an usable state yet.

## Usage

Prepare your known points as a vector of `Point2D`. For technical reasons in the employed package _VoronoiDelaunay.jl_, all points must be in the rectangle `[1, 2] x [1, 2]`.
```julia
pointlist = [0.45 0.6 0.4 0.6 0.4 0.6 0.55 0.7  0.4; 
             0.4  0.4 0.6 0.6 0.55 0.5 0.4  0.55 0.4] .+ 1
values = [pointlist[1,i] + pointlist[2,i] for i in 1:size(pointlist)[2]]
points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:size(pointlist)[2]]
```

Prepare your known values also a vector
```julia
values = [pointlist[1,i] + pointlist[2,i] for i in 1:size(pointlist)[2]]
```

Get the interpolator function
```
intp = NaturalNeighborsInterpolation.NaturalNeighborsInterpolator(points, values)
```

Use the interpolator to retrive interpolated values at any `Point2D`
```julia
interpolated_value = intp(Point2D(1.5, 1.6))
```
---
Note: Values outside the interpolation region are mapped to the value of their nearest neighbor.
