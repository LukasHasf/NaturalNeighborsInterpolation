# NaturalNeighborsInterpolation

A very work-in-progress Julia implementation of natural neighbor interpolation.

Currently not in an usable state yet.

## Usage

Prepare your known points as a vector of `Point2D`. For technical reasons in the employed package _VoronoiDelaunay.jl_, all points must be in the rectangle `[1, 2] x [1, 2]`.
```julia
using NaturalNeighborsInterpolation
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
```julia
intp = NaturalNeighborsInterpolator(points, values)
```

Use the interpolator to retrieve interpolated values at any `Point2D`
```julia
interpolated_value = intp(Point2D(1.5, 1.6))
```
![Alt text](interpolated_grid.png?raw=true "Example of using the interpolator at every gridpoint")

By default, extrapolation falls back to the nearest neighbor algorithm, but passing a `fallback` keyword argument to `NaturalNeighborsInterpolator` can change this behaviour:
```julia
intp = NaturalNeighborsInterpolator(points, values, fallback="natural")
```

![Alt text](interpolated_grid_natural.png?raw=true "Example of using the \"natural\" fallback")

```julia
intp = NaturalNeighborsInterpolator(points, values, fallback="nan")
```

![Alt text](interpolated_grid_nan.png?raw=true "Example of using the \"nan\" fallback")

