# TrainingPhantoms

![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build status](https://github.com/JuliaImageRecon/TrainingPhantoms.jl/workflows/CI/badge.svg)](https://github.com/JuliaImageRecon/TrainingPhantoms.jl/actions)
[![codecov.io](http://codecov.io/github/JuliaImageRecon/TrainingPhantoms.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaImageRecon/TrainingPhantoms.jl?branch=master)

## Purpose
This Julia package contains a collection of methods that generate different kinds of phantom with the aim of being used to synthesise machine learning datasets. Currently, there exist two generator methods
* [Ellipsoid Phantom Generator](src/Shape.jl): ellipsoids of different size, location, orientation and value should mimic the images seen, e.g., when contrast agents are administered in the form of boluses
* [Vessel Phantom Generator](src/Vessel.jl): randomly generated blood vessels

## Installation
```julia
using Pkg
Pkg.add("TrainingPhantoms")
```

## Examples
To help you get started, we have provided examples for each generator. The can be found in the [examples](examples/) folder.
<!--- ![ellipsoidPhantom2](https://github.com/JuliaImageRecon/TrainingPhantoms.jl/assets/115639115/863b0769-9643-4858-9201-3f94311ab2ba)
![vesselPhantom](https://github.com/JuliaImageRecon/TrainingPhantoms.jl/assets/115639115/79c95f1f-b284-4562-8464-9b12ac3edf7d) --->
![combined](https://github.com/JuliaImageRecon/TrainingPhantoms.jl/assets/115639115/5babfffa-c464-4b32-99b4-da8cd4e12e86)
