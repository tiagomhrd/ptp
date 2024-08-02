# Ptp (Polytope)

Ptp (Polytope) is a library for performing monomial integration over polytopes (so far polygons and polyhedra) using [Homogeneous Numerical Integration](https://www.sciencedirect.com/science/article/pii/S0167839620301011), and other functionalities associated with the implementation of the Virtual Element Method.

## Table of Contents

1. [Motivation](#motivation)
2. [Organization and Features](#organization-and-features)
3. [Installation Instruction](#installation-instructions)
4. [Usage and Examples](#usage-and-examples)

## Motivation

The motivation came from implementing the [Virtual Element Method](https://www.researchgate.net/publication/263884812_Basic_principles_of_Virtual_Element_Methods) (VEM), where integrals of monomials over polytopal (polygons and polyhedra) domains are required.

The representation of monomials follows the indexing notation proposed and implemented in [Mnl](https://github.com/tiagomhrd/mnl).

Apart from integrals, other information regarding polynomials might be required, such as the number of noncolinear lines/noncoplanar faces, etc.

This library consists in a lightweight way to compute these quantities.

## Organization and Features

The project contains a header (`ptp.h`) and a source (`ptp.cpp`) file, both located in `ptp/src/`.

A summary of the features can be found in the header file itself, reproduced below

```cpp
namespace ptp {
	namespace Polygon2D {
		const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector2d>& vertices, const int maxOrder);
		const double Diameter(const std::vector<Eigen::Vector2d>& vertices);
		const std::vector<Eigen::Vector3d> UniqueSides(const std::vector<Eigen::Vector2d>& vertices);
		const std::vector<Eigen::Vector3d> UniqueReentrantSides(const std::vector<Eigen::Vector2d>& vertices);
	}

	namespace Polygon3D {
		const std::vector<Eigen::Vector3d> GetVertices(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices);
		const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector3d>& vertices, int maxOrder);
		const double Diameter(const std::vector<Eigen::Vector3d>& vertices);
		const Eigen::Vector3d Normal(const std::vector<Eigen::Vector3d>& vertices);
		const Eigen::Vector4d PlaneEquation(const std::vector<Eigen::Vector3d>& vertices);
	}

	namespace Polyhedron {
		const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces, int maxOrder);
		const double Diameter(const std::vector<Eigen::Vector3d>& vertices);
		const std::vector<Eigen::Vector4d> UniquePlanes(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces);
	}

	namespace Transform {
		void TranslateVertices(std::vector<Eigen::Vector3d>& vertices, const Eigen::Vector3d& translation);
		void RotateVertices(std::vector<Eigen::Vector3d>& vertices, const Eigen::Matrix3d& Q);
	}
}
```

Points are explicitly treated as `Eigen::Vector2d` and `Eigen::Vector3d`, although in the future this dependency might be hidden or substituted.

All `MonomialIntegrals` functions return a vector of doubles, which are the integrals of monomials following the indexing notation mentioned before.
For `Polygon2D` these are the monomials over 2 variables, and for `Polygon3D` and `Polyhedron` those over 3 variables.
The `Diameter` functions return the largest distance between two vertices in the polytope.

`Polgon2D::UniqueSides` and `Polygon2D::UniqueReentrantSides` both return a list of `Eigen::Vector3d` representing the coefficients `{a, b, c}` for the line defined by `a + bx + cy = 0`.
Unique sides represent the support lines for the boundary of the `Polygon2D`, and Unique reentrant sides represents the support lines that define concavities.

The `Polygon3D::PlaneEquation` function returns an `Eigen::Vector4d` containing the scalar coefficients `{a,b,c,d}` for the plane equation in the format `a + bx + cy + dz = 0`.
This is, in turn, used by `Polyhedron::UniquePlanes` in similar manner to what was done for `Polygon2D`

This boundary information is necessary for high-order serendipity formulations of the virtual element method.

## Installation Instructions

The whole project is written in C++, and the core files are compatible with `C++17`.

Inside the `ptp` folder you can find premake (`premake.lua`) and CMake (`CMakeLists.txt`) files that can generate project files to compile this to a static library.

Feel free to clone the project or add as a submodule to your repository, but don't forget to do so recursively, as [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [Mnl](https://github.com/tiagomhrd/mnl) are submodules.

```bash
git clone --recursive https://github.com/tiagomhrd/ptp.git
```

### Dependencies

For external dependencies, this relies on `<algorithm>`, `Eigen` and `Mnl`.

### Development

For development, the project also supports two build systems: [CMake](https://cmake.org/) and [Premake](https://premake.github.io/).
The latter with its executable included in the project files inside the `third_party` directory along with its license.

These build systems will set the `tests` project as the main project, this uses the [Catch2](https://github.com/catchorg/Catch2) framework for testing.
The code for Catch2 is also included in its amalgamated version inside `third_party` (along with its License).
The version used is `v3.6.0`.

These tests can be useful for learning purposes.

To generate project files using Premake, just edit `GenerateProjects.bat` to your editor of choice using the [options provided](https://premake.github.io/docs/Using-Premake#using-premake-to-generate-project-files) and run it.
To avoid run a `.bat` file, you can also run the command directly on you command prompt:
```bash
\mnl> third_party\premake5.exe [action]
```

To generate them with CMake, do the following:
```bash
\mnl> mkdir build
\mnl> cd build
\mnl\build> cmake ..
```

## Usage and Examples

The use of the library is very direct.

### Two-dimensions
In 2D, points are represented by `Eigen::Vector2d`'s.

Polygons are represented by `std::vector<Eigen::Vector2d>` **ordered counterclockwise**, i.e., satisfying the [right-hand rule](https://en.wikipedia.org/wiki/Right-hand_rule).
This ordering is important, and is what guarantees you don't have to provide another list of indices.

### Three-dimensions

In 3D, points are represented by `Eigen::Vector3d`'s.

Polygons are represented by `std::vector<Eigen::Vector3d>` **ordered counterclockwise** relative to the external normal to the plane, following the [right-hand rule](https://en.wikipedia.org/wiki/Right-hand_rule).

For uses of `Polygon3D` namespace for faces of polyhedra, one can compose the function calls with the function
```cpp
const std::vector<Eigen::Vector3d> Polygon3D::GetVertices(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices);
```
For this function, you input a list of 3D points and a list of ordered indices corresponding to the face, and get the ordered list of vertices representing the face.

Polyhedra are represented by the list of 3D points and a list of list of ordered face indices.

There are also some additional utility functions for transforming list of points with translations (`Transform::TranslateVertices`) and rotation (`Transform::RotateVertices`), which are used internally but made available for the testing framework.

### Examples

The recommended examples to understand usage are present in `tests/src/Tests.cpp`.

There, analytical results for inertias of polygons in 2 and 3 dimensions and polyhedra are used to verify the monomial integrals.
The diameter functions are also verified.

As a representative example, let us show the test for the unit cube:

```cpp
struct Polyhedron {
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<size_t>> faces;
};

static Polyhedron cube() {
    Polyhedron poly;
    
    poly.vertices.reserve(8);
    poly.vertices.emplace_back(-0.5, -.5, -.5);
    poly.vertices.emplace_back(+0.5, -.5, -.5);
    poly.vertices.emplace_back(+0.5, +.5, -.5);
    poly.vertices.emplace_back(-0.5, +.5, -.5);
    poly.vertices.emplace_back(-0.5, -.5, +.5);
    poly.vertices.emplace_back(+0.5, -.5, +.5);
    poly.vertices.emplace_back(+0.5, +.5, +.5);
    poly.vertices.emplace_back(-0.5, +.5, +.5);

    poly.faces.reserve(6);
    poly.faces.push_back({0, 4, 7, 3}); // -x
    poly.faces.push_back({1, 2, 6, 5}); // +x
    poly.faces.push_back({0, 1, 5, 4}); // -y
    poly.faces.push_back({2, 3, 7, 6}); // +y
    poly.faces.push_back({0, 3, 2, 1}); // -z
    poly.faces.push_back({4, 5, 6, 7}); // +z

    return poly;
};
static Eigen::Matrix3d cubeInertiaTensor() {
    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    const double vol_12 = 1. / 12.;
    I(0, 0) = vol_12 * (_pow(dy,2) + _pow(dz,2));
    I(1, 1) = vol_12 * (_pow(dx,2) + _pow(dz,2));
    I(2, 2) = vol_12 * (_pow(dx,2) + _pow(dy,2));
    return I;
}

const double tol = 1e-10;
Polyhedron cube = cube();
Eigen::Matrix3d I = cubeInertiaTensor();
const double diameter = sqrt(3.);
const double vol = 1.;

REQUIRE_THAT(ptp::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));

const auto monInts = ptp::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(0.0, tol));
REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(0.0, tol));
REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(0.0, tol));
REQUIRE_THAT(monInts[5], Catch::Matchers::WithinAbs(0.0, tol));
REQUIRE_THAT(monInts[6], Catch::Matchers::WithinAbs(0.0, tol));
REQUIRE_THAT(monInts[8], Catch::Matchers::WithinAbs(0.0, tol));
REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));	
```
