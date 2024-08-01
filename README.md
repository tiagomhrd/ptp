# PTP (Polytope)

PTP (Polytope) is a library for performing monomial integration over polytopes (so far polygons and polyhedra) using [Homogeneous Numerical Integration](https://www.sciencedirect.com/science/article/pii/S0167839620301011), and other functionalities associated with the implementation of the Virtual Element Method.

## Table of Contents

1. [Motivation](#motivation)
2. [Organization and Features](#organization-and-features)
3. [Installation Instruction](#installation-instructions)
4. [Usage and Examples](#usage-and-examples)

## Motivation

The motivation came from implementing the [Virtual Element Method](https://www.researchgate.net/publication/263884812_Basic_principles_of_Virtual_Element_Methods) (VEM), where integrals of monomials over polytopal (polygons and polyhedra) domains are required.


DO NOT CONSIDER THIS FURTHER

Although the notation was used only to enumerate the monomials, it is here extended in a way one can find other useful information for the computation of the method while using only integers, such as the order of a given monomial, its derivatives and anti-derivatives, the monomial resulting from product of monomials, etc.

This approach has also been extended to arbitrary number of variabels, mostly because it was possible, although no use is proposed for them so far.

## Organization and Features

The project is organized into three files: `mnl.hpp`, `pnl.hpp` and `glq.hpp`.
These files can be found in the `include` directory.

The focus of this project is on treating monomials individually and sequentially, as this is how the matrices for the VEM projections are usually computed.
This is the core usage of the project, which is contained in the file `mnl.hpp`.
It allows certain information to be extracted based on these indices, e.g., for the monomial of index $\alpha$ ($m_\alpha$):
- Dimension of the polynomial space of given order;
- Order of the monomial;
- Exponent for each variable;
- Index for monomial resulting from the product of two other monomials;
- Index of the monomial which is the derivative of $m_\alpha$ with respect to one of the variables;
- Index of the monomial which is the antiderivative of $m_\alpha$ with respect to one of the variables.

There are some cases in which polynomials have to be employed.
For this, a simple framework for sparse polynomial representation based on these indices and hash tables (`std::unordered_map`) has been implemented in `pnl.hpp`.
This is not the best approach for more polynomial-heavy approaches, and I would suggest looking up other frameworks if this is the case for you.
Features for this framework include:
- Defining polynomials as collections of pairs of scalar and monomial index;
- Multiplying and adding polynomials.

As the use of monomials is usually associated with their integration over some domain, this project also contains Gauss-Legendre quadrature rules for the line up to order $k=61$, which are hardcoded compactly but can be retrieved using the functions in file `glq.hpp`.
- Rule retrieval is available for the rules based on the intervals $[-1,1]$ and $[0,1]$.

The choice to restrict the quadrature rules to the 1D ones is motivated by the use of [Homogeneous Numerical Integration](https://www.sciencedirect.com/science/article/pii/S0167839620301011) on general polytopes, which is currently being implemented in the sibling project [PTP](https://github.com/tiagomhrd/ptp) for two and three dimensions.
This integration scheme is aimed at the sequential integration of monomials by performing successive applications of the [Generalized Stokes Theorem](https://en.wikipedia.org/wiki/Generalized_Stokes_theorem) along with the leveraging [Euler's Theorem for Homogeneous Functions](https://en.wikipedia.org/wiki/Homogeneous_function#Euler's_theorem) to ultimately performs these integrations in the edges of the polytope without increase of the integrand's order, requiring only that the integrals of lesser orders are known.

The quadratures were obtained with the software [Mathematica](https://www.wolfram.com/mathematica/), using function `GaussianQuadratureWeights`, and compressing (using their symmetry) results into a lookup table.

## Installation Instructions

So far, the code is implemented only in C++, and uses only `.hpp` files which can be included directly into the user's code.

One can either copy them directly from github, or clone it.
```bash
git clone https://github.com/tiagomhrd/mnl.git
```

As a design decision, the main code has been implemented in `C++11`, to reach wider audiences.
The testing project has been implemented in `C++23`.

### Dependencies

Internal dependencies: `pnl.hpp` and `glq.hpp` include `mnl.hpp`.

External dependencies: 
- `mnl.hpp` includes `<array>`;
- `pnl.hpp` includes `mnl.hpp` (`<array>`) and `<unordered_map>`;
- `glq.hpp` includes `mnl.hpp` (`<array>`) and `<vector>`;

The dependency on `<array>` is associated with the use of lookup tables.
It can be substituted in the case of `mnl.hpp` and `pnl.hpp` by providing another way of computing factorials for the computation of Monomial Orders.
It is required for `glq.hpp` as the whole thing basically consists on a big lookup table.

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

This section introduces some notation and examples that can help understand the intended use and workflow of the project.

### Monomials

In a general sense, monomials can be thought of products of powers of variables, i.e., for $m_\alpha\in P_k(\mathbb{R}^d)$ we have

$$m_\alpha=\prod_{i=0}^{d-1} x_i^{e_i}$$

where $x_i$ is the variable and $e_i$ its corresponding exponent with

$$\sum_{i=0}^{d-1}e_i\leq k$$

because it is of order at most $k$.

Obs: Zero-based indexing for the variables is used throughout the code.

### Indexing

The core idea that motivated this approach is the enumeration of monomials and being able to decode information from these indices.

Two indices are common for all spaces contained here:
- $\alpha=-1$ corresponds to $0$, i.e., $m_{-1}=0$.
- $\alpha=0$ corresponds to the constant $1$, i.e., $m_0=1$.

The other indices depend on the number of variables, e.g. (adopting $x_0=x$, $x_1=y$, etc.)
- 1D: $m_1=x$, $m_2=x^2$, $\dots$ , $m_\alpha=x^\alpha$;
- 2D: $m_1=x$, $m_2=y$, $m_3=x^2$, $m_4=xy$, $m_5=y^2$, etc.
- 3D: $m_1=x$, $m_2=y$, $m_3=z$, etc.

In the examples, the macro `REQUIRE(statement);` is used to assert a statement.

All of the examples below are included inside the `Tests.cpp` file.

### Using the code in mnl.hpp

The design consists of classes templated for the number of variables $d$, containing static functions providing the functionality.
Briefly summarized as:

```cpp
namespace mnl{
    using monIndex = int;
    using monOrder = int;

    template<const int d>
    class Poly {
        constexpr static int        SpaceDim(const monOrder k);
        constexpr static monOrder   MonOrder(const monIndex alpha);
        constexpr static int        Exponent(const monIndex alpha, const int variable);
        constexpr static monIndex   Product (const monIndex alpha, const monIndex beta);
        constexpr static monIndex   D       (const monIndex alpha, const int variable);
        constexpr static monIndex   AD      (const monIndex alpha, const int variable);
    };
}
```

Aliases for $d\in [1,10]$ are provided as `PSpacedD`, i.e., `PSpace1D`, `PSpace2D`, `PSpace3D`, etc.

Let's take for example monomials in 2D.

One can check the dimension of the polynomial space of any order
```cpp
#include "mnl.hpp"
REQUIRE(mnl::PSpace2D::SpaceDim(0) == 1);  // {1}
REQUIRE(mnl::PSpace2D::SpaceDim(1) == 3);  // {1, x, y}
REQUIRE(mnl::PSpace2D::SpaceDim(2) == 6);  // {1, x, y, x^2, xy, y^2}
REQUIRE(mnl::PSpace2D::SpaceDim(3) == 10); // {1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3}
```

One can find the monomial order for any monomial based on its index
```cpp
#include "mnl.hpp"
REQUIRE(mnl::PSpace2D::MonOrder(0) == 0);  // m_0 = 1
REQUIRE(mnl::PSpace2D::MonOrder(1) == 1);  // m_1 = x
REQUIRE(mnl::PSpace2D::MonOrder(5) == 2);  // m_5 = y^2
REQUIRE(mnl::PSpace2D::MonOrder(7) == 3);  // m_7 = x^2y
```

One can find any of the exponents given the index
```cpp
#include "mnl.hpp"
// x is direction 0, y is direction 1
// m_0 = 1 = x^0 * y^0
REQUIRE(mnl::PSpace2D::Exponent(0, 0) == 0);
REQUIRE(mnl::PSpace2D::Exponent(0, 1) == 0);
// m_1 = x = x^1 * y^0
REQUIRE(mnl::PSpace2D::Exponent(1, 0) == 1);
REQUIRE(mnl::PSpace2D::Exponent(1, 1) == 0);
// m_7 = x^2y = x^2 * y^1
REQUIRE(mnl::PSpace2D::Exponent(7, 0) == 2);
REQUIRE(mnl::PSpace2D::Exponent(7, 1) == 1);
```

One can find the index of the monomial which is the product of two monomials
```cpp
#include "mnl.hpp"
// m_5 * m_7 = y^2 * x^2y = x^2y^3 = m_18
REQUIRE(mnl::PSpace2D::Product(5, 7) == 18);
```

One can find the index of the monomial which is the derivative with respect to some direction
```cpp
#include "mnl.hpp"
// x is direction 0, y is direction 1

// m_0 = 1 = x^0 * y^0
REQUIRE(mnl::PSpace2D::D(0, 0) == -1); // d/dx(1) = 0 = m_-1
REQUIRE(mnl::PSpace2D::D(0, 1) == -1); // d/dy(1) = 0 = m_-1

// m_1 = x = x^1 * y^0
REQUIRE(mnl::PSpace2D::D(1, 0) == 0);  // d/dx(x) = 1 = m_0
REQUIRE(mnl::PSpace2D::D(1, 1) == -1); // d/dy(x) = 0 = m_-1

// m_7 = x^2y = x^2 * y^1
REQUIRE(mnl::PSpace2D::D(7, 0) == 4);  // d/dx(x^2y) ~ xy = m_4
REQUIRE(mnl::PSpace2D::D(7, 1) == 3);  // d/dy(x^2y) ~ x^2 = m_3
```

One can find the index of the antiderivative with respect to some direction
```cpp
#include "mnl.hpp"
// x is direction 0, y is direction 1

// m_0 = 1 = x^0 * y^0
REQUIRE(mnl::PSpace2D::AD(0, 0) == 1); // ADx(1) = x = m_1
REQUIRE(mnl::PSpace2D::AD(0, 1) == 2); // ADy(1) = y = m_2

// m_1 = x = x^1 * y^0
REQUIRE(mnl::PSpace2D::AD(1, 0) == 3);  // ADx(x) ~ x^2 = m_3
REQUIRE(mnl::PSpace2D::AD(1, 1) == 4);  // ADy(x) ~ xy  = m_4

// m_7 = x^2y = x^2 * y^1
REQUIRE(mnl::PSpace2D::AD(7, 0) == 11);  // ADx(x^2y) ~ x^3y   = m_11
REQUIRE(mnl::PSpace2D::AD(7, 1) == 12);  // ADy(x^2y) ~ x^2y^2 = m_12
```

### Using the code in pnl.hpp

This code is not intended for heavy use, and therefore its interface is not exactly polished.
The initialization procedure is a bit rough for the polynomial classes, but the operations with them are friendlier.

A summary of the implementation of the polynomial structures is:
```cpp
#include "mnl.hpp"
#include <unordered_map>
namespace mnl{
    template <int d>
    struct Polynomial {
        std::unordered_map<monIndex, double> Terms;

        monOrder Order() const;
        Polynomial<d>& operator*=(const Polynomial<d>& p); // Checks Zeroes
        Polynomial<d>& operator+=(const Polynomial<d>& p); // Checks Zeroes
        void CheckZeroes(); // Checks for entries that should be zero (abs value of coefficient being less than 1e-10)
    };
    template<int d>
    Polynomial<d> operator*(const Polynomial<d>& p1, const Polynomial<d>& p2); // Checks Zeroes
}
```

Aliases are provided for $d\in[1,10]$ as Polynomial\<d\> being pnldD, i.e., `pnl2D`, `pnl3D`, etc.

Usage example in the 2D case
```cpp
#include "pnl.hpp"

mnl::pnl2D p, q;
// p is x^2 + y^2
p.Terms[3] = 1.0;
p.Terms[5] = 1.0;
// q is 2 * x^3 + 3 * y^3
q.Terms[6] = 2.0;
q.Terms[9] = 3.0;

mnl::pnl2D prod = p * q;
// prod = 2 * x^5 + 3 * x^2y^3 + 2 * x^3y^2 + 3 * y^5
REQUIRE(prod.Terms[15] == 2.0); // 2 * x^5     = 2 * m_15
REQUIRE(prod.Terms[18] == 3.0); // 3 * x^2y^3  = 3 * m_18
REQUIRE(prod.Terms[17] == 2.0); // 2 * x^3y^2  = 2 * m_17
REQUIRE(prod.Terms[20] == 3.0); // 3 * y^5     = 3 * m_20

// p^2 = x^4 + 2 * x^2y^2 + y^4
p *= p;
REQUIRE(p.Terms[10] == 1.0); // 1 * x^4     = 1 * m_10
REQUIRE(p.Terms[12] == 2.0); // 2 * x^2y^2  = 2 * m_12
REQUIRE(p.Terms[14] == 1.0); // 1 * y^4     = 1 * m_14

// p^2 + q = x^4 + 2 * x^2y^2 + y^4 + 2 * x^3 + 3 * y^3
p += q;
REQUIRE(p.Terms[10] == 1.0); // 1 * x^4     = 1 * m_10
REQUIRE(p.Terms[12] == 2.0); // 2 * x^2y^2  = 2 * m_12
REQUIRE(p.Terms[14] == 1.0); // 1 * y^4     = 1 * m_14
REQUIRE(p.Terms[6]  == 2.0); // 2 * x^3     = 2 * m_6
REQUIRE(p.Terms[9]  == 3.0); // 3 * y^3     = 3 * m_9
```

### Using the code in glq.hpp

The usage here is much more straightforward.
The code summary is
```cpp
#include "mnl.hpp"
#include <vector>
namespace mnl{
    const std::vector<std::array<double, 2>> GaussLegendre(const monOrder k);
    const std::vector<std::array<double, 2>> GaussLegendreR(const monOrder k);
}
```

`GaussLegendre(const monOrder k)` returns a vector with `{position, weight}` pairs in the local system for the interval $[-1,1]$.

`GaussLegendreR(const monOrder k)` returns a vector with `{position, weight}` pairs in the local system for the interval $[0,1]$.

