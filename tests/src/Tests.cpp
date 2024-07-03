#define _USE_MATH_DEFINES
#include <cmath>
#include "catch_amalgamated.hpp"
#include "src/pgn.h"

const double tol = std::numeric_limits<double>::epsilon() * 100.;

const std::vector<Eigen::Vector2d> regularPolygon(const int nv) {
	std::vector<Eigen::Vector2d> polygon;
	polygon.reserve(nv);
	const double halfTheta = M_PI / nv;
	const double r = pow(2. * sin(halfTheta), -1.);
	for (int i{}; i < nv; ++i)
		polygon.emplace_back(r * cos(2 * i * halfTheta), r * sin(2 * i * halfTheta));
	return polygon;
}

TEST_CASE("Polygon2D") {
	std::array<const std::vector<Eigen::Vector2d>, 4> regPolys = { regularPolygon(3), regularPolygon(4), regularPolygon(6), regularPolygon(8) };
	const auto& triangle = regPolys[0];
	const auto& square = regPolys[1];
	const auto& hexagon = regPolys[2];
	const auto& octagon = regPolys[3];
	SECTION("Triangle") {
		SECTION("Diameter") {
			REQUIRE(Polygon2D::Diameter(triangle) == Catch::Approx(1.));
		}
		SECTION("Monomials") {
			const auto monInts = Polygon2D::MonomialIntegrals(triangle, 2);
			REQUIRE(monInts[0] == Catch::Approx(sqrt(3.) / 4.));
			REQUIRE(abs(monInts[1]) < tol);
			REQUIRE(abs(monInts[2]) < tol);
			REQUIRE(monInts[3] == Catch::Approx(pow(32. * sqrt(3.), -1.)));
			REQUIRE(abs(monInts[4]) < tol);
			REQUIRE(monInts[5] == Catch::Approx(pow(32. * sqrt(3.), -1.)));
		}
	}
	SECTION("Square") {
		SECTION("Diameter") {
			REQUIRE(Polygon2D::Diameter(square) == Catch::Approx(sqrt(2.)));
		}
		SECTION("Monomials") {
			const auto monInts = Polygon2D::MonomialIntegrals(square, 2);
			REQUIRE(monInts[0] == Catch::Approx(1.));
			REQUIRE(abs(monInts[1]) < tol);
			REQUIRE(abs(monInts[2]) < tol);
			REQUIRE(monInts[3] == Catch::Approx(1. / 12.));
			REQUIRE(abs(monInts[4]) < tol);
			REQUIRE(monInts[5] == Catch::Approx(1. / 12.));
		}
	}
	SECTION("Hexagon") {
		SECTION("Diameter") {
			REQUIRE(Polygon2D::Diameter(hexagon) == Catch::Approx(2.));
		}
		SECTION("Monomials") {
			const auto monInts = Polygon2D::MonomialIntegrals(hexagon, 2);
			REQUIRE(monInts[0] == Catch::Approx(3.*sqrt(3)/2.));
			REQUIRE(abs(monInts[1]) < tol);
			REQUIRE(abs(monInts[2]) < tol);
			REQUIRE(monInts[3] == Catch::Approx(5. * sqrt(3.) / 16.));
			REQUIRE(abs(monInts[4]) < tol);
			REQUIRE(monInts[5] == Catch::Approx(5. * sqrt(3.) / 16.));
		}
	}
	SECTION("Octagon") {
		SECTION("Diameter") {
			REQUIRE(Polygon2D::Diameter(octagon) == Catch::Approx(sqrt(4. + 2. * sqrt(2.))));
		}
		SECTION("Monomials") {
			const auto monInts = Polygon2D::MonomialIntegrals(octagon, 2);
			REQUIRE(monInts[0] == Catch::Approx(2 * (sqrt(2.) +1.)));
			REQUIRE(abs(monInts[1]) < tol);
			REQUIRE(abs(monInts[2]) < tol);
			REQUIRE(monInts[3] == Catch::Approx((11. + 8. * sqrt(2.)) / 12.));
			REQUIRE(abs(monInts[4]) < tol);
			REQUIRE(monInts[5] == Catch::Approx((11. + 8. * sqrt(2.)) / 12.));
		}
	}
}

int main(int argc, char* argv[]) {
	int result = Catch::Session().run(argc, argv);
	return result;
}