#define _USE_MATH_DEFINES
#include <cmath>
#include "catch_amalgamated.hpp"
#include "src/pgn.h"
#include <iostream>

const double tol = std::numeric_limits<double>::epsilon() * 100.;
static double _pow(double base, int exp)
{
	double result = 1.;
	for (;;)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		if (!exp)
			break;
		base *= base;
	}

	return result;
}

static const std::vector<Eigen::Vector2d> regularPolygon(const int nv) {
	std::vector<Eigen::Vector2d> polygon;
	polygon.reserve(nv);
	const double halfTheta = M_PI / nv;
	const double r = pow(2. * sin(halfTheta), -1.);
	for (int i{}; i < nv; ++i)
		polygon.emplace_back(r * cos(2 * i * halfTheta), r * sin(2 * i * halfTheta));
	return polygon;
}

static double regPolyInertia(const int nv) {
	const double cot = cos(M_PI / nv) / sin(M_PI / nv);
	return nv / 96. * cot * (1 + 3 * cot * cot);
}

static Eigen::Matrix3d regPolyInertiaTensor(const int nv) {
	Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
	const double i = regPolyInertia(nv);
	I(0, 0) = I(1, 1) = .5 * i;
	I(2, 2) = i;
	return I;
}

static double regPolyArea(const int nv) {
	return nv * .25 * cos(M_PI / nv) / sin(M_PI / nv);
}

static double regPolyDiameter(const int nv) {
	const int halfway = nv / 2;
	const double halfTheta = M_PI / nv;
	const double r = pow(2. * sin(halfTheta), -1.);
	return r * sqrt(_pow(cos(halfway * halfTheta * 2) - 1, 2) + _pow(sin(halfway * halfTheta * 2), 2));
}



static const Eigen::Matrix3d RodriguesTensor(const Eigen::Vector3d& axis, const double angleInRadians) {
	Eigen::Matrix3d skewTheta;
	skewTheta << 0, -axis[2], axis[1], axis(2), 0, -axis(0), -axis(1), axis(0), 0.0;
	skewTheta *= angleInRadians;
	const double h1 = sin(angleInRadians) / angleInRadians;
	const double h2 = 0.5 * _pow(sin(0.5 * angleInRadians) / (.5 * angleInRadians), 2);
	return Eigen::Matrix3d::Identity() + h1 * skewTheta + h2 * skewTheta * skewTheta;
}

static std::vector<Eigen::Vector3d> regularPolygon3d(const int nv) {
	std::vector<Eigen::Vector3d> polygon;
	polygon.reserve(nv);
	const double halfTheta = M_PI / nv;
	const double r = pow(2. * sin(halfTheta), -1.);
	for (int i{}; i < nv; ++i)
		polygon.emplace_back(r * cos(2 * i * halfTheta), r * sin(2 * i * halfTheta), 0.0);
	return polygon;
}

static void rotatePolygon3D(std::vector<Eigen::Vector3d>& polygon, const Eigen::Matrix3d& rotMatrix) {
	std::transform(polygon.cbegin(), polygon.cend(), polygon.begin(), [&rotMatrix](const auto& pos) { return rotMatrix * pos; });
}

static void translatePolygon3D(std::vector<Eigen::Vector3d>& polygon, const Eigen::Vector3d& trnslVec) {
	std::transform(polygon.cbegin(), polygon.cend(), polygon.begin(), [&trnslVec](const auto& pos) { return pos + trnslVec; });
}

TEST_CASE("Polygon2D") {
	std::stringstream ss;
	for (int nv = 3, max = 15; nv <= max; ++nv) {
		ss << nv << "-gon"; 
		SECTION(ss.str()) {
			const auto poly = regularPolygon(nv);
			SECTION("Diameter") {
				REQUIRE(Polygon2D::Diameter(poly) == Catch::Approx(regPolyDiameter(nv)));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = Polygon2D::MonomialIntegrals(poly, 2); 
				REQUIRE(monInts[0] == Catch::Approx(regPolyArea(nv))); 
				REQUIRE(abs(monInts[1]) < tol); 
				REQUIRE(abs(monInts[2]) < tol); 
				const double sqInts = .5 * regPolyInertia(nv); 
				REQUIRE(monInts[3] == Catch::Approx(sqInts)); 
				REQUIRE(abs(monInts[4]) < tol); 
				REQUIRE(monInts[5] == Catch::Approx(sqInts)); 
			}
		}
		ss.str("");
	}
}

TEST_CASE("Polygon3D") {
	std::stringstream ss;
	const Eigen::Vector3d translation(1., 2., 3.);
	const Eigen::Matrix3d Q = RodriguesTensor(Eigen::Vector3d(sqrt(2.) * .5, sqrt(2.) * .5, 0.), M_PI_4);
	for (int nv = 3, max = 15; nv <= max; ++nv) {
		ss << nv << "-gon";
		SECTION(ss.str()) {
			auto poly = regularPolygon3d(nv);
			Eigen::Matrix3d I = regPolyInertiaTensor(nv);
			const double area = regPolyArea(nv);
			SECTION("No Transformation") {
				SECTION("Diameter") {
					REQUIRE(Polygon3D::Diameter(poly) == Catch::Approx(regPolyDiameter(nv)));
				}
				SECTION("Monomial Integrals") {
					const auto monInts = Polygon3D::MonomialIntegrals(poly, 2);
					REQUIRE(monInts[0] == Catch::Approx(area));
					REQUIRE(abs(monInts[1]) < tol);
					REQUIRE(abs(monInts[2]) < tol);
					REQUIRE(abs(monInts[3]) < tol);
					REQUIRE(abs(monInts[5]) < tol);
					REQUIRE(abs(monInts[6]) < tol);
					REQUIRE(abs(monInts[8]) < tol);
					REQUIRE(I(0, 0) == Catch::Approx(monInts[7] + monInts[9]));
					REQUIRE(I(1, 1) == Catch::Approx(monInts[4] + monInts[9]));
					REQUIRE(I(2, 2) == Catch::Approx(monInts[4] + monInts[7]));
				}
			}
			
			SECTION("Translation") {
				translatePolygon3D(poly, translation);
				I += area * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());

				SECTION("Diameter") {
					REQUIRE(Polygon3D::Diameter(poly) == Catch::Approx(regPolyDiameter(nv)));
				}				
				SECTION("Monomial Integrals") {
					const auto monInts = Polygon3D::MonomialIntegrals(poly, 2);
					REQUIRE(monInts[0] == Catch::Approx(area));
					REQUIRE(monInts[1] == Catch::Approx(translation(0) * area));
					REQUIRE(monInts[2] == Catch::Approx(translation(1) * area));
					REQUIRE(monInts[3] == Catch::Approx(translation(2) * area));
					REQUIRE(I(0, 1) == Catch::Approx(-monInts[5]));
					REQUIRE(I(0, 2) == Catch::Approx(-monInts[6]));
					REQUIRE(I(1, 2) == Catch::Approx(-monInts[8]));
					REQUIRE(I(0, 0) == Catch::Approx(monInts[7] + monInts[9]));
					REQUIRE(I(1, 1) == Catch::Approx(monInts[4] + monInts[9]));
					REQUIRE(I(2, 2) == Catch::Approx(monInts[4] + monInts[7]));
				}
			}
			SECTION("Rotation") {
				rotatePolygon3D(poly, Q);
				I = Q * I * Q.transpose();

				SECTION("Diameter") {
					REQUIRE(Polygon3D::Diameter(poly) == Catch::Approx(regPolyDiameter(nv)));
				}
				SECTION("Monomial Integrals") {
					const auto monInts = Polygon3D::MonomialIntegrals(poly, 2);
					REQUIRE(monInts[0] == Catch::Approx(area));
					REQUIRE(abs(monInts[1]) < tol);
					REQUIRE(abs(monInts[2]) < tol);
					REQUIRE(abs(monInts[3]) < tol);
					REQUIRE(I(0, 1) == Catch::Approx(-monInts[5]));
					REQUIRE(I(0, 2) == Catch::Approx(-monInts[6]));
					REQUIRE(I(1, 2) == Catch::Approx(-monInts[8]));
					REQUIRE(I(0, 0) == Catch::Approx(monInts[7] + monInts[9]));
					REQUIRE(I(1, 1) == Catch::Approx(monInts[4] + monInts[9]));
					REQUIRE(I(2, 2) == Catch::Approx(monInts[4] + monInts[7]));
				}

			}
			SECTION("Translation and Rotation") {
				translatePolygon3D(poly, translation);
				I += area * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
				rotatePolygon3D(poly, Q);
				I = Q * I * Q.transpose();
				const Eigen::Vector3d CG = Q * translation;

				SECTION("Diameter") {
					REQUIRE(Polygon3D::Diameter(poly) == Catch::Approx(regPolyDiameter(nv)));
				}
				SECTION("Monomial Integrals") {
					const auto monInts = Polygon3D::MonomialIntegrals(poly, 2);
					REQUIRE(monInts[0] == Catch::Approx(area));
					REQUIRE(monInts[1] == Catch::Approx(CG(0) * area));
					REQUIRE(monInts[2] == Catch::Approx(CG(1) * area));
					REQUIRE(monInts[3] == Catch::Approx(CG(2) * area));
					REQUIRE(I(0, 1) == Catch::Approx(-monInts[5]));
					REQUIRE(I(0, 2) == Catch::Approx(-monInts[6]));
					REQUIRE(I(1, 2) == Catch::Approx(-monInts[8]));
					REQUIRE(I(0, 0) == Catch::Approx(monInts[7] + monInts[9]));
					REQUIRE(I(1, 1) == Catch::Approx(monInts[4] + monInts[9]));
					REQUIRE(I(2, 2) == Catch::Approx(monInts[4] + monInts[7]));
				}
			}

		}
		ss.str("");
	}
	
	
}

int main(int argc, char* argv[]) {
	int result = Catch::Session().run(argc, argv);
	return result;
}