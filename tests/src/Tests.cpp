#define _USE_MATH_DEFINES
#include <cmath>
#include "catch_amalgamated.hpp"
#include "src/pgn.h"
#include <iostream>

constexpr double tol = 1e-10;

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

TEST_CASE("Polygon2D") {
	std::stringstream ss;
	for (int nv = 3, max = 15; nv <= max; ++nv) {
		ss << nv << "-gon"; 
		SECTION(ss.str()) {
			const auto poly = regularPolygon(nv);
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polygon2D::Diameter(poly), Catch::Matchers::WithinAbs(regPolyDiameter(nv), tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polygon2D::MonomialIntegrals(poly, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(regPolyArea(nv), tol)); 
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(0.0, tol));
				const double sqInts = .5 * regPolyInertia(nv); 
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(sqInts, tol)); 
				REQUIRE_THAT(monInts[4], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[5], Catch::Matchers::WithinAbs(sqInts, tol)); 
			}
		}
		ss.str("");
	}
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
			const double diameter = regPolyDiameter(nv);
			SECTION("No Transformation") {
				SECTION("Diameter") {
					REQUIRE_THAT(PTP::Polygon3D::Diameter(poly), Catch::Matchers::WithinAbs(diameter, tol));
				}
				SECTION("Monomial Integrals") {
					const auto monInts = PTP::Polygon3D::MonomialIntegrals(poly, 2);
					REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(area, tol));
					REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(monInts[5], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(monInts[6], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(monInts[8], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
					REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
					REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
				}
			}
			
			SECTION("Translation") {
				PTP::Transform::TranslateVertices(poly, translation);
				I += area * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());

				SECTION("Diameter") {
					REQUIRE_THAT(PTP::Polygon3D::Diameter(poly), Catch::Matchers::WithinAbs(diameter, tol));
				}				
				SECTION("Monomial Integrals") {
					const auto monInts = PTP::Polygon3D::MonomialIntegrals(poly, 2);
					REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(area, tol));
					REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(translation(0) * area, tol));
					REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(translation(1) * area, tol));
					REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(translation(2) * area, tol));
					REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
					REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
					REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
					REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
					REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
					REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
				}
			}
			SECTION("Rotation") {
				PTP::Transform::RotateVertices(poly, Q);
				I = Q * I * Q.transpose();

				SECTION("Diameter") {
					REQUIRE_THAT(PTP::Polygon3D::Diameter(poly), Catch::Matchers::WithinAbs(diameter, tol));
				}
				SECTION("Monomial Integrals") {
					const auto monInts = PTP::Polygon3D::MonomialIntegrals(poly, 2);
					REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(area, tol));
					REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(0.0, tol));
					REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
					REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
					REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
					REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
					REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
					REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
				}

			}
			SECTION("Translation and Rotation") {
				PTP::Transform::TranslateVertices(poly, translation);
				I += area * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
				PTP::Transform::RotateVertices(poly, Q);
				I = Q * I * Q.transpose();
				const Eigen::Vector3d CG = Q * translation;

				SECTION("Diameter") {
					REQUIRE_THAT(PTP::Polygon3D::Diameter(poly), Catch::Matchers::WithinAbs(diameter, tol));
				}
				SECTION("Monomial Integrals") {
					const auto monInts = PTP::Polygon3D::MonomialIntegrals(poly, 2);
					REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(area, tol));
					REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(CG(0) * area, tol));
					REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(CG(1) * area, tol));
					REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(CG(2) * area, tol));
					REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
					REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
					REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
					REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
					REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
					REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
				}
			}
		}
		ss.str("");
	}
}

struct Polyhedron {
	std::vector<Eigen::Vector3d> vertices;
	std::vector<std::vector<size_t>> faces;
};

static Polyhedron cuboid(const double dx, const double dy, const double dz) {
	const double hdx = .5 * dx, hdy = .5 * dy, hdz = .5 * dz;
	Polyhedron poly;
	
	poly.vertices.reserve(8);
	poly.vertices.emplace_back(-hdx, -hdy, -hdz);
	poly.vertices.emplace_back(+hdx, -hdy, -hdz);
	poly.vertices.emplace_back(+hdx, +hdy, -hdz);
	poly.vertices.emplace_back(-hdx, +hdy, -hdz);
	poly.vertices.emplace_back(-hdx, -hdy, +hdz);
	poly.vertices.emplace_back(+hdx, -hdy, +hdz);
	poly.vertices.emplace_back(+hdx, +hdy, +hdz);
	poly.vertices.emplace_back(-hdx, +hdy, +hdz);

	poly.faces.reserve(6);
	poly.faces.push_back({0, 4, 7, 3}); // -x
	poly.faces.push_back({1, 2, 6, 5}); // +x
	poly.faces.push_back({0, 1, 5, 4}); // -y
	poly.faces.push_back({2, 3, 7, 6}); // +y
	poly.faces.push_back({0, 3, 2, 1}); // -z
	poly.faces.push_back({4, 5, 6, 7}); // +z

	return poly;
};
static Eigen::Matrix3d cuboidInertiaTensor(const double dx, const double dy, const double dz) {
	Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
	const double vol_12 = dx * dy * dz / 12.;
	I(0, 0) = vol_12 * (_pow(dy,2) + _pow(dz,2));
	I(1, 1) = vol_12 * (_pow(dx,2) + _pow(dz,2));
	I(2, 2) = vol_12 * (_pow(dx,2) + _pow(dy,2));
	return I;
}

static Polyhedron icosahedron() {
	Polyhedron poly;

	poly.vertices.reserve(12);
	poly.vertices.emplace_back(0.288675134594813, 0.5, -0.75576131407617);
	poly.vertices.emplace_back(0.288675134594813, -0.5, -0.75576131407617);
	poly.vertices.emplace_back(0.934172358962715, 0, -0.178411044886545);
	poly.vertices.emplace_back(-0.577350269189625, 0, -0.75576131407617);
	poly.vertices.emplace_back(0.467086179481358, 0.809016994374945, 0.178411044886545);
	poly.vertices.emplace_back(-0.467086179481358, 0.809016994374945, -0.178411044886545);
	poly.vertices.emplace_back(0.467086179481358, -0.809016994374945, 0.178411044886545);
	poly.vertices.emplace_back(-0.467086179481358, -0.809016994374945, -0.178411044886545);
	poly.vertices.emplace_back(0.577350269189625, 0, 0.75576131407617);
	poly.vertices.emplace_back(-0.934172358962715, 0, 0.178411044886545);
	poly.vertices.emplace_back(-0.288675134594813, 0.5, 0.75576131407617);
	poly.vertices.emplace_back(-0.288675134594813, -0.5, 0.75576131407617);

	poly.faces.reserve(20);
	poly.faces.push_back({ 0, 5, 4 });
	poly.faces.push_back({ 0, 3, 5 });
	poly.faces.push_back({ 0, 4, 2 });
	poly.faces.push_back({ 0, 2, 1 });
	poly.faces.push_back({ 0, 1, 3 });
	poly.faces.push_back({ 1, 2, 6 });
	poly.faces.push_back({ 1, 7, 3 });
	poly.faces.push_back({ 1, 6, 7 });
	poly.faces.push_back({ 2, 4, 8 });
	poly.faces.push_back({ 2, 8, 6 });
	poly.faces.push_back({ 3, 9, 5 });
	poly.faces.push_back({ 3, 7, 9 });
	poly.faces.push_back({ 4, 5, 10 });
	poly.faces.push_back({ 4, 10, 8 });
	poly.faces.push_back({ 5, 9, 10 });
	poly.faces.push_back({ 6, 8, 11 });
	poly.faces.push_back({ 6, 11, 7 });
	poly.faces.push_back({ 7, 11, 9 });
	poly.faces.push_back({ 9, 11, 10 });
	poly.faces.push_back({ 8, 10, 11 });

	return poly;
}
static double icosahedronVolume() {
	return 5. * (3. + sqrt(5.)) / 12.;
}
static Eigen::Matrix3d icosahedronInertiaTensor() {
	return _pow((1. + sqrt(5.)) * .5, 2) * .1 * icosahedronVolume() * Eigen::Matrix3d::Identity();
}
static double icosahedronDiameter() {
	return sqrt((5. + sqrt(5.)) * .5);
}

static Polyhedron dodecahedron() {
	Polyhedron poly;
	
	poly.vertices.reserve(20);
	poly.vertices.emplace_back(0.5, 0.688190960235587, -1.1135163644116);
	poly.vertices.emplace_back(-0.5, 0.688190960235587, -1.1135163644116);
	poly.vertices.emplace_back(0.809016994374947, 1.1135163644116, -0.262865556059566);
	poly.vertices.emplace_back(0.809016994374947, -0.262865556059566, -1.1135163644116);
	poly.vertices.emplace_back(-0.809016994374947, 1.1135163644116, -0.262865556059566);
	poly.vertices.emplace_back(-0.809016994374947, -0.262865556059566, -1.1135163644116);
	poly.vertices.emplace_back(0, 1.37638192047117, 0.262865556059567);
	poly.vertices.emplace_back(0, -0.850650808352042, -1.1135163644116);
	poly.vertices.emplace_back(1.30901699437494, 0.425325404176019, 0.262865556059566);
	poly.vertices.emplace_back(1.30901699437494, -0.425325404176019, -0.262865556059566);
	poly.vertices.emplace_back(-1.30901699437494, 0.425325404176019, 0.262865556059566);
	poly.vertices.emplace_back(-1.30901699437494, -0.425325404176019, -0.262865556059566);
	poly.vertices.emplace_back(0, 0.850650808352042, 1.1135163644116);
	poly.vertices.emplace_back(0, -1.37638192047117, -0.262865556059567);
	poly.vertices.emplace_back(0.809016994374947, 0.262865556059566, 1.1135163644116);
	poly.vertices.emplace_back(-0.809016994374947, 0.262865556059566, 1.1135163644116);
	poly.vertices.emplace_back(0.809016994374947, -1.1135163644116, 0.262865556059566);
	poly.vertices.emplace_back(-0.809016994374947, -1.1135163644116, 0.262865556059566);
	poly.vertices.emplace_back(0.5, -0.688190960235587, 1.1135163644116);
	poly.vertices.emplace_back(-0.5, -0.688190960235587, 1.1135163644116);
	
	poly.faces.reserve(12);
	poly.faces.push_back({ 0, 3, 7, 5, 1 });
	poly.faces.push_back({ 0, 1, 4, 6, 2 });
	poly.faces.push_back({ 0, 2, 8, 9, 3 });
	poly.faces.push_back({ 1, 5, 11, 10, 4 });
	poly.faces.push_back({ 5, 7, 13, 17, 11 });
	poly.faces.push_back({ 3, 9, 16, 13, 7 });
	poly.faces.push_back({ 2, 6, 12, 14, 8 });
	poly.faces.push_back({ 4, 10, 15, 12, 6 });
	poly.faces.push_back({ 10, 11, 17, 19, 15 });
	poly.faces.push_back({ 8, 14, 18, 16, 9 });
	poly.faces.push_back({ 13, 16, 18, 19, 17 });
	poly.faces.push_back({ 12, 15, 19, 18, 14 });

	return poly;
}
static double dodecahedronVolume() {
	return (15. + 7. * sqrt(5.)) * .25;
}
static double dodecahedronDiameter() {
	return sqrt(3.) * (1. + sqrt(5.)) * .5;
}
static Eigen::Matrix3d dodecahedronInertiaTensor() {
	return (39. * (1. + sqrt(5.)) * .5 + 28.) / 150. * dodecahedronVolume() * Eigen::Matrix3d::Identity();
}

TEST_CASE("Polyhedra") {
	SECTION("Cube") {
		const double dx = 1., dy = 1., dz = 1.;
		Polyhedron cube = cuboid(dx, dy, dz);
		Eigen::Matrix3d I = cuboidInertiaTensor(dx, dy, dz);
		const Eigen::Vector3d translation(1., 2., 3.);
		const Eigen::Matrix3d Q = RodriguesTensor(Eigen::Vector3d(sqrt(2.) * .5, sqrt(2.) * .5, 0.), M_PI / 6);
		const double diameter = sqrt(_pow(dx, 2) + _pow(dy, 2) + _pow(dz, 2));
		const double vol = dx * dy * dz;
		SECTION("No Transformation") {
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
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
			}
		}
		SECTION("Translation") {
			PTP::Transform::TranslateVertices(cube.vertices, translation);
			I += vol * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(translation(0) * vol, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(translation(1) * vol, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(translation(2) * vol, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
		SECTION("Rotation") {
			PTP::Transform::RotateVertices(cube.vertices, Q);
			I = Q * I * Q.transpose();

			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
		SECTION("Translation and Rotation") {
			PTP::Transform::TranslateVertices(cube.vertices, translation);
			I += vol * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
			PTP::Transform::RotateVertices(cube.vertices, Q);
			I = Q * I * Q.transpose();
			const Eigen::Vector3d CG = Q * translation;

			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(CG(0) * vol, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(CG(1) * vol, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(CG(2) * vol, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
	}
	SECTION("Brick") {
		const double dx = 1., dy = 2., dz = 3.;
		Polyhedron cube = cuboid(dx, dy, dz);
		Eigen::Matrix3d I = cuboidInertiaTensor(dx, dy, dz);
		const Eigen::Vector3d translation(4., 5., 6.);
		const Eigen::Matrix3d Q = RodriguesTensor(Eigen::Vector3d(sqrt(2.) * .5, sqrt(2.) * .5, 0.), M_PI / 6);
		const double diameter = sqrt(_pow(dx, 2) + _pow(dy, 2) + _pow(dz, 2));
		const double vol = dx * dy * dz;
		SECTION("No Transformation") {
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
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
			}
		}
		SECTION("Translation") {
			PTP::Transform::TranslateVertices(cube.vertices, translation);
			I += vol * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(translation(0) * vol, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(translation(1) * vol, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(translation(2) * vol, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
		SECTION("Rotation") {
			PTP::Transform::RotateVertices(cube.vertices, Q);
			I = Q * I * Q.transpose();

			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
		SECTION("Translation and Rotation") {
			PTP::Transform::TranslateVertices(cube.vertices, translation);
			I += vol * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
			PTP::Transform::RotateVertices(cube.vertices, Q);
			I = Q * I * Q.transpose();
			const Eigen::Vector3d CG = Q * translation;

			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(cube.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(cube.vertices, cube.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(CG(0) * vol, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(CG(1) * vol, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(CG(2) * vol, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
	}
	SECTION("Icosahedron") {
		Polyhedron ico = icosahedron();
		Eigen::Matrix3d I = icosahedronInertiaTensor();
		const Eigen::Vector3d translation(4., 5., 6.);
		const Eigen::Matrix3d Q = RodriguesTensor(Eigen::Vector3d(sqrt(2.) * .5, sqrt(2.) * .5, 0.), M_PI / 6);
		const double diameter = icosahedronDiameter();
		const double vol = icosahedronVolume();
		SECTION("No Transformation") {
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(ico.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(ico.vertices, ico.faces, 2);
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
			}
		}
		SECTION("Translation") {
			PTP::Transform::TranslateVertices(ico.vertices, translation);
			I += vol * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(ico.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(ico.vertices, ico.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(translation(0) * vol, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(translation(1) * vol, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(translation(2) * vol, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
		SECTION("Rotation") {
			PTP::Transform::RotateVertices(ico.vertices, Q);
			I = Q * I * Q.transpose();

			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(ico.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(ico.vertices, ico.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
		SECTION("Translation and Rotation") {
			PTP::Transform::TranslateVertices(ico.vertices, translation);
			I += vol * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
			PTP::Transform::RotateVertices(ico.vertices, Q);
			I = Q * I * Q.transpose();
			const Eigen::Vector3d CG = Q * translation;

			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(ico.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(ico.vertices, ico.faces, 2);
				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(CG(0) * vol, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(CG(1) * vol, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(CG(2) * vol, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
	}

	SECTION("Dodecahedron") {
		Polyhedron dodeca = dodecahedron();
		Eigen::Matrix3d I = dodecahedronInertiaTensor();
		const Eigen::Vector3d translation(4., 5., 6.);
		const Eigen::Matrix3d Q = RodriguesTensor(Eigen::Vector3d(sqrt(2.) * .5, sqrt(2.) * .5, 0.), M_PI / 6);
		const double diameter = dodecahedronDiameter();
		const double vol = dodecahedronVolume();
		SECTION("No Transformation") {
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(dodeca.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(dodeca.vertices, dodeca.faces, 2);
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
			}
		}
		SECTION("Translation") {
			PTP::Transform::TranslateVertices(dodeca.vertices, translation);
			I += vol * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(dodeca.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(dodeca.vertices, dodeca.faces, 2);

				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(translation(0) * vol, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(translation(1) * vol, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(translation(2) * vol, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
		SECTION("Rotation") {
			PTP::Transform::RotateVertices(dodeca.vertices, Q);
			I = Q * I * Q.transpose();

			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(dodeca.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(dodeca.vertices, dodeca.faces, 2);

				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(0.0, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
		SECTION("Translation and Rotation") {
			PTP::Transform::TranslateVertices(dodeca.vertices, translation);
			I += vol * (translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose());
			PTP::Transform::RotateVertices(dodeca.vertices, Q);
			I = Q * I * Q.transpose();
			const Eigen::Vector3d CG = Q * translation;

			SECTION("Diameter") {
				REQUIRE_THAT(PTP::Polyhedron::Diameter(dodeca.vertices), Catch::Matchers::WithinAbs(diameter, tol));
			}
			SECTION("Monomial Integrals") {
				const auto monInts = PTP::Polyhedron::MonomialIntegrals(dodeca.vertices, dodeca.faces, 2);

				REQUIRE_THAT(monInts[0], Catch::Matchers::WithinAbs(vol, tol));
				REQUIRE_THAT(monInts[1], Catch::Matchers::WithinAbs(CG(0) * vol, tol));
				REQUIRE_THAT(monInts[2], Catch::Matchers::WithinAbs(CG(1) * vol, tol));
				REQUIRE_THAT(monInts[3], Catch::Matchers::WithinAbs(CG(2) * vol, tol));
				REQUIRE_THAT(I(0, 1), Catch::Matchers::WithinAbs(-monInts[5], tol));
				REQUIRE_THAT(I(0, 2), Catch::Matchers::WithinAbs(-monInts[6], tol));
				REQUIRE_THAT(I(1, 2), Catch::Matchers::WithinAbs(-monInts[8], tol));
				REQUIRE_THAT(I(0, 0), Catch::Matchers::WithinAbs(monInts[7] + monInts[9], tol));
				REQUIRE_THAT(I(1, 1), Catch::Matchers::WithinAbs(monInts[4] + monInts[9], tol));
				REQUIRE_THAT(I(2, 2), Catch::Matchers::WithinAbs(monInts[4] + monInts[7], tol));
			}
		}
	}
}

int main(int argc, char* argv[]) {
	int result = Catch::Session().run(argc, argv);
	return result;
}