#define _USE_MATH_DEFINES
#include <cmath>
#include "catch_amalgamated.hpp"
#include "src/pgn.h"
#include <iostream>

const double tol = std::numeric_limits<double>::epsilon() * 100.;

static const std::vector<Eigen::Vector2d> regularPolygon(const int nv) {
	std::vector<Eigen::Vector2d> polygon;
	polygon.reserve(nv);
	const double halfTheta = M_PI / nv;
	const double r = pow(2. * sin(halfTheta), -1.);
	for (int i{}; i < nv; ++i)
		polygon.emplace_back(r * cos(2 * i * halfTheta), r * sin(2 * i * halfTheta));
	return polygon;
}

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

static const Eigen::Matrix3d RodriguesTensor(const Eigen::Vector3d& axis, const double angleInRadians) {
	Eigen::Matrix3d skewTheta;
	skewTheta << 0, -axis[2], axis[1], axis(2), 0, -axis(0), -axis(1), axis(0), 0.0;
	skewTheta *= angleInRadians;
	const double h1 = sin(angleInRadians) / angleInRadians;
	const double h2 = 0.5 * _pow(sin(0.5 * angleInRadians) / (.5 * angleInRadians), 2);
	return Eigen::Matrix3d::Identity() + h1 * skewTheta + h2 * skewTheta * skewTheta;
}

static const std::vector<Eigen::Vector3d> regularPolygon3d(const int nv, const Eigen::Matrix3d& rotMatrix, const Eigen::Vector3d& translation) {
	std::vector<Eigen::Vector3d> polygon;
	polygon.reserve(nv);
	const double halfTheta = M_PI / nv;
	const double r = pow(2. * sin(halfTheta), -1.);
	for (int i{}; i < nv; ++i) {
		Eigen::Vector3d pos(r * cos(2 * i * halfTheta), r * sin(2 * i * halfTheta), 0.0);
		polygon.emplace_back(rotMatrix * (pos + translation));
	}
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

TEST_CASE("Polygon3D") {

	SECTION("No Transformation") {
		const Eigen::Matrix3d Q = Eigen::Matrix3d::Identity();
		const Eigen::Vector3d translation = Eigen::Vector3d::Zero();

		std::array<const std::vector<Eigen::Vector3d>, 4> regPolys = {
		regularPolygon3d(3, Q, translation),
		regularPolygon3d(4, Q, translation),
		regularPolygon3d(6, Q, translation),
		regularPolygon3d(8, Q, translation)
		};
		const auto& triangle = regPolys[0];
		const auto& square = regPolys[1];
		const auto& hexagon = regPolys[2];
		const auto& octagon = regPolys[3];

		SECTION("Triangle") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(triangle) == Catch::Approx(1.));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(triangle, 2);
				const double area = sqrt(3.) / 4.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = pow(32. * sqrt(3.), -1.);
				I(1, 1) = pow(32. * sqrt(3.), -1.);
				I(2, 2) = 2 * pow(32. * sqrt(3.), -1.);

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
		SECTION("Square") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(square) == Catch::Approx(sqrt(2.)));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(square, 2);
				const double area = 1.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = 1. / 12.;
				I(1, 1) = 1. / 12.;
				I(2, 2) = 2. / 12.;

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
		SECTION("Hexagon") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(hexagon) == Catch::Approx(2.));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(hexagon, 2);
				const double area = 3. * sqrt(3) / 2.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = 5. * sqrt(3.) / 16.;
				I(1, 1) = 5. * sqrt(3.) / 16.;
				I(2, 2) = 5. * sqrt(3.) / 8.;

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
		SECTION("Octagon") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(octagon) == Catch::Approx(sqrt(4. + 2. * sqrt(2.))));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(octagon, 2);
				const double area = 2 * (sqrt(2.) + 1.);

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = (11. + 8. * sqrt(2.)) / 12.;
				I(1, 1) = (11. + 8. * sqrt(2.)) / 12.;
				I(2, 2) = (11. + 8. * sqrt(2.)) / 6.;

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
	}
	SECTION("Translation Only") {
		const Eigen::Matrix3d Q = Eigen::Matrix3d::Identity();
		const Eigen::Vector3d translation(1., 2., 3.);
		const Eigen::Vector3d CG = Q * translation;

		std::array<const std::vector<Eigen::Vector3d>, 4> regPolys = {
		regularPolygon3d(3, Q, translation),
		regularPolygon3d(4, Q, translation),
		regularPolygon3d(6, Q, translation),
		regularPolygon3d(8, Q, translation)
		};
		const auto& triangle = regPolys[0];
		const auto& square = regPolys[1];
		const auto& hexagon = regPolys[2];
		const auto& octagon = regPolys[3];

		const Eigen::Matrix3d T = translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose();

		SECTION("Triangle") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(triangle) == Catch::Approx(1.));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(triangle, 2);
				const double area = sqrt(3.) / 4.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = pow(32. * sqrt(3.), -1.);
				I(1, 1) = pow(32. * sqrt(3.), -1.);
				I(2, 2) = 2 * pow(32. * sqrt(3.), -1.);

				// Translation
				I += area * T;

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
		SECTION("Square") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(square) == Catch::Approx(sqrt(2.)));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(square, 2);
				const double area = 1.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = 1. / 12.;
				I(1, 1) = 1. / 12.;
				I(2, 2) = 2. / 12.;

				// Translation
				I += area * T;

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
		SECTION("Hexagon") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(hexagon) == Catch::Approx(2.));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(hexagon, 2);
				const double area = 3. * sqrt(3) / 2.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = 5. * sqrt(3.) / 16.;
				I(1, 1) = 5. * sqrt(3.) / 16.;
				I(2, 2) = 5. * sqrt(3.) / 8.;

				// Translation
				I += area * T;

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
		SECTION("Octagon") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(octagon) == Catch::Approx(sqrt(4. + 2. * sqrt(2.))));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(octagon, 2);
				const double area = 2 * (sqrt(2.) + 1.);

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = (11. + 8. * sqrt(2.)) / 12.;
				I(1, 1) = (11. + 8. * sqrt(2.)) / 12.;
				I(2, 2) = (11. + 8. * sqrt(2.)) / 6.;

				// Translation
				I += area * T;

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
	SECTION("Rotation Only") {
		const Eigen::Matrix3d Q = RodriguesTensor(Eigen::Vector3d(sqrt(2.) * .5, sqrt(2.) * .5, 0.), M_PI_4);
		const Eigen::Vector3d translation(0., 0., 0.);

		std::array<const std::vector<Eigen::Vector3d>, 4> regPolys = {
		regularPolygon3d(3, Q, translation),
		regularPolygon3d(4, Q, translation),
		regularPolygon3d(6, Q, translation),
		regularPolygon3d(8, Q, translation)
		};
		const auto& triangle = regPolys[0];
		const auto& square = regPolys[1];
		const auto& hexagon = regPolys[2];
		const auto& octagon = regPolys[3];

		SECTION("Triangle") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(triangle) == Catch::Approx(1.));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(triangle, 2);
				const double area = sqrt(3.) / 4.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = pow(32. * sqrt(3.), -1.);
				I(1, 1) = pow(32. * sqrt(3.), -1.);
				I(2, 2) = 2 * pow(32. * sqrt(3.), -1.);

				// Rotation
				I = Q * I * Q.transpose();

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
		SECTION("Square") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(square) == Catch::Approx(sqrt(2.)));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(square, 2);
				const double area = 1.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = 1. / 12.;
				I(1, 1) = 1. / 12.;
				I(2, 2) = 2. / 12.;

				// Rotation
				I = Q * I * Q.transpose();

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
		SECTION("Hexagon") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(hexagon) == Catch::Approx(2.));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(hexagon, 2);
				const double area = 3. * sqrt(3) / 2.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = 5. * sqrt(3.) / 16.;
				I(1, 1) = 5. * sqrt(3.) / 16.;
				I(2, 2) = 5. * sqrt(3.) / 8.;

				// Rotation
				I = Q * I * Q.transpose();

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
		SECTION("Octagon") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(octagon) == Catch::Approx(sqrt(4. + 2. * sqrt(2.))));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(octagon, 2);
				const double area = 2 * (sqrt(2.) + 1.);

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = (11. + 8. * sqrt(2.)) / 12.;
				I(1, 1) = (11. + 8. * sqrt(2.)) / 12.;
				I(2, 2) = (11. + 8. * sqrt(2.)) / 6.;

				// Rotation
				I = Q * I * Q.transpose();

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
	}
	SECTION("Translation and Rotation") {
		const Eigen::Matrix3d Q = RodriguesTensor(Eigen::Vector3d(sqrt(2.) * .5, sqrt(2.) * .5, 0.), M_PI_4);
		const Eigen::Vector3d translation(1., 2., 3.);
		const Eigen::Vector3d CG = Q * translation;

		std::array<const std::vector<Eigen::Vector3d>, 4> regPolys = {
		regularPolygon3d(3, Q, translation),
		regularPolygon3d(4, Q, translation),
		regularPolygon3d(6, Q, translation),
		regularPolygon3d(8, Q, translation)
		};
		const auto& triangle = regPolys[0];
		const auto& square = regPolys[1];
		const auto& hexagon = regPolys[2];
		const auto& octagon = regPolys[3];

		const Eigen::Matrix3d T = translation.dot(translation) * Eigen::Matrix3d::Identity() - translation * translation.transpose();

		SECTION("Triangle") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(triangle) == Catch::Approx(1.));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(triangle, 2);
				const double area = sqrt(3.) / 4.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = pow(32. * sqrt(3.), -1.);
				I(1, 1) = pow(32. * sqrt(3.), -1.);
				I(2, 2) = 2 * pow(32. * sqrt(3.), -1.);

				// Rotation and Translation
				I = Q * (I + area * T) * Q.transpose();

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
		SECTION("Square") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(square) == Catch::Approx(sqrt(2.)));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(square, 2);
				const double area = 1.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = 1. / 12.;
				I(1, 1) = 1. / 12.;
				I(2, 2) = 2. / 12.;

				// Rotation and Translation
				I = Q * (I + area * T) * Q.transpose();

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
		SECTION("Hexagon") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(hexagon) == Catch::Approx(2.));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(hexagon, 2);
				const double area = 3. * sqrt(3) / 2.;

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = 5. * sqrt(3.) / 16.;
				I(1, 1) = 5. * sqrt(3.) / 16.;
				I(2, 2) = 5. * sqrt(3.) / 8.;

				// Rotation and Translation
				I = Q * (I + area * T) * Q.transpose();

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
		SECTION("Octagon") {
			SECTION("Diameter") {
				REQUIRE(Polygon3D::Diameter(octagon) == Catch::Approx(sqrt(4. + 2. * sqrt(2.))));
			}
			SECTION("Monomials") {
				const auto monInts = Polygon3D::MonomialIntegrals(octagon, 2);
				const double area = 2 * (sqrt(2.) + 1.);

				Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
				I(0, 0) = (11. + 8. * sqrt(2.)) / 12.;
				I(1, 1) = (11. + 8. * sqrt(2.)) / 12.;
				I(2, 2) = (11. + 8. * sqrt(2.)) / 6.;

				// Rotation and Translation
				I = Q * (I + area * T) * Q.transpose();

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
}

int main(int argc, char* argv[]) {
	int result = Catch::Session().run(argc, argv);
	return result;
}