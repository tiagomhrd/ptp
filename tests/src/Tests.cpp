#include "catch_amalgamated.hpp"
#include "src/pgn.h"

TEST_CASE("Polygon2D") {
	const std::vector<Eigen::Vector2d> square = {
		Eigen::Vector2d(-.5, -.5),
		Eigen::Vector2d(+.5, -.5),
		Eigen::Vector2d(+.5, +.5),
		Eigen::Vector2d(-.5, +.5)
	};
	SECTION("Monomials") {
		const auto monInts = Polygon2D::MonomialIntegrals(square, 2);
		REQUIRE(monInts[0] == Catch::Approx(1));
		REQUIRE(monInts[1] == Catch::Approx(0));
		REQUIRE(monInts[2] == Catch::Approx(0));
		REQUIRE(monInts[3] == Catch::Approx(1. / 12.));
		REQUIRE(monInts[4] == Catch::Approx(0));
		REQUIRE(monInts[5] == Catch::Approx(1. / 12.));
	}
}

int main(int argc, char* argv[]) {
	int result = Catch::Session().run(argc, argv);
	return result;
}