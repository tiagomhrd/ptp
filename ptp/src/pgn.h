#pragma once
#ifndef PTP_PGN
#define PTP_PGN

#include "Eigen/Eigen/Core"

class Polygon2D {
public:
	Polygon2D(const std::vector<Eigen::Vector2d>& vertices)
		: m_Vertices{ vertices }{}

	const std::vector<double> MonomialIntegrals(const int maxOrder) const { return MonomialIntegrals(m_Vertices, maxOrder); }
	static const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector2d>& vertices, const int maxOrder);
	static const double Diameter(const std::vector<Eigen::Vector2d>& vertices);

private:
	std::vector<Eigen::Vector2d> m_Vertices;
};

#endif