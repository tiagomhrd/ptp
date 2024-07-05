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

protected:
	std::vector<Eigen::Vector2d> m_Vertices;
};

class Polygon3D {
public:
	Polygon3D(const std::vector<Eigen::Vector3d>& vertices)
		: m_Vertices{vertices}{}
	Polygon3D(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices) { m_Vertices = GetVertices(polyhedronVertices, faceIndices); }

	// Static interface
	static const std::vector<Eigen::Vector3d> GetVertices(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices);
	static const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices, int maxOrder) { return MonomialIntegrals(GetVertices(polyhedronVertices, faceIndices), maxOrder); }
	static const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector3d>& vertices, int maxOrder);
	
	static const double Diameter(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices); // I think having this one specialized spares time
	static const double Diameter(const std::vector<Eigen::Vector3d>& vertices);
	
	static const Eigen::Vector3d Normal(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices) { return Normal(GetVertices(polyhedronVertices, faceIndices)); }
	static const Eigen::Vector3d Normal(const std::vector<Eigen::Vector3d>& vertices);

	const std::vector<double> MonomialIntegrals(int maxOrder) const { return MonomialIntegrals(m_Vertices, maxOrder); }
	const double Diameter() const { return Diameter(m_Vertices); }

protected:
	std::vector<Eigen::Vector3d> m_Vertices;
};

#endif