#pragma once
#ifndef PTP_PGN
#define PTP_PGN
#include "Eigen/Eigen/Core"

class Polygon2D {
public:
	Polygon2D(const std::vector<Eigen::Vector2d>& vertices)
		: m_Vertices{ vertices }{}

	static const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector2d>& vertices, const int maxOrder);
	static const double Diameter(const std::vector<Eigen::Vector2d>& vertices);
	static const std::vector<Eigen::Vector3d> UniqueSides(const std::vector<Eigen::Vector2d>& vertices);
	static const std::vector<Eigen::Vector3d> UniqueReentrantSides(const std::vector<Eigen::Vector2d>& vertices);
	
	const std::vector<double> MonomialIntegrals(const int maxOrder) const { return MonomialIntegrals(m_Vertices, maxOrder); }
	const std::vector<Eigen::Vector3d> UniqueSides() const { return UniqueSides(m_Vertices); }
	const std::vector<Eigen::Vector3d> UniqueReentrantSides() const { return UniqueReentrantSides(m_Vertices); }
	

protected:
	std::vector<Eigen::Vector2d> m_Vertices;
};

class Polygon3D {
public:
	Polygon3D(const std::vector<Eigen::Vector3d>& vertices)
		: m_Vertices{vertices}{}
	Polygon3D(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices) { m_Vertices = GetVertices(polyhedronVertices, faceIndices); }

	// Static interface
	
	static const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices, int maxOrder) { return MonomialIntegrals(GetVertices(polyhedronVertices, faceIndices), maxOrder); }
	static const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector3d>& vertices, int maxOrder);
	
	static const double Diameter(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices); // I think having this one specialized spares time
	static const double Diameter(const std::vector<Eigen::Vector3d>& vertices);
	
	static const Eigen::Vector3d Normal(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices) { return Normal(GetVertices(polyhedronVertices, faceIndices)); }
	static const Eigen::Vector3d Normal(const std::vector<Eigen::Vector3d>& vertices);

	static const Eigen::Vector4d PlaneEquation(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices) { return PlaneEquation(GetVertices(polyhedronVertices, faceIndices)); }
	static const Eigen::Vector4d PlaneEquation(const std::vector<Eigen::Vector3d>& vertices);

	const std::vector<double> MonomialIntegrals(int maxOrder) const { return MonomialIntegrals(m_Vertices, maxOrder); }
	const double Diameter() const { return Diameter(m_Vertices); }
	const Eigen::Vector3d Normal() const { return Normal(m_Vertices); }
	const Eigen::Vector4d PlaneEquation() const { return PlaneEquation(m_Vertices); }


protected:
	static const std::vector<Eigen::Vector3d> GetVertices(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices);
	std::vector<Eigen::Vector3d> m_Vertices;
};

class Polyhedron {
public:
	Polyhedron(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces)
		: m_Vertices{vertices}, m_Faces{ faces } {}

	static const std::vector<double> MonomialIntegrals(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces, int maxOrder);
	static const double Diameter(const std::vector<Eigen::Vector3d>& vertices);
	static const std::vector<Eigen::Vector4d> UniquePlanes(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces);


	const std::vector<double> MonomialIntegrals(int maxOrder) const { return MonomialIntegrals(m_Vertices, m_Faces, maxOrder); }
	const double Diameter() const { return Diameter(m_Vertices); }
	const std::vector<Eigen::Vector4d> UniquePlanes() const { return UniquePlanes(m_Vertices, m_Faces); }

	void TranslateVertices(const Eigen::Vector3d& translation);
	void RotateVertices(const Eigen::Matrix3d& Q);

protected:
	std::vector<Eigen::Vector3d> m_Vertices;
	std::vector<std::vector<size_t>> m_Faces;
};
#endif