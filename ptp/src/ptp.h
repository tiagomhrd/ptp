#pragma once
#ifndef PTP
#define PTP
#include <Eigen/Eigen/Core>
#include <vector>

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
#endif