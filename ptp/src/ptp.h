/**********************************************************************************************************************************
                                            Polytopes (PTP) Library by Tiago Fernandes Moherdaui

            This library provides polygon and polyhedra-related services aimed at implementations of the virtual element method.
			It relies on the monomial indexing proposed and implemented in the sibling library MNL (https://github.com/tiagomhrd/mnl)
            
			In summary:
			* Polygons are represented by a (couter-clockwise ordered) vector of points (Eigen::Vector2d).
			* Polyhedra are represented by an unordered list of points (Eigen::Vector3d) and 
				a vector of face connectivities (vector of ccw ordered indices for the face).

			The main services provided are:
			* Monomial Integrals: a vector of doubles with the integrals of monomials following the indexing in MNL
			* Diameter: largest distance between any two vertices
			* Unique Hyperplanes: vector of representation of unique support hyperplanes (lines in 2D, planes in 3D)
				that define the boundary of the polytope.
				This is usefull for implementation of Serendipity Virtual Element Methods.

			These operations are organized inside namespace "ptp", and then inside the namespace of the geometric entity:
			* ptp::Polygon2D
			* ptp::Polyhedron
			* ptp::Polygon3D (auxiliar for Polyhedron)
			* ptp::Transform (auxiliar for the Test framework)
			
			The code is produced in C++17.
			It includes only <Eigen/Core> and <vector> in the header.
			And includes <algorithm>, <Eigen/Dense>, <mnl.hpp> and <glq.hpp> in the source file.

**********************************************************************************************************************************/

#pragma once
#ifndef PTP
#define PTP
#include <Eigen/Eigen/Core>
#include <vector>

namespace ptp {
	namespace Polygon2D {
		const std::vector<double> 			MonomialIntegrals(const std::vector<Eigen::Vector2d>& vertices, const int maxOrder);
		const double 						Diameter(const std::vector<Eigen::Vector2d>& vertices);
		const std::vector<Eigen::Vector3d> 	UniqueSides(const std::vector<Eigen::Vector2d>& vertices);
		const std::vector<Eigen::Vector3d> 	UniqueReentrantSides(const std::vector<Eigen::Vector2d>& vertices);
	}

	namespace Polygon3D {
		const std::vector<Eigen::Vector3d> 	GetVertices(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices);
		const std::vector<double> 			MonomialIntegrals(const std::vector<Eigen::Vector3d>& vertices, int maxOrder);
		const double 						Diameter(const std::vector<Eigen::Vector3d>& vertices);
		const Eigen::Vector3d 				Normal(const std::vector<Eigen::Vector3d>& vertices);
		const Eigen::Vector4d 				PlaneEquation(const std::vector<Eigen::Vector3d>& vertices);
	}

	namespace Polyhedron {
		const std::vector<double> 			MonomialIntegrals(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces, int maxOrder);
		const double 						Diameter(const std::vector<Eigen::Vector3d>& vertices);
		const std::vector<Eigen::Vector4d> 	UniquePlanes(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces);
	}

	namespace Transform {
		void 								TranslateVertices(std::vector<Eigen::Vector3d>& vertices, const Eigen::Vector3d& translation);
		void 								RotateVertices(std::vector<Eigen::Vector3d>& vertices, const Eigen::Matrix3d& Q);
	}
}
#endif