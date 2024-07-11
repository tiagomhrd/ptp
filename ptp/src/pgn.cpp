#include "pgn.h"

#include "Eigen/Eigen/Dense"

#include "mnl/include/mnl.hpp"
#include "mnl/include/glq.hpp"

#ifndef M_PI2
#define M_PI2 1.570796326794897
#endif // !M_PI2

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

namespace PTP {

    const std::vector<double> Polygon2D::MonomialIntegrals(const std::vector<Eigen::Vector2d>& vertices, const int maxOrder)
    {
        std::vector<double> integrals(mnl::PSpace2D::SpaceDim(maxOrder));
        const Eigen::Vector2d& x0 = vertices[0];
        const size_t nv = vertices.size();

        // Create and initialize edge constants.
        const Eigen::VectorXd edgeConstants = [&x0, &vertices, &nv]() {
            Eigen::VectorXd res = Eigen::VectorXd::Zero(nv - 2);
            Eigen::Matrix2d Q;
            Q << 0., 1., -1., 0.;
            std::transform(vertices.cbegin() + 1, vertices.cend() - 1,
                vertices.cbegin() + 2, res.begin(),
                [&Q, &x0](const auto& xs, const auto& xe) {
                    return (xs - x0).dot(Q * (xe - xs));
                });
            return res;
            }();


            for (int n = 1, ngMax = ceil((maxOrder + 1.) * .5); n <= ngMax; ++n) {
                // Initialize quadrature structures
                const auto quadrature = mnl::GaussLegendreR(2 * n - 1);
                const size_t ng = quadrature.size();
                Eigen::VectorXd w = Eigen::VectorXd::Zero(ng);
                std::vector<std::vector<Eigen::Vector2d>> edgeGaussPos(nv - 2);
                for (auto& vec : edgeGaussPos) {
                    vec.reserve(ng);
                }
                for (size_t g{}; const auto & [xig, wg] : quadrature) {
                    w(g) = wg;
                    for (size_t e = 1; e <= nv - 2; ++e) {
                        edgeGaussPos[e - 1].push_back(vertices[e] * (1 - xig) + vertices[e + 1] * xig);
                    }
                    ++g;
                }

                // The order has to be iterated because of the homegeneous numerical integration
                for (int k = 2 * (n - 1), maxk = std::min(maxOrder, 2 * n - 1); k <= maxk; ++k) {
                    const double orderConst = pow((double)(k + 2), -1.);
                    for (int alpha = mnl::PSpace2D::SpaceDim(k - 1), maxAlpha = mnl::PSpace2D::SpaceDim(k); alpha < maxAlpha; ++alpha) {
                        const int ex = mnl::PSpace2D::Exponent(alpha, 0),
                            ey = mnl::PSpace2D::Exponent(alpha, 1),
                            dx = mnl::PSpace2D::D(alpha, 0),
                            dy = mnl::PSpace2D::D(alpha, 1);
                        Eigen::VectorXd partialIntegral = Eigen::VectorXd::Zero(nv - 2);
                        std::transform(edgeGaussPos.cbegin(), edgeGaussPos.cend(), partialIntegral.begin(),
                            [ng, ex, ey, &w](const auto& edgePosVector) {
                                Eigen::VectorXd values = Eigen::VectorXd::Zero(ng);
                                std::transform(edgePosVector.cbegin(), edgePosVector.cend(), values.begin(),
                                    [ex, ey](const auto& pos) { return _pow(pos[0], ex) * _pow(pos[1], ey); });
                                return values.dot(w);
                            });
                        integrals[alpha] = partialIntegral.dot(edgeConstants) + (dx >= 0 ? x0(0) * ex * integrals[dx] : 0.0) + (dy >= 0 ? x0(1) * ey * integrals[dy] : 0.0);
                        integrals[alpha] *= orderConst;
                    }
                }
            }

            return integrals;
    }

    const double Polygon2D::Diameter(const std::vector<Eigen::Vector2d>& vertices)
    {
        size_t nv = vertices.size();
        double diameter = 0.0;
        for (size_t i = 0; i < nv - 1; ++i)
            for (size_t j = i + 1; j < nv; ++j)
                diameter = std::max(diameter, (vertices[i] - vertices[j]).norm());
        return diameter;
    }

    static Eigen::Vector3d LineEquation(const Eigen::Vector2d& tau, const Eigen::Vector2d& point) {
        return Eigen::Vector3d(tau(0) * point(1) - tau(1) * point(0), -tau(1), tau(0)).normalized();
    };
    static void RemoveColinear(std::vector<Eigen::Vector3d>& lines) {
        if (lines.size() < 4)
            return;

        // Second step - remove repeated
        for (auto it = lines.begin(), begin = it, end = lines.end(); it != end; ++it) {
            for (auto it2 = it + 1; it2 != end; ++it2) {
                if ((*it).dot(*it2) >= 1. - 1e-9) {
                    lines.erase(it2);
                    end = lines.end();
                    break;
                }
            }
        }
    };
    const std::vector<Eigen::Vector3d> Polygon2D::UniqueSides(const std::vector<Eigen::Vector2d>& vertices)
    {
        std::vector<Eigen::Vector3d> out;
        const size_t nv = vertices.size();
        out.reserve(nv);

        for (size_t v = 0; v < nv; ++v) {
            out.push_back(LineEquation((vertices[(v + 1) % nv] - vertices[v]), vertices[v]));
        }

        if (nv == 3)
            return out;

        RemoveColinear(out);
        return out;
    }

    const std::vector<Eigen::Vector3d> Polygon2D::UniqueReentrantSides(const std::vector<Eigen::Vector2d>& vertices)
    {
        std::vector<Eigen::Vector3d> out;

        const size_t nv = vertices.size();
        out.reserve(nv);

        // First step - sliding window of 3 vertices
        for (size_t v = nv; v < 2 * nv; ++v) {
            const auto prev = vertices[(v - 1) % nv],
                curr = vertices[v % nv],
                next = vertices[(v + 1) % nv];
            const auto ni = Eigen::Vector2d(curr(1) - prev(1), prev(0) - curr(0));
            const Eigen::Vector2d taui1 = next - curr;
            if (ni.dot(taui1) > .0) {
                out.push_back(LineEquation(curr - prev, curr));
                out.push_back(LineEquation(taui1, curr));
            }
        }

        // Second step - remove repeated
        RemoveColinear(out);

        return out;
    }

    static const Eigen::Matrix3d RodriguesTensor(const Eigen::Vector3d& axis, const double angleInRadians) {
        Eigen::Matrix3d skewTheta;
        skewTheta << 0, -axis[2], axis[1], axis(2), 0, -axis(0), -axis(1), axis(0), 0.0;
        skewTheta *= angleInRadians;
        const double h1 = sin(angleInRadians) / angleInRadians;
        const double h2 = 0.5 * _pow(sin(0.5 * angleInRadians) / (.5 * angleInRadians), 2);
        return Eigen::Matrix3d::Identity() + h1 * skewTheta + h2 * skewTheta * skewTheta;
    }

    const std::vector<double> Polygon3D::MonomialIntegrals(const std::vector<Eigen::Vector3d>& vertices, int maxOrder)
    {
        std::vector<double> integrals(mnl::PSpace3D::SpaceDim(maxOrder));

        const Eigen::Vector3d& x0 = vertices[0];
        const size_t nv = vertices.size();
        const Eigen::Vector3d faceNormal = Normal(vertices);
        // Create and initialize edge constants.
        const Eigen::VectorXd edgeConstants = [&faceNormal, &x0, &vertices, &nv]() {
            Eigen::VectorXd res = Eigen::VectorXd::Zero(nv - 2);
            const Eigen::Matrix3d Q = RodriguesTensor(faceNormal, -M_PI2);
            std::transform(vertices.cbegin() + 1, vertices.cend() - 1,
                vertices.cbegin() + 2, res.begin(),
                [&Q, &x0](const auto& xs, const auto& xe) {
                    return (xs - x0).dot(Q * (xe - xs));
                });
            return res;
            }();


            for (int n = 1, ngMax = ceil((maxOrder + 1.) * .5); n <= ngMax; ++n) {
                // Initialize quadrature structures
                const auto quadrature = mnl::GaussLegendreR(2 * n - 1);
                const size_t ng = quadrature.size();
                Eigen::VectorXd w = Eigen::VectorXd::Zero(ng);
                std::vector<std::vector<Eigen::Vector3d>> edgeGaussPos(nv - 2);
                for (auto& vec : edgeGaussPos) {
                    vec.reserve(ng);
                }
                for (size_t g{}; const auto & [xig, wg] : quadrature) {
                    w(g) = wg;
                    for (size_t e = 1; e <= nv - 2; ++e) {
                        edgeGaussPos[e - 1].push_back(vertices[e] * (1 - xig) + vertices[e + 1] * xig);
                    }
                    ++g;
                }

                // The order has to be iterated because of the homegeneous numerical integration
                for (int k = 2 * (n - 1), maxk = std::min(maxOrder, 2 * n - 1); k <= maxk; ++k) {
                    const double orderConst = pow((double)(k + 2), -1.);
                    for (int alpha = mnl::PSpace3D::SpaceDim(k - 1), maxAlpha = mnl::PSpace3D::SpaceDim(k); alpha < maxAlpha; ++alpha) {
                        const int
                            ex = mnl::PSpace3D::Exponent(alpha, 0),
                            ey = mnl::PSpace3D::Exponent(alpha, 1),
                            ez = mnl::PSpace3D::Exponent(alpha, 2),
                            dx = mnl::PSpace3D::D(alpha, 0),
                            dy = mnl::PSpace3D::D(alpha, 1),
                            dz = mnl::PSpace3D::D(alpha, 2);
                        Eigen::VectorXd partialIntegral = Eigen::VectorXd::Zero(nv - 2);
                        std::transform(edgeGaussPos.cbegin(), edgeGaussPos.cend(), partialIntegral.begin(),
                            [ng, ex, ey, ez, &w](const auto& edgePosVector) {
                                Eigen::VectorXd values = Eigen::VectorXd::Zero(ng);
                                std::transform(edgePosVector.cbegin(), edgePosVector.cend(), values.begin(),
                                    [ex, ey, ez](const auto& pos) { return _pow(pos[0], ex) * _pow(pos[1], ey) * _pow(pos[2], ez); });
                                return values.dot(w);
                            });
                        integrals[alpha] = partialIntegral.dot(edgeConstants) + (dx >= 0 ? x0(0) * ex * integrals[dx] : 0.0) + (dy >= 0 ? x0(1) * ey * integrals[dy] : 0.0) + (dz >= 0 ? x0(2) * ez * integrals[dz] : 0.0);
                        integrals[alpha] *= orderConst;
                    }
                }
            }

            return integrals;
    }

    const double Polygon3D::Diameter(const std::vector<Eigen::Vector3d>& vertices)
    {
        size_t nv = vertices.size();
        double diameter = 0.0;
        for (size_t i = 0; i < nv - 1; ++i)
            for (size_t j = i + 1; j < nv; ++j)
                diameter = std::max(diameter, (vertices[i] - vertices[j]).norm());
        return diameter;
    }

    const Eigen::Vector3d Polygon3D::Normal(const std::vector<Eigen::Vector3d>& vertices)
    {
        return (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]).normalized();
    }

    const Eigen::Vector4d Polygon3D::PlaneEquation(const std::vector<Eigen::Vector3d>& vertices)
    {
        const Eigen::Vector3d normal = -Normal(vertices); // Minus sign to ensure this function takes positive values inside the polyhedron
        return Eigen::Vector4d(normal.dot(vertices[0]), normal(0), normal(1), normal(2)).normalized();
    }

    const std::vector<Eigen::Vector3d> Polygon3D::GetVertices(const std::vector<Eigen::Vector3d>& polyhedronVertices, const std::vector<size_t>& faceIndices)
    {
        std::vector<Eigen::Vector3d> vertices;
        vertices.reserve(faceIndices.size());
        std::for_each(faceIndices.cbegin(), faceIndices.cend(), [&vertices, &polyhedronVertices](const size_t index) { vertices.emplace_back(polyhedronVertices[index]); });
        return vertices;
    }

    const std::vector<double> Polyhedron::MonomialIntegrals(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces, int maxOrder)
    {
        std::vector<double> integrals(mnl::PSpace3D::SpaceDim(maxOrder));
        const Eigen::Vector3d& x0 = vertices[0];

        // Preprocess faces
        std::vector<std::vector<size_t>> effectiveFaces;
        std::vector<std::vector<double>> faceMonIntegrals;
        std::vector<Eigen::Vector3d> faceNormals;
        effectiveFaces.reserve(faces.size());
        faceMonIntegrals.reserve(faces.size());
        faceNormals.reserve(faces.size());
        for (const auto& face : faces) {
            if (std::find(face.cbegin(), face.cend(), 0) != face.cend())
                continue;
            effectiveFaces.emplace_back(face);
            faceNormals.emplace_back(Polygon3D::Normal(Polygon3D::GetVertices(vertices, face)));
            faceMonIntegrals.emplace_back(Polygon3D::MonomialIntegrals(Polygon3D::GetVertices(vertices, face), maxOrder));
        }
        const size_t nFaces = effectiveFaces.size();

        for (int k = 0; k <= maxOrder; ++k) {
            const double orderConst = pow((double)(k + 3), -1.);
            for (int alpha = mnl::PSpace3D::SpaceDim(k - 1), maxAlpha = mnl::PSpace3D::SpaceDim(k); alpha < maxAlpha; ++alpha) {
                const int ex = mnl::PSpace3D::Exponent(alpha, 0), ey = mnl::PSpace3D::Exponent(alpha, 1), ez = mnl::PSpace3D::Exponent(alpha, 2);
                const int dx = mnl::PSpace3D::D(alpha, 0), dy = mnl::PSpace3D::D(alpha, 1), dz = mnl::PSpace3D::D(alpha, 2);
                const Eigen::Vector3d intGradAlpha(
                    (dx == -1 ? 0.0 : ex * integrals[dx]),
                    (dy == -1 ? 0.0 : ey * integrals[dy]),
                    (dz == -1 ? 0.0 : ez * integrals[dz])
                );

                for (size_t f{}; f < nFaces; ++f)
                    integrals[alpha] += (vertices[effectiveFaces[f][0]] - x0).dot(faceNormals[f]) * faceMonIntegrals[f][alpha];

                integrals[alpha] += intGradAlpha.dot(x0);
                integrals[alpha] *= orderConst;
            }
        }

        return integrals;
    }

    const double Polyhedron::Diameter(const std::vector<Eigen::Vector3d>& vertices)
    {
        size_t nv = vertices.size();
        double diameter = 0.0;
        for (size_t i = 0; i < nv - 1; ++i)
            for (size_t j = i + 1; j < nv; ++j)
                diameter = std::max(diameter, (vertices[i] - vertices[j]).norm());
        return diameter;
    }

    static void RemoveCoplanar(std::vector<Eigen::Vector4d>& planes) {
        if (planes.size() < 5)
            return;

        // Second step - remove repeated
        for (auto it = planes.begin(), begin = it, end = planes.end(); it != end; ++it) {
            for (auto it2 = it + 1; it2 != end; ++it2) {
                if ((*it).dot(*it2) >= 1. - 1e-9) {
                    planes.erase(it2);
                    end = planes.end();
                    break;
                }
            }
        }
    }
    const std::vector<Eigen::Vector4d> Polyhedron::UniquePlanes(const std::vector<Eigen::Vector3d> vertices, const std::vector<std::vector<size_t>> faces)
    {
        std::vector<Eigen::Vector4d> out;
        out.reserve(faces.size());
        for (const auto& face : faces)
            out.push_back(Polygon3D::PlaneEquation(Polygon3D::GetVertices(vertices, face)));

        RemoveCoplanar(out);

        return out;
    }

    void Transform::TranslateVertices(std::vector<Eigen::Vector3d>& vertices, const Eigen::Vector3d& translation)
    {
        std::transform(vertices.cbegin(), vertices.cend(), vertices.begin(), [&translation](const auto& pos) { return pos + translation; });
    }
    void Transform::RotateVertices(std::vector<Eigen::Vector3d>& vertices, const Eigen::Matrix3d& Q)
    {
        std::transform(vertices.cbegin(), vertices.cend(), vertices.begin(), [&Q](const auto& pos) { return Q * pos; });
    }
};