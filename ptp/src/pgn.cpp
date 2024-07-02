#include "pgn.h"
#include "mnl/include/mnl.hpp"
#include "mnl/include/glq.hpp"


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

const std::vector<double> Polygon2D::MonomialIntegrals(const std::vector<Eigen::Vector2d>& vertices, const int maxOrder)
{
    std::vector<double> integrals(mnl::PSpace2D::SpaceDim(maxOrder));
    const Eigen::Vector2d& x0 = vertices[0];
    const size_t nv = vertices.size();
    Eigen::VectorXd edgeConstants = Eigen::VectorXd::Zero(nv - 2);
    // Initialize edgeConstants
    {
        Eigen::Matrix2d Q;
        Q << 0., 1., -1., 0.;
        std::transform(vertices.cbegin() + 1, vertices.cend() - 1,
            vertices.cbegin() + 2, edgeConstants.begin(),
            [&Q, &x0](const auto& xs, const auto& xe) {
                return (xs - x0).dot(Q * (xe - xs));
            });
    }

    for (int n = 1, ngMax = ceil((maxOrder + 1.) * .5); n <= ngMax; ++n) {
        // Initialize quadrature structures
        const auto quadrature = mnl::GaussLegendreR(2 * n - 1);
        
        const size_t ng = quadrature.size();
        
        Eigen::VectorXd w = Eigen::VectorXd::Zero(ng);
        
        std::vector<std::vector<Eigen::Vector2d>> edgeGaussPos(nv - 2);
        for (auto& vec : edgeGaussPos) {
            vec.reserve(ng);
        }
        for (size_t g{}; const auto& [xig, wg] : quadrature) {
            w(g) = wg;
            for (size_t e = 1; e <= nv - 2; ++e) {
                edgeGaussPos[e - 1].push_back(vertices[e] * (1 - xig) + vertices[e+1] * xig);
            }
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
                integrals[alpha] = partialIntegral.dot(edgeConstants) + (dx >= 0 ? x0(0) * integrals[dx] : 0.0) + (dy >= 0 ? x0(1) * integrals[dy] : 0.0);
                integrals[alpha] *= orderConst;
            }
        }
    }

    return integrals;
}
