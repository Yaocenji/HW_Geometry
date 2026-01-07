#include "MyNurbs.h"
#include <stdexcept>

namespace MyNurbs {

    bool RationalCurve::Interpolate(
        InterpolateMethod method,
        int degree,
        const std::vector<glm::dvec3>& points,
        RationalCurve*& curve
    ) {
        // 传入指针必须为空
        if (curve != nullptr) {
            printf_s("Curve pointer must be null!");
            return false;
        }

        // 插值点必须 大于等于 degree + 1
        if (points.size() < degree + 1) {
            printf_s("Points number must greater than or equals to degree + 1!");
            return false;
        }

        if (method == InterpolateMethod::Uniform) {
            return InterpolateUniform(degree, points, curve);
        }
    }

    bool RationalCurve::InterpolateUniform(
        int degree,
        const std::vector<glm::dvec3>& points,
        RationalCurve*& curve
    ) {
        int n = points.size();

        // 均匀地 构造每个插值点的参数向量
        std::vector<double> parameters;
        parameters.resize(n);
        for (int i = 0; i < n; i++) {
            parameters[i] = (double)i / ((double)1 / (double)(n - 1));
        }

        // 均匀地 构造节点向量
        std::vector<double> knots;
        knots.resize(n + degree + 1);
        for (int i = 0; i < knots.size(); i++) {
            if (i < degree + 1) {
                knots[i] = 0.0;
            }
            else if (i > n - 1) {
                knots[i] = 1.0;
            }
            else {
                knots[i] = (double)(i - degree) / (double)(n - degree);
            }
        }

        return true;
    }

    bool RationalCurve::Create(int deg,
		const std::vector<glm::dvec3>& ctrl_pts,
		const std::vector<double>& kns,
		const std::vector<double>& wts,
		RationalCurve*& curve) {

		// 传入指针必须为空
		if (curve != nullptr) {
			printf_s("Curve pointer must be null!");
			return false;
		}

		// 检查数量是否符合BSpline曲线规范
		if (deg < 1) {
            printf_s("Curve degree must greater than 0!");
			return false;
		}
		if (ctrl_pts.size() < deg + 1) {
            printf_s("Control points number must greater than degree!");
			return false;
		}
		if (kns.size() != ctrl_pts.size() + deg + 1) {
            printf_s("Invalid knot vector size.");
			return false;
		}

		// 创建曲线
		curve = new RationalCurve(deg, ctrl_pts, kns, wts);
		return true;
	}

	RationalCurve::RationalCurve(int deg,
		const std::vector<glm::dvec3>& ctrl_pts,
		const std::vector<double>& kns,
		const std::vector<double>& wts) : degree(deg), control_points(ctrl_pts), knots(kns), weights(wts){ }


    // 递归 Cox-de Boor 实现
    double RationalCurve::BasicCoxdeBoor(double u, int i, int p) {
        // p == 0
        // N_{i,0}(u) = 1 if u_i <= u < u_{i+1}, else 0
        if (p == 0) {
            // 左闭右开区间
            if ((knots[i] <= u) && (u < knots[i + 1])) {
                return 1.0;
            }
            else {
                return 0.0;
            }
        }

        // 递归步骤
        // N_{i,p}(u) = Term1 + Term2
        // Term1 = (u - u_i) / (u_{i+p} - u_i) * N_{i, p-1}(u)
        // Term2 = (u_{i+p+1} - u) / (u_{i+p+1} - u_{i+1}) * N_{i+1, p-1}(u)

        double leftTerm = 0.0;
        double rightTerm = 0.0;

        // 计算左项
        double denom1 = knots[i + p] - knots[i];
        if (denom1 != 0.0) { // 约定 0/0 = 0
            double coeff1 = (u - knots[i]) / denom1;
            leftTerm = coeff1 * BasicCoxdeBoor(u, i, p - 1);
        }

        // 计算右项
        double denom2 = knots[i + p + 1] - knots[i + 1];
        if (denom2 != 0.0) { // 约定 0/0 = 0
            double coeff2 = (knots[i + p + 1] - u) / denom2;
            rightTerm = coeff2 * BasicCoxdeBoor(u, i + 1, p - 1);
        }

        return leftTerm + rightTerm;
    }


    // 计算 NURBS 曲线上的点 C(u)
    glm::dvec3 RationalCurve::Evaluate(double u) {
        // 处理 u 的范围边界
        double min_u = knots[degree];
        double max_u = knots[knots.size() - degree - 1]; // 或者 knots[control_points.size()]

        // 简单的范围钳制
        if (u < min_u) u = min_u;
        if (u > max_u) u = max_u - 1e-9; // 防止取到最后一个节点导致 basis 为 0

        glm::dvec3 numerator(0.0);
        double denominator = 0.0;

        // NURBS 公式: C(u) = (Sum(N_{i,p} * w_i * P_i)) / (Sum(N_{i,p} * w_i))
        int n = control_points.size();

        for (int i = 0; i < n; ++i) {
            // 调用递归函数计算 N_{i,p}(u)
            double Nip = BasicCoxdeBoor(u, i, degree);

            if (Nip > 1e-9) { // 优化：只计算非零项
                double val = Nip * weights[i];
                numerator += val * control_points[i];
                denominator += val;
            }
        }

        if (denominator == 0.0) return glm::dvec3(0.0); // 避免除以0
        return numerator / denominator;
    }

    bool RationalCurve::isNurbs() {
        double maxWeight = -65536.0;
        double minWeight = +65536.0;
        for (int i = 0; i < weights.size(); i++) {
            if (weights[i] > maxWeight) {
                maxWeight = weights[i];
            }
            if (weights[i] < minWeight) {
                minWeight = weights[i];
            }
        }
        if (abs(maxWeight - minWeight) < Tolerance) {
            return false;// 说明是 B-Spline 曲线
        }
        else {
            return true;// 说明是 NURBS
        }
    }

    // 递归计算基函数导数
    // 公式: N'_{i,p} = p / (u_{i+p} - u_i) * N_{i,p-1} - p / (u_{i+p+1} - u_{i+1}) * N_{i+1,p-1}
    double RationalCurve::BasicDerivative(double u, int i, int p, int k) {
        // 0阶导数就是基函数本身
        if (k == 0) {
            return BasicCoxdeBoor(u, i, p);
        }

        // 如果求导阶数超过曲线阶数，导数为0 (常数的导数为0)
        if (k > p) {
            return 0.0;
        }

        // 递归计算
        // Term 1
        double leftTerm = 0.0;
        double denom1 = knots[i + p] - knots[i];
        if (denom1 != 0.0) {
            leftTerm = (p / denom1) * BasicDerivative(u, i, p - 1, k - 1);
        }

        // Term 2
        double rightTerm = 0.0;
        double denom2 = knots[i + p + 1] - knots[i + 1];
        if (denom2 != 0.0) {
            rightTerm = (p / denom2) * BasicDerivative(u, i + 1, p - 1, k - 1);
        }

        return leftTerm - rightTerm;
    }

    std::vector<glm::dvec3> RationalCurve::Derivatives(double u, int d) {
        // 限制 u 的范围
        double min_u = knots[degree];
        double max_u = knots[knots.size() - degree - 1];
        if (u < min_u) u = min_u;
        if (u > max_u) u = max_u - 1e-9;

        // 阶数 d 不能超过 min(degree, d) 其实可以，但高阶为0，这里不做硬性限制
        // 结果数组: index 0 = Point, 1 = 1st Deriv, 2 = 2nd Deriv ...
        std::vector<glm::dvec3> CK(d + 1);

        // 临时存储 A(u) 的各阶导数 (A是加权控制点 sum(Ni * wi * Pi))
        std::vector<glm::dvec3> Aders(d + 1, glm::dvec3(0.0));
        // 临时存储 w(u) 的各阶导数 (w是权重和 sum(Ni * wi))
        std::vector<double> wders(d + 1, 0.0);

        int n = control_points.size();

        // 1. 计算 A^{(k)} 和 w^{(k)}
        // 遍历每个求导阶数 k
        for (int k = 0; k <= d; ++k) {
            // 遍历所有控制点 (优化：实际只需要遍历受影响的区间，这里为了代码简单遍历所有)
            for (int i = 0; i < n; ++i) {
                // 计算第 i 个基函数的 k 阶导数
                double dN = BasicDerivative(u, i, degree, k);

                if (std::abs(dN) > 1e-9) {
                    glm::dvec3 P_weighted = control_points[i] * weights[i];

                    Aders[k] += dN * P_weighted; // A 的导数累加
                    wders[k] += dN * weights[i]; // w 的导数累加
                }
            }
        }

        // 2. 根据商求导法则计算 C^{(k)}(u)
        // 莱布尼茨公式变体: C^{(k)} = (1/w) * ( A^{(k)} - sum_{i=1 to k} ( C(k,i) * w^{(i)} * C^{(k-i)} ) )
        for (int k = 0; k <= d; ++k) {
            glm::dvec3 v = Aders[k];

            for (int i = 1; i <= k; ++i) {
                int binom = Utils::Binomial(k, i);
                v -= (double)binom * wders[i] * CK[k - i];
            }

            if (wders[0] != 0.0) { // wders[0] 就是 w(u)
                CK[k] = v / wders[0];
            }
            else {
                CK[k] = glm::dvec3(0.0);
            }
        }

        return CK;
    }




    // -----------------------------------------------------------------------
    // RationalSurface 实现
    // -----------------------------------------------------------------------

    bool RationalSurface::Create(int udeg, int vdeg,
        const std::vector<glm::dvec3>& ctrl_pts,
        const std::vector<double>& kns_u,
        const std::vector<double>& kns_v,
        const std::vector<double>& wts,
        RationalSurface*& surface) {

        if (surface != nullptr) {
            throw std::runtime_error("Surface pointer must be null!");
            return false;
        }
        if (udeg < 1 || vdeg < 1) {
            throw std::runtime_error("Degrees must be > 0");
            return false;
        }

        // 计算维度: KnotSize = ControlPointCount + Degree + 1
        // 所以: ControlPointCount = KnotSize - Degree - 1
        int n_u = kns_u.size() - udeg - 1;
        int n_v = kns_v.size() - vdeg - 1;

        if (n_u <= 0 || n_v <= 0) {
            throw std::runtime_error("Invalid knot vector sizes.");
            return false;
        }

        if (ctrl_pts.size() != n_u * n_v) {
            throw std::runtime_error("Control points size does not match knot vectors definition (Expect u_count * v_count).");
            return false;
        }
        if (wts.size() != ctrl_pts.size()) {
            throw std::runtime_error("Weights size mismatch.");
            return false;
        }

        surface = new RationalSurface(udeg, vdeg, ctrl_pts, kns_u, kns_v, wts);
        // 保存计算出的维度，方便后续使用
        surface->count_u = n_u;
        surface->count_v = n_v;

        return true;
    }

    RationalSurface::RationalSurface(int udeg, int vdeg,
        const std::vector<glm::dvec3>& ctrl_pts,
        const std::vector<double>& kns_u,
        const std::vector<double>& kns_v,
        const std::vector<double>& wts)
        : udegree(udeg), vdegree(vdeg),
        control_points(ctrl_pts), knots_u(kns_u), knots_v(kns_v), weights(wts) {
        // 构造函数里初始化一下维度，防止未通过Create调用的情况（虽然这里是private）
        count_u = kns_u.size() - udeg - 1;
        count_v = kns_v.size() - vdeg - 1;
    }

    // 通用的 Basis Function (复用逻辑，但需要传入特定的 knot 数组)
    double RationalSurface::BasicCoxdeBoor(double t, int i, int p, const std::vector<double>& knots) {
        if (p == 0) {
            if ((knots[i] <= t) && (t < knots[i + 1])) return 1.0;
            return 0.0;
        }

        double left = 0.0;
        double denom1 = knots[i + p] - knots[i];
        if (denom1 != 0.0) {
            left = ((t - knots[i]) / denom1) * BasicCoxdeBoor(t, i, p - 1, knots);
        }

        double right = 0.0;
        double denom2 = knots[i + p + 1] - knots[i + 1];
        if (denom2 != 0.0) {
            right = ((knots[i + p + 1] - t) / denom2) * BasicCoxdeBoor(t, i + 1, p - 1, knots);
        }

        return left + right;
    }

    glm::dvec3 RationalSurface::Evaluate(double u, double v) {
        // 范围钳制 (Clamp)
        double min_u = knots_u[udegree];
        double max_u = knots_u[knots_u.size() - udegree - 1];
        if (u < min_u) u = min_u;
        if (u > max_u) u = max_u - 1e-9;

        double min_v = knots_v[vdegree];
        double max_v = knots_v[knots_v.size() - vdegree - 1];
        if (v < min_v) v = min_v;
        if (v > max_v) v = max_v - 1e-9;

        glm::dvec3 numerator(0.0);
        double denominator = 0.0;

        // 遍历所有控制点 (O(N*M)) - 生产环境通常只遍历受影响的区间 (degree+1)*(degree+1)
        // 但为了递归实现的简洁性，这里全量遍历
        for (int i = 0; i < count_u; ++i) {
            // 优化：如果 u 方向基函数为 0，整行都不用算了
            double Nu = BasicCoxdeBoor(u, i, udegree, knots_u);
            if (std::abs(Nu) < 1e-9) continue;

            for (int j = 0; j < count_v; ++j) {
                double Nv = BasicCoxdeBoor(v, j, vdegree, knots_v);

                if (std::abs(Nv) > 1e-9) {
                    int idx = i * count_v + j; // 行优先索引
                    double w = weights[idx];
                    double val = Nu * Nv * w;

                    numerator += val * control_points[idx];
                    denominator += val;
                }
            }
        }

        if (denominator == 0.0) return glm::dvec3(0.0);
        return numerator / denominator;
    }

    // 通用导数基函数
    double RationalSurface::BasicDerivative(double t, int i, int p, int k, const std::vector<double>& knots) {
        if (k == 0) return BasicCoxdeBoor(t, i, p, knots);
        if (k > p) return 0.0;

        double leftTerm = 0.0;
        double denom1 = knots[i + p] - knots[i];
        if (denom1 != 0.0) {
            leftTerm = (p / denom1) * BasicDerivative(t, i, p - 1, k - 1, knots);
        }

        double rightTerm = 0.0;
        double denom2 = knots[i + p + 1] - knots[i + 1];
        if (denom2 != 0.0) {
            rightTerm = (p / denom2) * BasicDerivative(t, i + 1, p - 1, k - 1, knots);
        }
        return leftTerm - rightTerm;
    }

    // 计算曲面偏导数
    // 返回结果是平铺的，顺序：
    // d=0: [S]
    // d=1: [S, Su, Sv]
    // d=2: [S, Su, Sv, Suu, Suv, Svv]
    std::vector<glm::dvec3> RationalSurface::Derivatives(double u, double v, int d) {
        // Clamp u, v (同 Evaluate)
        double max_u = knots_u[knots_u.size() - udegree - 1];
        if (u > max_u) u = max_u - 1e-9;
        double max_v = knots_v[knots_v.size() - vdegree - 1];
        if (v > max_v) v = max_v - 1e-9;

        // 1. 计算分子 A(u,v) 和分母 w(u,v) 的所有混合偏导数
        // Aders[k][l] 表示 d^(k+l)A / du^k dv^l
        std::vector<std::vector<glm::dvec3>> Aders(d + 1, std::vector<glm::dvec3>(d + 1, glm::dvec3(0.0)));
        std::vector<std::vector<double>> wders(d + 1, std::vector<double>(d + 1, 0.0));

        // 预计算基函数导数以优化性能
        // U方向: [i][k] 表示第i个基函数的k阶导
        std::vector<std::vector<double>> Nu_ders(count_u, std::vector<double>(d + 1));
        std::vector<std::vector<double>> Nv_ders(count_v, std::vector<double>(d + 1));

        for (int i = 0; i < count_u; ++i)
            for (int k = 0; k <= d; ++k)
                Nu_ders[i][k] = BasicDerivative(u, i, udegree, k, knots_u);

        for (int j = 0; j < count_v; ++j)
            for (int l = 0; l <= d; ++l)
                Nv_ders[j][l] = BasicDerivative(v, j, vdegree, l, knots_v);

        // 累加计算 A 和 w 的偏导
        for (int i = 0; i < count_u; ++i) {
            for (int j = 0; j < count_v; ++j) {
                int idx = i * count_v + j;
                glm::dvec3 P_w = control_points[idx] * weights[idx];
                double W = weights[idx];

                // 仅当基函数非零时计算
                // 这里可以优化：只遍历 Nu_ders[i][0] != 0 的区间

                for (int k = 0; k <= d; ++k) {
                    double dNu = Nu_ders[i][k];
                    if (std::abs(dNu) < 1e-9) continue; // 剪枝

                    for (int l = 0; l <= d - k; ++l) { // 只需要计算 k+l <= d 的部分
                        double dNv = Nv_ders[j][l];

                        // Mixed Partial: d/du^k * d/dv^l
                        // For Tensor Product: Partial = N_i^(k)(u) * N_j^(l)(v)
                        double term = dNu * dNv;

                        Aders[k][l] += term * P_w;
                        wders[k][l] += term * W;
                    }
                }
            }
        }

        // 2. 使用广义莱布尼茨法则（商法则）求解 S 的偏导数 S^(k,l)
        // SKL[k][l] 存储结果
        std::vector<std::vector<glm::dvec3>> SKL(d + 1, std::vector<glm::dvec3>(d + 1, glm::dvec3(0.0)));
        std::vector<glm::dvec3> result;

        // 按阶数从小到大计算，因为高阶依赖低阶
        for (int order = 0; order <= d; ++order) {
            for (int k = 0; k <= order; ++k) {
                int l = order - k; // 当前计算 S_{k,l}

                glm::dvec3 v_num = Aders[k][l];

                // 减去所有 (i,j) < (k,l) 的 S_{i,j} * w_{k-i, l-j} * Binom
                // 公式: A^(k,l) = Sum_{i=0}^k Sum_{j=0}^l Binom(k,i)Binom(l,j) w^(i,j) S^(k-i, l-j)
                // 移项得: w^(0,0) S^(k,l) = A^(k,l) - [上述Sum 但排除 i=0,j=0 项]

                for (int i = 1; i <= k; ++i) {
                    double binom = (double)Utils::Binomial(k, i) * wders[i][0];
                    v_num -= binom * SKL[k - i][l];
                }
                for (int j = 1; j <= l; ++j) {
                    double binom = (double)Utils::Binomial(l, j) * wders[0][j];
                    v_num -= binom * SKL[k][l - j];
                }
                // 处理混合项 (上面的循环其实没有覆盖所有交叉项，需要完整的双重循环)
                // 正确的求和: Sum(i=0..k) Sum(j=0..l) of Binom * Binom * w * S
                // 排除 i=0, j=0 的情况

                // 重新实现减法部分：
                v_num = Aders[k][l];
                for (int i = 0; i <= k; ++i) {
                    for (int j = 0; j <= l; ++j) {
                        if (i == 0 && j == 0) continue; // 这一项包含我们要求的 S[k][l]，跳过

                        double bin = (double)Utils::Binomial(k, i) * (double)Utils::Binomial(l, j);
                        v_num -= bin * wders[i][j] * SKL[k - i][l - j];
                    }
                }

                if (wders[0][0] != 0.0) {
                    SKL[k][l] = v_num / wders[0][0];
                }

                result.push_back(SKL[k][l]);
            }
        }

        return result;
    }


    glm::dvec3 RationalSurface::GetNormal(double u, double v) {
        // 我们只需要一阶导数来计算切向量
        // Derivatives 返回结果顺序分析:
        // Index 0: Point
        // Order 1, k=0 (l=1): Sv -> Index 1
        // Order 1, k=1 (l=0): Su -> Index 2
        std::vector<glm::dvec3> derivs = Derivatives(u, v, 1);

        glm::dvec3 Sv = derivs[1];
        glm::dvec3 Su = derivs[2];

        // 法线 = Su x Sv (注意叉乘顺序，右手定则)
        glm::dvec3 normal = glm::cross(Su, Sv);

        // 归一化处理
        double len = glm::length(normal);

        // 处理奇异点（例如导数为0的情况），避免除以0
        if (len < 1e-9) {
            return glm::dvec3(0.0, 0.0, 1.0); // 或者返回 0 向量，视需求而定
        }

        return normal / len;
    }
    
}
