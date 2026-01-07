#pragma once

#include <vector>
#include "glm/glm.hpp"

#include <cmath>
#include <iostream>
#include <cassert>

namespace MyNurbs {

#define Tolerance 1e-6

    namespace Utils {

        // 辅助函数：计算组合数 C(n, k)
        int Binomial(int n, int k) {
            if (k < 0 || k > n) {
                return 0;
            }
            if (k == 0 || k == n) {
                return 1;
            }
            if (k > n / 2) {
                k = n - k;
            }

            long result = 1;
            for (int i = 1; i <= k; ++i) {
                result = result * (n - i + 1) / i;
            }
            return (int)result;
        }


        // 定义矩阵类型，简单的 vector<vector<double>>
        using Matrix = std::vector<std::vector<double>>;
        using Vector = std::vector<double>;
        // 打印矩阵辅助函数（调试用）
        void printMatrix(const Matrix& A) {
            for (const auto& row : A) {
                for (double val : row) printf("%.4f ", val);
                printf("\n");
            }
        }
        // 核心：高斯消元法求解 Ax = B
        // 带有列主元选择（Partial Pivoting），防止除以0，提高数值稳定性
        Vector solveLinearSystem(Matrix A, Vector B) {
            int n = A.size();
            if (n == 0) return {};
            assert(A[0].size() == n && B.size() == n); // 确保是方阵

            // 消元过程
            for (int i = 0; i < n; i++) {
                // 列主元选择
                int pivot = i;
                for (int j = i + 1; j < n; j++) {
                    if (std::abs(A[j][i]) > std::abs(A[pivot][i])) {
                        pivot = j;
                    }
                }
                // 交换行
                std::swap(A[i], A[pivot]);
                std::swap(B[i], B[pivot]);

                if (std::abs(A[i][i]) < 1e-9) {
                    std::cerr << "Error: Matrix is singular or nearly singular!" << std::endl;
                    return {}; // 矩阵奇异，无解或无穷多解
                }

                // 消元
                for (int j = i + 1; j < n; j++) {
                    double factor = A[j][i] / A[i][i];
                    B[j] -= factor * B[i];
                    for (int k = i; k < n; k++) {
                        A[j][k] -= factor * A[i][k];
                    }
                }
            }

            // 2. 回代过程 (Back Substitution)
            Vector x(n);
            for (int i = n - 1; i >= 0; i--) {
                double sum = 0;
                for (int j = i + 1; j < n; j++) {
                    sum += A[i][j] * x[j];
                }
                x[i] = (B[i] - sum) / A[i][i];
            }

            return x;
        }


        // 矩阵转置
        Matrix transpose(const Matrix& A) {
            if (A.empty()) return {};
            int rows = A.size();
            int cols = A[0].size();
            Matrix AT(cols, Vector(rows));
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    AT[j][i] = A[i][j];
                }
            }
            return AT;
        }

        // 矩阵乘法 C = A * B
        Matrix multiply(const Matrix& A, const Matrix& B) {
            if (A.empty() || B.empty()) return {};
            int r1 = A.size();
            int c1 = A[0].size();
            int r2 = B.size();
            int c2 = B[0].size();
            assert(c1 == r2);

            Matrix C(r1, Vector(c2, 0.0));
            for (int i = 0; i < r1; i++) {
                for (int k = 0; k < c1; k++) {
                    if (std::abs(A[i][k]) < 1e-9) continue; // 简单的稀疏优化
                    for (int j = 0; j < c2; j++) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            return C;
        }

        // 矩阵乘以向量 v_out = A * v_in
        Vector multiplyMV(const Matrix& A, const Vector& v) {
            int rows = A.size();
            int cols = A[0].size();
            assert(cols == v.size());
            Vector res(rows, 0.0);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    res[i] += A[i][j] * v[j];
                }
            }
            return res;
        }
    }

	class RationalCurve {
	public:
		// 曲线阶数
		int degree;
		// 控制点
		std::vector<glm::dvec3> control_points;
		// 节点
		std::vector<double> knots;
		// 权重
		std::vector<double> weights;

        enum InterpolateMethod {
            Uniform,
            ChordLength,
            Centripetal,
            Universal
        };

        /// <summary>
        /// 插值创建一个曲线
        /// </summary>
        /// <param name="degree"></param>
        /// <param name="points"></param>
        /// <param name="curve"></param>
        /// <returns></returns>
        static bool Interpolate(
            InterpolateMethod method,
            int degree,
            const std::vector<glm::dvec3>& points,
            RationalCurve*& curve
        );

		/// <summary>
		/// 静态创建一个曲线
		/// </summary>
		/// <param name="deg">曲线阶数</param>
		/// <param name="ctrl_pts">控制点数组</param>
		/// <param name="kns">节点数组</param>
		/// <param name="wts">权重</param>
		/// <param name="curve">返回结果</param>
		/// <returns></returns>
		static bool Create(int deg,
			const std::vector<glm::dvec3>& ctrl_pts,
			const std::vector<double>& kns,
			const std::vector<double>& wts,
			RationalCurve*& curve);



        /// <summary>
        /// 计算曲线上一点的位置
        /// </summary>
        /// <param name="u"></param>
        /// <returns></returns>
        glm::dvec3 Evaluate(double u);

        /// <summary>
        /// 计算曲线在 u 处的 n 阶导数
        /// </summary>
        /// <param name="u">参数</param>
        /// <param name="d">需要求导的最高阶数 (例如: 1求速度, 2求加速度)</param>
        /// <returns>返回一个数组，第0个是点坐标，第1个是一阶导，第2个是二阶导...</returns>
        std::vector<glm::dvec3> Derivatives(double u, int d);

        // 析构函数 (建议加上)
        ~RationalCurve() = default;

	private:
        /// <summary>
        /// 均匀参数化插值
        /// </summary>
        /// <param name="degree"></param>
        /// <param name="points"></param>
        /// <param name="curve"></param>
        /// <returns></returns>
        static bool InterpolateUniform(
            int degree,
            const std::vector<glm::dvec3>& points,
            RationalCurve*& curve
        );

		/// <summary>
		/// 构造函数
		/// </summary>
		/// <param name="deg">曲线阶数</param>
		/// <param name="ctrl_pts">控制点</param>
		/// <param name="kns">节点数组</param>
		/// <param name="wts">权重</param>
        /// <param name="curve">返回结果</param>
		RationalCurve(int deg,
			const std::vector<glm::dvec3>& ctrl_pts,
			const std::vector<double>& kns,
			const std::vector<double>& wts);

        /// <summary>
        /// 递归计算 B-Spline 基函数 N_{i,p}(u)
        /// </summary>
        /// <param name="u">参数值</param>
        /// <param name="i">节点索引</param>
        /// <param name="p">当前阶数</param>
        /// <returns>基函数的值</returns>
        double BasicCoxdeBoor(double u, int i, int p);

        /// <summary>
        /// 返回当前曲线是否为 NURBS 曲线，如果返回结果为 true，那么说明是 NURBS 曲线，否则是 B-Spline 曲线。
        /// </summary>
        /// <returns></returns>
        bool isNurbs();

        /// <summary>
        /// 递归计算基函数的 k 阶导数 N^{(k)}_{i,p}(u)
        /// </summary>
        double BasicDerivative(double u, int i, int p, int k);
	};

    class RationalSurface {
    public:
        // 曲面阶数
        int udegree;
        int vdegree;

        // 网格维度
        int count_u;
        int count_v;

        // 控制点（行优先）
        std::vector<glm::dvec3> control_points;
        // 节点
        // 节点
        std::vector<double> knots_u;
        std::vector<double> knots_v;
        // 权重（行优先）
        std::vector<double> weights;

        /// <summary>
        /// 创建曲面
        /// </summary>
        /// <param name="udeg">u方向阶数</param>
        /// <param name="vdeg">v方向阶数</param>
        /// <param name="ctrl_pts">控制点</param>
        /// <param name="kns">节点</param>
        /// <param name="wts">权重</param>
        /// <param name="surface">返回结果</param>
        /// <returns></returns>
        static bool Create(int udeg, int vdeg,
            const std::vector<glm::dvec3>& ctrl_pts,
            const std::vector<double>& kns_u,
            const std::vector<double>& kns_v,
            const std::vector<double>& wts,
            RationalSurface*& surface);

        /// <summary>
        /// 计算曲面上一点的位置
        /// </summary>
        /// <param name="u"></param>
        /// <returns></returns>
        glm::dvec3 Evaluate(double u, double v);


        /// <summary>
        /// 计算曲线在 u 处的 n 阶导数
        /// </summary>
        /// <param name="u">参数</param>
        /// <param name="d">需要求导的最高阶数 (例如: 1求速度, 2求加速度)</param>
        /// <returns>返回一个数组，第0个是点坐标，第1个是一阶导，第2个是二阶导...</returns>
        std::vector<glm::dvec3> Derivatives(double u, double v, int d);

        
        /// <summary>
        /// 计算曲面在 (u, v) 处的法线
        /// </summary>
        glm::dvec3 GetNormal(double u, double v);

        // 析构函数 (建议加上)
        ~RationalSurface() = default;

    private:
        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="udeg">u方向阶数</param>
        /// <param name="vdeg">v方向阶数</param>
        /// <param name="ctrl_pts">控制点</param>
        /// <param name="kns">节点数组</param>
        /// <param name="wts">权重</param>
        /// <param name="curve">返回结果</param>
        RationalSurface(int udeg, int vdeg,
            const std::vector<glm::dvec3>& ctrl_pts,
            const std::vector<double>& kns_u,
            const std::vector<double>& kns_v,
            const std::vector<double>& wts);


        /// <summary>
        /// 递归计算 B-Spline 基函数 N_{i,p}(u)
        /// </summary>
        /// <param name="u">u参数值</param
        /// <param name="v">v参数值</param>
        /// <param name="i">节点索引</param>
        /// <param name="p">当前阶数</param>
        /// <returns>基函数的值</returns>
        double BasicCoxdeBoor(double t, int i, int p, const std::vector<double>& knots);

        /// <summary>
        /// 返回当前曲面是否为 NURBS 曲面，如果返回结果为 true，那么说明是 NURBS 曲面，否则是 B-Spline 曲面。
        /// </summary>
        /// <returns></returns>
        bool isNurbs();

        /// <summary>
        /// 递归计算基函数的 k 阶导数 N^{(k)}_{i,p}(u)
        /// </summary>
        double BasicDerivative(double t, int i, int p, int k, const std::vector<double>& knots);
    };
}
