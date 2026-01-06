#pragma once

#include <vector>
#include "glm/glm.hpp"

namespace MyNurbs {

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

		bool Create(int deg,
			const std::vector<glm::dvec3>& ctrl_pts,
			const std::vector<double>& kns,
			const std::vector<double> wts,
			RationalCurve* curve);

	private:
		/// <summary>
		/// 构造函数
		/// </summary>
		/// <param name="deg">曲线阶数</param>
		/// <param name="ctrl_pts">控制点数组</param>
		/// <param name="kns">节点数组</param>
		/// <param name="wts">权重数组</param>
		RationalCurve(int deg,
			const std::vector<glm::dvec3>& ctrl_pts,
			const std::vector<double>& kns,
			const std::vector<double> wts);
	};

}
