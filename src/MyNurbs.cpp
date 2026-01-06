#include "MyNurbs.h"

namespace MyNurbs {
	
	bool RationalCurve::Create(int deg,
		const std::vector<glm::dvec3>& ctrl_pts,
		const std::vector<double>& kns,
		const std::vector<double> wts,
		RationalCurve* curve) {

		// 传入指针必须为空
		if (curve != nullptr) {
			std::throw_with_nested(std::runtime_error("Curve pointer must be null!");
			return false;
		}

		// 检查数量是否符合BSpline曲线规范
		if (deg < 1) {
			std::throw_with_nested(std::runtime_error("Curve degree must greater than 0!");
			return false;
		}
		if (ctrl_pts.size() < deg + 1) {
			std::throw_with_nested(std::runtime_error("Control points number must greater than degree!");
			return false;
		}
		if (kns.size() != ctrl_pts.size() + deg + 1) {
			std::throw_with_nested(std::runtime_error("Invalid knot vector size.");
			return false;
		}

		// 创建曲线
		curve = new RationalCurve(deg, ctrl_pts, kns, wts);
		return true;
	}

	RationalCurve::RationalCurve(int deg,
		const std::vector<glm::dvec3>& ctrl_pts,
		const std::vector<double>& kns,
		const std::vector<double> wts) : degree(deg), control_points(ctrl_pts), knots(kns), weights(wts);
}