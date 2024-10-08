#include <cmath>
#include <iostream>

#include <Math/Vector3D.h>
#include <TGraph.h>

constexpr double pi = 3.1415927;

/// @brief intersection of special plane and line 特殊平面和直线的交点
/// @param[in] line_point a point in line 直线上的任意一点
/// @param[in] line_direction line direction 直线方向
/// @param[in] plane_point a point in plane 平面上任意一点
/// @param[in] plane_normal_vector normal vector of plane, z should be zero
///		平面法向量, 要求 z 方向分量为 0
/// @param[out] intersection intersection 交点
/// @returns 0 if intersection is exists, -1 otherwise
///
int IntersectionOfPlaneAndLine(
	const ROOT::Math::XYZVector &line_point,
	const ROOT::Math::XYZVector &line_direction,
	const ROOT::Math::XYZVector &plane_point,
	const ROOT::Math::XYZVector &plane_normal_vector,
	ROOT::Math::XYZVector &intersection
) {
	double k1 = plane_normal_vector.Y() / plane_normal_vector.X();
	double k2 = line_direction.X() / line_direction.Y();
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	if (fabs(plane_normal_vector.X()) < 1e-7) {
		if (
			fabs(line_direction.Y()) < 1e-7
			&& fabs(line_point.Y() - plane_point.Y()) > 1e-7
		) {
			return -1;
		}
		y = plane_point.Y();
	} else if (fabs(line_direction.Y()) < 1e-7) {
		y = line_point.Y();
	} else if (fabs(k1+k2) < 1e-7) {
		return -1;
	} else {
		y = (
			plane_point.Y() * k1
			+ line_point.Y() * k2
			+ plane_point.X()
			- line_point.Y()
		) / (k1 + k2);
	}
	if (fabs(line_direction.Y()) < 1e-7) {
		if (fabs(line_direction.X()) < 1e-7) return -1;
		if (fabs(plane_normal_vector.X() < 1e-7)) return -1;
		x = k1 * (y - plane_point.Y()) + plane_point.X();
		z = line_point.Z()
			+ (x - line_point.X()) * line_direction.Z() / line_direction.X();
	} else {
		x = line_point.X()
			+ (y - line_point.Y()) * line_direction.X() / line_direction.Y();
		z = line_point.Z()
			+ (y - line_point.Y()) * line_direction.Z() / line_direction.Y();
	}
	intersection = ROOT::Math::XYZVector(x, y, z);
	return 0;
}


double WithPlane(
	const ROOT::Math::XYZVector &point,
	const ROOT::Math::XYZVector &plane_point,
	const ROOT::Math::XYZVector &plane_normal_vector
) {
	return (point - plane_point).Unit().Dot(plane_normal_vector.Unit());
}

int main() {
	// CsI geometry parameters
	// distance from most distant point to z-axis (mm)
	const double csi_r = 174.4;
	// bottom angle (alpha) of trapezoidal-shape surface
	const double bottom_angle = asin(104.0 / 107.35);
	// // sin(alpha)
	// constexpr double sin_alpha = 104.0 / 107.35;
	// // cos(alpha)
	// constexpr double cos_alpha = 26.6 / 107.35;
	// upper surface vertex
	ROOT::Math::XYZVector upper_surface_vertex[12];
	// upper surface normal vector
	ROOT::Math::XYZVector upper_surface_normal[12];
	// inner side surface (surface between two CsI) vertex
	ROOT::Math::XYZVector inner_surface_vertex[12];
	// inner side surface (surface between two CsI) normal vector
	ROOT::Math::XYZVector inner_surface_normal[12];
	// outer side surface vertex
	ROOT::Math::XYZVector outer_surface_vertex[12];
	// outer side surface normal vector
	ROOT::Math::XYZVector outer_surface_normal[12];
	// inner bottom vertex
	ROOT::Math::XYZVector inner_bottom_vertex[12];
	// construct plane surface for CsI
	for (size_t i = 0; i < 12; ++i) {
		// upper surface
		double upper_vertex_angle = pi * (1.0/2.0 - 1.0/3.0*(i/2));
		upper_surface_vertex[i] = ROOT::Math::XYZVector(
			csi_r * cos(upper_vertex_angle),
			csi_r * sin(upper_vertex_angle),
			151.8
		);
		double upper_normal_angle = i % 2 == 0
			? pi * (1.0 - (i/2)/3.0) - bottom_angle
			: bottom_angle - pi * (i/2) / 3.0;
		upper_surface_normal[i] = ROOT::Math::XYZVector(
			cos(upper_normal_angle),
			sin(upper_normal_angle),
			0.0
		);
		// inner surface
		inner_surface_vertex[i] = upper_surface_vertex[i];
		double inner_normal_angle = i % 2 == 0
			? -pi * (i/2) / 3.0
			: pi * (1.0 - (i/2)/3.0);
		inner_surface_normal[i] = ROOT::Math::XYZVector(
			cos(inner_normal_angle),
			sin(inner_normal_angle),
			0.0
		);
		// outer surface
		double outer_vertex_angle = i % 2 == 0
			? atan(64.33/157.94) + pi * (0.5 - (i/2)/3.0)
			: atan(157.94/64.33) - pi * (i/2)/3.0;
		double outer_vertex_r = sqrt(64.33*64.33 + 157.94*157.94);
		outer_surface_vertex[i] = ROOT::Math::XYZVector(
			outer_vertex_r * cos(outer_vertex_angle),
			outer_vertex_r * sin(outer_vertex_angle),
			151.8
		);
		double outer_normal_angle = i % 2 == 0
			? (2.0 - (i/2)/3.0) * pi - 2.0 * bottom_angle
			: 2.0 * bottom_angle - (1.0 + (i/2)/3.0) * pi;
		outer_surface_normal[i] = ROOT::Math::XYZVector(
			cos(outer_normal_angle),
			sin(outer_normal_angle),
			0.0
		);

		// bottom vertex
		inner_bottom_vertex[i] = ROOT::Math::XYZVector(
			(csi_r - 107.35) * cos(upper_vertex_angle),
			(csi_r - 107.35) * sin(upper_vertex_angle),
			151.8
		);
	}

	for (size_t i = 0; i < 12; ++i) {
		std::cout << "------------ " << i << " ----------------------\n";
		std::cout << "Outer upper point with upper surface "
			<< WithPlane(
				outer_surface_vertex[i],
				upper_surface_vertex[i],
				upper_surface_normal[i]
			) << "\n";

		ROOT::Math::XYZVector outer_bottom_vertex =
			inner_bottom_vertex[i]
			+ (outer_surface_vertex[i] - upper_surface_vertex[i]).Unit() * 13.2;

		std::cout << "Upper distance "
			<< (upper_surface_vertex[i] - outer_surface_vertex[i]).R() << "\n";
		std::cout << "Bottom distance "
			<< (inner_bottom_vertex[i] - outer_bottom_vertex).R() << "\n";
		std::cout << "Side distance "
			<< (upper_surface_vertex[i] - inner_bottom_vertex[i]).R() << ", "
			<< (outer_surface_vertex[i] - outer_bottom_vertex).R() << "\n";

		std::cout << "Inner bottom vertex with inner surface "
			<< WithPlane(
				inner_bottom_vertex[i],
				inner_surface_vertex[i],
				inner_surface_normal[i]
			) << "\n";
		std::cout << "Outer bottom vertex with outer surface "
			<< WithPlane(
				outer_bottom_vertex,
				outer_surface_vertex[i],
				outer_surface_normal[i]
			) << "\n";
	}

	// check intersection
	ROOT::Math::XYZVector intersection1;
	IntersectionOfPlaneAndLine(
		ROOT::Math::XYZVector(1.0, 1.0, 0.0),
		ROOT::Math::XYZVector(0.0, -1.0, 0.0),
		ROOT::Math::XYZVector(0.0, 0.0, 0.0),
		ROOT::Math::XYZVector(0.0, 1.0, 0.0),
		intersection1
	);
	std::cout << "Intersection1 "
		<< intersection1.X() << ", "
		<< intersection1.Y() << ", "
		<< intersection1.Z() << "\n";

	ROOT::Math::XYZVector intersection2;
	IntersectionOfPlaneAndLine(
		ROOT::Math::XYZVector(1.0, 1.0, 0.0),
		ROOT::Math::XYZVector(1.0, -1.0, 1.0),
		ROOT::Math::XYZVector(0.0, 0.0, 0.0),
		ROOT::Math::XYZVector(0.0, 1.0, 0.0),
		intersection2
	);
	std::cout << "Intersection2 "
		<< intersection2.X() << ", "
		<< intersection2.Y() << ", "
		<< intersection2.Z() << "\n";

	ROOT::Math::XYZVector intersection3;
	IntersectionOfPlaneAndLine(
		ROOT::Math::XYZVector(1.0, 1.0, 0.0),
		ROOT::Math::XYZVector(-1.0, -1.0, 1.0),
		ROOT::Math::XYZVector(0.0, 0.0, 0.0),
		ROOT::Math::XYZVector(1.0, 1.0, 0.0),
		intersection3
	);
	std::cout << "Intersection3 "
		<< intersection3.X() << ", "
		<< intersection3.Y() << ", "
		<< intersection3.Z() << "\n";

	ROOT::Math::XYZVector intersection4;
	IntersectionOfPlaneAndLine(
		ROOT::Math::XYZVector(1.0, 1.0, 0.0),
		ROOT::Math::XYZVector(-1.0, -2.0, 2.0),
		ROOT::Math::XYZVector(0.0, 0.0, 0.0),
		ROOT::Math::XYZVector(1.0, 3.0, 0.0),
		intersection4
	);
	std::cout << "Intersection4 "
		<< intersection4.X() << ", "
		<< intersection4.Y() << ", "
		<< intersection4.Z() << "\n";

	TGraph g1;
	g1.AddPoint(0.0, 0.0);
	g1.AddPoint(1.0, 0.0);
	g1.AddPoint(1.0, 1.0);
	g1.AddPoint(0.0, 1.0);
	std::cout << "(0.5, 0.3) in g1 " << g1.IsInside(0.5, 0.3) << "\n";

	return 0;
}