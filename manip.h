#ifndef MANIP_H
#define MANIP_H
#include <Eigen/Eigen>
#include <igl/copyleft/cgal/intersect_other.h>

Eigen::MatrixXd resize_mesh(Eigen::MatrixXd V, double scale);

bool meshes_intersect(
    const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA,
    const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB);

#endif // MANIP_H
