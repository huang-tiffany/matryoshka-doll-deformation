#ifndef MANIP_H
#define MANIP_H
#include <Eigen/Eigen>
#include <igl/copyleft/cgal/intersect_other.h>
#include <igl/per_vertex_normals.h>
#include <igl/copyleft/cgal/mesh_boolean.h>

Eigen::MatrixXd resize_mesh(Eigen::MatrixXd V, double scale);

bool meshes_intersect(
    const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA,
    const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB);

Eigen::MatrixXd offset_along_normals(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double offset);

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> create_shell(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double offset);

#endif // MANIP_H
