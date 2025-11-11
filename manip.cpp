#include "manip.h"

Eigen::MatrixXd resize_mesh(Eigen::MatrixXd V, double scale) {
    Eigen::RowVector3d centroid = V.colwise().mean();
    V.rowwise() -= centroid;
    V *= scale;
    V.rowwise() += centroid;
    return V;
}

bool meshes_intersect(
    const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA,
    const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB)
{
    Eigen::VectorXi IF;
    return igl::copyleft::cgal::intersect_other(VA, FA, VB, FB, true, IF);
}
