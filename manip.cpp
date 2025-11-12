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

Eigen::MatrixXd offset_along_normals(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double offset) {
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);
    return V + offset * N;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> create_shell(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double offset) {
    Eigen::MatrixXd V_outer = offset_along_normals(V, F, offset);

    Eigen::MatrixXd V_shell;
    Eigen::MatrixXi F_shell;

    igl::copyleft::cgal::mesh_boolean(V_outer, F, V, F, igl::MESH_BOOLEAN_TYPE_MINUS, V_shell, F_shell);

    return std::make_pair(V_shell, F_shell);
}
