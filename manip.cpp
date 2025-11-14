#include "manip.h"
#include <igl/copyleft/cgal/intersect_with_half_space.h>

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

bool split_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double y_coord) {
    Eigen::Vector3d point(0, y_coord, 0);
    Eigen::Vector3d normal_up(0, 1, 0);
    Eigen::Vector3d normal_down(0, -1, 0);

    Eigen::MatrixXd bottom_vertices;
    Eigen::MatrixXi bottom_faces;
    Eigen::MatrixXi bottom_indices;

    Eigen::MatrixXd top_vertices;
    Eigen::MatrixXi top_faces;
    Eigen::MatrixXi top_indices;

    // bool success_top = igl::copyleft::cgal::intersect_with_half_space(V, F, point, normal_up, top_vertices, top_faces, top_indices);
    // bool success_bottom = igl::copyleft::cgal::intersect_with_half_space(V, F, point, normal_up, bottom_vertices, bottom_faces, bottom_indices);

    // return success_top && success_bottom;
    return true;
}
