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


bool split_mesh_open(const Eigen::MatrixXd &V_in, const Eigen::MatrixXi &F_in,
                     Eigen::MatrixXd &top_V, Eigen::MatrixXi &top_F,
                     Eigen::MatrixXd &bottom_V, Eigen::MatrixXi &bottom_F,
                     double y_coord) {

    std::vector<Eigen::RowVector3d> top_vertices;
    std::vector<Eigen::RowVector3i> top_faces;
    std::vector<Eigen::RowVector3d> bottom_vertices;
    std::vector<Eigen::RowVector3i> bottom_faces;

    std::map<int, int> top_vertex_map;
    std::map<int, int> bottom_vertex_map;

    // Process each face
    for(int f = 0; f < F_in.rows(); f++) {
        Eigen::RowVector3i face = F_in.row(f);

        // Get vertices of the face
        Eigen::RowVector3d v0 = V_in.row(face(0));
        Eigen::RowVector3d v1 = V_in.row(face(1));
        Eigen::RowVector3d v2 = V_in.row(face(2));

        // Check which side of the plane each vertex is on
        bool above0 = v0(1) >= y_coord;
        bool above1 = v1(1) >= y_coord;
        bool above2 = v2(1) >= y_coord;

        int above_count = above0 + above1 + above2;

        // Face completely above - add to top
        if(above_count == 3) {
            Eigen::RowVector3i new_face;
            for(int i = 0; i < 3; i++) {
                int old_idx = face(i);
                if(top_vertex_map.find(old_idx) == top_vertex_map.end()) {
                    top_vertex_map[old_idx] = top_vertices.size();
                    top_vertices.push_back(V_in.row(old_idx));
                }
                new_face(i) = top_vertex_map[old_idx];
            }
            top_faces.push_back(new_face);
        }
        // Face completely below - add to bottom
        else if(above_count == 0) {
            Eigen::RowVector3i new_face;
            for(int i = 0; i < 3; i++) {
                int old_idx = face(i);
                if(bottom_vertex_map.find(old_idx) == bottom_vertex_map.end()) {
                    bottom_vertex_map[old_idx] = bottom_vertices.size();
                    bottom_vertices.push_back(V_in.row(old_idx));
                }
                new_face(i) = bottom_vertex_map[old_idx];
            }
            bottom_faces.push_back(new_face);
        }
        // Face crosses the plane - need to split it
        else {
            std::vector<Eigen::RowVector3d> vertices = {v0, v1, v2};
            std::vector<bool> above = {above0, above1, above2};
            std::vector<int> old_indices = {face(0), face(1), face(2)};

            std::vector<Eigen::RowVector3d> top_verts;
            std::vector<Eigen::RowVector3d> bottom_verts;
            std::vector<int> top_old_idx;
            std::vector<int> bottom_old_idx;

            // Separate vertices and compute intersections
            for(int i = 0; i < 3; i++) {
                int next = (i + 1) % 3;

                if(above[i]) {
                    top_verts.push_back(vertices[i]);
                    top_old_idx.push_back(old_indices[i]);
                } else {
                    bottom_verts.push_back(vertices[i]);
                    bottom_old_idx.push_back(old_indices[i]);
                }

                // Check if edge crosses plane
                if(above[i] != above[next]) {
                    double t = (y_coord - vertices[i](1)) / (vertices[next](1) - vertices[i](1));
                    Eigen::RowVector3d intersection = vertices[i] + t * (vertices[next] - vertices[i]);
                    top_verts.push_back(intersection);
                    bottom_verts.push_back(intersection);
                    top_old_idx.push_back(-1); // New vertex
                    bottom_old_idx.push_back(-1); // New vertex
                }
            }

            // Create top faces
            if(top_verts.size() >= 3) {
                Eigen::RowVector3i new_face;
                for(int i = 0; i < 3; i++) {
                    int old_idx = top_old_idx[i];
                    if(old_idx >= 0 && top_vertex_map.find(old_idx) != top_vertex_map.end()) {
                        new_face(i) = top_vertex_map[old_idx];
                    } else {
                        new_face(i) = top_vertices.size();
                        top_vertices.push_back(top_verts[i]);
                        if(old_idx >= 0) top_vertex_map[old_idx] = new_face(i);
                    }
                }
                top_faces.push_back(new_face);

                // If quad, add second triangle
                if(top_verts.size() == 4) {
                    Eigen::RowVector3i new_face2;
                    new_face2(0) = new_face(0);
                    new_face2(1) = new_face(2);
                    int old_idx = top_old_idx[3];
                    if(old_idx >= 0 && top_vertex_map.find(old_idx) != top_vertex_map.end()) {
                        new_face2(2) = top_vertex_map[old_idx];
                    } else {
                        new_face2(2) = top_vertices.size();
                        top_vertices.push_back(top_verts[3]);
                        if(old_idx >= 0) top_vertex_map[old_idx] = new_face2(2);
                    }
                    top_faces.push_back(new_face2);
                }
            }

            // Create bottom faces
            if(bottom_verts.size() >= 3) {
                Eigen::RowVector3i new_face;
                for(int i = 0; i < 3; i++) {
                    int old_idx = bottom_old_idx[i];
                    if(old_idx >= 0 && bottom_vertex_map.find(old_idx) != bottom_vertex_map.end()) {
                        new_face(i) = bottom_vertex_map[old_idx];
                    } else {
                        new_face(i) = bottom_vertices.size();
                        bottom_vertices.push_back(bottom_verts[i]);
                        if(old_idx >= 0) bottom_vertex_map[old_idx] = new_face(i);
                    }
                }
                bottom_faces.push_back(new_face);

                // If quad, add second triangle
                if(bottom_verts.size() == 4) {
                    Eigen::RowVector3i new_face2;
                    new_face2(0) = new_face(0);
                    new_face2(1) = new_face(2);
                    int old_idx = bottom_old_idx[3];
                    if(old_idx >= 0 && bottom_vertex_map.find(old_idx) != bottom_vertex_map.end()) {
                        new_face2(2) = bottom_vertex_map[old_idx];
                    } else {
                        new_face2(2) = bottom_vertices.size();
                        bottom_vertices.push_back(bottom_verts[3]);
                        if(old_idx >= 0) bottom_vertex_map[old_idx] = new_face2(2);
                    }
                    bottom_faces.push_back(new_face2);
                }
            }
        }
    }

    // Convert to Eigen matrices
    top_V.resize(top_vertices.size(), 3);
    for(size_t i = 0; i < top_vertices.size(); i++) {
        top_V.row(i) = top_vertices[i];
    }

    top_F.resize(top_faces.size(), 3);
    for(size_t i = 0; i < top_faces.size(); i++) {
        top_F.row(i) = top_faces[i];
    }

    bottom_V.resize(bottom_vertices.size(), 3);
    for(size_t i = 0; i < bottom_vertices.size(); i++) {
        bottom_V.row(i) = bottom_vertices[i];
    }

    bottom_F.resize(bottom_faces.size(), 3);
    for(size_t i = 0; i < bottom_faces.size(); i++) {
        bottom_F.row(i) = bottom_faces[i];
    }

    return true;
}
