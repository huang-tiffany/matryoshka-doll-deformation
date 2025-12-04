#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <imgui.h>
#include "portable-file-dialogs.h"
#include <igl/harmonic.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/swept_volume.h>
#include <igl/boundary_loop.h>
#include <iostream>
#include "manip.h"

enum class AppMode {
    EDIT_DEFORMATION,      // User is editing the base mesh
    NESTED_DOLL_DESIGN     // User is customizing nested dolls
};

int main(int argc, char *argv[])
{
    // ============ STATE MANAGEMENT ============
    AppMode current_mode = AppMode::EDIT_DEFORMATION;

    // Original mesh (never modified after loading)
    Eigen::MatrixXd V_original;
    Eigen::MatrixXi F_original;

    // Working mesh (used as base for deformation)
    Eigen::MatrixXd V_working;

    // Final deformed mesh (result after applying deformation)
    Eigen::MatrixXd V_deformed;

    // Mesh for nested doll design
    Eigen::MatrixXd V_base;  // Base template for nested mesh
    Eigen::MatrixXi F;       // Faces (shared across all meshes)

    Eigen::MatrixXd V_outer, V_outer_inner_shell;
    Eigen::MatrixXi F_inner_flipped;
    Eigen::MatrixXd nested_mesh;

    std::vector<Eigen::MatrixXd> V_dolls;
    std::vector<Eigen::MatrixXi> F_dolls;

    // ============ BIHARMONIC DEFORMATION STATE ============
    Eigen::VectorXi boundary_vertices;
    bool is_dragging = false;
    int selected_vertex = -1;
    Eigen::Vector3d drag_start_pos;
    float selection_radius = 0.05f;  // Radius for selecting vertices to drag
    std::vector<int> handle_vertices;  // Vertices in the handle region
    Eigen::MatrixXd handle_rest_positions;  // Rest positions of handle vertices

    // ============ UI PARAMETERS ============
    float slider_value = 1.f;
    float shell_thickness = 0.001f;
    int n_nests = 1;

    // Inner mesh translation
    float translate_x = 0.0f;
    float translate_y = 0.0f;
    float translate_z = 0.0f;

    // Deformation mode toggle
    int start_from_original = 0;  // 0 = start from last edited, 1 = start from original

    // ============ SWEPT VOLUME PARAMETERS ============
    const int grid_size = 75;
    const int time_steps = 200;
    const double isolevel = 0;

    // ============ CUTTING PLANE PARAMETERS ============
    float cutting_plane_y_coord = 0.5f;
    Eigen::MatrixXd V_plane;
    Eigen::MatrixXi F_plane;
    int plane_mesh_id = -1;
    double plane_size = 2.0;

    // ============ VIEWER SETUP ============
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    // ============ HELPER FUNCTIONS ============

    auto clear_all_meshes = [&]() {
        while(viewer.erase_mesh(viewer.selected_data_index)){};
        viewer.data().clear();
        plane_mesh_id = -1;
        cutting_plane_y_coord = 0.5f;
    };

    auto find_boundary_vertices = [&]() {
        // Find boundary vertices to use as fixed constraints
        std::vector<std::vector<int>> boundary_loops;
        igl::boundary_loop(F, boundary_loops);

        if(!boundary_loops.empty() && !boundary_loops[0].empty()) {
            boundary_vertices.resize(boundary_loops[0].size());
            for(int i = 0; i < boundary_loops[0].size(); i++) {
                boundary_vertices(i) = boundary_loops[0][i];
            }
            std::cout << "Found " << boundary_vertices.size() << " boundary vertices" << std::endl;
        } else {
            // If no boundary, select a few vertices as fixed anchors
            boundary_vertices.resize(4);
            int n = V_working.rows();
            boundary_vertices(0) = 0;
            boundary_vertices(1) = n/4;
            boundary_vertices(2) = n/2;
            boundary_vertices(3) = 3*n/4;
            std::cout << "No boundary found, using " << boundary_vertices.size() << " anchor vertices" << std::endl;
        }
    };

    auto find_vertices_in_radius = [&](const Eigen::Vector3d& center, double radius) -> std::vector<int> {
        std::vector<int> vertices;
        for(int i = 0; i < V_working.rows(); i++) {
            if((V_working.row(i).transpose() - center).norm() < radius) {
                vertices.push_back(i);
            }
        }
        return vertices;
    };

    auto apply_biharmonic_deformation = [&](const Eigen::Vector3d& target_pos) {
        if(handle_vertices.empty()) return;

        // Compute displacement for handle
        Eigen::Vector3d displacement = target_pos - drag_start_pos;

        // Prepare boundary conditions
        // b contains indices of constrained vertices (boundary + handle)
        // bc contains the DEFORMATION FIELD (displacement) for these vertices

        int n_boundary = boundary_vertices.size();
        int n_handle = handle_vertices.size();

        Eigen::VectorXi b(n_boundary + n_handle);
        Eigen::MatrixXd d_bc(n_boundary + n_handle, 3);

        // Boundary vertices have zero displacement (stay fixed)
        for(int i = 0; i < n_boundary; i++) {
            b(i) = boundary_vertices(i);
            d_bc.row(i) = Eigen::RowVector3d(0, 0, 0);
        }

        // Handle vertices have the computed displacement
        for(int i = 0; i < n_handle; i++) {
            b(n_boundary + i) = handle_vertices[i];
            d_bc.row(n_boundary + i) = displacement.transpose();
        }

        // Solve for biharmonic deformation field
        Eigen::MatrixXd d;
        igl::harmonic(V_working, F, b, d_bc, 2, d);  // k=2 for biharmonic

        // Apply deformation field to get new positions
        Eigen::MatrixXd V_new = V_working + d;

        // Update mesh
        viewer.data_list[0].set_vertices(V_new);
        V_working = V_new;
    };

    auto update_nested_mesh = [&]() {
        if(V_base.size() == 0) return;

        nested_mesh = resize_mesh(V_base, slider_value);
        Eigen::RowVector3d translation(translate_x, translate_y, translate_z);
        nested_mesh = nested_mesh.rowwise() + translation;

        if(viewer.data_list.size() >= 3) {
            viewer.data_list[2].set_vertices(nested_mesh);

            bool intersects_outer = meshes_intersect(nested_mesh, F, V_outer, F);
            bool intersects_inner = meshes_intersect(nested_mesh, F, V_outer_inner_shell, F_inner_flipped);

            if(intersects_outer || intersects_inner) {
                viewer.data_list[2].set_colors(Eigen::RowVector3d(1.0, 0.0, 0.0));
            } else {
                viewer.data_list[2].set_colors(Eigen::RowVector3d(0.0, 1.0, 0.0));
            }
        }
    };

    auto shell_changed = [&]() {
        if(V_base.size() == 0) return;

        V_outer = V_base;
        V_outer_inner_shell = offset_along_normals(V_outer, F, -shell_thickness);

        F_inner_flipped = F;
        F_inner_flipped.col(1).swap(F_inner_flipped.col(2));

        if(viewer.data_list.size() >= 3) {
            viewer.data_list[0].set_vertices(V_outer);
            viewer.data_list[0].show_faces = false;

            viewer.data_list[1].set_mesh(V_outer_inner_shell, F_inner_flipped);
            viewer.data_list[1].show_faces = false;
        }
    };

    auto update_plane = [&]() {
        Eigen::RowVector3d p0(0.0, cutting_plane_y_coord, 0.0);
        Eigen::RowVector3d u(1.0, 0.0, 0.0);
        Eigen::RowVector3d v(0.0, 0.0, 1.0);

        V_plane.resize(4,3);
        V_plane.row(0) = p0 + plane_size*u + plane_size*v;
        V_plane.row(1) = p0 - plane_size*u + plane_size*v;
        V_plane.row(2) = p0 - plane_size*u - plane_size*v;
        V_plane.row(3) = p0 + plane_size*u - plane_size*v;

        F_plane.resize(2,3);
        F_plane << 0,1,2, 0,2,3;

        if(plane_mesh_id == -1){
            viewer.append_mesh();
            plane_mesh_id = viewer.data_list.size()-1;
            viewer.data_list[plane_mesh_id].set_colors(Eigen::RowVector3d(0.5,0.5,1.0));
            viewer.data_list[plane_mesh_id].show_faces = true;
            viewer.data_list[plane_mesh_id].show_lines = true;
        }

        viewer.data_list[plane_mesh_id].set_mesh(V_plane, F_plane);
    };

    // ============ MODE SETUP FUNCTIONS ============

    auto setup_edit_mode = [&]() {
        clear_all_meshes();
        current_mode = AppMode::EDIT_DEFORMATION;
        is_dragging = false;
        selected_vertex = -1;
        handle_vertices.clear();

        // Determine which mesh to edit
        if(start_from_original == 1) {
            V_working = V_original;
        } else {
            // Start from last edited version (V_deformed if it exists, otherwise original)
            V_working = (V_deformed.size() > 0) ? V_deformed : V_original;
        }

        // Find boundary vertices for constraints
        find_boundary_vertices();

        // Display single mesh for editing
        viewer.append_mesh();
        viewer.data_list[0].set_mesh(V_working, F);
        viewer.data_list[0].show_faces = true;
        viewer.data_list[0].show_lines = true;
        viewer.data_list[0].set_colors(Eigen::RowVector3d(0.9, 0.9, 0.9));
    };

    auto apply_deformation = [&]() {
        // Save the current working mesh as the deformed result
        V_deformed = V_working;
        std::cout << "Deformation applied and saved!" << std::endl;
    };

    auto setup_nested_doll_mode = [&]() {
        clear_all_meshes();
        current_mode = AppMode::NESTED_DOLL_DESIGN;

        // Use the deformed mesh as base (or original if no deformation)
        V_base = (V_deformed.size() > 0) ? V_deformed : V_original;

        // Reset translation
        translate_x = translate_y = translate_z = 0.0f;

        // Create 3 meshes: outer shell, inner shell, nested mesh
        viewer.append_mesh();
        viewer.append_mesh();
        viewer.append_mesh();

        viewer.data_list[0].set_mesh(V_base, F);
        viewer.data_list[1].set_mesh(V_base, F);
        viewer.data_list[2].set_mesh(V_base, F);

        shell_changed();
        update_nested_mesh();
        update_plane();
    };

    auto visualize_shell_split = [&]() {
        if(V_base.size() == 0) return;

        if(plane_mesh_id != -1 && plane_mesh_id < viewer.data_list.size()) {
            viewer.erase_mesh(plane_mesh_id);
            plane_mesh_id = -1;
        }

        auto create_closed_shell_half = [](
                                            const Eigen::MatrixXd& outer_V, const Eigen::MatrixXi& outer_F,
                                            const Eigen::MatrixXd& inner_V, const Eigen::MatrixXi& inner_F,
                                            Eigen::MatrixXd& shell_V, Eigen::MatrixXi& shell_F) {

            int n_outer = outer_V.rows();
            int n_inner = inner_V.rows();

            shell_V.resize(n_outer + n_inner, 3);
            shell_V.block(0, 0, n_outer, 3) = outer_V;
            shell_V.block(n_outer, 0, n_inner, 3) = inner_V;

            int n_outer_faces = outer_F.rows();
            int n_inner_faces = inner_F.rows();

            shell_F.resize(n_outer_faces + n_inner_faces, 3);
            shell_F.block(0, 0, n_outer_faces, 3) = outer_F;

            for(int i = 0; i < n_inner_faces; i++) {
                shell_F.row(n_outer_faces + i) = inner_F.row(i).array() + n_outer;
            }
        };

        Eigen::MatrixXd outer_top_V, outer_bottom_V;
        Eigen::MatrixXi outer_top_F, outer_bottom_F;
        bool outer_success = split_mesh_open(V_outer, F, outer_top_V, outer_top_F,
                                             outer_bottom_V, outer_bottom_F, cutting_plane_y_coord);

        Eigen::MatrixXd inner_shell_top_V, inner_shell_bottom_V;
        Eigen::MatrixXi inner_shell_top_F, inner_shell_bottom_F;
        bool inner_success = split_mesh_open(V_outer_inner_shell, F,
                                             inner_shell_top_V, inner_shell_top_F,
                                             inner_shell_bottom_V, inner_shell_bottom_F, cutting_plane_y_coord);

        if(!outer_success || !inner_success) {
            std::cerr << "Failed to split shell meshes" << std::endl;
            return;
        }

        inner_shell_top_F.col(1).swap(inner_shell_top_F.col(2));
        inner_shell_bottom_F.col(1).swap(inner_shell_bottom_F.col(2));

        Eigen::MatrixXd top_shell_V, bottom_shell_V;
        Eigen::MatrixXi top_shell_F, bottom_shell_F;

        create_closed_shell_half(outer_top_V, outer_top_F,
                                 inner_shell_top_V, inner_shell_top_F,
                                 top_shell_V, top_shell_F);

        create_closed_shell_half(outer_bottom_V, outer_bottom_F,
                                 inner_shell_bottom_V, inner_shell_bottom_F,
                                 bottom_shell_V, bottom_shell_F);

        Eigen::RowVector3d translationUp(0, 0.15, 0);
        Eigen::RowVector3d translationDown(0, -0.15, 0);

        top_shell_V = top_shell_V.rowwise() + translationUp;
        bottom_shell_V = bottom_shell_V.rowwise() + translationDown;

        viewer.data_list[0].is_visible = false;
        viewer.data_list[1].is_visible = false;

        viewer.append_mesh();
        viewer.data_list[3].set_mesh(top_shell_V, top_shell_F);
        viewer.data_list[3].show_faces = true;
        viewer.data_list[3].show_lines = true;
        viewer.data_list[3].set_colors(Eigen::RowVector3d(0.7, 0.7, 0.9));
        viewer.data_list[3].face_based = false;

        viewer.append_mesh();
        viewer.data_list[4].set_mesh(bottom_shell_V, bottom_shell_F);
        viewer.data_list[4].show_faces = true;
        viewer.data_list[4].show_lines = true;
        viewer.data_list[4].set_colors(Eigen::RowVector3d(0.7, 0.7, 0.9));
        viewer.data_list[4].face_based = false;

        bool intersects_outer = meshes_intersect(nested_mesh, F, V_outer, F);
        bool intersects_inner = meshes_intersect(nested_mesh, F, V_outer_inner_shell, F_inner_flipped);

        if(intersects_outer || intersects_inner) {
            viewer.data_list[2].set_colors(Eigen::RowVector3d(1.0, 0.0, 0.0));
        } else {
            viewer.data_list[2].set_colors(Eigen::RowVector3d(0.0, 1.0, 0.0));
        }

        viewer.data_list[2].show_lines = true;
        viewer.data_list[2].face_based = false;
    };

    auto visualize_swept_volume = [&]() {
        if(V_base.size() == 0) return;

        Eigen::MatrixXd SV_top, SV_bottom;
        Eigen::MatrixXi SF_top, SF_bottom;

        const auto & transform_top = [](const double t)->Eigen::Affine3d {
            Eigen::Affine3d T = Eigen::Affine3d::Identity();
            T.translate(Eigen::Vector3d(-0.5,0,0));
            T.rotate(Eigen::AngleAxisd(t*0.5*igl::PI,Eigen::Vector3d(0,0,1)));
            T.translate(Eigen::Vector3d(0.5,0,0));
            return T;
        };

        const auto & transform_bottom = [](const double t)->Eigen::Affine3d {
            Eigen::Affine3d T = Eigen::Affine3d::Identity();
            T.translate(Eigen::Vector3d(0.5,0,0));
            T.rotate(Eigen::AngleAxisd(t*0.5*igl::PI,Eigen::Vector3d(0,0,1)));
            T.translate(Eigen::Vector3d(-0.5,0,0));
            return T;
        };

        if (viewer.data_list.size() < 5) {
            std::cerr << "Must visualize shell split first!" << std::endl;
            return;
        }

        igl::swept_volume(viewer.data_list[3].V, viewer.data_list[3].F,
                          transform_top, time_steps, grid_size, isolevel, SV_top, SF_top);

        igl::swept_volume(viewer.data_list[4].V, viewer.data_list[4].F,
                          transform_bottom, time_steps, grid_size, isolevel, SV_bottom, SF_bottom);

        Eigen::RowVector3d translationUp(0, 0.5, 0);
        Eigen::RowVector3d translationDown(0, -0.5, 0);

        SV_top = SV_top.rowwise() + translationUp;
        SV_bottom = SV_bottom.rowwise() + translationDown;

        viewer.append_mesh();
        viewer.data_list[viewer.data_list.size() - 1].set_mesh(SV_top, SF_top);
        viewer.data_list[viewer.data_list.size() - 1].show_faces = true;
        viewer.data_list[viewer.data_list.size() - 1].show_lines = true;
        viewer.data_list[viewer.data_list.size() - 1].set_colors(Eigen::RowVector3d(0.5, 0.8, 0.8));

        viewer.append_mesh();
        viewer.data_list[viewer.data_list.size() - 1].set_mesh(SV_bottom, SF_bottom);
        viewer.data_list[viewer.data_list.size() - 1].show_faces = true;
        viewer.data_list[viewer.data_list.size() - 1].show_lines = true;
        viewer.data_list[viewer.data_list.size() - 1].set_colors(Eigen::RowVector3d(0.5, 0.8, 0.8));
    };

    auto load_mesh = [&](const std::string& filename) {
        bool success = false;

        if (filename.find(".obj") != std::string::npos) {
            success = igl::readOBJ(filename, V_original, F_original);
        } else if (filename.find(".off") != std::string::npos) {
            success = igl::readOFF(filename, V_original, F_original);
        } else if (filename.find(".ply") != std::string::npos) {
            success = igl::readPLY(filename, V_original, F_original);
        }

        if (success) {
            F = F_original;
            V_working = V_original;
            V_deformed.resize(0, 0);  // Clear any previous deformation

            auto min_point = V_original.colwise().minCoeff();
            auto max_point = V_original.colwise().maxCoeff();
            auto bbox_size = (max_point - min_point).norm();
            selection_radius = bbox_size * 0.05f;  // 5% of bbox size

            setup_edit_mode();
        } else {
            std::cerr << "Failed to load: " << filename << std::endl;
        }

        return success;
    };

    // ============ UI CALLBACKS ============

    std::function<void()> ui_edit_mode, ui_nested_doll_mode;

    ui_edit_mode = [&]() {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.8f, 0.2f, 1.0f));
        ImGui::Text("=== EDIT DEFORMATION MODE ===");
        ImGui::PopStyleColor();
        ImGui::Spacing();

        if (ImGui::Button("Load mesh...", ImVec2(-1, 0))) {
            auto selection = pfd::open_file("Select a mesh file", ".",
                                            { "Mesh Files", "*.obj *.off *.ply",
                                             "All Files", "*" }).result();
            if (!selection.empty()) {
                load_mesh(selection[0]);
            }
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::Text("Deformation Source:");
        ImGui::RadioButton("Start from last edited version", &start_from_original, 0);
        ImGui::RadioButton("Start from original mesh", &start_from_original, 1);

        if (ImGui::Button("Reset to selected version", ImVec2(-1, 0))) {
            setup_edit_mode();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::SliderFloat("Selection Radius", &selection_radius, 0.01f, 0.5f)) {
            // Radius updated for next selection
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::TextWrapped("Biharmonic Deformation Instructions:");
        ImGui::BulletText("Click and drag on the mesh surface");
        ImGui::BulletText("Vertices within radius will be handles");
        ImGui::BulletText("Smooth, detail-preserving deformation");
        ImGui::BulletText("Boundary vertices stay fixed");
        ImGui::BulletText("Adjust 'Selection Radius' for area size");

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::Button("Apply Deformation & Continue", ImVec2(-1, 0))) {
            apply_deformation();
            setup_nested_doll_mode();
            menu.callback_draw_viewer_menu = ui_nested_doll_mode;
        }
    };

    ui_nested_doll_mode = [&]() {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.2f, 0.8f, 1.0f, 1.0f));
        ImGui::Text("=== NESTED DOLL DESIGN MODE ===");
        ImGui::PopStyleColor();
        ImGui::Spacing();

        if (ImGui::Button("Return to Edit Deformation", ImVec2(-1, 0))) {
            setup_edit_mode();
            menu.callback_draw_viewer_menu = ui_edit_mode;
            return;
        }

        ImGui::Spacing();

        if (ImGui::Button("Reset to Original Mesh", ImVec2(-1, 0))) {
            // Clear all deformations and customizations
            V_deformed.resize(0, 0);
            V_working = V_original;
            translate_x = translate_y = translate_z = 0.0f;
            slider_value = 1.0f;
            shell_thickness = 0.001f;
            cutting_plane_y_coord = 0.5f;

            // Return to edit mode with original mesh
            setup_edit_mode();
            menu.callback_draw_viewer_menu = ui_edit_mode;

            std::cout << "Reset to original mesh - all changes discarded" << std::endl;
            return;
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::SliderFloat("Nested Scale", &slider_value, 0.0f, 1.0f)) {
            update_nested_mesh();
        }

        if (ImGui::SliderFloat("Shell Thickness", &shell_thickness, 0.001f, 0.1f)) {
            shell_changed();
        }

        if(ImGui::SliderInt("Number of Nests", &n_nests, 1, 10)) {
            // Future implementation
        }

        if(ImGui::SliderFloat("Cutting Plane Y", &cutting_plane_y_coord, -1, 1)) {
            update_plane();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::Text("Translate Inner Mesh");

        bool translation_changed = false;
        if (ImGui::SliderFloat("Translate X", &translate_x, -1.0f, 1.0f)) {
            translation_changed = true;
        }
        if (ImGui::SliderFloat("Translate Y", &translate_y, -1.0f, 1.0f)) {
            translation_changed = true;
        }
        if (ImGui::SliderFloat("Translate Z", &translate_z, -1.0f, 1.0f)) {
            translation_changed = true;
        }

        if (translation_changed) {
            update_nested_mesh();
        }

        if (ImGui::Button("Reset Translation", ImVec2(-1, 0))) {
            translate_x = translate_y = translate_z = 0.0f;
            update_nested_mesh();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::Button("Visualize Shell Split", ImVec2(-1, 0))) {
            if(V_base.size() == 0) return;
            visualize_shell_split();
        }

        if (ImGui::Button("Visualize Swept Volume", ImVec2(-1, 0))) {
            if(V_base.size() == 0) return;
            visualize_swept_volume();
        }

        if (ImGui::Button("Generate dolls...", ImVec2(-1, 0))) {
            if(V_base.size() == 0) return;
            std::cout << "Generate dolls - Coming soon!" << std::endl;
        }
    };

    // ============ MOUSE CALLBACKS FOR BIHARMONIC DEFORMATION ============

    viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
        // Only allow deformation in EDIT mode
        if(current_mode != AppMode::EDIT_DEFORMATION) return false;

        Eigen::Vector3f bc_hit;
        int fid;
        auto x = viewer.current_mouse_x;
        auto y = viewer.core().viewport(3) - static_cast<float>(viewer.current_mouse_y);

        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y),
                                     viewer.core().view,
                                     viewer.core().proj,
                                     viewer.core().viewport,
                                     V_working,
                                     F,
                                     fid,
                                     bc_hit)) {
            // Get the 3D position of the click
            drag_start_pos = igl::unproject(Eigen::Vector3f(x, y, viewer.down_mouse_z),
                                            viewer.core().view,
                                            viewer.core().proj,
                                            viewer.core().viewport)
                                 .template cast<double>();

            // Find all vertices within selection radius
            handle_vertices = find_vertices_in_radius(drag_start_pos, selection_radius);

            if(!handle_vertices.empty()) {
                is_dragging = true;
                std::cout << "Selected " << handle_vertices.size() << " handle vertices" << std::endl;
                return true;
            }
        }
        return false;
    };

    viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
        if(current_mode != AppMode::EDIT_DEFORMATION || !is_dragging) return false;

        // Get current mouse position in 3D
        Eigen::Vector3d current_pos = igl::unproject(
                                          Eigen::Vector3f(viewer.current_mouse_x,
                                                          viewer.core().viewport[3] - static_cast<float>(viewer.current_mouse_y),
                                                          viewer.down_mouse_z),
                                          viewer.core().view,
                                          viewer.core().proj,
                                          viewer.core().viewport)
                                          .template cast<double>();

        // Apply biharmonic deformation
        apply_biharmonic_deformation(current_pos);

        return true;
    };

    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
        if(current_mode != AppMode::EDIT_DEFORMATION) return false;

        if (is_dragging) {
            is_dragging = false;
            handle_vertices.clear();
            std::cout << "Deformation complete" << std::endl;
            return true;
        }

        return false;
    };

    // ============ INITIALIZATION ============

    menu.callback_draw_viewer_menu = ui_edit_mode;
    viewer.launch();

    return 0;
}
