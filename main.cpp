#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <imgui.h>
#include "portable-file-dialogs.h"
#include <igl/kelvinlets.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/swept_volume.h>
#include <iostream>
#include "manip.h"

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V, V_base;  // V = current outer mesh, V_base = base template for nested mesh
    Eigen::MatrixXi F;

    Eigen::MatrixXd V_outer, V_outer_inner_shell;
    Eigen::MatrixXi F_inner_flipped;  // Flipped faces for inner shell - used consistently

    // Single nested mesh variable that persists across operations
    Eigen::MatrixXd nested_mesh;

    std::vector<Eigen::MatrixXd> V_dolls;
    std::vector<Eigen::MatrixXi> F_dolls;

    float slider_value = 1.f;
    float shell_thickness = 0.001f;
    int n_nests = 1;
    int mesh_selection = 0;
    // inner mesh
    float translate_x = 0.0f;
    float translate_y = 0.0f;
    float translate_z = 0.0f;

    const int grid_size = 75;
    const int time_steps = 200;
    const double isolevel = 0;

    Eigen::Vector3d posStart(0, 0, 0);
    Eigen::Vector3d posEnd;
    auto brush_strength = 1.;
    decltype(V) result;
    Eigen::Matrix3d mat;
    mat.setZero();

    auto brushRadius = .1;
    auto brushType = igl::BrushType::GRAB;
    auto scale = .05;

    float cutting_plane_y_coord = 0.5f;
    Eigen::MatrixXd V_plane;
    Eigen::MatrixXi F_plane;
    int plane_mesh_id = -1;
    double plane_size = 2.0;

    igl::opengl::glfw::Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    auto clear_all_meshes = [&]() {
        while(viewer.erase_mesh(viewer.selected_data_index)){};
        viewer.data().clear();

        plane_mesh_id = -1;
        cutting_plane_y_coord = 0.5f;
    };

    // Update nested_mesh based on scale and translation - THIS IS THE ONLY FUNCTION THAT MODIFIES nested_mesh
    auto update_nested_mesh = [&]() {
        if(V_base.size() == 0) return;

        // Scale based on V_base state
        nested_mesh = resize_mesh(V_base, slider_value);

        // Apply translation to nested mesh
        Eigen::RowVector3d translation(translate_x, translate_y, translate_z);
        nested_mesh = nested_mesh.rowwise() + translation;

        // Update mesh 2 (nested mesh)
        if(viewer.data_list.size() >= 3) {
            viewer.data_list[2].set_vertices(nested_mesh);

            // Check if nested mesh intersects with the outer shell surface
            bool intersects_outer = meshes_intersect(nested_mesh, F, V_outer, F);

            // Check if nested mesh intersects with the inner shell surface (using flipped faces)
            bool intersects_inner = meshes_intersect(nested_mesh, F, V_outer_inner_shell, F_inner_flipped);

            // The nested mesh should NOT intersect with outer surface
            // and SHOULD be completely inside (not intersecting inner surface either)
            // If it intersects either surface, it's in the shell wall - turn RED
            if(intersects_outer || intersects_inner) {
                viewer.data_list[2].set_colors(Eigen::RowVector3d(1.0, 0.0, 0.0));
            } else {
                viewer.data_list[2].set_colors(Eigen::RowVector3d(0.0, 1.0, 0.0));
            }
        }
    };

    auto shell_changed = [&]() {
        if(V.size() == 0) return;

        V_outer = V;
        V_outer_inner_shell = offset_along_normals(V_outer, F, -shell_thickness);

        // Create flipped faces for inner shell (used consistently for both visualization and collision)
        F_inner_flipped = F;
        F_inner_flipped.col(1).swap(F_inner_flipped.col(2));

        // Update visualization if we have 3 meshes
        if(viewer.data_list.size() >= 3) {
            // Update mesh 0 - outer shell surface
            viewer.data_list[0].set_vertices(V_outer);
            viewer.data_list[0].show_faces = false;

            // Update mesh 1 - inner shell surface (with flipped faces)
            viewer.data_list[1].set_mesh(V_outer_inner_shell, F_inner_flipped);  // ← CHANGED!
            viewer.data_list[1].show_faces = false;

        }
    };

    auto update_plane = [&]() {
        Eigen::RowVector3d p0(0.0, cutting_plane_y_coord, 0.0);
        Eigen::RowVector3d n(0.0, 1.0, 0.0);

        // Tangent vectors
        Eigen::RowVector3d u(1.0, 0.0, 0.0);
        Eigen::RowVector3d v(0.0, 0.0, 1.0);

        // Plane vertices
        V_plane.resize(4,3);
        V_plane.row(0) = p0 + plane_size*u + plane_size*v;
        V_plane.row(1) = p0 - plane_size*u + plane_size*v;
        V_plane.row(2) = p0 - plane_size*u - plane_size*v;
        V_plane.row(3) = p0 + plane_size*u - plane_size*v;

        // Plane faces
        F_plane.resize(2,3);
        F_plane << 0,1,2,
            0,2,3;

        if(plane_mesh_id == -1){
            viewer.append_mesh();
            plane_mesh_id = viewer.data_list.size()-1;
            viewer.data_list[plane_mesh_id].set_colors(Eigen::RowVector3d(0.5,0.5,1.0));
            viewer.data_list[plane_mesh_id].show_faces = true;
            viewer.data_list[plane_mesh_id].show_lines = true;
        }

        viewer.data_list[plane_mesh_id].set_mesh(V_plane, F_plane);
    };

    auto stage_1 = [&]() {
        clear_all_meshes();

        // Create 3 meshes:
        // Index 0: Outer shell surface (V_outer) - wireframe only
        // Index 1: Inner shell surface (V_outer_inner_shell with flipped normals) - wireframe only
        // Index 2: Nested mesh - visible with faces
        viewer.append_mesh();
        viewer.append_mesh();
        viewer.append_mesh();

        viewer.data_list[0].set_mesh(V, F);
        // viewer.data_list[0].show_faces = false;  // No faces, wireframe only
        // viewer.data_list[0].show_lines = true;
        // viewer.data_list[0].set_colors(Eigen::RowVector3d(0.8, 0.8, 0.9));

        viewer.data_list[1].set_mesh(V, F);
        // viewer.data_list[1].show_faces = false;  // No faces, wireframe only
        // viewer.data_list[1].show_lines = true;
        // viewer.data_list[1].set_colors(Eigen::RowVector3d(0.9, 0.8, 0.8));

        viewer.data_list[2].set_mesh(V, F);
        // viewer.data_list[2].show_faces = true;   // Nested mesh shows faces
        // viewer.data_list[2].show_lines = true;

        shell_changed();
        update_nested_mesh();  // Initialize nested_mesh
    };

    stage_1();

    auto visualize_shell_split = [&]() {
        if(V.size() == 0) return;

        // Permanently remove the plane mesh
        if(plane_mesh_id != -1 && plane_mesh_id < viewer.data_list.size()) {
            viewer.erase_mesh(plane_mesh_id);
            plane_mesh_id = -1;
        }

        // Helper function to combine outer and inner shell surfaces
        auto create_closed_shell_half = [](
                                            const Eigen::MatrixXd& outer_V, const Eigen::MatrixXi& outer_F,
                                            const Eigen::MatrixXd& inner_V, const Eigen::MatrixXi& inner_F,
                                            Eigen::MatrixXd& shell_V, Eigen::MatrixXi& shell_F) {

            // Combine outer and inner surfaces
            int n_outer = outer_V.rows();
            int n_inner = inner_V.rows();

            shell_V.resize(n_outer + n_inner, 3);
            shell_V.block(0, 0, n_outer, 3) = outer_V;
            shell_V.block(n_outer, 0, n_inner, 3) = inner_V;

            int n_outer_faces = outer_F.rows();
            int n_inner_faces = inner_F.rows();

            shell_F.resize(n_outer_faces + n_inner_faces, 3);
            shell_F.block(0, 0, n_outer_faces, 3) = outer_F;

            // Offset inner face indices
            for(int i = 0; i < n_inner_faces; i++) {
                shell_F.row(n_outer_faces + i) = inner_F.row(i).array() + n_outer;
            }
        };

        // Split outer shell surface
        Eigen::MatrixXd outer_top_V, outer_bottom_V;
        Eigen::MatrixXi outer_top_F, outer_bottom_F;

        bool outer_success = split_mesh_open(V_outer, F, outer_top_V, outer_top_F,
                                             outer_bottom_V, outer_bottom_F, cutting_plane_y_coord);

        // Split inner shell surface with original F
        Eigen::MatrixXd inner_shell_top_V, inner_shell_bottom_V;
        Eigen::MatrixXi inner_shell_top_F, inner_shell_bottom_F;

        bool inner_success = split_mesh_open(V_outer_inner_shell, F,
                                             inner_shell_top_V, inner_shell_top_F,
                                             inner_shell_bottom_V, inner_shell_bottom_F, cutting_plane_y_coord);

        if(!outer_success || !inner_success) {
            std::cerr << "Failed to split shell meshes" << std::endl;
            return;
        }

        // Flip normals AFTER splitting
        inner_shell_top_F.col(1).swap(inner_shell_top_F.col(2));
        inner_shell_bottom_F.col(1).swap(inner_shell_bottom_F.col(2));

        // Create combined shell halves using the helper function
        Eigen::MatrixXd top_shell_V, bottom_shell_V;
        Eigen::MatrixXi top_shell_F, bottom_shell_F;

        create_closed_shell_half(outer_top_V, outer_top_F,
                                 inner_shell_top_V, inner_shell_top_F,
                                 top_shell_V, top_shell_F);

        create_closed_shell_half(outer_bottom_V, outer_bottom_F,
                                 inner_shell_bottom_V, inner_shell_bottom_F,
                                 bottom_shell_V, bottom_shell_F);

        // Apply translations to separate the halves
        Eigen::RowVector3d translationUp(0, 0.15, 0);
        Eigen::RowVector3d translationDown(0, -0.15, 0);

        top_shell_V = top_shell_V.rowwise() + translationUp;
        bottom_shell_V = bottom_shell_V.rowwise() + translationDown;

        // Hide meshes 0 and 1 (outer and inner shell surfaces)
        viewer.data_list[0].is_visible = false;
        viewer.data_list[1].is_visible = false;

        // Add combined top shell - index 3
        viewer.append_mesh();
        viewer.data_list[3].set_mesh(top_shell_V, top_shell_F);
        viewer.data_list[3].show_faces = true;
        viewer.data_list[3].show_lines = true;
        viewer.data_list[3].set_colors(Eigen::RowVector3d(0.7, 0.7, 0.9));
        viewer.data_list[3].face_based = false;

        // Add combined bottom shell - index 4
        viewer.append_mesh();
        viewer.data_list[4].set_mesh(bottom_shell_V, bottom_shell_F);
        viewer.data_list[4].show_faces = true;
        viewer.data_list[4].show_lines = true;
        viewer.data_list[4].set_colors(Eigen::RowVector3d(0.7, 0.7, 0.9));
        viewer.data_list[4].face_based = false;

        // Nested mesh stays at index 2, update collision color
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
        if(V.size() == 0) return;

        Eigen::MatrixXd SV_top, SV_bottom;
        Eigen::MatrixXi SF_top, SF_bottom;

        const auto & transform_top = [](const double t)->Eigen::Affine3d
        {
            Eigen::Affine3d T = Eigen::Affine3d::Identity();
            T.translate(Eigen::Vector3d(-0.5,0,0));
            T.rotate(Eigen::AngleAxisd(t*0.5*igl::PI,Eigen::Vector3d(0,0,1)));
            T.translate(Eigen::Vector3d(0.5,0,0));
            return T;
        };

        const auto & transform_bottom = [](const double t)->Eigen::Affine3d
        {
            Eigen::Affine3d T = Eigen::Affine3d::Identity();
            T.translate(Eigen::Vector3d(0.5,0,0));
            T.rotate(Eigen::AngleAxisd(t*0.5*igl::PI,Eigen::Vector3d(0,0,1)));
            T.translate(Eigen::Vector3d(-0.5,0,0));
            return T;
        };

        // Changed from 7 to 5 (need indices 0,1,2,3,4)
        if (viewer.data_list.size() < 5) {
            std::cerr << "Must visualize shell split first!" << std::endl;
            return;
        }

        // Use index 3 (combined top shell)
        igl::swept_volume(
            viewer.data_list[3].V, viewer.data_list[3].F,
            transform_top, time_steps, grid_size, isolevel, SV_top, SF_top);

        // Use index 4 (combined bottom shell)
        igl::swept_volume(
            viewer.data_list[4].V, viewer.data_list[4].F,
            transform_bottom, time_steps, grid_size, isolevel, SV_bottom, SF_bottom);

        Eigen::RowVector3d translationUp(0, 0.5, 0);
        Eigen::RowVector3d translationDown(0, -0.5, 0);

        SV_top = SV_top.rowwise() + translationUp;
        SV_bottom = SV_bottom.rowwise() + translationDown;

        // Add swept volume meshes at indices 5 and 6
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
            success = igl::readOBJ(filename, V, F);
        } else if (filename.find(".off") != std::string::npos) {
            success = igl::readOFF(filename, V, F);
        } else if (filename.find(".ply") != std::string::npos) {
            success = igl::readPLY(filename, V, F);
        }

        if (success) {
            translate_x = translate_y = translate_z = 0.0f;  // Reset translations
            stage_1();

            V_base = V;  // Set base template for nested mesh
            auto min_point = V_base.colwise().minCoeff();
            auto max_point = V_base.colwise().maxCoeff();
            // to multiply brush force proportional to size of mesh
            brush_strength = (max_point - min_point).norm() * 2;
            update_plane();
        } else {
            std::cerr << "Failed to load: " << filename << std::endl;
        }

        return success;
    };

    std::function<void()> ui_1;

    ui_1 = [&]()
    {
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

        if (ImGui::SliderFloat("Nested Scale", &slider_value, 0.0f, 1.0f)) {
            update_nested_mesh();  // Only function that modifies nested_mesh
        }

        if (ImGui::SliderFloat("Shell Thickness", &shell_thickness, 0.001f, 0.1f)) {
            shell_changed();  // Does NOT affect nested_mesh
        }

        if(ImGui::SliderInt("Number of Nests", &n_nests, 1, 10)) {
            // Does NOT affect nested_mesh
        }

        if(ImGui::SliderFloat("Cutting Plane Y Translation", &cutting_plane_y_coord, -1, 1)) {
            update_plane();  // Does NOT affect nested_mesh
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        // Translation controls for inner mesh
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
            update_nested_mesh();  // Only function that modifies nested_mesh
        }

        if (ImGui::Button("Reset Translation", ImVec2(-1, 0))) {
            translate_x = translate_y = translate_z = 0.0f;
            update_nested_mesh();  // Only function that modifies nested_mesh
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::Button("Visualize Shell Split", ImVec2(-1, 0))) {
            if(V.size() == 0) return;
            visualize_shell_split();  // Uses nested_mesh but does NOT modify it
        }

        if (ImGui::Button("Visualize Swept Volume", ImVec2(-1, 0))) {
            if(V.size() == 0) return;
            visualize_swept_volume();  // Does NOT affect nested_mesh
        }

        if (ImGui::Button("Generate dolls...", ImVec2(-1, 0))) {
            if(V.size() == 0) return;
            // TODO: Implement doll generation logic
            std::cout << "Generate dolls - Coming soon!" << std::endl;
        }
    };

    viewer.callback_mouse_down =
        [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
        Eigen::Vector3f bc;
        int fid;
        auto x = viewer.current_mouse_x;
        auto y =
            viewer.core().viewport(3) - static_cast<float>(viewer.current_mouse_y);

        // Use V_base for mouse interaction base, or nested_mesh if available
        auto V_temp = V_base;
        if (nested_mesh.size() != 0) {
            V_temp = nested_mesh;
        }
        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y),
                                     viewer.core().view,
                                     viewer.core().proj,
                                     viewer.core().viewport,
                                     V_temp,
                                     F,
                                     fid,
                                     bc)) {
            posStart = igl::unproject(Eigen::Vector3f(x, y, viewer.down_mouse_z),
                                      viewer.core().view,
                                      viewer.core().proj,
                                      viewer.core().viewport)
                           .template cast<double>();
            return true;
        }
        return false;
    };

    viewer.callback_mouse_up =
        [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
        if (!posStart.isZero() && !posStart.hasNaN()) {
            Eigen::Vector3d oldPosEnd(posEnd.x(), posEnd.y(), posEnd.z());

            posEnd = igl::unproject(
                         Eigen::Vector3f(viewer.current_mouse_x,
                                         viewer.core().viewport[3] -
                                             static_cast<float>(viewer.current_mouse_y),
                                         viewer.down_mouse_z),
                         viewer.core().view,
                         viewer.core().proj,
                         viewer.core().viewport)
                         .template cast<double>();

            // exaggerate the force by a little bit
            Eigen::Vector3d forceVec = (posEnd - posStart) * brush_strength;

            int scaleFactor = forceVec.norm();
            if (posEnd.x() < posStart.x()) {
                // probably not the best way to determine direction.
                scaleFactor = -scaleFactor;
            }

            igl::kelvinlets(
                V_base,
                posStart,
                forceVec,
                mat,
                igl::KelvinletParams<double>(brushRadius, scale, brushType),
                result);

            viewer.data_list[2].set_vertices(result);
            V_base = result;
            posStart.setZero();

            update_nested_mesh();  // Update nested_mesh after brush deformation

            return true;
        }

        return false;
    };

    menu.callback_draw_viewer_menu = ui_1;

    viewer.launch();

    return 0;
}
