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
#include <iostream>
#include "manip.h"

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V, V1;
    Eigen::MatrixXi F;

    Eigen::MatrixXd V_outer, V_outer_inner_shell, V_nested;

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
    //update mesh_resized by add translation
    auto mesh_resized = [&]() {
        if(V1.size() == 0) return;
        V_nested = resize_mesh(V1, slider_value);

        // Apply translation to inner mesh
        Eigen::RowVector3d translation(translate_x, translate_y, translate_z);
        V_nested = V_nested.rowwise() + translation;

        viewer.data_list[1].set_vertices(V_nested);

        // Check if nested mesh intersects with the outer shell surface
        bool intersects_outer = meshes_intersect(V_nested, F, V_outer, F);

        // Check if nested mesh intersects with the inner shell surface
        bool intersects_inner = meshes_intersect(V_nested, F, V_outer_inner_shell, F);

        // The nested mesh should NOT intersect with outer surface
        // and SHOULD be completely inside (not intersecting inner surface either)
        // If it intersects either surface, it's in the shell wall - turn RED
        if(intersects_outer || intersects_inner) {
            viewer.data_list[1].set_colors(Eigen::RowVector3d(1.0, 0.0, 0.0));
        } else {
            viewer.data_list[1].set_colors(Eigen::RowVector3d(0.0, 1.0, 0.0));
        }
    };

    auto shell_changed = [&]() {
        if(V.size() == 0) return;

        V_outer = V;
        V_outer_inner_shell = offset_along_normals(V_outer, F, -shell_thickness);

        // Only update if we're in stage 1 (with 3 meshes)
        if(viewer.data_list.size() >= 3) {
            viewer.data_list[0].set_vertices(V_outer);
            viewer.data_list[0].show_faces = false;
            viewer.data_list[2].set_vertices(V_outer_inner_shell);
            viewer.data_list[2].show_faces = false;
        }
        mesh_resized();
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

        viewer.append_mesh();
        viewer.append_mesh();
        viewer.data_list[0].set_mesh(V, F);
        viewer.data_list[1].set_mesh(V, F);
        viewer.data_list[2].set_mesh(V, F);
        shell_changed();
        mesh_resized();
    };

    stage_1();

    auto visualize_shell_split = [&]() {
        if(V.size() == 0) return;

        // Split outer shell surface (open cut)
        Eigen::MatrixXd outer_top_V, outer_bottom_V;
        Eigen::MatrixXi outer_top_F, outer_bottom_F;

        bool outer_success = split_mesh_open(V_outer, F, outer_top_V, outer_top_F,
                                             outer_bottom_V, outer_bottom_F, cutting_plane_y_coord);

        // Flip normals for inner shell surface (swap columns 1 and 2 of face matrix)
        Eigen::MatrixXi F_inner_flipped = F;
        F_inner_flipped.col(1).swap(F_inner_flipped.col(2));

        // Split inner shell surface with flipped normals (open cut)
        Eigen::MatrixXd inner_shell_top_V, inner_shell_bottom_V;
        Eigen::MatrixXi inner_shell_top_F, inner_shell_bottom_F;

        bool inner_success = split_mesh_open(V_outer_inner_shell, F_inner_flipped,
                                             inner_shell_top_V, inner_shell_top_F,
                                             inner_shell_bottom_V, inner_shell_bottom_F, cutting_plane_y_coord);

        if(!outer_success || !inner_success) {
            std::cerr << "Failed to split shell meshes" << std::endl;
            return;
        }

        // Apply translations to separate the halves
        Eigen::RowVector3d translationUp(0, 0.15, 0);
        Eigen::RowVector3d translationDown(0, -0.15, 0);

        outer_top_V = outer_top_V.rowwise() + translationUp;
        inner_shell_top_V = inner_shell_top_V.rowwise() + translationUp;
        outer_bottom_V = outer_bottom_V.rowwise() + translationDown;
        inner_shell_bottom_V = inner_shell_bottom_V.rowwise() + translationDown;

        clear_all_meshes();

        // Add top half of shell (outer surface)
        viewer.append_mesh();
        viewer.data_list[0].set_mesh(outer_top_V, outer_top_F);
        viewer.data_list[0].show_faces = true;
        viewer.data_list[0].show_lines = true;
        viewer.data_list[0].set_colors(Eigen::RowVector3d(0.8, 0.8, 0.9));

        // Add top half of shell (inner surface - flipped normals)
        viewer.append_mesh();
        viewer.data_list[1].set_mesh(inner_shell_top_V, inner_shell_top_F);
        viewer.data_list[1].show_faces = true;
        viewer.data_list[1].show_lines = true;
        viewer.data_list[1].set_colors(Eigen::RowVector3d(0.9, 0.8, 0.8));

        // Add bottom half of shell (outer surface)
        viewer.append_mesh();
        viewer.data_list[2].set_mesh(outer_bottom_V, outer_bottom_F);
        viewer.data_list[2].show_faces = true;
        viewer.data_list[2].show_lines = true;
        viewer.data_list[2].set_colors(Eigen::RowVector3d(0.8, 0.8, 0.9));

        // Add bottom half of shell (inner surface - flipped normals)
        viewer.append_mesh();
        viewer.data_list[3].set_mesh(inner_shell_bottom_V, inner_shell_bottom_F);
        viewer.data_list[3].show_faces = true;
        viewer.data_list[3].show_lines = true;
        viewer.data_list[3].set_colors(Eigen::RowVector3d(0.9, 0.8, 0.8));

        // Add the nested doll mesh (V_nested) - keep it whole and centered
        if(V_nested.size() > 0) {
            viewer.append_mesh();
            viewer.data_list[4].set_mesh(V_nested, F);
            viewer.data_list[4].show_faces = true;
            viewer.data_list[4].show_lines = true;

            // Check if nested mesh intersects with shell
            bool intersects_outer = meshes_intersect(V_nested, F, V_outer, F);
            bool intersects_inner = meshes_intersect(V_nested, F, V_outer_inner_shell, F);

            if(intersects_outer || intersects_inner) {
                viewer.data_list[4].set_colors(Eigen::RowVector3d(1.0, 0.0, 0.0));
            } else {
                viewer.data_list[4].set_colors(Eigen::RowVector3d(0.0, 1.0, 0.0));
            }
        }
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


        if (viewer.data_list.size() < 4) {
            return;
        }

        igl::swept_volume(
            viewer.data_list[0].V,viewer.data_list[0].F,transform_top,time_steps,grid_size,isolevel,SV_top,SF_top);

        igl::swept_volume(
            viewer.data_list[2].V,viewer.data_list[2].F,transform_bottom,time_steps,grid_size,isolevel,SV_bottom,SF_bottom);


        // clear_all_meshes();

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

    auto mesh_split = [&]() {
        Eigen::MatrixXd top_vertices;
        Eigen::MatrixXi top_faces;
        Eigen::VectorXi top_indices;

        Eigen::MatrixXd bottom_vertices;
        Eigen::MatrixXi bottom_faces;
        Eigen::VectorXi bottom_indices;

        top_vertices.resize(0, 3);
        top_faces.resize(0, 3);
        top_indices.resize(0, 3);

        bottom_vertices.resize(0, 3);
        bottom_faces.resize(0, 3);
        bottom_indices.resize(0, 3);


        bool success = split_mesh(V, F, top_vertices, top_faces, top_indices, bottom_vertices, bottom_faces, bottom_indices, cutting_plane_y_coord);

        Eigen::RowVector3d translationUp(0, 0.05, 0);
        Eigen::RowVector3d translationDown(0, -0.05, 0);

        top_vertices = top_vertices.rowwise() + translationUp;
        bottom_vertices = bottom_vertices.rowwise() + translationDown;

        clear_all_meshes();
        viewer.append_mesh();
        viewer.append_mesh();

        viewer.data_list[0].set_mesh(top_vertices, top_faces);
        viewer.data_list[1].set_mesh(bottom_vertices, bottom_faces);

        viewer.data_list[0].show_faces = true;
        viewer.data_list[0].show_lines = true;

        viewer.data_list[1].show_faces = true;
        viewer.data_list[1].show_lines = true;
    };

    auto stage_2 = [&]() {
        mesh_split();
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

            V1 = V;
            auto min_point = V1.colwise().minCoeff();
            auto max_point = V1.colwise().maxCoeff();
            // to multiply brush force proportional to size of mesh
            brush_strength = (max_point - min_point).norm() * 2;
            update_plane();
        } else {
            std::cerr << "Failed to load: " << filename << std::endl;
        }

        return success;
    };

    std::function<void()> ui_1, ui_2;

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
            mesh_resized();
        }

        if (ImGui::SliderFloat("Shell Thickness", &shell_thickness, 0.001f, 0.1f)) {
            shell_changed();
        }

        if(ImGui::SliderInt("Number of Nests", &n_nests, 1, 10)) {

        }

        if(ImGui::SliderFloat("Cutting Plane Y Translation", &cutting_plane_y_coord, -1, 1)) {
            update_plane();
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
            mesh_resized();
        }

        if (ImGui::Button("Reset Translation", ImVec2(-1, 0))) {
            translate_x = translate_y = translate_z = 0.0f;
            mesh_resized();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::Button("Visualize Shell Split", ImVec2(-1, 0))) {
            if(V.size() == 0) return;
            visualize_shell_split();
        }

        if (ImGui::Button("Visualize Swept Volume", ImVec2(-1, 0))) {
            if(V.size() == 0) return;
            visualize_swept_volume();
        }

        if (ImGui::Button("Generate dolls...", ImVec2(-1, 0))) {
            if(V.size() == 0) return;
            // menu.callback_draw_viewer_menu = ui_2;
            stage_2();
        }
    };


    ui_2 = [&]()
    {
        ImGui::Spacing();

        if (ImGui::Button("Save meshes...", ImVec2(0, 0))) {

        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::Button("Go back...", ImVec2(0, 0))) {
            menu.callback_draw_viewer_menu = ui_1;
            stage_1();
        }
    };

    viewer.callback_mouse_down =
        [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
        Eigen::Vector3f bc;
        int fid;
        auto x = viewer.current_mouse_x;
        auto y =
            viewer.core().viewport(3) - static_cast<float>(viewer.current_mouse_y);

\
        auto V_temp = V1;
        if (V_nested.size() != 0) {
            V_temp = V_nested;
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

    // viewer.callback_mouse_move =
    //     [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
    //     if (!posStart.isZero() && !posStart.hasNaN()) {
    //         Eigen::Vector3d oldPosEnd(posEnd.x(), posEnd.y(), posEnd.z());

    //         posEnd = igl::unproject(
    //                      Eigen::Vector3f(viewer.current_mouse_x,
    //                                      viewer.core().viewport[3] -
    //                                          static_cast<float>(viewer.current_mouse_y),
    //                                      viewer.down_mouse_z),
    //                      viewer.core().view,
    //                      viewer.core().proj,
    //                      viewer.core().viewport)
    //                      .template cast<double>();

    //         std::cout << "move" << std::endl;

    //         // exaggerate the force by a little bit
    //         Eigen::Vector3d forceVec = (posEnd - posStart) * brush_strength;

    //         int scaleFactor = forceVec.norm();
    //         if (posEnd.x() < posStart.x()) {
    //             // probably not the best way to determine direction.
    //             scaleFactor = -scaleFactor;
    //         }

    //         igl::kelvinlets(
    //             V,
    //             posStart,
    //             forceVec,
    //             mat,
    //             igl::KelvinletParams<double>(brushRadius, scale, brushType),
    //             result);


    //         // viewer.data_list[0].set_vertices(result);


    //         return true;
    //     }
    //     return false;
    // };

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
                V1,
                posStart,
                forceVec,
                mat,
                igl::KelvinletParams<double>(brushRadius, scale, brushType),
                result);


            // viewer.data_list[0].set_vertices(result);
            viewer.data_list[1].set_vertices(result);
            V1 = result;
            posStart.setZero();

            mesh_resized();

            return true;
        }

        // if (!posStart.isZero()) {



        //      std::cout << "up" << std::endl;

        //     // clear_all_meshes();
        //     // viewer.append_mesh();

        //     // viewer.data_list[0].set_mesh(V, F);
        //     // viewer.data_list[0].show_faces = true;
        //     // viewer.data_list[0].show_lines = true;
        //     // viewer.data_list[1].set_colors(Eigen::RowVector3d(0.8, 0.8, 0.9));

        //     posStart.setZero();
        //     return true;
        // }
        return false;
    };


    menu.callback_draw_viewer_menu = ui_1;

    viewer.launch();

    return 0;
}
