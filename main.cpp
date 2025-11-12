#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <imgui.h>
#include "portable-file-dialogs.h"
#include <iostream>
#include "manip.h"

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Eigen::MatrixXd V_outer, V_outer_inner_shell, V_nested;

    std::vector<Eigen::MatrixXd> V_dolls;
    std::vector<Eigen::MatrixXi> F_dolls;


    float slider_value = 1.f;
    float shell_thickness = 0.01f;
    int n_nests = 1;
    int mesh_selection = 0;

    igl::opengl::glfw::Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    auto clear_all_meshes = [&]() {
        while(viewer.erase_mesh(viewer.selected_data_index)){};
        viewer.data().clear();
    };

    auto mesh_resized = [&]() {
        if(V.size() == 0) return;
        V_nested = resize_mesh(V, slider_value);

        viewer.data_list[1].set_vertices(V_nested);

        if(meshes_intersect(V_nested, F, V_outer_inner_shell, F)) {
            viewer.data_list[1].set_colors(Eigen::RowVector3d(1.0, 0.0, 0.0));
        } else {
            viewer.data_list[1].set_colors(Eigen::RowVector3d(0.0, 1.0, 0.0));
        }
    };

    auto shell_changed = [&]() {
        if(V.size() == 0) return;

        V_outer = V;
        V_outer_inner_shell = offset_along_normals(V_outer, F, -shell_thickness);

        viewer.data_list[0].set_vertices(V_outer);
        viewer.data_list[0].show_faces = false;
        viewer.data_list[2].set_vertices(V_outer_inner_shell);
        viewer.data_list[2].show_faces = false;
        mesh_resized();
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

    auto stage_2 = [&]() {
        clear_all_meshes();
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
            stage_1();
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

        if (ImGui::SliderFloat("Shell Thickness", &shell_thickness, 0.01f, 0.1f)) {
            shell_changed();
        }

        if(ImGui::SliderInt("Number of Nests", &n_nests, 1, 10)) {

        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::Button("Generate dolls...", ImVec2(-1, 0))) {
            if(V.size() == 0) return;
            menu.callback_draw_viewer_menu = ui_2;
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


    menu.callback_draw_viewer_menu = ui_1;

    viewer.launch();

    return 0;
}
