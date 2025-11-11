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

    Eigen::MatrixXd V_shell;
    Eigen::MatrixXi F_shell;

    Eigen::MatrixXd V_nested;
    Eigen::MatrixXi F_nested;

    float slider_value = 1.f;
    int mesh_selection = 0;

    igl::opengl::glfw::Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    viewer.append_mesh();

    auto mesh_resized = [&]() {
        if(V.size() == 0) return;
        V_nested = resize_mesh(V, slider_value);
        F_nested = F;

        viewer.data_list[1].set_vertices(V_nested);

        if(meshes_intersect(V_nested, F_nested, V_shell, F_shell)) {
            viewer.data_list[1].set_colors(Eigen::RowVector3d(1.0, 0.0, 0.0));
        } else {
            viewer.data_list[1].set_colors(Eigen::RowVector3d(0.0, 1.0, 0.0));
        }
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
            V_shell = V;
            F_shell = F;
            viewer.data_list[0].clear();
            viewer.data_list[0].set_mesh(V,F);
            viewer.data_list[1].clear();
            viewer.data_list[1].set_mesh(V,F);
            viewer.data_list[0].show_faces = false;
            mesh_resized();
        } else {
            std::cerr << "Failed to load: " << filename << std::endl;
        }

        return success;
    };

    menu.callback_draw_viewer_menu = [&]()
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

        if (ImGui::SliderFloat("Parameter", &slider_value, 0.0f, 1.0f)) {
            mesh_resized();
        }

    };

    viewer.launch();

    return 0;
}
