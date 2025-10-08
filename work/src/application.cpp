#include <iostream>
#include <string>
#include <chrono>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"
#include "SlimeSimulation/SlimeBlock.hpp" // Ensure this path is correct for your project

using namespace std;
using namespace cgra;
using namespace glm;

void basic_model::draw(const glm::mat4& view, const glm::mat4 proj) {
    mat4 modelview = view * modelTransform;
    glUseProgram(shader);
    glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, GL_FALSE, value_ptr(proj));
    glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, GL_FALSE, value_ptr(modelview));
    glUniform4fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));
    GLint lightPosLoc = glGetUniformLocation(shader, "uLightPos");
    if (lightPosLoc != -1) {
        glm::vec3 lightPosView = glm::vec3(view * glm::vec4(10, 10, 10, 1.0));
        glUniform3fv(lightPosLoc, 1, value_ptr(lightPosView));
    }
    mesh.draw();
}

Application::Application(GLFWwindow* window) : m_window(window) {
    m_slimeBlock = std::make_unique<SlimeBlock>();
    m_environment = std::make_unique<ProceduralEnvironment>();
    m_environment->initialize();
    cgra::mesh_builder plane_mb;
    plane_mb.push_vertex({ vec3(-10, 0, -10), vec3(0, 1, 0), vec2(0, 0) });
    plane_mb.push_vertex({ vec3(10, 0, -10), vec3(0, 1, 0), vec2(1, 0) });
    plane_mb.push_vertex({ vec3(10, 0, 10), vec3(0, 1, 0), vec2(1, 1) });
    plane_mb.push_vertex({ vec3(-10, 0, 10), vec3(0, 1, 0), vec2(0, 1) });
    plane_mb.push_indices({ 0, 1, 2, 0, 2, 3 });
    m_ground_model.mesh = plane_mb.build();
    shader_builder plane_sb;
    plane_sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("/res/shaders/color_vert.glsl"));
    plane_sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("/res/shaders/color_frag.glsl"));
    m_ground_model.shader = plane_sb.build();
    m_ground_model.color = vec4(0.4f, 0.4f, 0.4f, 1.0f);
    m_ground_model.modelTransform = glm::translate(glm::mat4(1.0f), glm::vec3(0, -2.0f, 0));
    shader_builder unlit_sb;
    unlit_sb.set_shader_source(GL_VERTEX_SHADER, R"RAW(
		#version 330 core
		layout (location = 0) in vec3 aPos;
		uniform mat4 uProjectionMatrix;
		uniform mat4 uModelViewMatrix;
		void main() {
			gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(aPos, 1.0);
		}
	)RAW");
    unlit_sb.set_shader_source(GL_FRAGMENT_SHADER, R"RAW(
		#version 330 core
		out vec4 FragColor;
		uniform vec4 uColor;
		void main() {
			FragColor = uColor;
		}
	)RAW");
    GLuint unlit_shader = unlit_sb.build();
    cgra::mesh_builder grab_cube_mb;
    {
        float s = 1.0f;
        std::vector<glm::vec3> positions = {
            {-s, -s,  s}, { s, -s,  s}, { s,  s,  s}, {-s,  s,  s},
            {-s, -s, -s}, { s, -s, -s}, { s,  s, -s}, {-s,  s, -s}
        };
        std::vector<unsigned int> indices = {
            0, 1, 2, 0, 2, 3, 5, 4, 7, 5, 7, 6, 1, 5, 6, 1, 6, 2,
            4, 0, 3, 4, 3, 7, 3, 2, 6, 3, 6, 7, 4, 5, 1, 4, 1, 0
        };
        for (const auto& pos : positions) {
            grab_cube_mb.push_vertex({ pos, glm::vec3(0,0,0), glm::vec2(0,0) });
        }
        grab_cube_mb.indices = indices;
    }
    m_grab_sphere_model.mesh = grab_cube_mb.build();
    m_grab_sphere_model.shader = unlit_shader;
    m_grab_sphere_model.color = glm::vec4(1.0, 1.0, 0.0, 1.0);
    m_last_frame_time = std::chrono::high_resolution_clock::now();
}

Application::~Application() = default;

void Application::render() {
    auto now = std::chrono::high_resolution_clock::now();
    m_delta_time = std::chrono::duration_cast<std::chrono::duration<float>>(now - m_last_frame_time).count();
    m_last_frame_time = now;
    m_time_accumulator += m_delta_time;
    while (m_time_accumulator >= m_fixed_time_step) {
        m_slimeBlock->update(m_fixed_time_step, *m_environment);
        m_time_accumulator -= m_fixed_time_step;
    }

    // --- Camera Movement ---
    float current_speed = m_camera_speed * m_delta_time;
    if (m_w_down) m_camera_pos += current_speed * m_camera_front;
    if (m_s_down) m_camera_pos -= current_speed * m_camera_front;
    if (m_a_down) m_camera_pos -= glm::normalize(glm::cross(m_camera_front, m_camera_up)) * current_speed;
    if (m_d_down) m_camera_pos += glm::normalize(glm::cross(m_camera_front, m_camera_up)) * current_speed;
    if (m_space_down) m_camera_pos += current_speed * m_camera_up;
    if (m_shift_down) m_camera_pos -= current_speed * m_camera_up;


    int width, height;
    glfwGetFramebufferSize(m_window, &width, &height);
    m_windowsize = vec2(width, height);
    glViewport(0, 0, width, height);
    glClearColor(0.3f, 0.3f, 0.4f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // --- View and Projection Matrix Setup ---
    mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);
    mat4 view = glm::lookAt(m_camera_pos, m_camera_pos + m_camera_front, m_camera_up);

    vec3 camera_pos = m_camera_pos;
    if (m_show_grid) drawGrid(view, proj);
    if (m_show_axis) drawAxis(view, proj);
    glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);
    m_environment->render(view, proj);
    m_ground_model.draw(view, proj);
    m_slimeBlock->draw(view, proj, camera_pos, m_light_pos);
    if (m_surface_grab.is_active) {
        float time = glfwGetTime();
        float pulse = 0.5f * (1.0f + sin(time * 8.0f));
        float scale_pulse = 0.15f + pulse * 0.1f;
        m_grab_sphere_model.color = glm::vec4(1.0, 1.0, 0.0, 0.5f + pulse * 0.5f);
        const auto& p0 = m_slimeBlock->getParticle(m_surface_grab.p_indices[0]);
        const auto& p1 = m_slimeBlock->getParticle(m_surface_grab.p_indices[1]);
        const auto& p2 = m_slimeBlock->getParticle(m_surface_grab.p_indices[2]);
        glm::vec3 current_grab_pos = p0.position * m_surface_grab.bary_coords.x +
            p1.position * m_surface_grab.bary_coords.y +
            p2.position * m_surface_grab.bary_coords.z;
        m_grab_sphere_model.modelTransform = glm::translate(glm::mat4(1.0f), current_grab_pos) * glm::scale(glm::mat4(1.0f), glm::vec3(scale_pulse));
        glDisable(GL_DEPTH_TEST);
        m_grab_sphere_model.draw(view, proj);
        glEnable(GL_DEPTH_TEST);
    }
    glDisable(GL_BLEND);
}


void Application::renderGUI() {
    ImGui::SetNextWindowPos(ImVec2(5, 5), ImGuiSetCond_Once);
    ImGui::SetNextWindowSize(ImVec2(350, 480), ImGuiSetCond_Once);
    ImGui::Begin("Options", 0);

    ImGui::Text("Application %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    ImGui::Checkbox("Show axis", &m_show_axis); ImGui::SameLine();
    ImGui::Checkbox("Show grid", &m_show_grid);
    ImGui::Checkbox("Wireframe", &m_showWireframe); ImGui::SameLine();
    if (ImGui::Button("Screenshot")) rgba_image::screenshot(true);

    ImGui::Separator();
    ImGui::Text("Camera Control");
    ImGui::SliderFloat("Pitch", &m_pitch, -pi<float>() / 2, pi<float>() / 2, "%.2f");
    ImGui::SliderFloat("Yaw", &m_yaw, -pi<float>(), pi<float>(), "%.2f");
    ImGui::SliderFloat("Camera Speed", &m_camera_speed, 1.0f, 100.0f, "%.1f");
    ImGui::SliderFloat3("Light Position", value_ptr(m_light_pos), -20.0f, 20.0f, "%.1f");

    ImGui::Separator();
    ImGui::Text("Environment Parameters");
    static float grassHeight = m_environment->getGrassHeight();
    static float grassRoughness = m_environment->getGrassRoughness();
    static float mountainHeight = m_environment->getMountainHeight();
    static float mountainRoughness = m_environment->getMountainRoughness();
    static float blendDistance = m_environment->getBlendDistance();
    static int seed = static_cast<int>(m_environment->getSeed());

    bool environmentChanged = false;
    if (ImGui::SliderFloat("Grass Height", &grassHeight, 0.1f, 10.0f, "%.1f")) environmentChanged = true;
    if (ImGui::SliderFloat("Grass Roughness", &grassRoughness, 0.1f, 2.0f, "%.2f")) environmentChanged = true;
    if (ImGui::SliderFloat("Mountain Height", &mountainHeight, 10.0f, 200.0f, "%.1f")) environmentChanged = true;
    if (ImGui::SliderFloat("Mountain Roughness", &mountainRoughness, 0.1f, 3.0f, "%.2f")) environmentChanged = true;
    if (ImGui::SliderFloat("Blend Distance", &blendDistance, 10.0f, 200.0f, "%.1f")) environmentChanged = true;

    ImGui::SliderFloat("Grass Contrast", &m_grassContrast, 0.5f, 4.0f, "%.2f");

    if (ImGui::SliderInt("Seed", &seed, 0, 1000)) environmentChanged = true;

    m_environment->setGrassContrast(m_grassContrast);

    if (environmentChanged) {
        m_environment->updateParameters(grassHeight, grassRoughness, mountainHeight,
            mountainRoughness, blendDistance, static_cast<unsigned int>(seed));
    }
    if (ImGui::Button("Regenerate Terrain")) {
        m_environment->regenerate();
    }

    ImGui::Separator();
    ImGui::Text("Slime Simulation");
    if (ImGui::Button("Reset Slime")) { m_slimeBlock->reset(); }
    ImGui::SameLine();
    ImGui::Checkbox("Rainbow Mode", m_slimeBlock->getRainbowModePtr());
    if (ImGui::SliderInt("Subdivisions", m_slimeBlock->getSubdivisionsPtr(), 0, 9)) { m_slimeBlock->reset(); }
    ImGui::SliderFloat3("Gravity", value_ptr(*m_slimeBlock->getGravityPtr()), -20.0f, 20.0f, "%.2f");
    ImGui::SliderInt("Solver Iterations", m_slimeBlock->getSolverIterationsPtr(), 1, 20);
    ImGui::SliderFloat("Stiffness", m_slimeBlock->getStiffnessPtr(), 0.0f, 1.0f, "%.3f");
    ImGui::SliderFloat("Damping", m_slimeBlock->getDampingPtr(), 0.0f, 0.1f, "%.3f");

    ImGui::Separator();
    ImGui::Text("Debug Stats");
    float initial_volume = m_slimeBlock->getInitialVolume();
    float current_volume = m_slimeBlock->getCurrentVolume();
    float volume_error = (initial_volume != 0) ? ((current_volume - initial_volume) / initial_volume) * 100.0f : 0.0f;
    ImGui::Text("Initial Volume: %.2f", initial_volume);
    ImGui::Text("Current Volume: %.2f", current_volume);
    ImGui::Text("Volume Error: %.2f %%", volume_error);
    ImGui::Text("Kinetic Energy: %.2f", m_slimeBlock->getKineticEnergy());

    ImGui::End();
}

void Application::cursorPosCallback(double xpos, double ypos) {
    if (m_leftMouseDown) {
        float xoffset = xpos - m_mousePosition.x;
        float yoffset = m_mousePosition.y - ypos; // reversed since y-coordinates go from top to bottom

        float sensitivity = 0.005f;
        xoffset *= sensitivity;
        yoffset *= sensitivity;

        m_yaw += xoffset;
        m_pitch += yoffset;

        // Clamp pitch to avoid flipping
        m_pitch = glm::clamp(m_pitch, -glm::pi<float>() / 2.0f + 0.01f, glm::pi<float>() / 2.0f - 0.01f);

        // Recalculate camera front vector
        glm::vec3 front;
        front.x = cos(m_yaw) * cos(m_pitch);
        front.y = sin(m_pitch);
        front.z = sin(m_yaw) * cos(m_pitch);
        m_camera_front = glm::normalize(front);
    }
    else if (m_rightMouseDown && m_surface_grab.is_active) {
        // Get mouse movement delta for depth adjustment
        float yoffset = m_mousePosition.y - ypos;

        // If holding Left Ctrl, adjust the grab depth
        if (m_ctrl_down) {
            float depth_sensitivity = 0.1f;
            m_surface_grab.grab_depth += yoffset * depth_sensitivity;
            if (m_surface_grab.grab_depth < 1.0f) {
                m_surface_grab.grab_depth = 1.0f;
            }
        }

        // Calculate new target position using ray-plane intersection
        int width, height;
        glfwGetFramebufferSize(m_window, &width, &height);
        mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);
        mat4 view = glm::lookAt(m_camera_pos, m_camera_pos + m_camera_front, m_camera_up);
        glm::vec4 viewport = glm::vec4(0, 0, width, height);

        // 1. Create a ray from the current mouse position
        glm::vec3 screen_pos_near(xpos, m_windowsize.y - ypos, 0.0f);
        glm::vec3 screen_pos_far(xpos, m_windowsize.y - ypos, 1.0f);
        glm::vec3 ray_origin = glm::unProject(screen_pos_near, view, proj, viewport);
        glm::vec3 ray_far = glm::unProject(screen_pos_far, view, proj, viewport);
        glm::vec3 ray_dir = glm::normalize(ray_far - ray_origin);

        // 2. Define the dragging plane (perpendicular to camera's view)
        glm::vec3 plane_normal = -m_camera_front;
        glm::vec3 plane_point = m_camera_pos + m_camera_front * m_surface_grab.grab_depth;

        // 3. Calculate intersection of the ray with the plane
        float denom = dot(ray_dir, plane_normal);
        if (abs(denom) > 1e-6) { // Avoid division by zero
            float t = dot(plane_point - ray_origin, plane_normal) / denom;
            glm::vec3 new_target_pos = ray_origin + ray_dir * t;
            m_slimeBlock->dragSurface(new_target_pos);
        }
    }
    m_mousePosition = vec2(xpos, ypos);
}
void Application::mouseButtonCallback(int button, int action, int mods) {
    (void)mods;
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        m_leftMouseDown = (action == GLFW_PRESS);
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        m_rightMouseDown = (action == GLFW_PRESS);
        if (m_rightMouseDown) {
            int width, height;
            glfwGetFramebufferSize(m_window, &width, &height);
            mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);
            mat4 view = glm::lookAt(m_camera_pos, m_camera_pos + m_camera_front, m_camera_up);
            glm::vec3 screen_pos_near(m_mousePosition.x, m_windowsize.y - m_mousePosition.y, 0.0f);
            glm::vec3 screen_pos_far(m_mousePosition.x, m_windowsize.y - m_mousePosition.y, 1.0f);
            glm::vec4 viewport = glm::vec4(0, 0, width, height);
            glm::vec3 world_pos_near = glm::unProject(screen_pos_near, view, proj, viewport);
            glm::vec3 world_pos_far = glm::unProject(screen_pos_far, view, proj, viewport);
            glm::vec3 ray_dir = glm::normalize(world_pos_far - world_pos_near);
            glm::vec3 ray_origin = world_pos_near;
            int tri_idx;
            glm::vec3 hit_pos, bary_coords;
            if (m_slimeBlock->findClosestTriangle(ray_origin, ray_dir, tri_idx, hit_pos, bary_coords)) {
                const auto& indices = m_slimeBlock->getMeshIndices();
                m_surface_grab.is_active = true;
                m_surface_grab.p_indices[0] = indices[tri_idx * 3];
                m_surface_grab.p_indices[1] = indices[tri_idx * 3 + 1];
                m_surface_grab.p_indices[2] = indices[tri_idx * 3 + 2];
                m_surface_grab.bary_coords = bary_coords;
                m_surface_grab.initial_hit_pos = hit_pos;

                // Calculate and store initial grab depth
                m_surface_grab.grab_depth = glm::distance(m_camera_pos, hit_pos);

                m_slimeBlock->grabSurface(m_surface_grab.p_indices[0], m_surface_grab.p_indices[1], m_surface_grab.p_indices[2], bary_coords);
            }
        }
        else {
            m_slimeBlock->releaseSurface();
            m_surface_grab.is_active = false;
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS) {
        m_slimeBlock->unstuck();
    }
}
void Application::scrollCallback(double xoffset, double yoffset) {
    (void)xoffset;
    m_camera_speed += yoffset * 2.0f;
    if (m_camera_speed < 1.0f) m_camera_speed = 1.0f;
}
void Application::keyCallback(int key, int scancode, int action, int mods) {
    (void)scancode, (void)mods;
    if (key == GLFW_KEY_W) m_w_down = (action != GLFW_RELEASE);
    if (key == GLFW_KEY_A) m_a_down = (action != GLFW_RELEASE);
    if (key == GLFW_KEY_S) m_s_down = (action != GLFW_RELEASE);
    if (key == GLFW_KEY_D) m_d_down = (action != GLFW_RELEASE);
    if (key == GLFW_KEY_SPACE) m_space_down = (action != GLFW_RELEASE);
    if (key == GLFW_KEY_LEFT_SHIFT) m_shift_down = (action != GLFW_RELEASE);
    if (key == GLFW_KEY_LEFT_CONTROL) m_ctrl_down = (action != GLFW_RELEASE);
}
void Application::charCallback(unsigned int c) {
    (void)c;
}
