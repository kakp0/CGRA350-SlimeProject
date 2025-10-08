#pragma once

// std
#include <memory>
#include <chrono>

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "ProceduralEnvironment.hpp"

class SlimeBlock;

// Basic model that holds the shader, mesh and transform for drawing.
struct basic_model {
    GLuint shader = 0;
    cgra::gl_mesh mesh;
    glm::vec4 color{ 0.7, 0.7, 0.7, 1.0 };
    glm::mat4 modelTransform{ 1.0 };
    GLuint texture;

    void draw(const glm::mat4& view, const glm::mat4 proj);
};

// Struct to hold information about a surface grab
struct SurfaceGrab {
    bool is_active = false;
    int p_indices[3]{ -1, -1, -1 };
    glm::vec3 bary_coords{ 0.0f };
    glm::vec3 initial_hit_pos{ 0.0f };
    glm::vec3 grab_plane_normal{ 0.0f };
    float grab_depth = 0.0f; // Distance from camera to grabbed point
};

// Main application class
class Application {
private:
    // window
    glm::vec2 m_windowsize;
    GLFWwindow* m_window;

    // Free-look camera parameters
    glm::vec3 m_camera_pos = glm::vec3(0.0f, 15.0f, 40.0f);
    glm::vec3 m_camera_front = glm::vec3(0.0f, 0.0f, -1.0f);
    glm::vec3 m_camera_up = glm::vec3(0.0f, 1.0f, 0.0f);
    float m_camera_speed = 20.0f;

    // Camera orientation
    float m_pitch = 0.0f;
    float m_yaw = -glm::pi<float>() / 2.0f; // Start looking down the -Z axis

    // last input
    bool m_leftMouseDown = false;
    bool m_rightMouseDown = false;
    glm::vec2 m_mousePosition;
    SurfaceGrab m_surface_grab;

    // Key states for movement
    bool m_w_down = false;
    bool m_a_down = false;
    bool m_s_down = false;
    bool m_d_down = false;
    bool m_space_down = false;
    bool m_shift_down = false;
    bool m_ctrl_down = false; // To toggle depth control


    // drawing flags
    bool m_show_axis = false;
    bool m_show_grid = false;
    bool m_showWireframe = false;

    // Scene lighting
    glm::vec3 m_light_pos{ 10, 10, 10 };

    // geometry
    std::unique_ptr<SlimeBlock> m_slimeBlock;
    std::unique_ptr<ProceduralEnvironment> m_environment;
    basic_model m_ground_model;
    basic_model m_grab_sphere_model;

    // Physics simulation timing
    std::chrono::high_resolution_clock::time_point m_last_frame_time;
    float m_time_accumulator = 0.0f;
    const float m_fixed_time_step = 1.0f / 60.0f;
    float m_delta_time = 0;
    bool m_paused = false;

    float m_grassContrast = 2.0f;

public:
    // setup
    Application(GLFWwindow*);
    ~Application();

    // disable copy constructors (for safety)
    Application(const Application&) = delete;
    Application& operator=(const Application&) = delete;

    // rendering callbacks (every frame)
    void render();
    void renderGUI();

    // input callbacks
    void cursorPosCallback(double xpos, double ypos);
    void mouseButtonCallback(int button, int action, int mods);
    void scrollCallback(double xoffset, double yoffset);
    void keyCallback(int key, int scocode, int action, int mods);
    void charCallback(unsigned int c);
};