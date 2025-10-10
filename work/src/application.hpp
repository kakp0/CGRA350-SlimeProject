#pragma once

#include <memory>
#include <chrono>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "ProceduralEnvironment.hpp"

class SlimeBlock; // Forward-declaration to avoid circular includes.

// A simple aggregate of data for a renderable object.
struct basic_model {
	GLuint shader = 0;
	cgra::gl_mesh mesh;
	glm::vec4 color{ 0.7, 0.7, 0.7, 1.0 };
	glm::mat4 modelTransform{ 1.0 };
	GLuint texture;

	void draw(const glm::mat4& view, const glm::mat4 proj);
};

// State for interactively grabbing and manipulating a point on a mesh surface.
struct SurfaceGrab {
	bool is_active = false;
	int p_indices[3]{ -1, -1, -1 };      // Vertex indices of the grabbed triangle.
	glm::vec3 bary_coords{ 0.0f };       // Barycentric coords of the grab point within the triangle.
	glm::vec3 initial_hit_pos{ 0.0f };   // Initial intersection point in world space.
	glm::vec3 grab_plane_normal{ 0.0f }; // The normal of the plane used for dragging calculations.
	float grab_depth = 0.0f;             // Distance from camera to the grab point.
};

// Main application class orchestrating the window, rendering, and user input.
class Application {
private:
	glm::vec2 m_windowsize;
	GLFWwindow* m_window;

	// Camera state using a free-look model.
	glm::vec3 m_camera_pos = glm::vec3(0.0f, 15.0f, 40.0f);
	glm::vec3 m_camera_front = glm::vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 m_camera_up = glm::vec3(0.0f, 1.0f, 0.0f);
	float m_camera_speed = 20.0f;
	float m_pitch = 0.0f;                      // Camera rotation around the x-axis.
	float m_yaw = -glm::pi<float>() / 2.0f;    // Camera rotation around the y-axis.

	// Input state tracking.
	bool m_leftMouseDown = false;
	bool m_rightMouseDown = false;
	glm::vec2 m_mousePosition;
	SurfaceGrab m_surface_grab;

	// Key press states for smooth movement.
	bool m_w_down = false;
	bool m_a_down = false;
	bool m_s_down = false;
	bool m_d_down = false;
	bool m_space_down = false;
	bool m_shift_down = false;
	bool m_ctrl_down = false; // Toggles depth control for surface grab.

	// Debug and rendering flags.
	bool m_show_axis = false;
	bool m_show_grid = false;
	bool m_showWireframe = false;

	glm::vec3 m_light_pos{ 10, 10, 10 }; // A single point light for the scene.

	// Scene object management.
	std::unique_ptr<SlimeBlock> m_slimeBlock;
	std::unique_ptr<ProceduralEnvironment> m_environment;
	basic_model m_ground_model;
	basic_model m_grab_sphere_model; // Visual indicator for the grab point.

	// Timing for a fixed-step physics simulation.
	std::chrono::high_resolution_clock::time_point m_last_frame_time;
	float m_time_accumulator = 0.0f;            // Accumulates frame time to be consumed by fixed steps.
	const float m_fixed_time_step = 1.0f / 60.0f; // Ensures physics is framerate-independent.
	float m_delta_time = 0;
	bool m_paused = false;

	float m_grassContrast = 2.0f; // A shader-specific parameter.

public:
	Application(GLFWwindow*);
	~Application();

	// Disallow copying to prevent accidental duplication of resources.
	Application(const Application&) = delete;
	Application& operator=(const Application&) = delete;

	void render();
	void renderGUI();

	// GLFW input callbacks.
	void cursorPosCallback(double xpos, double ypos);
	void mouseButtonCallback(int button, int action, int mods);
	void scrollCallback(double xoffset, double yoffset);
	void keyCallback(int key, int scocode, int action, int mods);
	void charCallback(unsigned int c);
};
