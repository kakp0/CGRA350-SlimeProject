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

// *** ADDED *** Struct to hold information about a surface grab
struct SurfaceGrab {
	bool is_active = false;
	int p_indices[3]{ -1, -1, -1 };
	glm::vec3 bary_coords{ 0.0f };
	glm::vec3 initial_hit_pos{ 0.0f };
	glm::vec3 grab_plane_normal{ 0.0f };
};

// Main application class
class Application {
private:
	// window
	glm::vec2 m_windowsize;
	GLFWwindow* m_window;

	// orbital camera
	float m_pitch = .86;
	float m_yaw = -.86;
	float m_distance = 20;

	// last input
	bool m_leftMouseDown = false;
	glm::vec2 m_mousePosition;

	// *** MODIFIED *** Input state for new surface interaction
	bool m_rightMouseDown = false;
	SurfaceGrab m_surface_grab;

	// drawing flags
	bool m_show_axis = false;
	bool m_show_grid = false;
	bool m_showWireframe = false;

	// Scene lighting
	glm::vec3 m_light_pos{ 10, 10, 10 };

	// geometry
	std::unique_ptr<SlimeBlock> m_slimeBlock;
	basic_model m_ground_model;
	basic_model m_grab_sphere_model;

	// Physics simulation timing
	std::chrono::high_resolution_clock::time_point m_last_frame_time;
	float m_time_accumulator = 0.0f;
	const float m_fixed_time_step = 1.0f / 60.0f;
	float m_delta_time = 0;
	bool m_paused = false;

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
	void keyCallback(int key, int scancode, int action, int mods);
	void charCallback(unsigned int c);
};

