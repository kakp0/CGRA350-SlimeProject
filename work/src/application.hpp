#pragma once

// std
#include <memory> // ADDED: For std::unique_ptr

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"

// FORWARD DECLARATION: This tells the compiler that a class named "SlimeBlock" exists,
// without needing to include the full header. This is key to breaking the dependency.
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

	// drawing flags
	bool m_show_axis = false;
	bool m_show_grid = false;
	bool m_showWireframe = false;

	// geometry
	// CHANGED: We now use a smart pointer to our forward-declared class.
	std::unique_ptr<SlimeBlock> m_slimeBlock;

public:
	// setup
	Application(GLFWwindow*);
	// ADDED: A destructor is now required.
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

