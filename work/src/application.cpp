// std
#include <iostream>
#include <string>
#include <chrono>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"
#include "SlimeSimulation/SlimeBlock.hpp" 

using namespace std;
using namespace cgra;
using namespace glm;

// basic_model's draw function, which was missing
void basic_model::draw(const glm::mat4& view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform4fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));

	// Pass lighting info if the shader needs it
	glm::vec3 lightPosView = glm::vec3(view * glm::vec4(10, 10, 10, 1.0)); // Example light pos
	glUniform3fv(glGetUniformLocation(shader, "uLightPos"), 1, value_ptr(lightPosView));

	mesh.draw(); // draw
}


Application::Application(GLFWwindow* window) : m_window(window) {
	m_slimeBlock = std::make_unique<SlimeBlock>();

	// Build the ground plane mesh
	cgra::mesh_builder plane_mb;
	plane_mb.push_vertex({ vec3(-10, 0, -10), vec3(0, 1, 0), vec2(0, 0) });
	plane_mb.push_vertex({ vec3(10, 0, -10), vec3(0, 1, 0), vec2(1, 0) });
	plane_mb.push_vertex({ vec3(10, 0, 10), vec3(0, 1, 0), vec2(1, 1) });
	plane_mb.push_vertex({ vec3(-10, 0, 10), vec3(0, 1, 0), vec2(0, 1) });
	plane_mb.push_indices({ 0, 1, 2, 0, 2, 3 });
	m_ground_model.mesh = plane_mb.build();

	// Build the ground plane shader
	shader_builder plane_sb;
	plane_sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("/res/shaders/color_vert.glsl"));
	plane_sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("/res/shaders/color_frag.glsl"));
	m_ground_model.shader = plane_sb.build();
	m_ground_model.color = vec4(0.4f, 0.4f, 0.4f, 1.0f);
	m_ground_model.modelTransform = glm::translate(glm::mat4(1.0f), glm::vec3(0, -2.0f, 0));

	m_last_frame_time = std::chrono::high_resolution_clock::now();
}

Application::~Application() = default;

void Application::render() {
	auto now = std::chrono::high_resolution_clock::now();
	m_delta_time = std::chrono::duration_cast<std::chrono::duration<float>>(now - m_last_frame_time).count();
	m_last_frame_time = now;
	m_time_accumulator += m_delta_time;

	while (m_time_accumulator >= m_fixed_time_step) {
		m_slimeBlock->update(m_fixed_time_step);
		m_time_accumulator -= m_fixed_time_step;
	}

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

	mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);
	mat4 view = translate(mat4(1), vec3(0, 0, -m_distance))
		* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
		* rotate(mat4(1), m_yaw, vec3(0, 1, 0));

	vec3 camera_pos = vec3(inverse(view)[3]);

	if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);
	glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);

	m_slimeBlock->draw(view, proj, camera_pos, m_light_pos);

	// Draw the ground plane
	m_ground_model.draw(view, proj);


	glDisable(GL_BLEND);
}

void Application::renderGUI() {
	ImGui::SetNextWindowPos(ImVec2(5, 5), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(350, 450), ImGuiSetCond_Once);
	ImGui::Begin("Options", 0);

	ImGui::Text("Application %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::Checkbox("Show axis", &m_show_axis);
	ImGui::SameLine();
	ImGui::Checkbox("Show grid", &m_show_grid);
	ImGui::Checkbox("Wireframe", &m_showWireframe);
	ImGui::SameLine();
	if (ImGui::Button("Screenshot")) rgba_image::screenshot(true);

	ImGui::Separator();
	ImGui::Text("Camera Control");
	ImGui::SliderFloat("Pitch", &m_pitch, -pi<float>() / 2, pi<float>() / 2, "%.2f");
	ImGui::SliderFloat("Yaw", &m_yaw, -pi<float>(), pi<float>(), "%.2f");
	ImGui::SliderFloat("Distance", &m_distance, 0, 100, "%.2f", 2.0f);
	ImGui::SliderFloat3("Light Position", value_ptr(m_light_pos), -20.0f, 20.0f, "%.1f");

	ImGui::Separator();
	ImGui::Text("Slime Simulation");
	if (ImGui::Button("Reset Slime")) { m_slimeBlock->reset(); }
	ImGui::SameLine();
	ImGui::Checkbox("Rainbow Mode", m_slimeBlock->getRainbowModePtr());

	// CORRECTED: Subdivision limit changed to 9
	if (ImGui::SliderInt("Subdivisions", m_slimeBlock->getSubdivisionsPtr(), 0, 9)) { m_slimeBlock->reset(); }
	ImGui::SliderFloat3("Gravity", value_ptr(*m_slimeBlock->getGravityPtr()), -20.0f, 20.0f, "%.2f");
	ImGui::SliderInt("Solver Iterations", m_slimeBlock->getSolverIterationsPtr(), 1, 20);

	ImGui::SliderFloat("Stiffness", m_slimeBlock->getStiffnessPtr(), 0.0f, 1.0f, "%.3f");

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
		vec2 whsize = m_windowsize / 2.0f;
		m_pitch += float(acos(glm::clamp((m_mousePosition.y - whsize.y) / whsize.y, -1.0f, 1.0f)) - acos(glm::clamp((float(ypos) - whsize.y) / whsize.y, -1.0f, 1.0f)));
		m_pitch = float(glm::clamp(m_pitch, -pi<float>() / 2, pi<float>() / 2));
		m_yaw += float(acos(glm::clamp((m_mousePosition.x - whsize.x) / whsize.x, -1.0f, 1.0f)) - acos(glm::clamp((float(xpos) - whsize.x) / whsize.x, -1.0f, 1.0f)));
		if (m_yaw > pi<float>()) m_yaw -= float(2 * pi<float>());
		else if (m_yaw < -pi<float>()) m_yaw += float(2 * pi<float>());
	}
	m_mousePosition = vec2(xpos, ypos);
}

void Application::mouseButtonCallback(int button, int action, int mods) {
	(void)mods;
	if (button == GLFW_MOUSE_BUTTON_LEFT) m_leftMouseDown = (action == GLFW_PRESS);
}

void Application::scrollCallback(double xoffset, double yoffset) {
	(void)xoffset;
	m_distance *= pow(1.1f, -yoffset);
}

void Application::keyCallback(int key, int scancode, int action, int mods) {
	(void)key, (void)scancode, (void)action, (void)mods;
}

void Application::charCallback(unsigned int c) {
	(void)c;
}

