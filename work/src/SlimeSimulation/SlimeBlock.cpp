// src/SlimeSimulation/SlimeBlock.cpp

#include "SlimeBlock.hpp"
#include "cgra/cgra_shader.hpp"
#include <glm/gtc/matrix_transform.hpp>

SlimeBlock::SlimeBlock() {
    // Use the shader builder from the original application code
    cgra::shader_builder sb;
    sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
    sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_frag.glsl"));
    m_model.shader = sb.build();

    // MOVED: The cube-building logic is now directly inside SlimeBlock.
    // This makes the class self-contained.
    cgra::mesh_builder mb;

    // Add vertices for a cube.
    // We need 24 vertices because each face has unique normals.
    // Front face
    mb.push_vertex({ glm::vec3(-0.5, -0.5, 0.5), glm::vec3(0, 0, 1), glm::vec2(0, 0) });
    mb.push_vertex({ glm::vec3(0.5, -0.5, 0.5), glm::vec3(0, 0, 1), glm::vec2(1, 0) });
    mb.push_vertex({ glm::vec3(0.5, 0.5, 0.5), glm::vec3(0, 0, 1), glm::vec2(1, 1) });
    mb.push_vertex({ glm::vec3(-0.5, 0.5, 0.5), glm::vec3(0, 0, 1), glm::vec2(0, 1) });
    // Back face
    mb.push_vertex({ glm::vec3(0.5, -0.5, -0.5), glm::vec3(0, 0, -1), glm::vec2(0, 0) });
    mb.push_vertex({ glm::vec3(-0.5, -0.5, -0.5), glm::vec3(0, 0, -1), glm::vec2(1, 0) });
    mb.push_vertex({ glm::vec3(-0.5, 0.5, -0.5), glm::vec3(0, 0, -1), glm::vec2(1, 1) });
    mb.push_vertex({ glm::vec3(0.5, 0.5, -0.5), glm::vec3(0, 0, -1), glm::vec2(0, 1) });
    // Right face
    mb.push_vertex({ glm::vec3(0.5, -0.5, 0.5), glm::vec3(1, 0, 0), glm::vec2(0, 0) });
    mb.push_vertex({ glm::vec3(0.5, -0.5, -0.5), glm::vec3(1, 0, 0), glm::vec2(1, 0) });
    mb.push_vertex({ glm::vec3(0.5, 0.5, -0.5), glm::vec3(1, 0, 0), glm::vec2(1, 1) });
    mb.push_vertex({ glm::vec3(0.5, 0.5, 0.5), glm::vec3(1, 0, 0), glm::vec2(0, 1) });
    // Left face
    mb.push_vertex({ glm::vec3(-0.5, -0.5, -0.5), glm::vec3(-1, 0, 0), glm::vec2(0, 0) });
    mb.push_vertex({ glm::vec3(-0.5, -0.5, 0.5), glm::vec3(-1, 0, 0), glm::vec2(1, 0) });
    mb.push_vertex({ glm::vec3(-0.5, 0.5, 0.5), glm::vec3(-1, 0, 0), glm::vec2(1, 1) });
    mb.push_vertex({ glm::vec3(-0.5, 0.5, -0.5), glm::vec3(-1, 0, 0), glm::vec2(0, 1) });
    // Top face
    mb.push_vertex({ glm::vec3(-0.5, 0.5, 0.5), glm::vec3(0, 1, 0), glm::vec2(0, 0) });
    mb.push_vertex({ glm::vec3(0.5, 0.5, 0.5), glm::vec3(0, 1, 0), glm::vec2(1, 0) });
    mb.push_vertex({ glm::vec3(0.5, 0.5, -0.5), glm::vec3(0, 1, 0), glm::vec2(1, 1) });
    mb.push_vertex({ glm::vec3(-0.5, 0.5, -0.5), glm::vec3(0, 1, 0), glm::vec2(0, 1) });
    // Bottom face
    mb.push_vertex({ glm::vec3(-0.5, -0.5, -0.5), glm::vec3(0, -1, 0), glm::vec2(0, 0) });
    mb.push_vertex({ glm::vec3(0.5, -0.5, -0.5), glm::vec3(0, -1, 0), glm::vec2(1, 0) });
    mb.push_vertex({ glm::vec3(0.5, -0.5, 0.5), glm::vec3(0, -1, 0), glm::vec2(1, 1) });
    mb.push_vertex({ glm::vec3(-0.5, -0.5, 0.5), glm::vec3(0, -1, 0), glm::vec2(0, 1) });

    // Add indices for the 12 triangles of the cube
    mb.push_indices({ 0, 1, 2, 0, 2, 3 });       // Front
    mb.push_indices({ 4, 5, 6, 4, 6, 7 });       // Back
    mb.push_indices({ 8, 9, 10, 8, 10, 11 });    // Right
    mb.push_indices({ 12, 13, 14, 12, 14, 15 }); // Left
    mb.push_indices({ 16, 17, 18, 16, 18, 19 }); // Top
    mb.push_indices({ 20, 21, 22, 20, 22, 23 }); // Bottom

    m_model.mesh = mb.build();

    // Set the color to a translucent green (Red, Green, Blue, Alpha)
    m_model.color = glm::vec4(0.1, 0.8, 0.2, 0.5); // 50% transparent

    // Optional: Scale the cube to make it bigger
    m_model.modelTransform = glm::scale(glm::mat4(1.0), glm::vec3(4.0));
}

void SlimeBlock::draw(const glm::mat4& view, const glm::mat4& proj) {
    // Use the basic_model's own draw method
    m_model.draw(view, proj);
}

