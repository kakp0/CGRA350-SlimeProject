#include "SlimeBlock.hpp"
#include "ProceduralEnvironment.hpp" // Include the environment for collision checks
#include "cgra/cgra_shader.hpp"
#include "BlockTypes.hpp" // Include to get BlockID enum
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <numeric>
#include <cmath>
#include <limits>

SlimeBlock::SlimeBlock() {
    cgra::shader_builder sb;
    sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("/res/shaders/color_vert.glsl"));
    sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("/res/shaders/color_frag.glsl"));
    m_shader = sb.build();

    reset();
}

void SlimeBlock::reset() {
    m_particles.clear();
    m_constraints.clear();
    m_mesh_vertices.clear();
    m_mesh_indices.clear();

    int dim = m_subdivisions + 2;
    float step = 1.0f / (dim - 1);

    // Create particles
    for (int k = 0; k < dim; ++k) {
        for (int j = 0; j < dim; ++j) {
            for (int i = 0; i < dim; ++i) {
                Particle p;
                glm::vec3 local_pos = glm::vec3(i * step, j * step, k * step) - glm::vec3(0.5f);
                p.position = local_pos * m_scale + m_spawn_position;
                p.predicted_position = p.position;
                p.velocity = glm::vec3(0.0f);
                p.inverse_mass = 1.0f;
                m_particles.push_back(p);
            }
        }
    }

    auto get_index = [&](int i, int j, int k) {
        return i + j * dim + k * dim * dim;
        };

    // Create constraints
    for (int k = 0; k < dim; ++k) {
        for (int j = 0; j < dim; ++j) {
            for (int i = 0; i < dim; ++i) {
                int p0_idx = get_index(i, j, k);
                if (i < dim - 1) m_constraints.push_back({ p0_idx, get_index(i + 1, j, k), step * m_scale });
                if (j < dim - 1) m_constraints.push_back({ p0_idx, get_index(i, j + 1, k), step * m_scale });
                if (k < dim - 1) m_constraints.push_back({ p0_idx, get_index(i, j, k + 1), step * m_scale });
                if (i < dim - 1 && j < dim - 1) m_constraints.push_back({ p0_idx, get_index(i + 1, j + 1, k), glm::length(glm::vec3(step, step, 0.0f)) * m_scale });
                if (i < dim - 1 && k < dim - 1) m_constraints.push_back({ p0_idx, get_index(i + 1, j, k + 1), glm::length(glm::vec3(step, 0.0f, step)) * m_scale });
                if (j < dim - 1 && k < dim - 1) m_constraints.push_back({ p0_idx, get_index(i, j + 1, k + 1), glm::length(glm::vec3(0.0f, step, step)) * m_scale });
            }
        }
    }

    float main_diag_dist = glm::length(glm::vec3(1.0f, 1.0f, 1.0f)) * m_scale;
    m_constraints.push_back({ get_index(0, 0, 0), get_index(dim - 1, dim - 1, dim - 1), main_diag_dist });
    m_constraints.push_back({ get_index(dim - 1, 0, 0), get_index(0, dim - 1, dim - 1), main_diag_dist });
    m_constraints.push_back({ get_index(0, dim - 1, 0), get_index(dim - 1, 0, dim - 1), main_diag_dist });
    m_constraints.push_back({ get_index(0, 0, dim - 1), get_index(dim - 1, dim - 1, 0), main_diag_dist });

    // Build mesh for rendering
    for (int j = 0; j < dim - 1; ++j) {
        for (int i = 0; i < dim - 1; ++i) {
            m_mesh_indices.push_back(get_index(i, j, 0)); m_mesh_indices.push_back(get_index(i + 1, j, 0)); m_mesh_indices.push_back(get_index(i + 1, j + 1, 0));
            m_mesh_indices.push_back(get_index(i, j, 0)); m_mesh_indices.push_back(get_index(i + 1, j + 1, 0)); m_mesh_indices.push_back(get_index(i, j + 1, 0));
            m_mesh_indices.push_back(get_index(i, j, dim - 1)); m_mesh_indices.push_back(get_index(i + 1, j + 1, dim - 1)); m_mesh_indices.push_back(get_index(i + 1, j, dim - 1));
            m_mesh_indices.push_back(get_index(i, j, dim - 1)); m_mesh_indices.push_back(get_index(i, j + 1, dim - 1)); m_mesh_indices.push_back(get_index(i + 1, j + 1, dim - 1));
        }
    }
    for (int k = 0; k < dim - 1; ++k) {
        for (int i = 0; i < dim - 1; ++i) {
            m_mesh_indices.push_back(get_index(i, dim - 1, k)); m_mesh_indices.push_back(get_index(i + 1, dim - 1, k)); m_mesh_indices.push_back(get_index(i + 1, dim - 1, k + 1));
            m_mesh_indices.push_back(get_index(i, dim - 1, k)); m_mesh_indices.push_back(get_index(i + 1, dim - 1, k + 1)); m_mesh_indices.push_back(get_index(i, dim - 1, k + 1));
            m_mesh_indices.push_back(get_index(i, 0, k)); m_mesh_indices.push_back(get_index(i + 1, 0, k + 1)); m_mesh_indices.push_back(get_index(i + 1, 0, k));
            m_mesh_indices.push_back(get_index(i, 0, k)); m_mesh_indices.push_back(get_index(i, 0, k + 1)); m_mesh_indices.push_back(get_index(i + 1, 0, k + 1));
        }
    }
    for (int k = 0; k < dim - 1; ++k) {
        for (int j = 0; j < dim - 1; ++j) {
            m_mesh_indices.push_back(get_index(0, j, k)); m_mesh_indices.push_back(get_index(0, j + 1, k + 1)); m_mesh_indices.push_back(get_index(0, j, k + 1));
            m_mesh_indices.push_back(get_index(0, j, k)); m_mesh_indices.push_back(get_index(0, j + 1, k)); m_mesh_indices.push_back(get_index(0, j + 1, k + 1));
            m_mesh_indices.push_back(get_index(dim - 1, j, k)); m_mesh_indices.push_back(get_index(dim - 1, j, k + 1)); m_mesh_indices.push_back(get_index(dim - 1, j + 1, k + 1));
            m_mesh_indices.push_back(get_index(dim - 1, j, k)); m_mesh_indices.push_back(get_index(dim - 1, j + 1, k + 1)); m_mesh_indices.push_back(get_index(dim - 1, j + 1, k));
        }
    }
    for (const auto& p : m_particles) {
        m_mesh_vertices.push_back({ p.position, glm::vec3(0, 1, 0), glm::vec2(0) });
    }

    cgra::mesh_builder mb;
    mb.vertices = m_mesh_vertices;
    mb.indices = m_mesh_indices;
    m_mesh = mb.build();
    glBindBuffer(GL_ARRAY_BUFFER, m_mesh.vbo);
    glBufferData(GL_ARRAY_BUFFER, m_mesh_vertices.size() * sizeof(cgra::mesh_vertex), m_mesh_vertices.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    m_initial_volume = std::abs(calculateCurrentVolume(false));
}

void SlimeBlock::update(float deltaTime, ProceduralEnvironment& environment) {
    if (m_rainbow_mode) {
        m_rainbow_time += deltaTime;
        m_color = glm::vec4(sin(m_rainbow_time * 0.5f) * 0.5f + 0.5f, sin(m_rainbow_time * 0.5f + 2.0f) * 0.5f + 0.5f, sin(m_rainbow_time * 0.5f + 4.0f) * 0.5f + 0.5f, 0.5f);
    }
    else {
        m_color = m_base_color;
    }

    if (m_particles.empty()) return;

    for (auto& p : m_particles) {
        if (p.inverse_mass > 0) {
            p.velocity += m_gravity * deltaTime;
            p.predicted_position = p.position + p.velocity * deltaTime;
        }
    }

    projectConstraints(environment);

    for (auto& p : m_particles) {
        if (m_surface_grab_data.is_active && (
            &p == &m_particles[m_surface_grab_data.p_indices[0]] ||
            &p == &m_particles[m_surface_grab_data.p_indices[1]] ||
            &p == &m_particles[m_surface_grab_data.p_indices[2]]))
        {
            p.position = p.predicted_position;
        }
        else if (p.inverse_mass > 0) {
            p.velocity = (p.predicted_position - p.position) / deltaTime;
            p.position = p.predicted_position;
        }
    }
    dampVelocities();
    m_current_volume = calculateCurrentVolume(false);
    m_kinetic_energy = 0.0f;
    for (auto& p : m_particles) {
        if (p.inverse_mass > 0) m_kinetic_energy += 0.5f * (1.0f / p.inverse_mass) * glm::dot(p.velocity, p.velocity);
    }
    updateMesh();
}

void SlimeBlock::dampVelocities() {
    if (m_particles.size() <= 1 || m_damping == 0.0f) return;
    glm::vec3 x_cm(0.0f), v_cm(0.0f); float total_mass = 0.0f;
    for (const auto& p : m_particles) {
        if (p.inverse_mass > 0) { float mass = 1.0f / p.inverse_mass; x_cm += p.position * mass; v_cm += p.velocity * mass; total_mass += mass; }
    }
    if (total_mass == 0.0f) return;
    x_cm /= total_mass; v_cm /= total_mass;
    glm::vec3 L(0.0f); glm::mat3 I(0.0f);
    for (const auto& p : m_particles) {
        if (p.inverse_mass > 0) {
            float mass = 1.0f / p.inverse_mass; glm::vec3 r = p.position - x_cm; L += glm::cross(r, p.velocity * mass);
            glm::mat3 r_tilde(0, r.z, -r.y, -r.z, 0, r.x, r.y, -r.x, 0); I += r_tilde * glm::transpose(r_tilde) * mass;
        }
    }
    if (abs(glm::determinant(I)) < 1e-9) return;
    glm::vec3 omega = glm::inverse(I) * L;
    omega.y *= 0.1f;
    for (auto& p : m_particles) {
        if (p.inverse_mass > 0) { glm::vec3 r = p.position - x_cm; glm::vec3 dv = (v_cm + glm::cross(omega, r)) - p.velocity; p.velocity += m_damping * dv; }
    }
}

float SlimeBlock::calculateCurrentVolume(bool use_predicted_pos) const {
    if (m_particles.size() < 4) return 0.0f;
    float volume = 0.0f;
    for (size_t i = 0; i < m_mesh_indices.size(); i += 3) {
        const auto& p1 = use_predicted_pos ? m_particles[m_mesh_indices[i]].predicted_position : m_particles[m_mesh_indices[i]].position;
        const auto& p2 = use_predicted_pos ? m_particles[m_mesh_indices[i + 1]].predicted_position : m_particles[m_mesh_indices[i + 1]].position;
        const auto& p3 = use_predicted_pos ? m_particles[m_mesh_indices[i + 2]].predicted_position : m_particles[m_mesh_indices[i + 2]].position;
        volume += glm::dot(p1, glm::cross(p2, p3));
    }
    return volume / 6.0f;
}

void SlimeBlock::projectConstraints(ProceduralEnvironment& environment) {
    float effective_stiffness = glm::clamp(m_stiffness * (m_subdivisions + 1.0f), 0.0f, 1.0f);
    float k_prime = 1.0f - pow(1.0f - effective_stiffness, 1.0f / m_solver_iterations);

    float volume_k_prime = 1.0f - pow(1.0f - m_volume_stiffness, 1.0f / m_solver_iterations);

    for (int iter = 0; iter < m_solver_iterations; ++iter) {
        // Apply grab constraints
        if (m_surface_grab_data.is_active) {
            Particle& p0 = m_particles[m_surface_grab_data.p_indices[0]];
            Particle& p1 = m_particles[m_surface_grab_data.p_indices[1]];
            Particle& p2 = m_particles[m_surface_grab_data.p_indices[2]];

            float b0 = m_surface_grab_data.bary_coords.x;
            float b1 = m_surface_grab_data.bary_coords.y;
            float b2 = m_surface_grab_data.bary_coords.z;

            glm::vec3 current_grab_point = p0.predicted_position * b0 + p1.predicted_position * b1 + p2.predicted_position * b2;
            glm::vec3 correction_vec = m_surface_grab_data.target_pos - current_grab_point;

            float w0 = p0.inverse_mass, w1 = p1.inverse_mass, w2 = p2.inverse_mass;
            float denominator = w0 * b0 * b0 + w1 * b1 * b1 + w2 * b2 * b2;

            if (denominator > 1e-6) {
                float k_grab = 1.0f; // Grab is always stiff
                glm::vec3 correction = (correction_vec / denominator) * k_grab;
                p0.predicted_position += correction * w0 * b0;
                p1.predicted_position += correction * w1 * b1;
                p2.predicted_position += correction * w2 * b2;
            }
        }

        // Distance constraints
        for (const auto& c : m_constraints) {
            Particle& p1 = m_particles[c.p1]; Particle& p2 = m_particles[c.p2];
            glm::vec3 p1_p2 = p2.predicted_position - p1.predicted_position;
            float current_dist = glm::length(p1_p2);
            if (current_dist > 1e-6) {
                float error = current_dist - c.rest_distance;
                glm::vec3 correction = p1_p2 * (error / current_dist);
                float w1 = p1.inverse_mass; float w2 = p2.inverse_mass;
                if (w1 + w2 > 0) {
                    p1.predicted_position += (w1 / (w1 + w2)) * correction * k_prime;
                    p2.predicted_position -= (w2 / (w1 + w2)) * correction * k_prime;
                }
            }
        }

        // Volume preservation
        float current_volume = calculateCurrentVolume(true);
        float volume_error = current_volume - m_initial_volume;
        if (std::abs(volume_error) > 1e-6) {
            std::vector<glm::vec3> gradients(m_particles.size(), glm::vec3(0.0f));
            float sum_gradient_sq = 0.0f;
            for (size_t j = 0; j < m_mesh_indices.size(); j += 3) {
                int i1 = m_mesh_indices[j], i2 = m_mesh_indices[j + 1], i3 = m_mesh_indices[j + 2];
                const auto& p1 = m_particles[i1].predicted_position, & p2 = m_particles[i2].predicted_position, & p3 = m_particles[i3].predicted_position;
                gradients[i1] += glm::cross(p2, p3); gradients[i2] += glm::cross(p3, p1); gradients[i3] += glm::cross(p1, p2);
            }
            for (size_t i = 0; i < gradients.size(); ++i) {
                gradients[i] /= 6.0f;
                sum_gradient_sq += glm::dot(gradients[i], gradients[i]) * m_particles[i].inverse_mass;
            }
            if (sum_gradient_sq > 1e-9) {
                float lambda = -volume_error / sum_gradient_sq;
                for (size_t i = 0; i < m_particles.size(); ++i) {
                    if (m_particles[i].inverse_mass > 0.0f) {
                        m_particles[i].predicted_position += lambda * m_particles[i].inverse_mass * gradients[i] * volume_k_prime;
                    }
                }
            }
        }

        // --- REVISED COLLISION LOGIC ---
        for (auto& p : m_particles) {
            // Check for collision with the environment
            BlockID block = environment.getBlockAtWorldPos(p.predicted_position);

            // If the block is solid (not Air or Glass)
            if (block != BlockID::Air && block != BlockID::Glass) {
                // This is a simple but effective response:
                // Stop the particle by resetting its predicted position to its last known good position.
                p.predicted_position = p.position;
                continue; // Skip the ground plane check if we already collided
            }

            // Fallback Ground collision
            if (p.predicted_position.y < m_ground_level) {
                p.predicted_position.y = m_ground_level;
            }
        }
    }
}

void SlimeBlock::updateMesh() {
    if (m_particles.empty() || m_mesh_vertices.empty()) return;
    for (size_t i = 0; i < m_particles.size(); ++i) {
        m_mesh_vertices[i].pos = m_particles[i].position;
    }
    std::vector<glm::vec3> temp_normals(m_particles.size(), glm::vec3(0.0f));
    for (size_t i = 0; i < m_mesh_indices.size(); i += 3) {
        unsigned int i0 = m_mesh_indices[i], i1 = m_mesh_indices[i + 1], i2 = m_mesh_indices[i + 2];
        glm::vec3 p0 = m_particles[i0].position, p1 = m_particles[i1].position, p2 = m_particles[i2].position;
        glm::vec3 face_normal = glm::cross(p1 - p0, p2 - p0);
        temp_normals[i0] += face_normal; temp_normals[i1] += face_normal; temp_normals[i2] += face_normal;
    }
    for (size_t i = 0; i < m_particles.size(); ++i) {
        if (glm::length(temp_normals[i]) > 0.0f) {
            m_mesh_vertices[i].norm = glm::normalize(temp_normals[i]);
        }
    }
    glBindBuffer(GL_ARRAY_BUFFER, m_mesh.vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, m_mesh_vertices.size() * sizeof(cgra::mesh_vertex), m_mesh_vertices.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void SlimeBlock::draw(const glm::mat4& view, const glm::mat4& proj, const glm::vec3& view_pos, const glm::vec3& light_pos) {
    glm::mat4 model_matrix = glm::mat4(1.0f);
    glm::mat4 modelview = view * model_matrix;
    glUseProgram(m_shader);
    glUniformMatrix4fv(glGetUniformLocation(m_shader, "uProjectionMatrix"), 1, false, glm::value_ptr(proj));
    glUniformMatrix4fv(glGetUniformLocation(m_shader, "uModelViewMatrix"), 1, false, glm::value_ptr(modelview));
    glUniform4fv(glGetUniformLocation(m_shader, "uColor"), 1, glm::value_ptr(m_color));
    glm::vec3 lightPosView = glm::vec3(view * glm::vec4(light_pos, 1.0));
    glUniform3fv(glGetUniformLocation(m_shader, "uLightPos"), 1, glm::value_ptr(lightPosView));
    m_mesh.draw();
}

bool SlimeBlock::findClosestTriangle(const glm::vec3& ray_origin, const glm::vec3& ray_dir, int& closest_tri_idx, glm::vec3& hit_pos, glm::vec3& bary_coords) {
    float min_dist = std::numeric_limits<float>::max();
    bool found_intersection = false;

    for (size_t i = 0; i < m_mesh_indices.size(); i += 3) {
        const auto& p0 = m_particles[m_mesh_indices[i]].position;
        const auto& p1 = m_particles[m_mesh_indices[i + 1]].position;
        const auto& p2 = m_particles[m_mesh_indices[i + 2]].position;

        glm::vec3 edge1 = p1 - p0;
        glm::vec3 edge2 = p2 - p0;
        glm::vec3 h = glm::cross(ray_dir, edge2);
        float a = glm::dot(edge1, h);

        if (a > -1e-6 && a < 1e-6) continue;

        float f = 1.0 / a;
        glm::vec3 s = ray_origin - p0;
        float u = f * glm::dot(s, h);

        if (u < 0.0 || u > 1.0) continue;

        glm::vec3 q = glm::cross(s, edge1);
        float v = f * glm::dot(ray_dir, q);

        if (v < 0.0 || u + v > 1.0) continue;

        float t = f * glm::dot(edge2, q);

        if (t > 1e-6 && t < min_dist) {
            min_dist = t;
            closest_tri_idx = i / 3;
            hit_pos = ray_origin + ray_dir * t;
            bary_coords = { 1.0f - u - v, u, v };
            found_intersection = true;
        }
    }
    return found_intersection;
}

void SlimeBlock::grabSurface(int p_idx0, int p_idx1, int p_idx2, const glm::vec3& barycentric_coords) {
    m_surface_grab_data.is_active = true;
    m_surface_grab_data.p_indices[0] = p_idx0;
    m_surface_grab_data.p_indices[1] = p_idx1;
    m_surface_grab_data.p_indices[2] = p_idx2;
    m_surface_grab_data.bary_coords = barycentric_coords;
}

void SlimeBlock::dragSurface(const glm::vec3& new_target_pos) {
    if (m_surface_grab_data.is_active) {
        m_surface_grab_data.target_pos = new_target_pos;
    }
}

void SlimeBlock::releaseSurface() {
    if (m_surface_grab_data.is_active) {
        m_surface_grab_data.is_active = false;
    }
}

const Particle& SlimeBlock::getParticle(int index) const {
    return m_particles[index];
}

const std::vector<unsigned int>& SlimeBlock::getMeshIndices() const {
    return m_mesh_indices;
}