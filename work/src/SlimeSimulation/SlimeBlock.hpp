#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "cgra/cgra_mesh.hpp"

// A single particle in our simulation
struct Particle {
    glm::vec3 position;
    glm::vec3 predicted_position;
    glm::vec3 velocity;
    float inverse_mass;
};

// A constraint between two particles
struct DistanceConstraint {
    int p1, p2;
    float rest_distance;
};

// *** MODIFIED *** Simplified the grab data; we no longer need to store original mass
struct SurfaceGrabData {
    bool is_active = false;
    int p_indices[3]{ -1, -1, -1 };
    glm::vec3 bary_coords{ 0.0f };
    glm::vec3 target_pos{ 0.0f };
};


class SlimeBlock {
private:
    // Simulation data
    std::vector<Particle> m_particles;
    std::vector<DistanceConstraint> m_constraints;

    // Rendering data
    cgra::gl_mesh m_mesh;
    GLuint m_shader = 0;
    std::vector<cgra::mesh_vertex> m_mesh_vertices;
    std::vector<unsigned int> m_mesh_indices;

    // Simulation parameters
    int m_subdivisions = 1;
    glm::vec3 m_gravity{ 0.0f, -9.81f, 0.0f };
    int m_solver_iterations = 5;
    float m_stiffness = 0.5f;
    float m_volume_stiffness = 1.0f;
    float m_damping = 0.02f;
    float m_initial_volume = 0.0f;
    float m_current_volume = 0.0f;
    float m_kinetic_energy = 0.0f;

    float m_ground_level = -2.0f;

    // World placement
    float m_scale = 4.0f;
    glm::vec3 m_spawn_position{ 0.0f, 10.0f, 0.0f };

    // Appearance
    glm::vec4 m_color{ 0.1, 0.8, 0.2, 0.5 };
    glm::vec4 m_base_color{ 0.1, 0.8, 0.2, 0.5 };
    bool m_rainbow_mode = false;
    float m_rainbow_time = 0.0f;

    // Interaction state
    SurfaceGrabData m_surface_grab_data;

    // Private methods
    void projectConstraints();
    void updateMesh();
    float calculateCurrentVolume(bool use_predicted_pos) const;
    void dampVelocities();

public:
    SlimeBlock();
    void update(float deltaTime);
    void draw(const glm::mat4& view, const glm::mat4& proj, const glm::vec3& view_pos, const glm::vec3& light_pos);
    void reset();

    // Methods for new surface interaction
    bool findClosestTriangle(const glm::vec3& ray_origin, const glm::vec3& ray_dir, int& tri_idx, glm::vec3& hit_pos, glm::vec3& bary_coords);
    void grabSurface(int p_idx0, int p_idx1, int p_idx2, const glm::vec3& barycentric_coords);
    void dragSurface(const glm::vec3& new_target_pos);
    void releaseSurface();

    // Public getters for particles and indices
    const Particle& getParticle(int index) const;
    const std::vector<unsigned int>& getMeshIndices() const;


    // Getters for GUI
    int* getSubdivisionsPtr() { return &m_subdivisions; }
    glm::vec3* getGravityPtr() { return &m_gravity; }
    int* getSolverIterationsPtr() { return &m_solver_iterations; }
    float* getStiffnessPtr() { return &m_stiffness; }
    float* getDampingPtr() { return &m_damping; }
    float getInitialVolume() const { return m_initial_volume; }
    float getCurrentVolume() const { return m_current_volume; }
    float getKineticEnergy() const { return m_kinetic_energy; }
    bool* getRainbowModePtr() { return &m_rainbow_mode; }
};

