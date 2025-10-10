#pragma once
// Made with the assistance of AI (Gemini 2.5 Pro)

#include "cgra/cgra_mesh.hpp"
#include "opengl.hpp"
#include <glm/glm.hpp>
#include <vector>
#include <functional>

class ProceduralEnvironment;

// A simple struct representing a point mass in our simulation.

struct Particle {
    glm::vec3 position{ 0.0f };
    glm::vec3 predicted_position{ 0.0f };
    glm::vec3 velocity{ 0.0f };
    float inverse_mass = 1.0f; // Using inverse mass simplifies projections, especially for static points (inverse_mass = 0).
};

// Defines a simple distance constraint, the most basic building block for cloth or deformable bodies.
struct Constraint {
    int p1 = 0; // Index of the first particle.
    int p2 = 0; // Index of the second particle.
    float rest_distance = 0.0f; // The distance this constraint will try to maintain.
};

// Holds data for user interaction, allowing direct manipulation of the particle positions.
struct SurfaceGrabData {
    bool is_active = false;
    int p_indices[3]{ 0, 0, 0 }; // The three vertices of the grabbed triangle.
    glm::vec3 bary_coords{ 0.0f }; // Barycentric coordinates to pinpoint the grab location on the triangle.
    glm::vec3 target_pos{ 0.0f }; // Where the user is dragging the point to.
    glm::vec3 initial_hit_pos{ 0.0f };
    glm::vec3 grab_plane_normal{ 0.0f };
};


class SlimeBlock {
public:
    SlimeBlock();

    void update(float deltaTime, ProceduralEnvironment& environment);
    void draw(const glm::mat4& view, const glm::mat4& proj, const glm::vec3& view_pos, const glm::vec3& light_pos);
    void reset();

    // Mouse interaction for grabbing and manipulating the simulated object.
    bool findClosestTriangle(const glm::vec3& ray_origin, const glm::vec3& ray_dir, int& closest_tri_idx, glm::vec3& hit_pos, glm::vec3& bary_coords);
    void grabSurface(int p_idx0, int p_idx1, int p_idx2, const glm::vec3& barycentric_coords);
    void dragSurface(const glm::vec3& new_target_pos);
    void releaseSurface();
    void unstuck();

    // Basic accessors for particle and mesh data.
    const Particle& getParticle(int index) const;
    const std::vector<unsigned int>& getMeshIndices() const;

    // Accessors for GUI controls to tweak simulation parameters at runtime.
    bool* getRainbowModePtr() { return &m_rainbow_mode; }
    int* getSubdivisionsPtr() { return &m_subdivisions; }
    glm::vec3* getGravityPtr() { return &m_gravity; }
    int* getSolverIterationsPtr() { return &m_solver_iterations; } // More iterations = stiffer, but more expensive.
    float* getStiffnessPtr() { return &m_stiffness; }
    float* getDampingPtr() { return &m_damping; }
    float getInitialVolume() const { return m_initial_volume; }
    float getCurrentVolume() const { return m_current_volume; }
    float getKineticEnergy() const { return m_kinetic_energy; }

private:
    // The core of the PBD solver loop where constraints are enforced.
    void projectConstraints(ProceduralEnvironment& environment);

    void updateMesh(); // Helper to sync the visual mesh with the simulation particles.
    void dampVelocities(); // Applies damping to reduce energy and increase stability. 
    float calculateCurrentVolume(bool use_predicted_pos) const; // Used for the volume preservation constraint.

    // Simulation state
    std::vector<Particle> m_particles;
    std::vector<Constraint> m_constraints;
    cgra::gl_mesh m_mesh;
    GLuint m_shader;
    std::vector<cgra::mesh_vertex> m_mesh_vertices;
    std::vector<unsigned int> m_mesh_indices;

    // Tunable simulation parameters
    glm::vec3 m_spawn_position{ 0.0f, 50.0f, 0.0f };
    float m_scale = 4.0f;
    int m_subdivisions = 4;
    int m_solver_iterations = 10; 
    glm::vec3 m_gravity{ 0.0f, -9.81f, 0.0f }; 
    float m_stiffness = 0.07f; 
    float m_volume_stiffness = 0.8f;
    float m_damping = 0.01f;
    float m_ground_level = -2.0f;

    // Rendering properties
    glm::vec4 m_color{ 0.2f, 0.8f, 0.3f, 0.7f };
    glm::vec4 m_base_color{ 0.2f, 0.8f, 0.3f, 0.7f };
    bool m_rainbow_mode = false;
    float m_rainbow_time = 0.0f;

    // Diagnostic and state variables
    float m_initial_volume = 0.0f;
    float m_current_volume = 0.0f;
    float m_kinetic_energy = 0.0f;
    SurfaceGrabData m_surface_grab_data;
};
