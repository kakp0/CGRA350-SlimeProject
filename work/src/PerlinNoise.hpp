#pragma once
#include <vector>
#include <glm/glm.hpp>

class PerlinNoise {
public:
    // Constructor
    PerlinNoise(unsigned int seed);

    // Re-initialize with a new seed
    void initialize(unsigned int seed);

    // Noise functions
    float noise(float x, float y) const;
    float noise(float x, float y, float z) const;
    float fractalNoise(float x, float y, int octaves, float persistence, float scale) const;

    // Terrain specific functions
    float terrainHeight(float x, float y, float heightScale, float roughness) const;
    float mountainHeight(float x, float y, float heightScale, float roughness) const;

    // NEW: Getters to allow the GPU to use the same noise data
    const std::vector<int>& getPermutationTable() const { return p; }
    const std::vector<glm::vec3>& getGradients() const { return gradients; }

private:
    std::vector<int> p;
    std::vector<glm::vec3> gradients;

    // Helper methods
    float fade(float t) const;
    float lerp(float a, float b, float t) const;
    float dotGridGradient(int ix, int iy, float x, float y) const;
    float dotGridGradient(int ix, int iy, int iz, float x, float y, float z) const;
};