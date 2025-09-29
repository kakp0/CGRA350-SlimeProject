#pragma once

#include <vector>
#include <glm/glm.hpp>

class PerlinNoise {
public:
    // Constructor
    PerlinNoise(unsigned int seed = 0);

    // Generate 2D Perlin noise
    float noise(float x, float y) const;
    // ADDED: Generate 3D Perlin noise
    float noise(float x, float y, float z) const;

    // Generate fractal noise by combining multiple octaves
    float fractalNoise(float x, float y, int octaves, float persistence, float scale) const;

    // Generate 2D terrain height with configurable parameters
    float terrainHeight(float x, float y, float heightScale, float roughness) const;

    // Generate mountain height with more dramatic variations
    float mountainHeight(float x, float y, float heightScale, float roughness) const;

private:
    // Permutation table for improved randomness
    std::vector<int> p;

    // CHANGED: Gradient vectors are now 3D
    std::vector<glm::vec3> gradients;

    // Helper functions
    float fade(float t) const;
    float lerp(float a, float b, float t) const;
    float dotGridGradient(int ix, int iy, float x, float y) const;
    // ADDED: 3D version of dotGridGradient
    float dotGridGradient(int ix, int iy, int iz, float x, float y, float z) const;


    // Initialize the permutation table and gradients
    void initialize(unsigned int seed);
};