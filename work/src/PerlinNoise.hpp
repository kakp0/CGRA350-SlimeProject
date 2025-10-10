#pragma once // Prevents the header from being included multiple times in a single compilation unit.

#include <vector>
#include <glm/glm.hpp> // Using the GLM library for vector types, common in graphics programming.

// Implements Ken Perlin's "classic" noise algorithm, a procedural texture primitive.
class PerlinNoise {
public:
    PerlinNoise(unsigned int seed); // Constructor to initialize the noise generator with a specific seed.
    void initialize(unsigned int seed); // Re-seeds the generator, allowing for different noise patterns without creating a new object.

    // Generates 2D Perlin noise for a given coordinate.
    float noise(float x, float y) const;
    // Generates 3D Perlin noise for a given coordinate.
    float noise(float x, float y, float z) const;
    // Generates fractal noise by summing multiple "octaves" of Perlin noise.
    float fractalNoise(float x, float y, int octaves, float persistence, float scale) const;

    // Domain-specific application of fractal noise, often used for procedural terrain generation.
    float terrainHeight(float x, float y, float heightScale, float roughness) const;
    // Another domain-specific variant, possibly tuned for more mountainous features.
    float mountainHeight(float x, float y, float heightScale, float roughness) const;

    // Accessor for the permutation table, useful for GPU-based noise generation where this data is needed in a shader.
    const std::vector<int>& getPermutationTable() const { return p; }
    // Accessor for the gradient vectors, also for GPU/shader use.
    const std::vector<glm::vec3>& getGradients() const { return gradients; }

private:
    std::vector<int> p; // The permutation table, effectively a pre-shuffled array of 0-255 used for hashing grid coordinates. It's doubled to 512 elements to avoid modulo operations for wrap-around.
    std::vector<glm::vec3> gradients; // A lookup table of pre-defined, pseudo-random gradient vectors.

    // The quintic easing function (6t^5 - 15t^4 + 10t^3). Its first and second derivatives are zero at t=0 and t=1, which eliminates lattice-based artifacts seen with simpler interpolation.
    float fade(float t) const;
    // Standard linear interpolation between two values.
    float lerp(float a, float b, float t) const;
    // Calculates the dot product between the gradient vector at a grid point and the vector from that grid point to the sample point. This determines the gradient's influence. (2D version)
    float dotGridGradient(int ix, int iy, float x, float y) const;
    // 3D version of the dot product calculation.
    float dotGridGradient(int ix, int iy, int iz, float x, float y, float z) const;
};
