#include "PerlinNoise.hpp"
#include <cmath>
#include <algorithm>
#include <random>

PerlinNoise::PerlinNoise(unsigned int seed) {
    initialize(seed);
}

void PerlinNoise::initialize(unsigned int seed) {
    p.resize(512);
    std::vector<int> permutation(256);

    for (int i = 0; i < 256; ++i) { // Create a sequence of 0-255.
        permutation[i] = i;
    }

    std::mt19937 generator(seed); // Seed the Mersenne Twister engine.
    std::shuffle(permutation.begin(), permutation.end(), generator); // Shuffle the sequence to create a pseudo-random permutation.

    for (int i = 0; i < 256; ++i) {
        p[i] = permutation[i];
        p[i + 256] = permutation[i]; // Duplicate the permutation table to avoid expensive modulo operations during gradient selection.
    }

    gradients.resize(16); // Pre-defined gradient vectors. These are directions from the center towards the edges/vertices of a cube.
    gradients[0] = glm::vec3(1, 1, 0);
    gradients[1] = glm::vec3(-1, 1, 0);
    gradients[2] = glm::vec3(1, -1, 0);
    gradients[3] = glm::vec3(-1, -1, 0);
    gradients[4] = glm::vec3(1, 0, 1);
    gradients[5] = glm::vec3(-1, 0, 1);
    gradients[6] = glm::vec3(1, 0, -1);
    gradients[7] = glm::vec3(-1, 0, -1);
    gradients[8] = glm::vec3(0, 1, 1);
    gradients[9] = glm::vec3(0, -1, 1);
    gradients[10] = glm::vec3(0, 1, -1);
    gradients[11] = glm::vec3(0, -1, -1);
    gradients[12] = glm::vec3(1, 1, 0); // Some duplicates are acceptable for wrapping.
    gradients[13] = glm::vec3(-1, 1, 0);
    gradients[14] = glm::vec3(0, -1, 1);
    gradients[15] = glm::vec3(0, -1, -1);
}

float PerlinNoise::fade(float t) const {
    // Use a quintic interpolation curve: 6t^5 - 15t^4 + 10t^3.
    // This provides C2 continuity (zero first and second derivatives at t=0 and t=1), removing grid-line artifacts.
    return t * t * t * (t * (t * 6 - 15) + 10);
}

float PerlinNoise::lerp(float a, float b, float t) const {
    return a + t * (b - a); // Standard linear interpolation.
}

float PerlinNoise::dotGridGradient(int ix, int iy, float x, float y) const {
    int ix_wrapped = ix & 255; // Use bitwise AND for efficient modulo 256.
    int iy_wrapped = iy & 255;

    // Use the permutation table to hash the grid coordinates and select a pseudo-random gradient vector.
    int gradientIndex = p[p[ix_wrapped] + iy_wrapped] & 7; // Modulo 8 for 2D gradients.
    glm::vec2 gradient = glm::vec2(gradients[gradientIndex].x, gradients[gradientIndex].y);

    float dx = x - static_cast<float>(ix); // Compute the distance vector from the grid point to the input point.
    float dy = y - static_cast<float>(iy);

    return dx * gradient.x + dy * gradient.y; // Calculate the dot product.
}

float PerlinNoise::dotGridGradient(int ix, int iy, int iz, float x, float y, float z) const {
    int ix_wrapped = ix & 255;
    int iy_wrapped = iy & 255;
    int iz_wrapped = iz & 255;

    // The hashing process is extended for 3D by another table lookup.
    int gradientIndex = p[p[p[ix_wrapped] + iy_wrapped] + iz_wrapped] & 15; // Modulo 16 for 3D gradients.
    glm::vec3 gradient = gradients[gradientIndex];

    float dx = x - static_cast<float>(ix);
    float dy = y - static_cast<float>(iy);
    float dz = z - static_cast<float>(iz);

    return dx * gradient.x + dy * gradient.y + dz * gradient.z; // Dot product in 3D.
}

float PerlinNoise::noise(float x, float y) const {
    int x0 = static_cast<int>(std::floor(x)); // Find the integer coordinates of the grid cell containing the point.
    int x1 = x0 + 1;
    int y0 = static_cast<int>(std::floor(y));
    int y1 = y0 + 1;

    float sx = fade(x - static_cast<float>(x0)); // Calculate eased interpolation weights.
    float sy = fade(y - static_cast<float>(y0));

    float n0 = dotGridGradient(x0, y0, x, y); // Compute dot products for the four corners of the grid cell.
    float n1 = dotGridGradient(x1, y0, x, y);
    float ix0 = lerp(n0, n1, sx); // Interpolate along the x-axis.

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);
    float ix1 = lerp(n0, n1, sx);

    return lerp(ix0, ix1, sy); // Interpolate along the y-axis to get the final noise value.
}

float PerlinNoise::noise(float x, float y, float z) const {
    int x0 = static_cast<int>(std::floor(x)); // Integer coordinates for the containing cube.
    int x1 = x0 + 1;
    int y0 = static_cast<int>(std::floor(y));
    int y1 = y0 + 1;
    int z0 = static_cast<int>(std::floor(z));
    int z1 = z0 + 1;

    float sx = fade(x - static_cast<float>(x0)); // Eased interpolation weights for each axis.
    float sy = fade(y - static_cast<float>(y0));
    float sz = fade(z - static_cast<float>(z0));

    // The process extends to trilinear interpolation across the 8 vertices of the cube.
    float n0 = dotGridGradient(x0, y0, z0, x, y, z);
    float n1 = dotGridGradient(x1, y0, z0, x, y, z);
    float ix0 = lerp(n0, n1, sx); // Interpolate along x for the "bottom" face.

    n0 = dotGridGradient(x0, y1, z0, x, y, z);
    n1 = dotGridGradient(x1, y1, z0, x, y, z);
    float ix1 = lerp(n0, n1, sx);

    float iy0 = lerp(ix0, ix1, sy); // Interpolate along y for the "bottom" face.

    n0 = dotGridGradient(x0, y0, z1, x, y, z);
    n1 = dotGridGradient(x1, y0, z1, x, y, z);
    ix0 = lerp(n0, n1, sx); // Interpolate along x for the "top" face.

    n0 = dotGridGradient(x0, y1, z1, x, y, z);
    n1 = dotGridGradient(x1, y1, z1, x, y, z);
    ix1 = lerp(n0, n1, sx);

    float iy1 = lerp(ix0, ix1, sy); // Interpolate along y for the "top" face.

    return lerp(iy0, iy1, sz); // Finally, interpolate along z between the two faces.
}


float PerlinNoise::fractalNoise(float x, float y, int octaves, float persistence, float scale) const {
    // Implements Fractal Brownian Motion (fBm) by summing layers (octaves) of noise.
    float total = 0.0f;
    float frequency = scale; // Frequency doubles each octave.
    float amplitude = 1.0f; // Amplitude is modulated by persistence each octave.
    float maxValue = 0.0f; // Used for normalization to keep the output range predictable.

    for (int i = 0; i < octaves; ++i) {
        total += noise(x * frequency, y * frequency) * amplitude;
        maxValue += amplitude;
        amplitude *= persistence;
        frequency *= 2.0f;
    }

    return (total / maxValue + 1.0f) * 0.5f; // Normalize the result to the [0, 1] range.
}

float PerlinNoise::terrainHeight(float x, float y, float heightScale, float roughness) const {
    // A domain-specific application of fractal noise, combining octaves to simulate terrain.
    float height = fractalNoise(x, y, 4, 0.5f, 0.01f) * heightScale; // Base terrain
    height += fractalNoise(x, y, 6, 0.3f, 0.05f) * heightScale * 0.3f; // Medium details
    height += fractalNoise(x, y, 8, 0.2f, 0.1f) * heightScale * 0.1f; // Fine details
    return height * roughness;
}

float PerlinNoise::mountainHeight(float x, float y, float heightScale, float roughness) const {
    // Another variation of fBm, tuned to create more mountainous, dramatic features.
    float height = fractalNoise(x, y, 3, 0.6f, 0.005f) * heightScale;
    height += fractalNoise(x, y, 4, 0.4f, 0.02f) * heightScale * 0.4f;
    height += std::abs(fractalNoise(x, y, 5, 0.3f, 0.04f)) * heightScale * 0.6f; // `abs` creates sharper ridges.
    height = std::pow(height, 0.8f); // Power function can flatten valleys and sharpen peaks.
    return height * roughness;
}
