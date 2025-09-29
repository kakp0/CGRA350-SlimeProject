#include "PerlinNoise.hpp"
#include <cmath>
#include <algorithm>
#include <random>

PerlinNoise::PerlinNoise(unsigned int seed) {
    initialize(seed);
}

void PerlinNoise::initialize(unsigned int seed) {
    // Initialize permutation table with values 0-255
    p.resize(512);
    std::vector<int> permutation(256);

    for (int i = 0; i < 256; ++i) {
        permutation[i] = i;
    }

    std::mt19937 generator(seed);
    std::shuffle(permutation.begin(), permutation.end(), generator);

    for (int i = 0; i < 256; ++i) {
        p[i] = permutation[i];
        p[i + 256] = permutation[i];
    }

    // CHANGED: Initialize 3D gradient vectors
    gradients.resize(16);
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
    gradients[12] = glm::vec3(1, 1, 0); // Duplicates for wrapping
    gradients[13] = glm::vec3(-1, 1, 0);
    gradients[14] = glm::vec3(0, -1, 1);
    gradients[15] = glm::vec3(0, -1, -1);
}

float PerlinNoise::fade(float t) const {
    return t * t * t * (t * (t * 6 - 15) + 10);
}

float PerlinNoise::lerp(float a, float b, float t) const {
    return a + t * (b - a);
}

float PerlinNoise::dotGridGradient(int ix, int iy, float x, float y) const {
    int ix_wrapped = ix & 255;
    int iy_wrapped = iy & 255;

    int gradientIndex = p[p[ix_wrapped] + iy_wrapped] & 7;
    glm::vec2 gradient = glm::vec2(gradients[gradientIndex].x, gradients[gradientIndex].y);

    float dx = x - static_cast<float>(ix);
    float dy = y - static_cast<float>(iy);

    return dx * gradient.x + dy * gradient.y;
}

// ADDED: 3D dotGridGradient implementation
float PerlinNoise::dotGridGradient(int ix, int iy, int iz, float x, float y, float z) const {
    int ix_wrapped = ix & 255;
    int iy_wrapped = iy & 255;
    int iz_wrapped = iz & 255;

    int gradientIndex = p[p[p[ix_wrapped] + iy_wrapped] + iz_wrapped] & 15;
    glm::vec3 gradient = gradients[gradientIndex];

    float dx = x - static_cast<float>(ix);
    float dy = y - static_cast<float>(iy);
    float dz = z - static_cast<float>(iz);

    return dx * gradient.x + dy * gradient.y + dz * gradient.z;
}

float PerlinNoise::noise(float x, float y) const {
    int x0 = static_cast<int>(std::floor(x));
    int x1 = x0 + 1;
    int y0 = static_cast<int>(std::floor(y));
    int y1 = y0 + 1;

    float sx = fade(x - static_cast<float>(x0));
    float sy = fade(y - static_cast<float>(y0));

    float n0 = dotGridGradient(x0, y0, x, y);
    float n1 = dotGridGradient(x1, y0, x, y);
    float ix0 = lerp(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);
    float ix1 = lerp(n0, n1, sx);

    return lerp(ix0, ix1, sy);
}

// ADDED: 3D noise implementation
float PerlinNoise::noise(float x, float y, float z) const {
    int x0 = static_cast<int>(std::floor(x));
    int x1 = x0 + 1;
    int y0 = static_cast<int>(std::floor(y));
    int y1 = y0 + 1;
    int z0 = static_cast<int>(std::floor(z));
    int z1 = z0 + 1;

    float sx = fade(x - static_cast<float>(x0));
    float sy = fade(y - static_cast<float>(y0));
    float sz = fade(z - static_cast<float>(z0));

    float n0 = dotGridGradient(x0, y0, z0, x, y, z);
    float n1 = dotGridGradient(x1, y0, z0, x, y, z);
    float ix0 = lerp(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, z0, x, y, z);
    n1 = dotGridGradient(x1, y1, z0, x, y, z);
    float ix1 = lerp(n0, n1, sx);

    float iy0 = lerp(ix0, ix1, sy);

    n0 = dotGridGradient(x0, y0, z1, x, y, z);
    n1 = dotGridGradient(x1, y0, z1, x, y, z);
    ix0 = lerp(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, z1, x, y, z);
    n1 = dotGridGradient(x1, y1, z1, x, y, z);
    ix1 = lerp(n0, n1, sx);

    float iy1 = lerp(ix0, ix1, sy);

    return lerp(iy0, iy1, sz);
}

// ... (fractalNoise, terrainHeight, mountainHeight remain the same) ...

float PerlinNoise::fractalNoise(float x, float y, int octaves, float persistence, float scale) const {
    float total = 0.0f;
    float frequency = scale;
    float amplitude = 1.0f;
    float maxValue = 0.0f;

    for (int i = 0; i < octaves; ++i) {
        total += noise(x * frequency, y * frequency) * amplitude;
        maxValue += amplitude;
        amplitude *= persistence;
        frequency *= 2.0f;
    }

    return (total / maxValue + 1.0f) * 0.5f;
}

float PerlinNoise::terrainHeight(float x, float y, float heightScale, float roughness) const {
    float height = fractalNoise(x, y, 4, 0.5f, 0.01f) * heightScale;
    height += fractalNoise(x, y, 6, 0.3f, 0.05f) * heightScale * 0.3f;
    height += fractalNoise(x, y, 8, 0.2f, 0.1f) * heightScale * 0.1f;
    return height * roughness;
}

float PerlinNoise::mountainHeight(float x, float y, float heightScale, float roughness) const {
    float height = fractalNoise(x, y, 3, 0.6f, 0.005f) * heightScale;
    height += fractalNoise(x, y, 4, 0.4f, 0.02f) * heightScale * 0.4f;
    height += std::abs(fractalNoise(x, y, 5, 0.3f, 0.04f)) * heightScale * 0.6f;
    height = std::pow(height, 0.8f);
    return height * roughness;
}