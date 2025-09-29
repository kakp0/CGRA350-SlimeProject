#pragma once

#include "PerlinNoise.hpp"
#include "cgra/cgra_mesh.hpp"
#include "cgra/cgra_shader.hpp"
#include <glm/glm.hpp>
#include <vector>
#include <memory>
#include <string>

// A simple structure to hold data loaded from a structure file
struct Structure {
    int width = 0, depth = 0, height = 0;
    std::vector<int> blockTypes;

    bool loadFromFile(const std::string& filepath);
    int getBlock(int x, int y, int z) const;
};

class ProceduralEnvironment {
public:
    ProceduralEnvironment();
    ~ProceduralEnvironment() = default;

    void initialize();
    void render(const glm::mat4& view, const glm::mat4& projection);
    void updateParameters(float grassHeight, float grassRoughness, float mountainHeight, float mountainRoughness, float blendDistance, unsigned int seed);
    void regenerate();

    float getGrassHeight() const { return m_grassHeight; }
    float getGrassRoughness() const { return m_grassRoughness; }
    float getMountainHeight() const { return m_mountainHeight; }
    float getMountainRoughness() const { return m_mountainRoughness; }
    float getBlendDistance() const { return m_blendDistance; }
    unsigned int getSeed() const { return m_seed; }

private:
    float m_grassHeight, m_grassRoughness, m_mountainHeight, m_mountainRoughness, m_blendDistance;
    unsigned int m_seed;

    // --- CHANGE APPLIED HERE ---
    // Dimensions increased from 96 to 128 to create a larger area
    static const int WORLD_SIZE_X = 128;
    static const int WORLD_SIZE_Z = 128;
    static const int WORLD_HEIGHT = 64;

    std::unique_ptr<PerlinNoise> m_perlinNoise;
    GLuint m_voxelShader;

    struct VoxelMesh {
        GLuint vao;
        GLuint vbo[3];
        GLuint ibo;
        GLsizei indexCount;
    };
    VoxelMesh m_voxelMesh;

    std::vector<glm::vec3> m_voxelVertices;
    std::vector<glm::vec3> m_voxelNormals;
    std::vector<float> m_voxelBlockType;
    std::vector<unsigned int> m_voxelIndices;

    // 3D Grid for World Data
    std::vector<int> m_blockGrid;
    int getBlock(int x, int y, int z) const;
    void setBlock(int x, int y, int z, int type);

    // Restructured Generation Functions
    void populateBlockGrid();
    void placeStructures();
    void buildMeshFromGrid();
    void addFaceToMesh(const glm::vec3& pos, float blockType, const std::string& face);
    void createVoxelMesh();

    void createShaders();
    GLuint createShaderProgram(const char* vertexSource, const char* fragmentSource);
};