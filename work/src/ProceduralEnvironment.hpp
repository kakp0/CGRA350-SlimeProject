#pragma once
#include <vector>
#include <string>
#include <memory>
#include "opengl.hpp"
#include <glm/glm.hpp>
#include "PerlinNoise.hpp"
#include "BlockTypes.hpp"

// Forward-declare SlimeBlock to prevent circular header dependencies
class SlimeBlock;

struct VoxelMesh {
    GLuint vao = 0;
    GLuint vbo[3]{ 0, 0, 0 };
    GLuint ibo = 0;
    GLsizei indexCount = 0;
};

struct Structure {
    int width = 0, depth = 0, height = 0;
    std::vector<BlockID> blockTypes;
    bool loadFromFile(const std::string& filepath);
    BlockID getBlock(int x, int y, int z) const;
};

class ProceduralEnvironment {
public:
    ProceduralEnvironment();
    void initialize();
    void updateParameters(float gh, float gr, float mh, float mr, float bd, unsigned int s);
    void regenerate();
    void setGrassContrast(float contrast);
    void render(const glm::mat4& view, const glm::mat4& projection);

    // NEW: Public method to check for blocks at a specific world position.
    BlockID getBlockAtWorldPos(const glm::vec3& worldPos) const;

    // Getters for GUI
    float getGrassHeight() const { return m_grassHeight; }
    float getGrassRoughness() const { return m_grassRoughness; }
    float getMountainHeight() const { return m_mountainHeight; }
    float getMountainRoughness() const { return m_mountainRoughness; }
    float getBlendDistance() const { return m_blendDistance; }
    unsigned int getSeed() const { return m_seed; }

private:
    // World Constants
    static const int WORLD_SIZE_X = 128;
    static const int WORLD_SIZE_Z = 128;
    static const int WORLD_HEIGHT = 64;

    // Member variables
    std::vector<BlockID> m_blockGrid;
    std::vector<BlockMaterial> m_blockMaterials;
    VoxelMesh m_voxelMesh;
    GLuint m_voxelShader;
    std::unique_ptr<PerlinNoise> m_perlinNoise;
    float m_grassHeight, m_grassRoughness, m_mountainHeight, m_mountainRoughness, m_blendDistance;
    unsigned int m_seed;
    std::vector<glm::vec3> m_voxelVertices;
    std::vector<glm::vec3> m_voxelNormals;
    std::vector<float> m_voxelBlockType;
    std::vector<unsigned int> m_voxelIndices;

    // Member variables to hold noise data for the GPU
    std::vector<int> m_noisePermutationTable;
    std::vector<glm::vec3> m_noiseGradients;
    float m_grassContrast = 2.0f;
    // Private Methods
    void populateBlockGrid();
    void placeStructures();
    void buildMeshFromGrid();
    void createVoxelMesh();
    void addFaceToMesh(const glm::vec3& pos, BlockID blockType, const std::string& face);
    void createShaders();
    void uploadBlockMaterials();
    void uploadNoiseUniforms();
    GLuint createShaderProgram(const char* vertexSource, const char* fragmentSource);

    BlockID getBlock(int x, int y, int z) const;
    void setBlock(int x, int y, int z, BlockID type);
};