#pragma once
#include <vector>
#include <string>
#include <memory>
#include "opengl.hpp"
#include <glm/glm.hpp>
#include "PerlinNoise.hpp"
#include "BlockTypes.hpp"

class SlimeBlock; // Forward-declaration to prevent circular header dependencies with SlimeBlock.hpp

// Encapsulates all OpenGL objects needed to render the voxel world mesh.
struct VoxelMesh {
    GLuint vao = 0; // Vertex Array Object: Stores the state of the mesh's vertex configuration.
    GLuint vbo[3]{ 0, 0, 0 }; // Vertex Buffer Objects: For vertex positions, normals, and block type data.
    GLuint ibo = 0; // Index Buffer Object: For indexed drawing to reduce vertex duplication.
    GLsizei indexCount = 0; // The number of indices to draw.
};

// Represents a schematic for a structure that can be placed in the world.
struct Structure {
    int width = 0, depth = 0, height = 0; // Dimensions of the structure's grid.
    std::vector<BlockID> blockTypes; // 1D array storing the structure's block data, mapped to 3D space.
    bool loadFromFile(const std::string& filepath); // Deserializes structure data from a file.
    BlockID getBlock(int x, int y, int z) const; // Accessor for block data using 3D coordinates.
};

// Manages the procedural generation, meshing, and rendering of the voxel environment.
class ProceduralEnvironment {
public:
    ProceduralEnvironment();
    void initialize(); // Sets up initial state, shaders, and generates the first world.
    void updateParameters(float gh, float gr, float mh, float mr, float bd, unsigned int s); // Updates terrain generation variables.
    void regenerate(); // Re-runs the procedural generation and meshing algorithms.
    void setGrassContrast(float contrast); // Adjusts a shader uniform for visual effect.
    void render(const glm::mat4& view, const glm::mat4& projection); // Renders the environment from a given camera perspective.

    BlockID getBlockAtWorldPos(const glm::vec3& worldPos) const; // Converts world-space coordinates to grid coordinates to identify a block.

    // Getters for exposing generation parameters, often to a GUI.
    float getGrassHeight() const { return m_grassHeight; }
    float getGrassRoughness() const { return m_grassRoughness; }
    float getMountainHeight() const { return m_mountainHeight; }
    float getMountainRoughness() const { return m_mountainRoughness; }
    float getBlendDistance() const { return m_blendDistance; }
    unsigned int getSeed() const { return m_seed; }

private:
    // Defines the static dimensions of the world grid.
    static const int WORLD_SIZE_X = 128;
    static const int WORLD_SIZE_Z = 128;
    static const int WORLD_HEIGHT = 64;

    // A flattened 1D vector representing the 3D grid of blocks.
    std::vector<BlockID> m_blockGrid;
    std::vector<BlockMaterial> m_blockMaterials; // Material properties (e.g., color) for each block type.
    VoxelMesh m_voxelMesh; // The OpenGL mesh representation of the world.
    GLuint m_voxelShader; // The shader program used to render the voxels.
    std::unique_ptr<PerlinNoise> m_perlinNoise; // Noise generator for terrain.

    // Terrain generation parameters.
    float m_grassHeight, m_grassRoughness, m_mountainHeight, m_mountainRoughness, m_blendDistance;
    unsigned int m_seed;

    // CPU-side buffers for accumulating mesh data before uploading to the GPU.
    std::vector<glm::vec3> m_voxelVertices;
    std::vector<glm::vec3> m_voxelNormals;
    std::vector<float> m_voxelBlockType;
    std::vector<unsigned int> m_voxelIndices;

    // Perlin noise data duplicated on the CPU for potential GPU-side procedural texturing/effects.
    std::vector<int> m_noisePermutationTable;
    std::vector<glm::vec3> m_noiseGradients;
    float m_grassContrast = 2.0f;

    // Core private methods for world generation and rendering setup.
    void populateBlockGrid(); // Fills m_blockGrid using layered Perlin noise.
    void placeStructures(); // Places pre-defined structures into the block grid.
    void buildMeshFromGrid(); // Traverses m_blockGrid to generate a renderable mesh using greedy meshing or a similar algorithm.
    void createVoxelMesh(); // Allocates and uploads mesh data from CPU vectors to GPU VBOs/IBO.
    void addFaceToMesh(const glm::vec3& pos, BlockID blockType, const std::string& face); // Helper to add a single quad face to the mesh data vectors.
    void createShaders(); // Compiles vertex and fragment shaders and links them into a shader program.
    void uploadBlockMaterials(); // Uploads the m_blockMaterials vector to a shader uniform buffer or texture.
    void uploadNoiseUniforms(); // Uploads noise data to the GPU for shader use.
    GLuint createShaderProgram(const char* vertexSource, const char* fragmentSource); // Utility to compile/link a shader program.

    BlockID getBlock(int x, int y, int z) const; // Helper to access m_blockGrid with bounds checking.
    void setBlock(int x, int y, int z, BlockID type); // Helper to modify m_blockGrid with bounds checking.
};
