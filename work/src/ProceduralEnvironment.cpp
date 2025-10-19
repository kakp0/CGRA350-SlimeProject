#include "ProceduralEnvironment.hpp"
#include "BlockMaterials.hpp"
#include <iostream>
#include <fstream>
#include <map>
#include <random>
#include <glm/gtc/type_ptr.hpp>

// Made with the assistance of AI (Gemini 2.5 Pro)

bool Structure::loadFromFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open structure file: " << filepath << std::endl;
        return false;
    }

    file >> width >> depth >> height;
    if (width == 0 || depth == 0 || height == 0) return false;

    blockTypes.assign(width * depth * height, BlockID::Air); // Pre-allocate a 1D vector to store 3D structure data.

    // Map characters from the file to corresponding block type enums for parsing.
    std::map<char, BlockID> blockMap = {
        {'.', BlockID::Air}, {'c', BlockID::Cobblestone}, {'w', BlockID::Wood}, {'g', BlockID::Glass}
    };

    // Iterate through the file, reading character by character to populate the block data.
    for (int y = 0; y < height; ++y) {
        for (int z = 0; z < depth; ++z) {
            for (int x = 0; x < width; ++x) {
                char blockChar;
                file >> blockChar;
                if (blockMap.count(blockChar)) {
                    // Convert 3D coordinates (x, y, z) to a 1D array index.
                    blockTypes[y * (width * depth) + z * width + x] = blockMap[blockChar];
                }
            }
        }
    }
    return true;
}

BlockID Structure::getBlock(int x, int y, int z) const {
    // Perform bounds checking to prevent out-of-range access.
    if (x < 0 || x >= width || y < 0 || y >= height || z < 0 || z >= depth) {
        return BlockID::Air; // Treat out-of-bounds as empty space.
    }
    // Convert 3D coordinates to a 1D index to access the flat vector.
    return blockTypes[y * (width * depth) + z * width + x];
}

ProceduralEnvironment::ProceduralEnvironment()
    : m_grassHeight(5.0f), m_grassRoughness(0.5f), m_mountainHeight(30.0f), m_mountainRoughness(1.0f), m_blendDistance(30.0f), m_seed(42) {
    m_perlinNoise = std::make_unique<PerlinNoise>(m_seed); // Initialize Perlin noise generator with a default seed.
}

void ProceduralEnvironment::initialize() {
    m_noisePermutationTable = m_perlinNoise->getPermutationTable(); // Cache noise data for performance.
    m_noiseGradients = m_perlinNoise->getGradients();
    BlockMaterialRegistry::createAllMaterials(m_blockMaterials); // Define the properties (colors, etc.) for each block type.
    createShaders(); // Compile and link the GLSL shaders.
    regenerate(); // Perform the initial world generation.
}

void ProceduralEnvironment::updateParameters(float gh, float gr, float mh, float mr, float bd, unsigned int s) {
    m_grassHeight = gh;
    m_grassRoughness = gr;
    m_mountainHeight = mh;
    m_mountainRoughness = mr;
    m_blendDistance = bd;

    // If the seed has changed, the entire noise function must be re-initialized.
    if (m_seed != s) {
        m_seed = s;
        m_perlinNoise = std::make_unique<PerlinNoise>(m_seed);
        m_noisePermutationTable = m_perlinNoise->getPermutationTable();
        m_noiseGradients = m_perlinNoise->getGradients();
    }
}

void ProceduralEnvironment::setGrassContrast(float contrast) {
    m_grassContrast = contrast;
}

void ProceduralEnvironment::regenerate() {
    populateBlockGrid(); // Generate the raw voxel data.
    placeStructures(); // Add pre-designed structures to the voxel grid.
    buildMeshFromGrid(); // Convert the voxel grid into a renderable mesh, culling hidden faces.
    createVoxelMesh(); // Upload the generated mesh data to GPU buffers (VAO, VBOs).
}

BlockID ProceduralEnvironment::getBlock(int x, int y, int z) const {
    if (x < 0 || x >= WORLD_SIZE_X || y < 0 || y >= WORLD_HEIGHT || z < 0 || z >= WORLD_SIZE_Z) {
        return BlockID::Air;
    }
    return m_blockGrid[y * (WORLD_SIZE_X * WORLD_SIZE_Z) + z * WORLD_SIZE_X + x];
}

void ProceduralEnvironment::setBlock(int x, int y, int z, BlockID type) {
    if (x < 0 || x >= WORLD_SIZE_X || y < 0 || y >= WORLD_HEIGHT || z < 0 || z >= WORLD_SIZE_Z) {
        return;
    }
    m_blockGrid[y * (WORLD_SIZE_X * WORLD_SIZE_Z) + z * WORLD_SIZE_X + x] = type;
}

void ProceduralEnvironment::populateBlockGrid() {
    m_blockGrid.assign(WORLD_SIZE_X * WORLD_SIZE_Z * WORLD_HEIGHT, BlockID::Air);
    float halfSizeX = WORLD_SIZE_X * 0.5f;
    float halfSizeZ = WORLD_SIZE_Z * 0.5f;

    const float caveFrequency = 0.08f;
    const float caveThreshold = 0.6f;

    for (int x = 0; x < WORLD_SIZE_X; ++x) {
        for (int z = 0; z < WORLD_SIZE_Z; ++z) {
            float worldX = static_cast<float>(x) - halfSizeX;
            float worldZ = static_cast<float>(z) - halfSizeZ;

            // Combine multiple noise functions (octaves) for more interesting terrain.
            float baseHeight = m_perlinNoise->terrainHeight(worldX, worldZ, m_grassHeight, m_grassRoughness);
            float surfaceDetail = m_perlinNoise->noise(worldX * 8.0f, worldZ * 8.0f) * 1.5f; // Higher frequency noise for small details.
            float mountainVal = m_perlinNoise->mountainHeight(worldX, worldZ, m_mountainHeight, m_mountainRoughness);

            float distanceFromCenter = glm::length(glm::vec2(worldX, worldZ));
            // Use smoothstep to create a gradual blend from plains to mountains near the world edge.
            float blendFactor = glm::smoothstep(halfSizeX - m_blendDistance, halfSizeX, distanceFromCenter);

            int surfaceHeight = static_cast<int>(baseHeight + surfaceDetail + mountainVal * blendFactor);
            if (surfaceHeight >= WORLD_HEIGHT) surfaceHeight = WORLD_HEIGHT - 1;

            for (int y = 0; y < surfaceHeight; ++y) {
                if (y == 0) { setBlock(x, y, z, BlockID::Bedrock); continue; }

                // Use 3D Perlin noise to carve out caves.
                float caveNoise = m_perlinNoise->noise(worldX * caveFrequency, y * caveFrequency, worldZ * caveFrequency);
                caveNoise += ((10.0f - y) / 10.0f) * 0.1f; // Bias noise to make caves more common near the surface.
                if (caveNoise > caveThreshold) continue; // If noise value is above threshold, carve an air block.

                int depth = surfaceHeight - 1 - y;
                BlockID blockType = BlockID::Stone;
                // Layer block types based on depth from the surface.
                if (depth == 0) { // Top layer
                    if (surfaceHeight > 25) blockType = BlockID::Snow; else blockType = BlockID::Grass;
                }
                else if (depth <= 3) { // Sub-surface
                    blockType = BlockID::Dirt;
                }
                setBlock(x, y, z, blockType);
            }
        }
    }
}

void ProceduralEnvironment::placeStructures() {
    std::string currentFilePath = __FILE__;
    std::string currentDir = currentFilePath.substr(0, currentFilePath.find_last_of("/\\"));
    std::string housePath = currentDir + "/small_house.txt";

    Structure house;
    if (!house.loadFromFile(housePath)) return;

    // Use a 2D grid to track occupied areas, preventing structures from spawning too close together.
    std::vector<bool> occupied(WORLD_SIZE_X * WORLD_SIZE_Z, false);
    const int spacingBuffer = 45; // Minimum distance between structures.
    std::mt19937 rng(m_seed);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    const float spawnChance = 0.02f;

    for (int z = 10; z < WORLD_SIZE_Z - 10 - house.depth; ++z) {
        for (int x = 10; x < WORLD_SIZE_X - 10 - house.width; ++x) {
            if (occupied[z * WORLD_SIZE_X + x]) continue;

            // Find the ground height across the structure's footprint to check for flatness.
            int minY = WORLD_HEIGHT; int maxY = 0;
            for (int dz = 0; dz < house.depth; ++dz) {
                for (int dx = 0; dx < house.width; ++dx) {
                    int groundY = 0;
                    for (int y = WORLD_HEIGHT - 1; y >= 0; --y) {
                        if (getBlock(x + dx, y, z + dz) != BlockID::Air) { groundY = y; break; }
                    }
                    if (groundY < minY) minY = groundY;
                    if (groundY > maxY) maxY = groundY;
                }
            }

            // If the area is reasonably flat and passes a random chance check, place the structure.
            if (maxY - minY <= 2 && dist(rng) < spawnChance) {
                // Mark a buffered area around the new structure as occupied.
                for (int dz = -spacingBuffer; dz < house.depth + spacingBuffer; ++dz) {
                    for (int dx = -spacingBuffer; dx < house.width + spacingBuffer; ++dx) {
                        int markX = x + dx; int markZ = z + dz;
                        if (markX >= 0 && markX < WORLD_SIZE_X && markZ >= 0 && markZ < WORLD_SIZE_Z) {
                            occupied[markZ * WORLD_SIZE_X + markX] = true;
                        }
                    }
                }

                // Flatten the ground to create a foundation.
                for (int dz = 0; dz < house.depth; ++dz) {
                    for (int dx = 0; dx < house.width; ++dx) {
                        for (int y = minY - 2; y < minY; ++y) setBlock(x + dx, y, z + dz, BlockID::Dirt);
                        setBlock(x + dx, minY, z + dz, BlockID::Grass);
                    }
                }

                // Copy the block data from the Structure object into the world grid.
                for (int y = 0; y < house.height; ++y) {
                    for (int z_local = 0; z_local < house.depth; ++z_local) {
                        for (int x_local = 0; x_local < house.width; ++x_local) {
                            BlockID blockType = house.getBlock(x_local, y, z_local);
                            if (blockType != BlockID::Air) {
                                setBlock(x + x_local, minY + 1 + y, z + z_local, blockType);
                            }
                        }
                    }
                }
            }
        }
    }
}

void ProceduralEnvironment::buildMeshFromGrid() {
    m_voxelVertices.clear();
    m_voxelNormals.clear();
    m_voxelBlockType.clear();
    m_voxelIndices.clear();

    for (int y = 0; y < WORLD_HEIGHT; ++y) {
        for (int z = 0; z < WORLD_SIZE_Z; ++z) {
            for (int x = 0; x < WORLD_SIZE_X; ++x) {
                BlockID currentBlock = getBlock(x, y, z);
                if (currentBlock == BlockID::Air) continue;

                // Center the world origin at (0, y, 0) for easier camera manipulation.
                glm::vec3 worldPos = { static_cast<float>(x) - WORLD_SIZE_X * 0.5f, static_cast<float>(y), static_cast<float>(z) - WORLD_SIZE_Z * 0.5f };

                // This is a critical optimization: only generate a face if it's exposed to air.
                if (getBlock(x, y, z + 1) == BlockID::Air) addFaceToMesh(worldPos, currentBlock, "front");
                if (getBlock(x, y, z - 1) == BlockID::Air) addFaceToMesh(worldPos, currentBlock, "back");
                if (getBlock(x + 1, y, z) == BlockID::Air) addFaceToMesh(worldPos, currentBlock, "right");
                if (getBlock(x - 1, y, z) == BlockID::Air) addFaceToMesh(worldPos, currentBlock, "left");
                if (getBlock(x, y + 1, z) == BlockID::Air) addFaceToMesh(worldPos, currentBlock, "top");
                if (getBlock(x, y - 1, z) == BlockID::Air) addFaceToMesh(worldPos, currentBlock, "bottom");
            }
        }
    }
}

void ProceduralEnvironment::addFaceToMesh(const glm::vec3& pos, BlockID blockType, const std::string& face) {
    unsigned int baseIndex = m_voxelVertices.size();
    glm::vec3 n;
    std::vector<glm::vec3> v(4);

    // Define the 8 vertices of a unit cube centered at the origin.
    glm::vec3 p0 = { pos.x - 0.5f, pos.y - 0.5f, pos.z + 0.5f };
    glm::vec3 p1 = { pos.x + 0.5f, pos.y - 0.5f, pos.z + 0.5f };
    glm::vec3 p2 = { pos.x + 0.5f, pos.y + 0.5f, pos.z + 0.5f };
    glm::vec3 p3 = { pos.x - 0.5f, pos.y + 0.5f, pos.z + 0.5f };
    glm::vec3 p4 = { pos.x - 0.5f, pos.y - 0.5f, pos.z - 0.5f };
    glm::vec3 p5 = { pos.x + 0.5f, pos.y - 0.5f, pos.z - 0.5f };
    glm::vec3 p6 = { pos.x + 0.5f, pos.y + 0.5f, pos.z - 0.5f };
    glm::vec3 p7 = { pos.x - 0.5f, pos.y + 0.5f, pos.z - 0.5f };

    // Select the four correct vertices and the normal vector based on the face direction.
    if (face == "front") { v = { p0, p1, p2, p3 }; n = { 0,0,1 }; }
    if (face == "back") { v = { p5, p4, p7, p6 }; n = { 0,0,-1 }; }
    if (face == "right") { v = { p1, p5, p6, p2 }; n = { 1,0,0 }; }
    if (face == "left") { v = { p4, p0, p3, p7 }; n = { -1,0,0 }; }
    if (face == "top") { v = { p3, p2, p6, p7 }; n = { 0,1,0 }; }
    if (face == "bottom") { v = { p4, p5, p1, p0 }; n = { 0,-1,0 }; }

    m_voxelVertices.insert(m_voxelVertices.end(), v.begin(), v.end());
    for (int i = 0; i < 4; ++i) {
        m_voxelNormals.push_back(n);
        // Pass block type as a float vertex attribute to be used by the fragment shader.
        m_voxelBlockType.push_back(static_cast<float>(blockType));
    }

    // A quad is made of two triangles, defined by 6 indices.
    m_voxelIndices.insert(m_voxelIndices.end(), { baseIndex, baseIndex + 1, baseIndex + 2, baseIndex, baseIndex + 2, baseIndex + 3 });
}

void ProceduralEnvironment::createVoxelMesh() {
    // If a mesh already exists, delete old GPU buffers to prevent memory leaks.
    if (m_voxelMesh.vao) {
        glDeleteVertexArrays(1, &m_voxelMesh.vao);
        glDeleteBuffers(3, m_voxelMesh.vbo);
        glDeleteBuffers(1, &m_voxelMesh.ibo);
    }
    glGenVertexArrays(1, &m_voxelMesh.vao);
    glBindVertexArray(m_voxelMesh.vao);
    glGenBuffers(3, m_voxelMesh.vbo);

    // VBO for vertex positions
    glBindBuffer(GL_ARRAY_BUFFER, m_voxelMesh.vbo[0]);
    glBufferData(GL_ARRAY_BUFFER, m_voxelVertices.size() * sizeof(glm::vec3), m_voxelVertices.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr); // layout (location = 0)

    // VBO for normal vectors
    glBindBuffer(GL_ARRAY_BUFFER, m_voxelMesh.vbo[1]);
    glBufferData(GL_ARRAY_BUFFER, m_voxelNormals.size() * sizeof(glm::vec3), m_voxelNormals.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr); // layout (location = 1)

    // VBO for block type ID
    glBindBuffer(GL_ARRAY_BUFFER, m_voxelMesh.vbo[2]);
    glBufferData(GL_ARRAY_BUFFER, m_voxelBlockType.size() * sizeof(float), m_voxelBlockType.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, nullptr); // layout (location = 2)

    // Index Buffer Object (IBO) for drawing order
    glGenBuffers(1, &m_voxelMesh.ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_voxelMesh.ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_voxelIndices.size() * sizeof(unsigned int), m_voxelIndices.data(), GL_STATIC_DRAW);

    m_voxelMesh.indexCount = static_cast<GLsizei>(m_voxelIndices.size());
    glBindVertexArray(0); // Unbind VAO.
}

void ProceduralEnvironment::createShaders() {
    // The vertex shader's main job is to apply model-view-projection transformations.
    // It also passes per-vertex data like position, normal, and block type to the fragment shader.
    const char* vertexSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;
        layout (location = 2) in float aBlockType;

        out vec3 v_worldPos;
        out vec3 v_worldNormal;
        out float v_blockType;

        uniform mat4 uProjectionMatrix;
        uniform mat4 uModelViewMatrix;

        void main() {
            v_worldPos = aPos; // Pass world-space position.
            v_worldNormal = aNormal; // Pass normal.
            v_blockType = aBlockType; // Pass block type.
            gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(aPos, 1.0);
        }
    )";

    // The fragment shader calculates the final color for each pixel.
    // It uses the block type to select a material and then runs a procedural texturing algorithm.
    const char* fragmentSource = R"(
        #version 330 core
        out vec4 FragColor;

        in vec3 v_worldPos;
        in vec3 v_worldNormal;
        in float v_blockType;

        struct Material {
            vec3 baseColor;
            float alpha;
            vec3 detailColor;
            vec4 materialParams; // Generic parameters for shader functions (e.g., frequency, density).
        };

        const int MAX_BLOCK_TYPES = 16;
        uniform Material uMaterials[MAX_BLOCK_TYPES];
        uniform vec3 uLightDirection;

        // Noise data uploaded from the CPU.
        uniform int p[512];
        uniform vec3 gradients[16];
        uniform float uGrassContrast;

        // --- GLSL Noise Function Implementations ---
        // Simple hash-based random number generator.
        float valueRand(vec2 co){ return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453); }

        // Value noise: generates random values and interpolates between them.
        float valueNoise(vec2 st) {
            vec2 i = floor(st);
            vec2 f = fract(st);
            float a = valueRand(i); float b = valueRand(i + vec2(1.0, 0.0));
            float c = valueRand(i + vec2(0.0, 1.0)); float d = valueRand(i + vec2(1.0, 1.0));
            vec2 u = f * f * (3.0 - 2.0 * f); // Smoothstep interpolation.
            return mix(a, b, u.x) + (c - a)* u.y * (1.0 - u.x) + (d - b) * u.y * u.x;
        }

        // --- Perlin Noise Functions (from Ken Perlin) ---
        float fade(float t) { return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); }
        vec2 fade(vec2 t) { return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); }
        int hash(int x, int y) { return p[p[x & 255] + (y & 255)] & 15; }
        float dotGridGradient(ivec2 grid, vec2 pos) {
            vec2 gradient = gradients[hash(grid.x, grid.y)].xy;
            vec2 dist = pos - vec2(grid);
            return dot(dist, gradient);
        }
        float perlinNoise(vec2 pos) {
            ivec2 i = ivec2(floor(pos));
            vec2 f = fract(pos);
            vec2 u = fade(f);
            float g00 = dotGridGradient(i, pos);
            float g10 = dotGridGradient(i + ivec2(1,0), pos);
            float g01 = dotGridGradient(i + ivec2(0,1), pos);
            float g11 = dotGridGradient(i + ivec2(1,1), pos);
            return mix(mix(g00, g10, u.x), mix(g01, g11, u.x), u.y);
        }

        // --- Block-Specific Procedural Coloring Functions ---
        vec4 getSimpleColor(Material mat, vec3 worldPos) {
            float n = valueNoise(worldPos.xz * mat.materialParams.x);
            vec3 color = mix(mat.baseColor, mat.detailColor, n * 0.4);
            return vec4(color, mat.alpha);
        }
        
        vec4 getGrassColor(Material mat, vec3 worldPos) {
            // Use multiple layers of noise for varied grass color.
            float regionValue = (perlinNoise(worldPos.xz * 0.015) + 1.0) * 0.5; // Low-frequency noise for large regions.
            vec3 darkRegionColor = vec3(0.08, 0.35, 0.1);
            vec3 regionalBaseColor = mix(mat.baseColor, darkRegionColor, pow(regionValue, uGrassContrast));
            float patchNoise = valueNoise(worldPos.xz * mat.materialParams.x); // High-frequency noise for small patches.
            vec3 finalColor = regionalBaseColor;
            if (patchNoise > 0.65) {
                float patchIntensity = (patchNoise - 0.65) / (1.0 - 0.65);
                finalColor = mix(regionalBaseColor, mat.detailColor, patchIntensity * 0.8);
            }
            return vec4(finalColor, mat.alpha);
        }

        vec4 getDirtColor(Material mat, vec3 worldPos) {
            float n = valueRand(worldPos.xz * mat.materialParams.x);
            vec3 color = mat.baseColor;
            if (n > 0.85) { color = mix(color, mat.detailColor, 0.5); } // Add random specks.
            return vec4(color, mat.alpha);
        }

        vec4 getWoodColor(Material mat, vec3 worldPos, vec3 normal) {
            float grainFrequency = mat.materialParams.x;
            float plankHeight = mat.materialParams.y;
            float gapSharpness = mat.materialParams.z;
            vec2 uv;
            // Project texture coordinates based on face normal to avoid stretching (tri-planar mapping).
            if (abs(normal.y) > 0.5) uv = worldPos.xz;
            else if (abs(normal.x) > 0.5) uv = worldPos.zy;
            else uv = worldPos.xy;
            float plankLine = abs(fract(uv.y / plankHeight - 0.5) * 2.0 - 1.0);
            float plankGap = smoothstep(gapSharpness, 1.0, plankLine); // Create soft gaps between planks.
            float grain = valueNoise(uv * vec2(1.0, grainFrequency)); // Wood grain effect.
            grain = pow(grain, 2.0) * 0.3;
            vec3 color = mat.baseColor - grain;
            color = mix(mat.detailColor, color, plankGap); // Darken the gaps.
            return vec4(color, mat.alpha);
        }

        vec4 getCobblestoneColor(Material mat, vec3 worldPos, vec3 normal) {
            float cellDensity = mat.materialParams.x;
            float groutWidth = mat.materialParams.z;
            vec2 p;
            if(abs(normal.y) > 0.5) p = worldPos.xz; else if(abs(normal.x) > 0.5) p = worldPos.zy; else p = worldPos.xy;

            // Use Worley (cellular) noise to generate cobblestone shapes.
            vec2 i = floor(p * cellDensity);
            vec2 f = fract(p * cellDensity);
            float min_dist = 1.0;
            // Find distance to the closest point in a 3x3 grid of cells.
            for (int y = -1; y <= 1; y++) {
                for (int x = -1; x <= 1; x++) {
                    vec2 neighbor_cell = vec2(float(x), float(y));
                    vec2 point_in_cell = vec2(valueRand(i + neighbor_cell), valueRand(i + neighbor_cell + vec2(5.2, 1.3)));
                    vec2 diff = neighbor_cell + point_in_cell - f;
                    min_dist = min(min_dist, length(diff));
                }
            }
            vec3 color = mat.baseColor + valueRand(i) * 0.05; // Vary stone color per-cell.
            // Use smoothstep with the distance field to draw the grout lines.
            color = mix(mat.detailColor, color, smoothstep(groutWidth, groutWidth + 0.1, min_dist));
            return vec4(color, mat.alpha);
        }

        void main() {
            int blockIndex = int(round(v_blockType));
            if (blockIndex <= 0 || blockIndex >= MAX_BLOCK_TYPES) { discard; }

            Material mat = uMaterials[blockIndex];
            vec4 surfaceColor;

            // Select the appropriate coloring function based on the block type.
            if (blockIndex == 1) surfaceColor = getGrassColor(mat, v_worldPos);
            else if (blockIndex == 2) surfaceColor = getDirtColor(mat, v_worldPos);
            else if (blockIndex == 3) surfaceColor = getSimpleColor(mat, v_worldPos);
            else if (blockIndex == 5) surfaceColor = getSimpleColor(mat, v_worldPos);
            else if (blockIndex == 6) surfaceColor = getWoodColor(mat, v_worldPos, v_worldNormal);
            else if (blockIndex == 7) surfaceColor = getCobblestoneColor(mat, v_worldPos, v_worldNormal);
            else surfaceColor = vec4(mat.baseColor, mat.alpha); // Fallback for unhandled types.

            // --- Basic Lighting Model ---
            vec3 normal = normalize(v_worldNormal);
            // Fake ambient occlusion: top-facing surfaces are brighter.
            float ao = 0.7 + (max(0.0, normal.y) * 0.3);
            vec3 lightDir = normalize(uLightDirection);
            // Lambertian diffuse reflection.
            float diff = max(dot(normal, lightDir), 0.0);
            vec3 lighting = (vec3(0.4) + diff * vec3(1.0)) * ao; // Ambient + Diffuse
            
            vec3 finalColor = surfaceColor.rgb * lighting;
            FragColor = vec4(finalColor, surfaceColor.a);
        }
    )";
    m_voxelShader = createShaderProgram(vertexSource, fragmentSource);
}

void ProceduralEnvironment::uploadBlockMaterials() {
    glUseProgram(m_voxelShader);
    // Loop through all defined materials and upload their properties to the shader's uniform array.
    for (size_t i = 0; i < m_blockMaterials.size(); ++i) {
        std::string base = "uMaterials[" + std::to_string(i) + "].";
        glUniform3fv(glGetUniformLocation(m_voxelShader, (base + "baseColor").c_str()), 1, &m_blockMaterials[i].baseColor[0]);
        glUniform1f(glGetUniformLocation(m_voxelShader, (base + "alpha").c_str()), m_blockMaterials[i].alpha);
        glUniform3fv(glGetUniformLocation(m_voxelShader, (base + "detailColor").c_str()), 1, &m_blockMaterials[i].detailColor[0]);
        glUniform4fv(glGetUniformLocation(m_voxelShader, (base + "materialParams").c_str()), 1, &m_blockMaterials[i].materialParams[0]);
    }
}

void ProceduralEnvironment::uploadNoiseUniforms() {
    glUseProgram(m_voxelShader);
    // The permutation table and gradient vectors for Perlin noise must be available to the GPU.
    GLint pLoc = glGetUniformLocation(m_voxelShader, "p");
    if (pLoc != -1) glUniform1iv(pLoc, m_noisePermutationTable.size(), m_noisePermutationTable.data());
    GLint gradLoc = glGetUniformLocation(m_voxelShader, "gradients");
    if (gradLoc != -1) glUniform3fv(gradLoc, m_noiseGradients.size(), glm::value_ptr(m_noiseGradients[0]));
}

void ProceduralEnvironment::render(const glm::mat4& view, const glm::mat4& projection) {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE); // Cull back-faces for performance.
    glUseProgram(m_voxelShader);
    uploadBlockMaterials(); // Ensure material data is up-to-date.
    uploadNoiseUniforms();

    glUniform1f(glGetUniformLocation(m_voxelShader, "uGrassContrast"), m_grassContrast);

    // Set up matrices and other global uniforms.
    glm::vec3 lightDirection = glm::normalize(glm::vec3(0.5f, 1.0f, 0.3f));
    glUniformMatrix4fv(glGetUniformLocation(m_voxelShader, "uProjectionMatrix"), 1, GL_FALSE, &projection[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(m_voxelShader, "uModelViewMatrix"), 1, GL_FALSE, &view[0][0]);
    glUniform3fv(glGetUniformLocation(m_voxelShader, "uLightDirection"), 1, &lightDirection[0]);

    glBindVertexArray(m_voxelMesh.vao);
    if (m_voxelMesh.indexCount > 0) {
        glDrawElements(GL_TRIANGLES, m_voxelMesh.indexCount, GL_UNSIGNED_INT, 0); // Issue the draw call.
    }
    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_CULL_FACE);
}

GLuint ProceduralEnvironment::createShaderProgram(const char* vertexSource, const char* fragmentSource) {
    // Standard OpenGL shader compilation and linking boilerplate.
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, nullptr);
    glCompileShader(vertexShader);
    GLint success;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) { GLchar infoLog[512]; glGetShaderInfoLog(vertexShader, 512, nullptr, infoLog); std::cerr << "Vertex shader compilation failed: " << infoLog << std::endl; }

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, nullptr);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) { GLchar infoLog[512]; glGetShaderInfoLog(fragmentShader, 512, nullptr, infoLog); std::cerr << "Fragment shader compilation failed: " << infoLog << std::endl; }

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) { GLchar infoLog[512]; glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog); std::cerr << "Shader program linking failed: " << infoLog << std::endl; }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return shaderProgram;
}

BlockID ProceduralEnvironment::getBlockAtWorldPos(const glm::vec3& worldPos) const {
    // Convert from continuous world-space coordinates back to discrete grid indices.
    // This is the inverse of the transformation used in buildMeshFromGrid().
    int gridX = static_cast<int>(floor(worldPos.x + WORLD_SIZE_X * 0.5f));
    int gridY = static_cast<int>(floor(worldPos.y));
    int gridZ = static_cast<int>(floor(worldPos.z + WORLD_SIZE_Z * 0.5f));

    return getBlock(gridX, gridY, gridZ);
}
