#include "ProceduralEnvironment.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <string> 
#include <random> // Required for random placement

// Block type IDs
enum BlockID {
    Air = 0, Grass = 1, Dirt = 2, Stone = 3, Snow = 4, Bedrock = 5,
    Wood = 6, Cobblestone = 7, Glass = 8
};

// --- Structure Loading ---
bool Structure::loadFromFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open structure file: " << filepath << std::endl;
        return false;
    }

    file >> width >> depth >> height;
    if (width == 0 || depth == 0 || height == 0) return false;

    blockTypes.assign(width * depth * height, Air);

    std::map<char, int> blockMap = {
        {'.', Air}, {'c', Cobblestone}, {'w', Wood}, {'g', Glass}
    };

    for (int y = 0; y < height; ++y) {
        for (int z = 0; z < depth; ++z) {
            for (int x = 0; x < width; ++x) {
                char blockChar;
                file >> blockChar;
                if (blockMap.count(blockChar)) {
                    blockTypes[y * (width * depth) + z * width + x] = blockMap[blockChar];
                }
            }
        }
    }
    return true;
}

int Structure::getBlock(int x, int y, int z) const {
    if (x < 0 || x >= width || y < 0 || y >= height || z < 0 || z >= depth) return Air;
    return blockTypes[y * (width * depth) + z * width + x];
}

ProceduralEnvironment::ProceduralEnvironment()
    : m_grassHeight(5.0f), m_grassRoughness(0.5f),
    m_mountainHeight(30.0f), m_mountainRoughness(1.0f),
    m_blendDistance(30.0f), m_seed(42) {
    m_perlinNoise = std::make_unique<PerlinNoise>(m_seed);
}

void ProceduralEnvironment::initialize() {
    createShaders();
    regenerate();
}

void ProceduralEnvironment::updateParameters(float gh, float gr, float mh, float mr, float bd, unsigned int s) {
    m_grassHeight = gh; m_grassRoughness = gr; m_mountainHeight = mh; m_mountainRoughness = mr; m_blendDistance = bd;
    if (m_seed != s) {
        m_seed = s;
        m_perlinNoise = std::make_unique<PerlinNoise>(m_seed);
    }
}

void ProceduralEnvironment::regenerate() {
    populateBlockGrid();
    placeStructures();
    buildMeshFromGrid();
    createVoxelMesh();
}

int ProceduralEnvironment::getBlock(int x, int y, int z) const {
    if (x < 0 || x >= WORLD_SIZE_X || y < 0 || y >= WORLD_HEIGHT || z < 0 || z >= WORLD_SIZE_Z) {
        return Air;
    }
    return m_blockGrid[y * (WORLD_SIZE_X * WORLD_SIZE_Z) + z * WORLD_SIZE_X + x];
}

void ProceduralEnvironment::setBlock(int x, int y, int z, int type) {
    if (x < 0 || x >= WORLD_SIZE_X || y < 0 || y >= WORLD_HEIGHT || z < 0 || z >= WORLD_SIZE_Z) {
        return;
    }
    m_blockGrid[y * (WORLD_SIZE_X * WORLD_SIZE_Z) + z * WORLD_SIZE_X + x] = type;
}

void ProceduralEnvironment::populateBlockGrid() {
    m_blockGrid.assign(WORLD_SIZE_X * WORLD_SIZE_Z * WORLD_HEIGHT, Air);

    float halfSizeX = WORLD_SIZE_X * 0.5f;
    float halfSizeZ = WORLD_SIZE_Z * 0.5f;

    const float caveFrequency = 0.08f;
    const float caveThreshold = 0.6f;

    for (int x = 0; x < WORLD_SIZE_X; ++x) {
        for (int z = 0; z < WORLD_SIZE_Z; ++z) {
            float worldX = static_cast<float>(x) - halfSizeX;
            float worldZ = static_cast<float>(z) - halfSizeZ;

            float baseHeight = m_perlinNoise->terrainHeight(worldX, worldZ, m_grassHeight, m_grassRoughness);
            float surfaceDetail = m_perlinNoise->noise(worldX * 8.0f, worldZ * 8.0f) * 1.5f;
            float mountainVal = m_perlinNoise->mountainHeight(worldX, worldZ, m_mountainHeight, m_mountainRoughness);
            float distanceFromCenter = glm::length(glm::vec2(worldX, worldZ));
            float blendFactor = glm::smoothstep(halfSizeX - m_blendDistance, halfSizeX, distanceFromCenter);
            int surfaceHeight = static_cast<int>(baseHeight + surfaceDetail + mountainVal * blendFactor);
            if (surfaceHeight >= WORLD_HEIGHT) surfaceHeight = WORLD_HEIGHT - 1;

            for (int y = 0; y < surfaceHeight; ++y) {
                if (y == 0) {
                    setBlock(x, y, z, Bedrock);
                    continue;
                }

                float caveNoise = m_perlinNoise->noise(worldX * caveFrequency, y * caveFrequency, worldZ * caveFrequency);
                caveNoise += ((10.0f - y) / 10.0f) * 0.1f;
                if (caveNoise > caveThreshold) {
                    continue;
                }

                int depth = surfaceHeight - 1 - y;
                int blockType = Stone;
                if (depth == 0) {
                    if (surfaceHeight > 25) blockType = Snow;
                    else blockType = Grass;
                }
                else if (depth <= 3) {
                    blockType = Dirt;
                }
                setBlock(x, y, z, blockType);
            }
        }
    }
}


void ProceduralEnvironment::placeStructures() {
    // Keep the __FILE__ method for finding the asset file
    std::string currentFilePath = __FILE__;
    std::string currentDir = currentFilePath.substr(0, currentFilePath.find_last_of("/\\"));
    std::string housePath = currentDir + "/small_house.txt";

    Structure house;
    if (!house.loadFromFile(housePath)) {
        return;
    }

    // --- MERGED LOGIC FOR RANDOM AND SPACED PLACEMENT ---

    std::vector<bool> occupied(WORLD_SIZE_X * WORLD_SIZE_Z, false);
    const int spacingBuffer = 45;

    std::mt19937 rng(m_seed);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    const float spawnChance = 0.02f;

    for (int z = 10; z < WORLD_SIZE_Z - 10 - house.depth; ++z) {
        for (int x = 10; x < WORLD_SIZE_X - 10 - house.width; ++x) {

            if (occupied[z * WORLD_SIZE_X + x]) {
                continue;
            }

            int minY = WORLD_HEIGHT;
            int maxY = 0;

            for (int dz = 0; dz < house.depth; ++dz) {
                for (int dx = 0; dx < house.width; ++dx) {
                    int groundY = 0;
                    for (int y = WORLD_HEIGHT - 1; y >= 0; --y) {
                        if (getBlock(x + dx, y, z + dz) != Air) {
                            groundY = y;
                            break;
                        }
                    }
                    if (groundY < minY) minY = groundY;
                    if (groundY > maxY) maxY = groundY;
                }
            }

            if (maxY - minY <= 2 && dist(rng) < spawnChance) {

                std::cout << "Found a flat spot. Placing structure at grid location ("
                    << x << ", " << minY + 1 << ", " << z << ")." << std::endl;

                // Mark the area as occupied BEFORE placing the house
                for (int dz = -spacingBuffer; dz < house.depth + spacingBuffer; ++dz) {
                    for (int dx = -spacingBuffer; dx < house.width + spacingBuffer; ++dx) {
                        int markX = x + dx;
                        int markZ = z + dz;
                        if (markX >= 0 && markX < WORLD_SIZE_X && markZ >= 0 && markZ < WORLD_SIZE_Z) {
                            occupied[markZ * WORLD_SIZE_X + markX] = true;
                        }
                    }
                }

                // Build foundation and place structure
                for (int dz = 0; dz < house.depth; ++dz) {
                    for (int dx = 0; dx < house.width; ++dx) {
                        for (int y = minY - 2; y < minY; ++y) {
                            setBlock(x + dx, y, z + dz, Dirt);
                        }
                        setBlock(x + dx, minY, z + dz, Grass);
                    }
                }

                for (int y = 0; y < house.height; ++y) {
                    for (int z_local = 0; z_local < house.depth; ++z_local) {
                        for (int x_local = 0; x_local < house.width; ++x_local) {
                            int blockType = house.getBlock(x_local, y, z_local);
                            if (blockType != Air) {
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
                int currentBlock = getBlock(x, y, z);
                if (currentBlock == Air) continue;

                glm::vec3 worldPos = {
                    static_cast<float>(x) - WORLD_SIZE_X * 0.5f,
                    static_cast<float>(y),
                    static_cast<float>(z) - WORLD_SIZE_Z * 0.5f
                };

                if (getBlock(x, y, z + 1) == Air) addFaceToMesh(worldPos, currentBlock, "front");
                if (getBlock(x, y, z - 1) == Air) addFaceToMesh(worldPos, currentBlock, "back");
                if (getBlock(x + 1, y, z) == Air) addFaceToMesh(worldPos, currentBlock, "right");
                if (getBlock(x - 1, y, z) == Air) addFaceToMesh(worldPos, currentBlock, "left");
                if (getBlock(x, y + 1, z) == Air) addFaceToMesh(worldPos, currentBlock, "top");
                if (getBlock(x, y - 1, z) == Air) addFaceToMesh(worldPos, currentBlock, "bottom");
            }
        }
    }
}

void ProceduralEnvironment::addFaceToMesh(const glm::vec3& pos, float blockType, const std::string& face) {
    unsigned int baseIndex = m_voxelVertices.size();
    glm::vec3 n;
    std::vector<glm::vec3> v(4);

    glm::vec3 p0 = { pos.x - 0.5f, pos.y - 0.5f, pos.z + 0.5f }; glm::vec3 p1 = { pos.x + 0.5f, pos.y - 0.5f, pos.z + 0.5f };
    glm::vec3 p2 = { pos.x + 0.5f, pos.y + 0.5f, pos.z + 0.5f }; glm::vec3 p3 = { pos.x - 0.5f, pos.y + 0.5f, pos.z + 0.5f };
    glm::vec3 p4 = { pos.x - 0.5f, pos.y - 0.5f, pos.z - 0.5f }; glm::vec3 p5 = { pos.x + 0.5f, pos.y - 0.5f, pos.z - 0.5f };
    glm::vec3 p6 = { pos.x + 0.5f, pos.y + 0.5f, pos.z - 0.5f }; glm::vec3 p7 = { pos.x - 0.5f, pos.y + 0.5f, pos.z - 0.5f };

    if (face == "front") { v = { p0, p1, p2, p3 }; n = { 0,0,1 }; }
    if (face == "back") { v = { p5, p4, p7, p6 }; n = { 0,0,-1 }; }
    if (face == "right") { v = { p1, p5, p6, p2 }; n = { 1,0,0 }; }
    if (face == "left") { v = { p4, p0, p3, p7 }; n = { -1,0,0 }; }
    if (face == "top") { v = { p3, p2, p6, p7 }; n = { 0,1,0 }; }
    if (face == "bottom") { v = { p4, p5, p1, p0 }; n = { 0,-1,0 }; }

    m_voxelVertices.insert(m_voxelVertices.end(), v.begin(), v.end());
    for (int i = 0; i < 4; ++i) { m_voxelNormals.push_back(n); m_voxelBlockType.push_back(blockType); }
    m_voxelIndices.insert(m_voxelIndices.end(), { baseIndex, baseIndex + 1, baseIndex + 2, baseIndex, baseIndex + 2, baseIndex + 3 });
}

void ProceduralEnvironment::createVoxelMesh() {
    if (m_voxelMesh.vao) {
        glDeleteVertexArrays(1, &m_voxelMesh.vao);
        glDeleteBuffers(3, m_voxelMesh.vbo);
        glDeleteBuffers(1, &m_voxelMesh.ibo);
    }

    glGenVertexArrays(1, &m_voxelMesh.vao);
    glBindVertexArray(m_voxelMesh.vao);
    glGenBuffers(3, m_voxelMesh.vbo);

    glBindBuffer(GL_ARRAY_BUFFER, m_voxelMesh.vbo[0]);
    glBufferData(GL_ARRAY_BUFFER, m_voxelVertices.size() * sizeof(glm::vec3), m_voxelVertices.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    glBindBuffer(GL_ARRAY_BUFFER, m_voxelMesh.vbo[1]);
    glBufferData(GL_ARRAY_BUFFER, m_voxelNormals.size() * sizeof(glm::vec3), m_voxelNormals.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    glBindBuffer(GL_ARRAY_BUFFER, m_voxelMesh.vbo[2]);
    glBufferData(GL_ARRAY_BUFFER, m_voxelBlockType.size() * sizeof(float), m_voxelBlockType.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, nullptr);

    glGenBuffers(1, &m_voxelMesh.ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_voxelMesh.ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_voxelIndices.size() * sizeof(unsigned int), m_voxelIndices.data(), GL_STATIC_DRAW);

    m_voxelMesh.indexCount = static_cast<GLsizei>(m_voxelIndices.size());
    glBindVertexArray(0);
}

void ProceduralEnvironment::createShaders() {
    const char* vertexSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;
        layout (location = 2) in float aBlockType;

        out vec3 Normal;
        out float BlockType;
        
        uniform mat4 uProjectionMatrix;
        uniform mat4 uModelViewMatrix;
        
        void main() {
            Normal = aNormal;
            BlockType = aBlockType;
            gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(aPos, 1.0);
        }
    )";

    const char* fragmentSource = R"(
        #version 330 core
        out vec4 FragColor;
        
        in vec3 Normal;
        in float BlockType;
        
        uniform vec3 uLightDirection;
        
        void main() {
            vec3 grassColor = vec3(0.2, 0.6, 0.2);
            vec3 dirtColor = vec3(0.5, 0.35, 0.2);
            vec3 stoneColor = vec3(0.5, 0.5, 0.5);
            vec3 snowColor = vec3(0.9, 0.9, 0.95);
            vec3 bedrockColor = vec3(0.2, 0.2, 0.25);
            vec3 woodColor = vec3(0.6, 0.4, 0.2);
            vec3 cobblestoneColor = vec3(0.45, 0.45, 0.45);
            vec3 glassColor = vec3(0.8, 0.9, 1.0);
            
            vec3 baseColor;
            float alpha = 1.0;
            
            float roundedBlockType = round(BlockType);
            
            if (roundedBlockType == 1.0) baseColor = grassColor;
            else if (roundedBlockType == 2.0) baseColor = dirtColor;
            else if (roundedBlockType == 3.0) baseColor = stoneColor;
            else if (roundedBlockType == 4.0) baseColor = snowColor;
            else if (roundedBlockType == 5.0) baseColor = bedrockColor;
            else if (roundedBlockType == 6.0) baseColor = woodColor;
            else if (roundedBlockType == 7.0) baseColor = cobblestoneColor;
            else if (roundedBlockType == 8.0) {
                baseColor = glassColor;
                alpha = 0.4;
            }
            else baseColor = vec3(1.0, 0.0, 1.0);

            float ao = 0.7 + (max(0.0, Normal.y) * 0.3);
            vec3 lightDir = normalize(uLightDirection);
            float diff = max(dot(Normal, lightDir), 0.0);
            vec3 diffuse = diff * vec3(1.0);
            vec3 ambient = vec3(0.4);
            
            vec3 result = (ambient + diffuse) * baseColor * ao;
            FragColor = vec4(result, alpha);
        }
    )";
    m_voxelShader = createShaderProgram(vertexSource, fragmentSource);
}

GLuint ProceduralEnvironment::createShaderProgram(const char* vertexSource, const char* fragmentSource) {
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, nullptr);
    glCompileShader(vertexShader);
    GLint success;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        GLchar infoLog[512];
        glGetShaderInfoLog(vertexShader, 512, nullptr, infoLog);
        std::cerr << "Vertex shader compilation failed: " << infoLog << std::endl;
    }

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, nullptr);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        GLchar infoLog[512];
        glGetShaderInfoLog(fragmentShader, 512, nullptr, infoLog);
        std::cerr << "Fragment shader compilation failed: " << infoLog << std::endl;
    }

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        GLchar infoLog[512];
        glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
        std::cerr << "Shader program linking failed: " << infoLog << std::endl;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return shaderProgram;
}

void ProceduralEnvironment::render(const glm::mat4& view, const glm::mat4& projection) {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    glm::vec3 lightDirection = glm::normalize(glm::vec3(0.5f, 1.0f, 0.3f));

    glUseProgram(m_voxelShader);
    glUniformMatrix4fv(glGetUniformLocation(m_voxelShader, "uProjectionMatrix"), 1, GL_FALSE, &projection[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(m_voxelShader, "uModelViewMatrix"), 1, GL_FALSE, &view[0][0]);
    glUniform3fv(glGetUniformLocation(m_voxelShader, "uLightDirection"), 1, &lightDirection[0]);

    glBindVertexArray(m_voxelMesh.vao);
    if (m_voxelMesh.indexCount > 0) {
        glDrawElements(GL_TRIANGLES, m_voxelMesh.indexCount, GL_UNSIGNED_INT, 0);
    }

    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_CULL_FACE);
}