#pragma once
#include <glm/glm.hpp>

// Defines all the unique block identifiers 
enum class BlockID : int {
    Air = 0,
    Grass = 1,
    Dirt = 2,
    Stone = 3,
    Snow = 4,
    Bedrock = 5,
    Wood = 6,
    Cobblestone = 7,
    Glass = 8,
    Count // A special value to get the total number of block types
};

// Represents the material properties of a block, including parameters for procedural generation.
struct BlockMaterial {
    glm::vec3 baseColor{ 1.0f, 0.0f, 1.0f };
    float alpha = 1.0f;

    // Parameters for procedural texturing
    glm::vec3 detailColor{ 0.0f, 0.0f, 0.0f };
    glm::vec4 materialParams{ 0.0f, 0.0f, 0.0f, 0.0f };
    // x: primary frequency/density
    // y: secondary frequency/size
    // z: blend factor/sharpness
    // w: pattern selector/type
};