#include "BlockMaterials.hpp"
// Made with the assistance of AI (Gemini 2.5 Pro)

// The following functions are simple factories for creating specific material types.
// Each function instantiates a BlockMaterial struct and populates its fields with
// hardcoded values that define the visual properties of a block in the renderer.

BlockMaterial BlockMaterialRegistry::createGrassMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.2f, 0.6f, 0.2f);    // A vibrant green base color.
    mat.detailColor = glm::vec3(0.1f, 0.4f, 0.15f); // A darker green for procedural texturing (e.g., noise).
    mat.materialParams.x = 12.0f;                  // Likely a frequency/tiling factor for a noise function in the shader.
    return mat;
}

BlockMaterial BlockMaterialRegistry::createDirtMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.5f, 0.35f, 0.2f);
    mat.detailColor = glm::vec3(0.35f, 0.25f, 0.15f);
    mat.materialParams.x = 20.0f; // A higher frequency suggests a finer, more detailed texture pattern than grass.
    return mat;
}

BlockMaterial BlockMaterialRegistry::createStoneMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.5f, 0.5f, 0.5f);       // Mid-gray for the main stone color.
    mat.detailColor = glm::vec3(0.45f, 0.45f, 0.45f);  // Slightly darker gray for texture details.
    mat.materialParams.x = 10.0f;                     // Controls the scale of the stone's procedural texture.
    return mat;
}

BlockMaterial BlockMaterialRegistry::createSnowMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.9f, 0.9f, 0.95f); // A bright, slightly blueish-white for snow.
    // No detail color or params; snow is likely rendered as a uniform color.
    return mat;
}

BlockMaterial BlockMaterialRegistry::createBedrockMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.2f, 0.2f, 0.25f);    // Very dark, slightly blue base.
    mat.detailColor = glm::vec3(0.15f, 0.15f, 0.2f); // Even darker detail color.
    mat.materialParams.x = 8.0f;                    // A lower frequency for larger, less repetitive patterns.
    return mat;
}

BlockMaterial BlockMaterialRegistry::createWoodMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.6f, 0.4f, 0.2f);
    mat.detailColor = glm::vec3(0.4f, 0.2f, 0.05f);
    // Uses more shader parameters than simpler materials.
    // x: Primary noise frequency.
    // y: Secondary frequency or distortion factor for wood grain.
    // z: Controls the sharpness or contrast of the grain pattern.
    mat.materialParams = glm::vec4(18.0f, 4.0f, 0.95f, 0.0f);
    return mat;
}

BlockMaterial BlockMaterialRegistry::createCobblestoneMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.45f, 0.45f, 0.45f);
    mat.detailColor = glm::vec3(0.2f, 0.2f, 0.2f);
    // Parameters likely for a Voronoi or cellular noise shader to create the stone shapes.
    // x: Cell density/scale.
    // z: Mortar thickness or color blend factor.
    mat.materialParams = glm::vec4(1.5f, 0.0f, 0.08f, 0.0f);
    return mat;
}

BlockMaterial BlockMaterialRegistry::createGlassMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.8f, 0.9f, 1.0f); // A bright, translucent color.
    mat.alpha = 0.4f;                           // Sets transparency; requires the renderer to handle blending.
    return mat;
}

// This function orchestrates the creation of all block materials and populates a lookup table.
// Using a vector indexed by BlockID provides an O(1) lookup for a block's material.
void BlockMaterialRegistry::createAllMaterials(std::vector<BlockMaterial>& materials) {
    // Pre-allocate the vector to the exact size needed, preventing reallocations.
    materials.resize(static_cast<size_t>(BlockID::Count));

    // Populate the vector, using the enum's integer value as the index for each material.
    // This establishes a direct mapping from a block's ID to its rendering properties.
    materials[static_cast<int>(BlockID::Air)] = BlockMaterial(); // Air has a default, empty material.
    materials[static_cast<int>(BlockID::Grass)] = createGrassMaterial();
    materials[static_cast<int>(BlockID::Dirt)] = createDirtMaterial();
    materials[static_cast<int>(BlockID::Stone)] = createStoneMaterial();
    materials[static_cast<int>(BlockID::Snow)] = createSnowMaterial();
    materials[static_cast<int>(BlockID::Bedrock)] = createBedrockMaterial();
    materials[static_cast<int>(BlockID::Wood)] = createWoodMaterial();
    materials[static_cast<int>(BlockID::Cobblestone)] = createCobblestoneMaterial();
    materials[static_cast<int>(BlockID::Glass)] = createGlassMaterial();
}
