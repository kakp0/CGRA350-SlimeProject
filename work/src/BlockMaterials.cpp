#include "BlockMaterials.hpp"

BlockMaterial BlockMaterialRegistry::createGrassMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.2f, 0.6f, 0.2f);
    mat.detailColor = glm::vec3(0.1f, 0.4f, 0.15f);
    mat.materialParams.x = 12.0f;
    return mat;
}

BlockMaterial BlockMaterialRegistry::createDirtMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.5f, 0.35f, 0.2f);
    mat.detailColor = glm::vec3(0.35f, 0.25f, 0.15f);
    mat.materialParams.x = 20.0f;
    return mat;
}

BlockMaterial BlockMaterialRegistry::createStoneMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.5f, 0.5f, 0.5f);
    mat.detailColor = glm::vec3(0.45f, 0.45f, 0.45f);
    mat.materialParams.x = 10.0f;
    return mat;
}

BlockMaterial BlockMaterialRegistry::createSnowMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.9f, 0.9f, 0.95f);
    return mat;
}

BlockMaterial BlockMaterialRegistry::createBedrockMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.2f, 0.2f, 0.25f);
    mat.detailColor = glm::vec3(0.15f, 0.15f, 0.2f);
    mat.materialParams.x = 8.0f;
    return mat;
}

BlockMaterial BlockMaterialRegistry::createWoodMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.6f, 0.4f, 0.2f);
    mat.detailColor = glm::vec3(0.4f, 0.2f, 0.05f);
    mat.materialParams = glm::vec4(18.0f, 4.0f, 0.95f, 0.0f);
    return mat;
}

BlockMaterial BlockMaterialRegistry::createCobblestoneMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.45f, 0.45f, 0.45f);
    mat.detailColor = glm::vec3(0.2f, 0.2f, 0.2f);
    mat.materialParams = glm::vec4(1.5f, 0.0f, 0.08f, 0.0f);
    return mat;
}

BlockMaterial BlockMaterialRegistry::createGlassMaterial() {
    BlockMaterial mat;
    mat.baseColor = glm::vec3(0.8f, 0.9f, 1.0f);
    mat.alpha = 0.4f;
    return mat;
}

void BlockMaterialRegistry::createAllMaterials(std::vector<BlockMaterial>& materials) {
    materials.resize(static_cast<size_t>(BlockID::Count));
    materials[static_cast<int>(BlockID::Air)] = BlockMaterial();
    materials[static_cast<int>(BlockID::Grass)] = createGrassMaterial();
    materials[static_cast<int>(BlockID::Dirt)] = createDirtMaterial();
    materials[static_cast<int>(BlockID::Stone)] = createStoneMaterial();
    materials[static_cast<int>(BlockID::Snow)] = createSnowMaterial();
    materials[static_cast<int>(BlockID::Bedrock)] = createBedrockMaterial();
    materials[static_cast<int>(BlockID::Wood)] = createWoodMaterial();
    materials[static_cast<int>(BlockID::Cobblestone)] = createCobblestoneMaterial();
    materials[static_cast<int>(BlockID::Glass)] = createGlassMaterial();
}