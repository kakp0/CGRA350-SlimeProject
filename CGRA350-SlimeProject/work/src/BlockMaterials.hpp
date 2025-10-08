#pragma once
#include <vector>
#include "BlockTypes.hpp"

// A utility class to create and manage all block materials.
class BlockMaterialRegistry {
public:
    // Fills a vector with the material properties for all defined block types.
    static void createAllMaterials(std::vector<BlockMaterial>& materials);

private:
    // These helpers define the complex details for each block.
    static BlockMaterial createGrassMaterial();
    static BlockMaterial createDirtMaterial();
    static BlockMaterial createStoneMaterial();
    static BlockMaterial createSnowMaterial();
    static BlockMaterial createBedrockMaterial();
    static BlockMaterial createWoodMaterial();
    static BlockMaterial createCobblestoneMaterial();
    static BlockMaterial createGlassMaterial();
};