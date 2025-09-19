// src/SlimeSimulation/SlimeBlock.hpp

#pragma once
#include "application.hpp" // We need the basic_model struct

class SlimeBlock {
private:
    basic_model m_model;

public:
    // Constructor will create the cube mesh
    SlimeBlock();

    // Draw method
    void draw(const glm::mat4& view, const glm::mat4& proj);
};
