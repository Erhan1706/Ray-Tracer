#pragma once

#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include "common.h"
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <map>

// Forward declarations.
struct Image;

struct MipMap {
    int level;
    int width, height;
    std::vector<glm::vec3> pixels;
};
;
// Given an image and a texture coordinate, return the corresponding texel.
glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features);

MipMap initialiseMipMap(const Image& image);

MipMap generateMipMapPyramid(int height, int width, std::vector<glm::vec3> pixels, int level);

std::map<int, MipMap> getMipMapPyramid();

glm::vec3 acquireTexelMipMap(int level, const glm::vec2& texCoord, const Features& features);