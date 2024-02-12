#pragma once
#include "framework/image.h"
#include <framework/ray.h>
#include "common.h"
#include <glm/glm.hpp>

extern Image environment_map;

glm::vec2 convert_xyz_to_uv(const float& x, const float& y, const float& z);

glm::vec3 getColorEnvironment(Ray& ray, const Features& features);
