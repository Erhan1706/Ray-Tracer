#include "cube_map.h"
#include "texture.h"
#include <glm/glm.hpp>

Image environment_map = Image("../../../data/environment.jpg");

glm::vec2 convert_xyz_to_uv(const float& x, const float& y, const float& z) {
    float X, Y, Z, U, V, maxSide, offsetU, offsetV;

    U = V = 0;
    offsetU = offsetV = 0;

    X = std::abs(x);
    Y = std::abs(y);
    Z = std::abs(z);

    maxSide = std::max({ X, Y, Z });

    if (x <= 0 && X == maxSide) {
        U = z;
        V = y; 

        offsetV = 1.0f / 3.0f;
    }

    if (x > 0 && X == maxSide) {
        U = -z;
        V = y;

        offsetV = 1.0f / 3.0f;
        offsetU = 2.0f / 4.0f;
    }

    if (y <= 0 && Y == maxSide) {
        U = x;
        V = z;

        offsetU = 1.0f / 4.0f;
    }

    if (y > 0 && Y == maxSide) {
        U = x;
        V = -z;

        offsetV = 2.0f / 3.0f;
        offsetU = 1.0f / 4.0f;
    }

    if (z <= 0 && Z == maxSide) {
        U = -x;
        V = y;

        offsetV = 1.0f / 3.0f;
        offsetU = 3.0f / 4.0f;
    }

    if (z > 0 && Z == maxSide) {
        U = x;
        V = y;

        offsetV = 1.0f / 3.0f;
        offsetU = 1.0f / 4.0f;
    }

    U = 0.5f * (U + 1);
    V = 0.5f * (V + 1);

    return glm::vec2 { U / 4.0f + offsetU, V / 3.0f + offsetV };
}

glm::vec3 getColorEnvironment(Ray& ray, const Features& features)
{
    glm::vec3 pointOnCube = ray.origin + ray.t * ray.direction;

    pointOnCube /= std::max({ std::abs(pointOnCube.x), std::abs(pointOnCube.y), std::abs(pointOnCube.z) });

    glm::vec2 point = convert_xyz_to_uv(pointOnCube.x, pointOnCube.y, pointOnCube.z);

    return acquireTexel(environment_map, point, features);
}
