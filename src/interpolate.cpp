#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 normal = glm::cross((v1 - v0), (v2 - v0));
    float normalLength = sqrt(pow(normal.x, 2) + pow(normal.y, 2) + pow(normal.z, 2));

    glm::vec3 normalA = glm::cross((v2 - v1), (p - v1));
    glm::vec3 normalB = glm::cross((v0 - v2), (p - v2));
    glm::vec3 normalC = glm::cross((v1 - v0), (p - v0));

    float alfa = (glm::dot(normal, normalA)) / (pow(normalLength, 2));
    float beta = (glm::dot(normal, normalB)) / (pow(normalLength, 2));
    float gamma = (glm::dot(normal, normalC)) / (pow(normalLength, 2));

    return glm::vec3(alfa, beta, gamma);
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    glm::vec3 interpolatedNormal = (barycentricCoord.x * n0) + (barycentricCoord.y * n1) + (barycentricCoord.z * n2);
    interpolatedNormal = glm::normalize(interpolatedNormal);
    return interpolatedNormal;
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return barycentricCoord.x * t0 + barycentricCoord.y * t1 + barycentricCoord.z * t2;
}
