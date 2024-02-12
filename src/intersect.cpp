#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    glm::vec3 vec0 = v1 - v0;
    glm::vec3 vec1 = v2 - v0;
    glm::vec3 vec2 = p - v0;

    float d00 = glm::dot(vec0, vec0);
    float d01 = glm::dot(vec0, vec1);
    float d11 = glm::dot(vec1, vec1);
    float d20 = glm::dot(vec2, vec0);
    float d21 = glm::dot(vec2, vec1);
    float div = d00 * d11 - d01 * d01;

    if (div == 0)
        return false;

    float a = (d11 * d20 - d01 * d21) / div;
    float b = (d00 * d21 - d01 * d20) / div;

    if (a < 0 || b < 0)
        return false;
    if (a + b > 1)
        return false;

    return true;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    if (abs(glm::dot(ray.origin, plane.normal) - plane.D) < 0.01)
        return false;

    float div = glm::dot(plane.normal, ray.direction);
    if (abs(div) > 0.00001f) {
        float t = (plane.D - glm::dot(ray.origin, plane.normal)) / div;
        if (t > 0) {
            if (ray.t > t) {
                ray.t = t;
                return true;
            }
        }
    }
    return false;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 vertice1 = v2 - v0;
    glm::vec3 vertice2 = v2 - v1;
    glm::vec3 newNormal = glm::cross(vertice1, vertice2);
    plane.normal = glm::normalize(newNormal);
    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    Plane triPlane = trianglePlane(v0, v1, v2);
    float oldT = ray.t;
    if (intersectRayWithPlane(triPlane, ray)) {
        glm::vec3 p = ray.origin + ray.t * ray.direction;
        if (pointInTriangle(v0, v1, v2, triPlane.normal, p)) {         
            // This dot product detects points outside of the Cornell box, that should be completely disprovided of light.
            // By inverting the normal we make sure that these points will get no shading and be dark.
            hitInfo.normal = glm::dot(triPlane.normal, glm::normalize(ray.direction)) < 0 ? triPlane.normal : -triPlane.normal;
            return true;
        }
    }
    ray.t = oldT;
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    glm::vec3 dirNormal = glm::normalize(ray.direction);
    float deltaPart1 = std::pow(glm::dot(dirNormal, ray.origin - sphere.center), 2);
    float deltaPart2 = std::pow(glm::length(ray.origin - sphere.center), 2) - std::pow(sphere.radius, 2);
    float delta = deltaPart1 - deltaPart2;

    if (delta < 0)
        return false;

    float distance = -(glm::dot(dirNormal, ray.origin - sphere.center)) - std::sqrt(delta);

    if (distance < 0) {
        distance = -(glm::dot(dirNormal, ray.origin - sphere.center)) + std::sqrt(delta);
    }
    if (ray.t > distance && distance > 0) {
        ray.t = distance;
        hitInfo.material = sphere.material;
        glm::vec3 surfacePoint = ray.origin + ray.t * ray.direction;
        hitInfo.normal = glm::normalize(surfacePoint - sphere.center);
        return true;
    }
    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    // TODO: implement this function.
    if (ray.direction.x == 0 || ray.direction.y == 0 || ray.direction.z == 0)
        return false;
    float txmin = (box.lower.x - ray.origin.x) / ray.direction.x;
    float txmax = (box.upper.x - ray.origin.x) / ray.direction.x;
    float tymin = (box.lower.y - ray.origin.y) / ray.direction.y;
    float tymax = (box.upper.y - ray.origin.y) / ray.direction.y;
    float tzmin = (box.lower.z - ray.origin.z) / ray.direction.z;
    float tzmax = (box.upper.z - ray.origin.z) / ray.direction.z;

    float tinX = std::min(txmin, txmax);
    float toutX = std::max(txmin, txmax);
    float tinY = std::min(tymin, tymax);
    float toutY = std::max(tymin, tymax);
    float tinZ = std::min(tzmin, tzmax);
    float toutZ = std::max(tzmin, tzmax);

    float tin = std::max(tinX, tinY);
    tin = std::max(tin, tinZ);
    float tout = std::min(toutX, toutY);
    tout = std::min(tout, toutZ);

    if (((tin > 0) && (tin > tout)) || tout < 0)
        return false;

    if (tin < 0 && tout > 0) {
        ray.t = tout;
        return true;
    } else if (ray.t > tin) {
        ray.t = tin;
        return true;
    }

    return false;
}

