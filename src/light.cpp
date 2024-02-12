#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <vector>
#include <iostream>
#include <random>

int numSamples = 10;

void getNumSamples(int numOfSamples)
{
    numSamples = numOfSamples;
}

void softShadowsDebug(const Scene& scene, Ray ray, const BvhInterface& bvh, const Features& features)
{
    std::vector<std::variant<SegmentLight, ParallelogramLight>> areaLights = getAreaLights(scene);
    if (areaLights.empty())
        return;
    glm::vec3 surfacePoint = ray.origin + ray.t * ray.direction;
    for (std::variant<SegmentLight, ParallelogramLight> light : areaLights) {
        if (std::holds_alternative<SegmentLight>(light)) {
            const SegmentLight segmentLight = std::get<SegmentLight>(light);
            for (int i = 0; i < numSamples; i++) {
                PointLight sample;
                sampleSegmentLight(segmentLight, sample.position, sample.color);
                drawShadowRay(sample, surfacePoint, bvh, features, ray);
            }
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
            for (int i = 0; i < numSamples; i++) {
                PointLight sample;
                sampleParallelogramLight(parallelogramLight, sample.position, sample.color);
                drawShadowRay(sample, surfacePoint, bvh, features, ray);
            }
        }
    }
    drawRay(ray, glm::vec3 { 0.0f, 0, 1.0f });
}
// Auxiliary method that returns a list of all Segment and Parallelogram lights on the scene.
std::vector<std::variant<SegmentLight, ParallelogramLight>> getAreaLights(const Scene& scene) {
    std::vector<std::variant<SegmentLight, ParallelogramLight>> result;
    for (const auto& light : scene.lights) {
        if (std::holds_alternative<SegmentLight>(light)) {
            const SegmentLight segmentLight = std::get<SegmentLight>(light);
            result.push_back(segmentLight);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
            result.push_back(parallelogramLight);
        }
    }
    return result;
}

std::vector<PointLight> getPointLights(const Scene& scene) {
    std::vector<PointLight> lights;
    for (const auto& light : scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            const PointLight pointLight = std::get<PointLight>(light);
            lights.push_back(pointLight);
        }
    }
    return lights;
}

void hardShadowsDebug(const Scene& scene, Ray ray, const BvhInterface& bvh, const Features& features) {
    for (PointLight light: getPointLights(scene)) {
        drawShadowRay(light, ray.origin + ray.direction * ray.t, bvh, features, ray);
    }
}

//Auxiliary method to draw each shadow debug ray for soft shadows
void drawShadowRay(PointLight sample, glm::vec3 surfacePoint, const BvhInterface& bvh, const Features& features, Ray ray)
{
    Ray shadowRay;
    shadowRay.origin = surfacePoint + 0.0001f;
    shadowRay.direction = sample.position - surfacePoint;
    HitInfo shadowRayHitInfo;
    shadowRay.t = 1.0f;
    bvh.intersect(shadowRay, shadowRayHitInfo, features);
    testVisibilityLightSample(sample.position, glm::vec3 { 0, 0, 0 }, bvh, features, ray, shadowRayHitInfo) == 1.0f ? 
        drawRay(shadowRay, sample.color) : drawRay(shadowRay, glm::vec3 {1.0f ,0 , 0});
}

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    // Generate a random uniformly distributed number between 0 and 1.
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0, 1);
    float x1 = dist(gen);
    glm::vec3 segmentVector = segmentLight.endpoint0 - segmentLight.endpoint1;
    position = (segmentVector * x1) + segmentLight.endpoint1;
    // Interpolate the color
    color = x1 * segmentLight.color0 + (1.0f - x1) * segmentLight.color1;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    // Generate a random uniformly distributed number between 0 and 1.
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0, 1);
    float x1 = dist(gen);
    float x2 = dist(gen);
    position = parallelogramLight.v0 + x1 * parallelogramLight.edge01 + x2 * parallelogramLight.edge02;
    //Compute the color by bilinear interpolation
    color = ((1.0f - x1) * parallelogramLight.color0 + x1 * parallelogramLight.color1) * (1.0f - x2) +
        ((1.0f - x1) * parallelogramLight.color2 + x1 * parallelogramLight.color3) * x2;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 surfacePoint = ray.origin + ray.direction * ray.t;
    // Vector from the light point to the intersection position.
    Ray shadowRay;
    shadowRay.direction = samplePos - surfacePoint;
    // Offset to reduce noise a bit.
    shadowRay.origin = surfacePoint + 0.0001f;
    HitInfo shadowRayHitInfo;
    // The t parameter for each light source is 1.
    shadowRay.t = 1.0f;
    bvh.intersect(shadowRay, shadowRayHitInfo, features);
    // The intersection function will modify the t value if hits anything closer than the light source.
    // See if ray.t is smaller than the previously defined value (1).
    if (shadowRay.t < 1.0f)
        return 0.f;
    else
        return 1.0f;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading || features.enableTextureMapping) {
        // If shading is enabled, compute the contribution from all lights.
        glm::vec3 color;
        scene.lights.empty() ? color = hitInfo.material.kd : color = { 0, 0, 0 };
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                // Test if point is in the shadow
                if (features.enableHardShadow) {
                    if (testVisibilityLightSample(pointLight.position, glm::vec3 { 1.0 }, bvh, features, ray, hitInfo) == 1.0f)
                        color += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                } else 
                    color += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                // Perform your calculations for a segment light.
                for (int i = 0; i < numSamples; i++) {
                    PointLight sample;
                    HitInfo dummy;
                    sampleSegmentLight(segmentLight, sample.position, sample.color);
                    if (features.enableSoftShadow) {
                        if (testVisibilityLightSample(sample.position, glm::vec3 { 0, 0, 0 }, bvh, features, ray, dummy) == 1.0f) {
                            color += computeShading(sample.position, sample.color, features, ray, hitInfo);
                        }
                    } else
                        color += computeShading(sample.position, sample.color, features, ray, hitInfo);
                }
                color = color / (float) numSamples;
            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                // Perform your calculations for a parallelogram light.
                for (int i = 0; i < numSamples; i++) {
                    PointLight sample;
                    HitInfo dummy;
                    sampleParallelogramLight(parallelogramLight, sample.position, sample.color);
                    if (features.enableSoftShadow) {
                        if (testVisibilityLightSample(sample.position, glm::vec3 { 0, 0, 0 }, bvh, features, ray, dummy) == 1.0f) {
                            color += computeShading(sample.position, sample.color, features, ray, hitInfo);
                        }
                    } else
                        color += computeShading(sample.position, sample.color, features, ray, hitInfo);
                }
                color = color / (float) numSamples;
            }
        }
        return color;
    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}