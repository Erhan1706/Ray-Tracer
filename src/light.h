#pragma once
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "shading.h"

void sampleSegmentLight (const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color);

void sampleParallelogramLight (const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color);

float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo);

glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo);

void getNumSamples(int numOfSamples);

void softShadowsDebug(const Scene& scene, Ray ray, const BvhInterface& bvh, const Features& features);

void hardShadowsDebug(const Scene& scene, Ray ray, const BvhInterface& bvh, const Features& features);

std::vector<std::variant<SegmentLight, ParallelogramLight>> getAreaLights(const Scene& scene);

void drawShadowRay(PointLight sample, glm::vec3 surfacePoint, const BvhInterface& bvh, const Features& features, Ray ray);
