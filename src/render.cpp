#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "texture.h"
#include "cube_map.h"
#include <framework/trackball.h>
#include <random>
#ifdef NDEBUG
#include <omp.h>
#endif

const bool DEBUG = true;

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    Ray reflection;
    HitInfo hitInfo;
    bool intersect = false;

    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 color = computeLightContribution(scene, bvh, features, ray, hitInfo);

        // Draw a white debug ray if the ray hits and shading is disabled.
        features.enableShading || features.extra.enableTransparency ? drawRay(ray, color) : drawRay(ray, glm::vec3(1.0f));

        if (features.enableSoftShadow)
            softShadowsDebug(scene, ray, bvh, features);

        if (features.enableHardShadow) {
            hardShadowsDebug(scene, ray, bvh, features);
        }
       

        if (features.extra.enableTransparency) {
            if (hitInfo.material.transparency != 1.0f) {
                ray.origin = (ray.origin + ray.direction * ray.t) + ray.direction * 0.0001f;
                ray.t = std::numeric_limits<float>::max();
                color = color * hitInfo.material.transparency + (1.0f - hitInfo.material.transparency) * getFinalColor(scene, bvh, ray, features);
            }
        }
        
        if (features.enableRecursive) {
            reflection = computeReflectionRay(ray, hitInfo);
            reflection.origin += (reflection.direction * 0.001f);
            reflection.t = std::numeric_limits<float>::max();
            rayDepth = rayDepth + 1;
            if (rayDepth < 4 && hitInfo.material.ks != glm::vec3{ 0.0f }) {
                color += hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth);
            }
            //option 2: (is same as example)
            /*
            if (rayDepth < 300 && hitInfo.material.ks != glm::vec3 { 0.0f }) {
                if (rayDepth <= 1) {
                    color = color * (glm::vec3(1.0f) - hitInfo.material.ks);
                }
                glm::vec3 tempColor = hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth);
                color += tempColor;
            }
            */
        } else if (features.extra.enableEnvironmentMapping) {
            reflection = computeReflectionRay(ray, hitInfo);

            color = getColorEnvironment(reflection, features);
        }

        // Set the color of the pixel to white if the ray hits.
        return color;
    } else {
        // Draw a red debug ray if the ray missed.

        if (features.extra.enableEnvironmentMapping) {
            return getColorEnvironment(ray, features);
        }

        if (DEBUG) {
            drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        }
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    std::mt19937 rng;
    std::uniform_real_distribution<float> dist(0, 1);
    const int n = 4;
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            glm::vec3 color;
            std::vector<glm::vec2> points;

            // NOTE: (-1, -1) at the bottom left of the screen, (-1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };

            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            color = getFinalColor(scene, bvh, cameraRay, features);

            if (features.extra.enableDepthOfField) {
                Trackball temp_camera = camera;
                glm::vec3 focal_point = cameraRay.origin + cameraRay.direction * camera.m_focal_distance;
                Ray temp_ray;
                int m = 80;

                for (int i = 0; i < m; ++i) {
                    glm::vec3 rand_dir = glm::vec3(dist(rng) * 2 - 1, dist(rng) * 2 - 1, dist(rng) * 2 - 1);
                    temp_camera.setCamera(camera.lookAt(), camera.rotationEulerAngles() + camera.m_aperture * rand_dir * 0.01f, camera.distanceFromLookAt());

                    temp_ray.origin = temp_camera.position();
                    temp_ray.direction = glm::normalize(focal_point - temp_camera.position());
                    temp_ray.t = std::numeric_limits<float>::max();

                    color += getFinalColor(scene,
                        bvh,
                        temp_ray,
                        features);
                }

                color = color / (m * 1.0f);
            }

            if (features.extra.enableMultipleRaysPerPixel) {
                for (int i = 1; i < n; ++i) {
                    for (int j = 1; j < n; ++j) { 
                        color += getFinalColor(scene,
                            bvh,
                            camera.generateRay(
                                { 
                                    float(x + 1.0 * (i + dist(rng)) / n) / float(windowResolution.x) * 2.0f - 1.0f,
                                    float(y + 1.0 * (j + dist(rng)) / n) / float(windowResolution.y) * 2.0f - 1.0f
                                }), 
                            features
                        );
                    }
                }

                color /= (n - 1) * (n - 1);
            }

            screen.setPixel(x, y, color);
        }
    }
}

// Special render method that takes in some extra parameters to set the motion.
void renderRayTracingMotion(const Scene& scene, Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features, float time, bool debug, float translationX, float translationY, float translationZ)
{
    glm::ivec2 windowResolution = screen.resolution();
    std::mt19937 rng;
    std::uniform_real_distribution<float> dist(0, 1);
    const int n = 4;
    glm::vec3 initialPosition = camera.lookAt();
    // This is for the visual debug of motion blur where all the tranlated mesh is shown.
    if (debug) {
        glm::vec3 newCameraPosition = { camera.lookAt().x + translationX * time, camera.lookAt().y + translationY * time, camera.lookAt().z + translationZ * time };
        camera.setCamera(newCameraPosition, camera.rotationEulerAngles() , camera.distanceFromLookAt());
    }
    
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            glm::vec3 color;
            std::vector<glm::vec2> points;

            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            Ray cameraRay;
            if (DEBUG)
                cameraRay = camera.generateRay(normalizedPixelPos);
            else 
                // Special ray method that takes random time between 0 and 1 to set where the ray originates from.
                cameraRay = camera.generateRayMotion(normalizedPixelPos, dist(rng), translationX, translationY, translationZ);

            color = getFinalColor(scene, bvh, cameraRay, features);

            if (features.extra.enableMultipleRaysPerPixel) {
                for (int i = 1; i < n; ++i) {
                    for (int j = 1; j < n; ++j) {
                        color += getFinalColor(scene,
                            bvh,
                            camera.generateRay(
                                { float(x + 1.0 * (i + dist(rng)) / n) / float(windowResolution.x) * 2.0f - 1.0f,
                                    float(y + 1.0 * (j + dist(rng)) / n) / float(windowResolution.y) * 2.0f - 1.0f }),
                            features);
                    }
                }

                color /= n * n;
            }

            screen.setPixel(x, y, color);
        }
    }
    // To avoid the camera moving away we set it at the initial position at the end of the method.
    if (debug)
        camera.setCamera(initialPosition, camera.rotationEulerAngles(), camera.distanceFromLookAt());
}

std::vector<Ray> debug_toDraw;

void debug_dof(const Scene& scene, Ray ray, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    if (!DEBUG) {
        return;
    }

    if (camera.m_aperture == 0) {
        debug_toDraw.clear();

        glm::ivec2 windowResolution = screen.resolution();
        std::mt19937 rng(0);
        std::uniform_real_distribution<float> dist(0, 1);

        HitInfo dummy {};

        const Ray cameraRay = ray;

        Trackball temp_camera = camera;

        temp_camera.m_aperture = 1.0f;    

        glm::vec3 focal_point = cameraRay.origin + cameraRay.direction * camera.m_focal_distance;
        Ray temp_ray;
        int m = 16;
        for (int i = 0; i < m; ++i) {
            glm::vec3 rand_dir = glm::vec3(dist(rng) * 2 - 1, dist(rng) * 2 - 1, dist(rng) * 2 - 1);
            temp_camera.setCamera(camera.lookAt(), camera.rotationEulerAngles() + temp_camera.m_aperture * rand_dir * 0.01f, camera.distanceFromLookAt());

            temp_ray.origin = temp_camera.position();
            temp_ray.direction = glm::normalize(focal_point - temp_camera.position());
            temp_ray.t = std::numeric_limits<float>::max();

            bvh.intersect(temp_ray, dummy, features);

            debug_toDraw.push_back(temp_ray);
        }
    }

    for (auto it : debug_toDraw) {
        drawRay(it, glm::vec3(0.05f, 0.8f, 0.05f));
    }
}
