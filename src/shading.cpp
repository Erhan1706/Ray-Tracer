#include "texture.h"
#include "draw.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

int mipMapLevel = 0;
bool debugMipMap = false;

void getMipMapLevel(int level) {
    mipMapLevel = level;
}

void setDebugMipMap(bool debug) {
    debugMipMap = debug;
}

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading || features.enableTextureMapping) {
        //Diffuse
        glm::vec3 surfacePoint = ray.origin + ray.direction * ray.t;
        glm::vec3 lightVector = glm::normalize(lightPosition - surfacePoint);
        glm::vec3 normal = glm::normalize(hitInfo.normal);
        float dotDiffuse = glm::dot(normal, lightVector);
        dotDiffuse = dotDiffuse < 0 ? 0 : dotDiffuse;
        // Specular
        glm::vec3 viewVector = glm::normalize(ray.origin - surfacePoint);
        glm::vec3 reflectedVector = 2 * glm::dot(lightVector, normal) * normal - lightVector;
        float dotSpecular = glm::dot(glm::normalize(reflectedVector), viewVector);
        dotSpecular = dotSpecular < 0 ? 0 : dotSpecular;
        glm::vec3 kd = hitInfo.material.kd;
        if (features.enableTextureMapping) {
            if (hitInfo.material.kdTexture) {
                Image* img = hitInfo.material.kdTexture.get();
                kd = acquireTexel(*img, hitInfo.texCoord, features);
                // If mipmap is enabling we have to generate the entire pyramid.
                if (features.extra.enableMipmapTextureFiltering) {
                    if (hitInfo.material.kdTexture) {
                        // Fetching the texel
                        // If debug is on we do not compute the footprint of a pixel, we just use a mipmap level as the 
                        // the entire texture.
                        if (debugMipMap) { 
                            kd = acquireTexelMipMap(mipMapLevel, hitInfo.texCoord, features);
                            if (kd == glm::vec3(NULL))
                                kd = hitInfo.material.kd;
                        } else {
                        // Calculate the pixel footprint to know where to apply the mipmap.
                            float d = glm::distance(ray.origin, surfacePoint);                     
                            float mipMapel = log2(d);
                            float weight = mipMapel - (int) mipMapel;
                            // Trilinear interpolation
                            mipMapLevel = std::min(6,(int)floor(mipMapel));
                            kd = acquireTexelMipMap(mipMapLevel, hitInfo.texCoord, features) * (1- weight) 
                                + weight * acquireTexelMipMap(mipMapLevel+1, hitInfo.texCoord, features);
                        }
                    }
                }
            } else {
                // To prevent shading being computed when only the texture checkbox is checked for a mesh with no textures.
                return hitInfo.material.kd;
            }
        }
        return (lightColor * kd * dotDiffuse) + (lightColor * hitInfo.material.ks 
            * pow(dotSpecular, hitInfo.material.shininess));
    } 
    return hitInfo.material.kd;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    glm::vec3 dRay = { ray.direction.x, ray.direction.y, ray.direction.z };
    dRay = normalize(dRay);
    glm::vec3 reflectionRayDirection = dRay - 2 * dot(dRay, hitInfo.normal) * hitInfo.normal;

    Ray reflectionRay = {{ ray.origin.x + ray.direction.x * ray.t, ray.origin.y + ray.direction.y * ray.t, ray.origin.z + ray.direction.z * ray.t }, reflectionRayDirection, ray.t };
    
    return reflectionRay;
}