#include "texture.h"
#include <framework/image.h>
#include <map>

std::map<int, MipMap> mipMapPyramid;

std::map<int, MipMap> getMipMapPyramid() {
    return mipMapPyramid;
}

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    if (features.extra.enableBilinearTextureFiltering) {
        // Source: Fundamentals of Computer Graphics by Steve Marschner and Peter Shirley
        float x = texCoord.x;
        float y = 1 - texCoord.y;
        // Transform texture coordinates between 0 and 1 to coordinates between the width and height of the texture.
        // Also subtract by 0.5 because the center is on the image coordinates (0.5), (0.5).
        float u_p = x * image.width - 0.5;
        float v_p = y * image.height - 0.5;
        // Get the floored index as well as index + 1 to use in the interpolation for both width and height.
        int iu0 = floor(u_p);
        int iu1 = iu0 + 1;
        int iv0 = floor(v_p);
        int iv1 = iv0 + 1;
        // Here we take the decimal parts of the tex coordinates to serve as weights in the interpolation.
        float a_u = (iu1 - u_p);
        float b_u = 1 - a_u;
        float a_v = (iv1 - v_p);
        float b_v = 1 - a_v;
        return a_u * a_v * image.pixels[iv0 * image.width + iu0]
            + a_u * b_v * image.pixels[iv1 * image.width + iu0]
            + b_u * a_v * image.pixels[iv0 * image.width + iu1]
            + b_u * b_v * image.pixels[iv1 * image.width + iu1];
    }
    // Because texCoords is between 0 and 1, have to multiply it by the width/height.
    int x = texCoord.x * image.width;
    // Apparently the image is inverted and because of that we have to take y = 1- y
    int y = (1 - texCoord.y) * image.height;
    
    int index = y * image.width + x;
    glm::vec3 result;
    if (index < 0 || index >= image.pixels.size())
        result = image.pixels[0];
    else
        result = image.pixels[index];
    return result;
}
// Method that puts the original texture as the first level of the mipMap (level 0)
MipMap initialiseMipMap(const Image& image) {
    MipMap mipMap;
    mipMap.level = 0;
    mipMap.height = image.height;
    mipMap.width = image.width;
    mipMap.pixels = image.pixels;
    mipMapPyramid[mipMap.level] = mipMap;
    return mipMap;
}

// Auxiliary method for acquiring the texels when using mimaps, instead of passing the image we just pass the level,
// to know which mipmap level to fetch for in the pyramid.
glm::vec3 acquireTexelMipMap(int level, const glm::vec2& texCoord, const Features& features)
{
    auto mipMapLevel = mipMapPyramid.find(level);
    if (mipMapLevel == mipMapPyramid.end()) {
        // Element not found.
        return glm::vec3(NULL);
    } else {
        MipMap mipmap = mipMapLevel->second;
        int x = texCoord.x * mipmap.width;
        // Apparently the image is inverted and because of that we have to take y = 1- y
        int y = (1 - texCoord.y) * mipmap.height;

        int index = y * mipmap.width + x;
        glm::vec3 result;
        if (index < 0 || index >= mipmap.pixels.size())
            result = mipmap.pixels[0];
        else
            result = mipmap.pixels[index];
        return result;
    }

}


// Auxiliary method that will generate the subsequent mipmap levels.
MipMap generateMipMapPyramid(int height, int width, std::vector<glm::vec3> pixels, int level)
{
    MipMap mipMap;
    mipMap.height = height / 2;
    mipMap.width = width / 2;
    mipMap.level = level;
    // Make a 2x2 box around the pixel and average its value.
    for (int i = 0; i < height; i += 2) {
        for (int j = 0; j < width; j += 2) {
            int indexCol1 = j * width + i;
            int indexCol2 = j * width + (i + 1);
            mipMap.pixels.push_back((pixels[indexCol1] + pixels[indexCol1 + 1] + pixels[indexCol2] + pixels[indexCol2 + 1]) / 4.0f);
        }
    }
    mipMapPyramid[mipMap.level] = mipMap;
    return mipMap;
}