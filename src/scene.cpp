#include "scene.h"
#include "texture.h"
#include <cmath>
#include <iostream>

Scene loadScenePrebuilt(SceneType type, const std::filesystem::path& dataDir)
{
    Scene scene;
    scene.type = type;
    switch (type) {
    case SingleTriangle: {
        // Load a 3D model with a single triangle
        auto subMeshes = loadMesh(dataDir / "triangle.obj");
        subMeshes[0].material.kd = glm::vec3(1.0f);
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
    } break;
    case Cube: {
        // Load a 3D model of a cube with 12 triangles
        auto subMeshes = loadMesh(dataDir / "cube.obj");
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        // scene.lights.push_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        scene.lights.emplace_back(SegmentLight {
            .endpoint0 = glm::vec3(1.5f, 0.5f, -0.6f),
            .endpoint1 = glm::vec3(-1, 0.5f, -0.5f),
            .color0 = glm::vec3(0.9f, 0.2f, 0.1f), // Red-ish
            .color1 = glm::vec3(0.2f, 1, 0.3f) // Green-ish
        });
    } break;
    case CubeTextured: {
        auto subMeshes = loadMesh(dataDir / "cube-textured.obj");
        for (Mesh& m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1.0, 1.5, -1.0), glm::vec3(1) });
    } break;
    case CornellBox: {
        // Load a 3D model of a Cornell Box
        auto subMeshes = loadMesh(dataDir / "CornellBox-Mirror-Rotated.obj", true);
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(0, 0.58f, 0), glm::vec3(1) }); // Light at the top of the box
    } break;
    case CornellBoxParallelogramLight: {
        // Load a 3D model of a Cornell Box
        auto subMeshes = loadMesh(dataDir / "CornellBox-Mirror-Rotated.obj", true);
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        // Light at the top of the box.
        scene.lights.emplace_back(ParallelogramLight {
            .v0 = glm::vec3(-0.2f, 0.5f, 0),
            .edge01 = glm::vec3(0.4f, 0, 0),
            .edge02 = glm::vec3(0.0f, 0.0f, 0.4f),
            .color0 = glm::vec3(1, 0, 0), // Red
            .color1 = glm::vec3(0, 1, 0), // Green
            .color2 = glm::vec3(0, 0, 1), // Blue
            .color3 = glm::vec3(0, 1, 1), // Purple
        });
    } break;
    case Monkey: {
        // Load a 3D model of a Monkey
        auto subMeshes = loadMesh(dataDir / "monkey.obj", true);
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        scene.lights.emplace_back(PointLight { glm::vec3(1, -1, -1), glm::vec3(1) });
    } break;
    case Teapot: {
        // Load a 3D model of a Teapot
        auto subMeshes = loadMesh(dataDir / "teapot.obj", true);
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
    } break;
    case Dragon: {
        // Load a 3D model of a Dragon
        auto subMeshes = loadMesh(dataDir / "dragon.obj", true);
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
    } break;
    case Spheres: {
        scene.spheres.push_back(Sphere { glm::vec3(3.0f, -2.0f, 10.2f), 1.0f, Material { glm::vec3(0.8f, 0.2f, 0.2f) } });
        scene.spheres.push_back(Sphere { glm::vec3(-2.0f, 2.0f, 4.0f), 2.0f, Material { glm::vec3(0.6f, 0.8f, 0.2f) } });
        scene.spheres.push_back(Sphere { glm::vec3(0.0f, 0.0f, 6.0f), 0.75f, Material { glm::vec3(0.2f, 0.2f, 0.8f) } });
        scene.lights.emplace_back(PointLight { glm::vec3(3, 0, 3), glm::vec3(15) });
    } break;
    case Custom: {
        // === Replace custom.obj by your own 3D model (or call your 3D model custom.obj) ===
        auto subMeshes = loadMesh(dataDir / "custom.obj");
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        // === CHANGE THE LIGHTING IF DESIRED ===
        scene.lights.emplace_back(PointLight { glm::vec3(0, 1, 0), glm::vec3(1) });
        // Spherical light: position, radius, color
        // scene.lights.push_back(SphericalLight{ glm::vec3(0, 1.5f, 0), 0.2f, glm::vec3(1) });
    } break;
    case TextureDebug: {
        auto subMeshes = loadMesh(dataDir / "tdebug.obj");
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        // === CHANGE THE LIGHTING IF DESIRED ===
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        // Spherical light: position, radius, color
        // scene.lights.push_back(SphericalLight{ glm::vec3(0, 1.5f, 0), 0.2f, glm::vec3(1) });
    } break;
    case TransparentDebug: {
        auto subMeshes = loadMesh(dataDir / "transparent.obj", true);
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        //scene.lights.emplace_back(PointLight { glm::vec3(1, -1, -1), glm::vec3(1) });
    } break;
    case Mipmap: {
        auto subMeshes = loadMesh(dataDir / "alias.obj");
        for (Mesh m : subMeshes) {
            createMipMap(m);
        }
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        // scene.lights.emplace_back(PointLight { glm::vec3(1, -1, -1), glm::vec3(1) });
    } break;
    };

    return scene;
}

Scene loadSceneFromFile(const std::filesystem::path& path, const std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight>>& lights)
{
    Scene scene;
    scene.lights = std::move(lights);

    auto subMeshes = loadMesh(path);
    std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));

    return scene;
}
// Auxiliary method that creates the mipmap chain for all the textures, to be used later.
void createMipMap(Mesh mesh)
{
    if (mesh.material.kdTexture) {
        // Generating the mipmap
        Image img = *mesh.material.kdTexture.get();
        int level = 1;
        MipMap mipmap = initialiseMipMap(img);
        // We are only working with square textures, either height/width can be cheked here.
        float i = img.width;
        // Check until we have a mipMap with width/height = 1, the last level.
        while (i > 1) {
            mipmap = generateMipMapPyramid(mipmap.height, mipmap.width, mipmap.pixels, level);
            level++;
            i /= 2;
        }
    }
}
