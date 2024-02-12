#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;

struct Node {
    bool isLeaf; 
    AxisAlignedBox box;
    /*
        #1 leaf node: (meshId, triangleId) of contained triangles
        #2 internal node: (left Children, right Children)
    */
    std::vector<std::pair<int, int>> children;
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const Features& features);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

    // Recursive funcion to generate the AABB's
    void build(int node, int depth, std::vector<std::pair<int, int>> meshes, const glm::vec3 splitDirection);

private:
    float getPositionMax(const int meshId, const int triangleId, int coord);
    float getPositionMin(const int meshId, const int triangleId, int coord);
    void medianSplit(
        std::vector<std::pair<int, int>>& leftTriangles, 
        std::vector<std::pair<int, int>>& rightTriangles,
        std::vector<std::pair<float, int>>& centroids,
        const std::vector<std::pair<int, int>>& triangles
    );
    void SAH_BinningSplit(
        std::vector<std::pair<int, int>>& leftTriangles,
        std::vector<std::pair<int, int>>& rightTriangles,
        const std::vector<std::pair<float, int>>& centroids,
        const std::vector<std::pair<int, int>>& triangles,
        const int& node,
        const int& axis
    );
    void buildAxisAlignedBox(AxisAlignedBox& box, std::vector<std::pair<int, int>>& triangles);
    void buildCentroids(
        std::vector<std::pair<float, int>>& centroids,
        const glm::vec3 direction,
        const std::vector<std::pair<int, int>>& triangles
    );
    void build(bool sahBinning);
    void updateBox(AxisAlignedBox& box, const int& meshId, const int& triangleId);

    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Node> m_tree;
    std::vector<int> m_LeafNodes;
};