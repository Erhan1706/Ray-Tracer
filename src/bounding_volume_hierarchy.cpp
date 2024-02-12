#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <queue>

// KEEP ON FALSE, causes LAG otherwise
const bool DEBUG = false;

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene) {
    build(features.extra.enableBvhSahBinning);
}

void BoundingVolumeHierarchy::build(bool sahBinning)
{
    // build the root node
    // build triangles vector
    std::vector<std::pair<int, int>> triangles;

    m_tree.resize(1);

    for (int i = 0; i < m_pScene->meshes.size(); ++i) {
        for (int j = 0; j < m_pScene->meshes[i].triangles.size(); ++j) {
            triangles.push_back({ i, j });
        }
    }

    buildAxisAlignedBox(m_tree[0].box, triangles);

    // meta data
    m_numLevels = std::min(18.0, ceil(log2(triangles.size()))) + 1;

    m_numLeaves = 0;

    m_tree.resize((1 << (m_numLevels + 1)) - 1);

    glm::vec3 oX(1, 0, 0);

    build(0, m_numLevels, triangles, sahBinning ? glm::vec3(0) : oX);
}

float getArea(AxisAlignedBox& box) {
    float area = 1;

    for (int i = 0; i < 3; ++i) {
        area *= (box.upper[i] - box.lower[i]);
    }

    return area;
}

void BoundingVolumeHierarchy::buildAxisAlignedBox(AxisAlignedBox& box, std::vector<std::pair<int, int>>& triangles) {
    box.lower = glm::vec3(std::numeric_limits<float>::max());
    box.upper = glm::vec3(-std::numeric_limits<float>::max());

    for (const auto& it : triangles) {
        updateBox(box, it.first, it.second);
    }
}

/**
 * Recursively generate every node of the tree
 * 
 * param: triangles: {meshId, triangleId}
 *        splitDirection: if <0, 0, 0> then do SAH binning criterion
 * 
 */
void BoundingVolumeHierarchy::build(int node, int depth, std::vector<std::pair<int, int>> triangles, const glm::vec3 splitDirection)
{
    // termination condition
    if (depth == 0 || triangles.size() <= 1) {   // leaf node
        m_tree[node].isLeaf = true;

        m_tree[node].children = triangles;

        ++m_numLeaves;
        m_LeafNodes.push_back(node);

        return;
    }

    std::vector <std::pair<float, int>> centroids; 
    glm::vec3 splitAxis;

    int leftNode = node * 2 + 1;
    int rightNode = node * 2 + 2;

    std::pair<float, int> maxAxis;
    // choose splitAxis
    if (splitDirection == glm::vec3(0)) {
        maxAxis = std::max({
            std::make_pair(
                m_tree[node].box.upper[0] - m_tree[node].box.lower[0],
                0
            ),
            std::make_pair(
                m_tree[node].box.upper[1] - m_tree[node].box.lower[1],
                1
            ),
            std::make_pair(
                m_tree[node].box.upper[2] - m_tree[node].box.lower[2],
                2
            )
        });
        --maxAxis.second;
        float bin_width;
        
        do {
            maxAxis.second = (maxAxis.second + 1) % 3;

            splitAxis = glm::vec3(0);

            bin_width = (m_tree[node].box.upper[maxAxis.second] - m_tree[node].box.lower[maxAxis.second]) / 16.0f;
        } while (bin_width == 0);
        
        
        splitAxis[maxAxis.second] = 1;
    } else {
        splitAxis = splitDirection;
    }

    buildCentroids(centroids, splitAxis, triangles);

    // children nodes
    AxisAlignedBox leftBox, rightBox;
    std::vector<std::pair<int, int>> leftTriangles, rightTriangles;

    if (splitDirection != glm::vec3(0)) {
        medianSplit(leftTriangles, rightTriangles, centroids, triangles);
    } else {
        SAH_BinningSplit(leftTriangles, rightTriangles, centroids, triangles, node, maxAxis.second);
    }

    buildAxisAlignedBox(leftBox, leftTriangles);
    buildAxisAlignedBox(rightBox, rightTriangles);
    
    //--- finished computation ---

    // free up the memory
    triangles.resize(0);
    triangles.shrink_to_fit();

    centroids.resize(0);
    centroids.shrink_to_fit();
    
    // update current node

    m_tree[node].children.push_back({ leftNode, rightNode });

    m_tree[leftNode].box = leftBox;
    m_tree[rightNode].box = rightBox;

    m_tree[node].isLeaf = false;

    // next split direction
    // x -> y -> z -> x
    glm::vec3 nextSplitDirection = glm::vec3(0);

    if (splitDirection.x)
        nextSplitDirection.y = 1;
    else if (splitDirection.y)
        nextSplitDirection.z = 1;
    else if (splitDirection.z)
        nextSplitDirection.x = 1;

    // recursively call generate on the two new nodes

    // left Node
    build(leftNode, depth - 1, leftTriangles, nextSplitDirection);

    //right Node
    build(rightNode, depth - 1, rightTriangles, nextSplitDirection);
}

void BoundingVolumeHierarchy::updateBox(AxisAlignedBox& box, const int& meshId, const int& triangleId) {
    for (int p = 0; p < 3; ++p) {
        // upper
        box.upper[p] = std::max(box.upper[p], getPositionMax(meshId, triangleId, p));

        // lower
        box.lower[p] = std::min(box.lower[p], getPositionMin(meshId, triangleId, p));
    }
}

// merge B into A
void mergeBox(AxisAlignedBox& A, const AxisAlignedBox& B)
{
    for (int p = 0; p < 3; ++p) {
        A.lower[p] = std::min(A.lower[p], B.lower[p]);
        A.upper[p] = std::max(A.upper[p], B.upper[p]);
    }
}

void BoundingVolumeHierarchy::SAH_BinningSplit(
    std::vector<std::pair<int, int>>& leftTriangles,
    std::vector<std::pair<int, int>>& rightTriangles,
    const std::vector<std::pair<float, int>>& centroids,
    const std::vector<std::pair<int, int>>& triangles,
    const int& node,
    const int& axis
) {
    auto getCost = [](std::pair<AxisAlignedBox, int> X) {
        if (!X.second)
            return 0.0f;
        return getArea(X.first) * X.second;
    };

    const float bin_width = (m_tree[node].box.upper[axis] - m_tree[node].box.lower[axis]) / 16.0f;
    int splitBin;   // bin where we split
    float costMin = std::numeric_limits<float>::max();

    std::pair<AxisAlignedBox, int> bins[18], right[18], left;

    // reset data
    for (int i = 0; i < 18; ++i) {
        bins[i].first = { glm::vec3(std::numeric_limits<float>::max()), glm::vec3(-std::numeric_limits<float>::max()) };
        bins[i].second = 0;

        right[i].first = { glm::vec3(std::numeric_limits<float>::max()), glm::vec3(-std::numeric_limits<float>::max()) };
        right[i].second = 0;
    }

    left.first = { glm::vec3(std::numeric_limits<float>::max()), glm::vec3(-std::numeric_limits<float>::max()) };
    left.second = 0;

    // put triangles in bins
    for (const auto& it : centroids) {
        int index = (it.first - m_tree[node].box.lower[axis]) / bin_width;

        ++bins[index].second;

        updateBox(bins[index].first, triangles[it.second].first, triangles[it.second].second);
    }

    for (int i = 16; i >= 0; --i) {
        right[i] = right[i + 1];

        mergeBox(right[i].first, bins[i].first);
        right[i].second += bins[i].second;
    }

    for (int i = 0; i < 16; ++i) {
        float cost = 0;

        // update left
        
        mergeBox(left.first, bins[i].first);
        left.second += bins[i].second;

        cost = getCost(left) + getCost(right[i+1]);

        if (cost < costMin) {
            costMin = cost;
            splitBin = i + 1;
        }
    }

    // build left & right triangles
    for (const auto& it : centroids) {
        if ((it.first - m_tree[node].box.lower[axis]) / bin_width < splitBin) {
            leftTriangles.push_back(triangles[it.second]);
        } else {
            rightTriangles.push_back(triangles[it.second]);
        }
    }
}

void BoundingVolumeHierarchy::buildCentroids(
    std::vector<std::pair<float, int>>& centroids,
    const glm::vec3 direction,
    const std::vector<std::pair<int, int>>& triangles
) {
    for (int i = 0; i < triangles.size(); ++i) {
        Mesh* mesh = &(m_pScene->meshes[triangles[i].first]);

        glm::uvec3 triangle = mesh->triangles[triangles[i].second];

        glm::vec3 centroid = (mesh->vertices[triangle.x].position
                                 + mesh->vertices[triangle.y].position
                                 + mesh->vertices[triangle.z].position)
            / 3.0f;

        centroids.push_back({ glm::dot(centroid * direction, glm::vec3(1)), i });
    }
}

void BoundingVolumeHierarchy::medianSplit(
    std::vector<std::pair<int, int>>& leftTriangles,
    std::vector<std::pair<int, int>>& rightTriangles,
    std::vector<std::pair<float, int>>& centroids,
    const std::vector<std::pair<int, int>>& triangles
) {
    std::sort(centroids.begin(), centroids.end());

    for (int i = 0; i < centroids.size(); ++i) {
        int meshId = triangles[centroids[i].second].first;
        int triangleId = triangles[centroids[i].second].second;

        if (i < centroids.size() / 2) {
            leftTriangles.push_back({ meshId, triangleId });
        } else {
            rightTriangles.push_back({ meshId, triangleId });
        }
    }
}

float BoundingVolumeHierarchy::getPositionMax(const int meshId, const int triangleId, int coord)
{
    return std::max({ 
            m_pScene->meshes[meshId].vertices[m_pScene->meshes[meshId].triangles[triangleId][0]].position[coord],
            m_pScene->meshes[meshId].vertices[m_pScene->meshes[meshId].triangles[triangleId][1]].position[coord],
            m_pScene->meshes[meshId].vertices[m_pScene->meshes[meshId].triangles[triangleId][2]].position[coord] 
        });
}

float BoundingVolumeHierarchy::getPositionMin(const int meshId, const int triangleId, int coord)
{
    return std::min({ m_pScene->meshes[meshId].vertices[m_pScene->meshes[meshId].triangles[triangleId][0]].position[coord],
        m_pScene->meshes[meshId].vertices[m_pScene->meshes[meshId].triangles[triangleId][1]].position[coord],
        m_pScene->meshes[meshId].vertices[m_pScene->meshes[meshId].triangles[triangleId][2]].position[coord] });
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    // - 1 because else the last leaf in the interface goes outside vector range and crashes the program
    return m_numLeaves;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    level = std::min(level, m_numLevels);

    int startingIndex = 0;
    for (int i = 0; i < level; i++) {
        startingIndex += pow(2, i);
    }

    int cap = startingIndex + pow(2, level);

    AxisAlignedBox aabb;

    for (startingIndex; startingIndex < cap; startingIndex++) {
        aabb = { m_tree[startingIndex].box.lower, m_tree[startingIndex].box.upper };

        if (m_tree[startingIndex].children.size())
            drawAABB(aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.8f);
    }
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    leafIdx = std::min(leafIdx, (int)m_LeafNodes.size() - 1);

    drawAABB(m_tree[m_LeafNodes[leafIdx]].box, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 1.0f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
// !!! DOESNT update HitInfo.barycentric coord
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    bool hit = false;

    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hit = true;
                    glm::vec3 surfacePoint = ray.origin + ray.t * ray.direction;
                    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, surfacePoint);
                    hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                    
                    if (features.enableNormalInterpDebug) {
                        Ray normal0 { v0.position, v0.normal, 1.0f };
                        Ray normal1 { v1.position, v1.normal, 1.0f };
                        Ray normal2 { v2.position, v2.normal, 1.0f };
                        glm::vec3 interpolatedNormal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
                        Ray interpolatedNormalRay { surfacePoint, interpolatedNormal, 1.0f };

                        drawRay(normal0, { 0.0f, 1.0f, 0.0f });
                        drawRay(normal1, { 0.0f, 1.0f, 0.0f });
                        drawRay(normal2, { 0.0f, 1.0f, 0.0f });
                        drawRay(interpolatedNormalRay, { 1.0f, 0.0f, 0.0f });
                    }
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else { 
        if (!intersectRayWithShape(m_tree[0].box, ray)) {
            return false;
        }
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> pq;
        int triangleNode, node, meshId;
        glm::vec3 surfacePoint;
        std::pair<int, int> triangleHit;
        Vertex *A, *B, *C;
        float old_t, t;

        // ray.t, node
        ray.t = std::numeric_limits<float>::max();
        pq.push({ ray.t, 0 });

        while (!pq.empty()) {
            t = pq.top().first;
            node = pq.top().second;

            pq.pop();
            
            if (t > ray.t) {
                break;
            }

            if (m_tree[node].isLeaf) {

                for (const auto& it : m_tree[node].children) {
                    A = &m_pScene->meshes[it.first].vertices[m_pScene->meshes[it.first].triangles[it.second][0]];
                    B = &m_pScene->meshes[it.first].vertices[m_pScene->meshes[it.first].triangles[it.second][1]];
                    C = &m_pScene->meshes[it.first].vertices[m_pScene->meshes[it.first].triangles[it.second][2]];

                    old_t = ray.t;
                    ray.t = std::numeric_limits<float>::max();

                    if (intersectRayWithTriangle(A->position, B->position, C->position, ray, hitInfo)) {
                        hit = true;

                        if (ray.t < old_t) {
                            triangleHit = it;
                            triangleNode = node;
                                             
                            if (features.enableNormalInterp && DEBUG) {
                                Ray normal0 { A->position, A->normal, 1.0f };
                                Ray normal1 { B->position, B->normal, 1.0f };
                                Ray normal2 { C->position, C->normal, 1.0f };
                                glm::vec3 interpolatedNormal = interpolateNormal(A->normal, B->normal, C->normal, hitInfo.barycentricCoord);
                                Ray interpolatedNormalRay { surfacePoint, interpolatedNormal, 1.0f };

                                drawRay(normal0, { 0.0f, 1.0f, 0.0f });
                                drawRay(normal1, { 0.0f, 1.0f, 0.0f });
                                drawRay(normal2, { 0.0f, 1.0f, 0.0f });
                                drawRay(interpolatedNormalRay, { 1.0f, 0.0f, 0.0f });
                            }

                            old_t = ray.t;
                        } 
                    }

                    ray.t = old_t;
                }

                continue;
            }

            old_t = ray.t;

            ray.t = std::numeric_limits<float>::max();
            if (intersectRayWithShape(m_tree[node * 2 + 1].box, ray)) {
                if (ray.origin.x >= m_tree[node * 2 + 1].box.lower.x && ray.origin.x <= m_tree[node * 2 + 1].box.upper.x
                    && ray.origin.y >= m_tree[node * 2 + 1].box.lower.y && ray.origin.y <= m_tree[node * 2 + 1].box.upper.y
                    && ray.origin.z >= m_tree[node * 2 + 1].box.lower.z && ray.origin.z <= m_tree[node * 2 + 1].box.upper.z) 
                {    
                    ray.t = 0;
                }

                pq.push({ ray.t, node * 2 + 1 });
            }

            ray.t = std::numeric_limits<float>::max();
            if (intersectRayWithShape(m_tree[node * 2 + 2].box, ray)) {
                if (ray.origin.x >= m_tree[node * 2 + 2].box.lower.x && ray.origin.x <= m_tree[node * 2 + 2].box.upper.x
                    && ray.origin.y >= m_tree[node * 2 + 2].box.lower.y && ray.origin.y <= m_tree[node * 2 + 2].box.upper.y
                    && ray.origin.z >= m_tree[node * 2 + 2].box.lower.z && ray.origin.z <= m_tree[node * 2 + 2].box.upper.z) {
                    ray.t = 0;
                }

                pq.push({ ray.t, node * 2 + 2 });
            }

            ray.t = old_t;
        }

        // update hitInfo
        if (hit) {
            A = &m_pScene->meshes[triangleHit.first].vertices[m_pScene->meshes[triangleHit.first].triangles[triangleHit.second][0]];
            B = &m_pScene->meshes[triangleHit.first].vertices[m_pScene->meshes[triangleHit.first].triangles[triangleHit.second][1]];
            C = &m_pScene->meshes[triangleHit.first].vertices[m_pScene->meshes[triangleHit.first].triangles[triangleHit.second][2]];

            ray.t = std::numeric_limits<float>::max();

            // hit the triangle again
            intersectRayWithTriangle(A->position, B->position, C->position, ray, hitInfo);

            surfacePoint = ray.origin + ray.t * ray.direction;

            hitInfo.barycentricCoord = computeBarycentricCoord(A->position, B->position, C->position, surfacePoint);
            hitInfo.texCoord = interpolateTexCoord(A->texCoord, B->texCoord, C->texCoord, hitInfo.barycentricCoord);

            hitInfo.material = m_pScene->meshes[triangleHit.first].material;

            if (features.enableNormalInterp) {
                hitInfo.normal = interpolateNormal(A->normal, B->normal, C->normal, hitInfo.barycentricCoord);
            }
        }
        
        // press 'R' under RASTERIZATION to see the ray visual debug
        // make sure that BVH flag is 'ON'
        if (hit && DEBUG) {
            float eps = 0.05;
 
            while (true) {
                drawAABB(m_tree[triangleNode].box, DrawMode::Wireframe, glm::vec3(1.0f - eps, 0.05f, eps), 1.0f - eps);

                if (!triangleNode)
                    break;

                eps += 0.04;
                triangleNode = (triangleNode - 1) / 2;
            }

            A = &m_pScene->meshes[triangleHit.first].vertices[m_pScene->meshes[triangleHit.first].triangles[triangleHit.second][0]];
            B = &m_pScene->meshes[triangleHit.first].vertices[m_pScene->meshes[triangleHit.first].triangles[triangleHit.second][1]];
            C = &m_pScene->meshes[triangleHit.first].vertices[m_pScene->meshes[triangleHit.first].triangles[triangleHit.second][2]];

            drawTriangle(*A, *B, *C);
        }

        return hit;
    }
}