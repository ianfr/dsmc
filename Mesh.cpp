//
// Mesh.cpp
// 3D Mesh loader implementation for DSMC
//

#include "Mesh.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <map>

bool Mesh::loadOBJ(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open OBJ file: " << filename << std::endl;
        return false;
    }
    
    m_vertices.clear();
    m_normals.clear();
    m_triangles.clear();
    
    std::vector<std::array<float, 3>> tempNormals;
    std::vector<std::array<int, 3>> faceVertices;
    std::vector<std::array<int, 3>> faceNormals;
    
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        
        if (prefix == "v") {
            // Vertex position
            float x, y, z;
            iss >> x >> y >> z;
            m_vertices.push_back({x, y, z});
        }
        else if (prefix == "vn") {
            // Vertex normal
            float nx, ny, nz;
            iss >> nx >> ny >> nz;
            tempNormals.push_back({nx, ny, nz});
        }
        else if (prefix == "f") {
            // Face (triangle or polygon)
            std::vector<int> vIndices;
            std::vector<int> nIndices;
            
            std::string vertexData;
            while (iss >> vertexData) {
                int vIdx = 0, tIdx = 0, nIdx = 0;
                
                // Parse "v", "v/t", "v/t/n", or "v//n" formats
                size_t firstSlash = vertexData.find('/');
                if (firstSlash == std::string::npos) {
                    vIdx = std::stoi(vertexData);
                } else {
                    vIdx = std::stoi(vertexData.substr(0, firstSlash));
                    size_t secondSlash = vertexData.find('/', firstSlash + 1);
                    if (secondSlash != std::string::npos) {
                        if (secondSlash > firstSlash + 1) {
                            tIdx = std::stoi(vertexData.substr(firstSlash + 1, secondSlash - firstSlash - 1));
                        }
                        nIdx = std::stoi(vertexData.substr(secondSlash + 1));
                    } else {
                        tIdx = std::stoi(vertexData.substr(firstSlash + 1));
                    }
                }
                
                // OBJ indices are 1-based
                vIndices.push_back(vIdx - 1);
                if (nIdx > 0) {
                    nIndices.push_back(nIdx - 1);
                }
            }
            
            // Triangulate polygon (fan triangulation)
            for (size_t i = 1; i + 1 < vIndices.size(); i++) {
                faceVertices.push_back({vIndices[0], vIndices[i], vIndices[i + 1]});
                if (!nIndices.empty()) {
                    faceNormals.push_back({nIndices[0], nIndices[i], nIndices[i + 1]});
                }
            }
        }
    }
    
    file.close();
    
    // Build triangles
    m_triangles.resize(faceVertices.size());
    for (size_t i = 0; i < faceVertices.size(); i++) {
        const auto& v0 = m_vertices[faceVertices[i][0]];
        const auto& v1 = m_vertices[faceVertices[i][1]];
        const auto& v2 = m_vertices[faceVertices[i][2]];
        
        m_triangles[i].v0[0] = v0[0];
        m_triangles[i].v0[1] = v0[1];
        m_triangles[i].v0[2] = v0[2];
        
        m_triangles[i].v1[0] = v1[0];
        m_triangles[i].v1[1] = v1[1];
        m_triangles[i].v1[2] = v1[2];
        
        m_triangles[i].v2[0] = v2[0];
        m_triangles[i].v2[1] = v2[1];
        m_triangles[i].v2[2] = v2[2];
    }
    
    // Compute normals
    computeTriangleNormals();
    
    std::cout << "Loaded OBJ: " << m_vertices.size() << " vertices, " 
              << m_triangles.size() << " triangles" << std::endl;
    
    return true;
}

bool Mesh::loadSTL(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Failed to open STL file: " << filename << std::endl;
        return false;
    }
    
    m_vertices.clear();
    m_normals.clear();
    m_triangles.clear();
    
    // Check if binary or ASCII
    char header[80];
    file.read(header, 80);
    
    // Read number of triangles
    uint32_t numTriangles;
    file.read(reinterpret_cast<char*>(&numTriangles), 4);
    
    // Check if this looks like a valid binary STL
    std::streampos currentPos = file.tellg();
    file.seekg(0, std::ios::end);
    std::streampos fileSize = file.tellg();
    file.seekg(currentPos);
    
    // Binary STL: 80 + 4 + (50 * numTriangles) bytes
    bool isBinary = (fileSize == std::streampos(84 + 50 * numTriangles));
    
    if (isBinary) {
        // Binary STL
        m_triangles.resize(numTriangles);
        
        for (uint32_t i = 0; i < numTriangles; i++) {
            float normal[3], v0[3], v1[3], v2[3];
            uint16_t attrib;
            
            file.read(reinterpret_cast<char*>(normal), 12);
            file.read(reinterpret_cast<char*>(v0), 12);
            file.read(reinterpret_cast<char*>(v1), 12);
            file.read(reinterpret_cast<char*>(v2), 12);
            file.read(reinterpret_cast<char*>(&attrib), 2);
            
            m_triangles[i].v0[0] = v0[0];
            m_triangles[i].v0[1] = v0[1];
            m_triangles[i].v0[2] = v0[2];
            
            m_triangles[i].v1[0] = v1[0];
            m_triangles[i].v1[1] = v1[1];
            m_triangles[i].v1[2] = v1[2];
            
            m_triangles[i].v2[0] = v2[0];
            m_triangles[i].v2[1] = v2[1];
            m_triangles[i].v2[2] = v2[2];
            
            m_triangles[i].normal[0] = normal[0];
            m_triangles[i].normal[1] = normal[1];
            m_triangles[i].normal[2] = normal[2];
        }
    } else {
        // ASCII STL - reopen as text
        file.close();
        std::ifstream asciiFile(filename);
        
        std::string line;
        Triangle currentTriangle;
        int vertexIndex = 0;
        
        while (std::getline(asciiFile, line)) {
            std::istringstream iss(line);
            std::string keyword;
            iss >> keyword;
            
            if (keyword == "facet") {
                std::string normal;
                iss >> normal;
                iss >> currentTriangle.normal[0] 
                    >> currentTriangle.normal[1] 
                    >> currentTriangle.normal[2];
                vertexIndex = 0;
            }
            else if (keyword == "vertex") {
                float x, y, z;
                iss >> x >> y >> z;
                
                if (vertexIndex == 0) {
                    currentTriangle.v0[0] = x;
                    currentTriangle.v0[1] = y;
                    currentTriangle.v0[2] = z;
                } else if (vertexIndex == 1) {
                    currentTriangle.v1[0] = x;
                    currentTriangle.v1[1] = y;
                    currentTriangle.v1[2] = z;
                } else if (vertexIndex == 2) {
                    currentTriangle.v2[0] = x;
                    currentTriangle.v2[1] = y;
                    currentTriangle.v2[2] = z;
                }
                vertexIndex++;
            }
            else if (keyword == "endfacet") {
                m_triangles.push_back(currentTriangle);
            }
        }
        
        asciiFile.close();
    }
    
    file.close();
    
    // Ensure normals are valid
    computeTriangleNormals();
    
    std::cout << "Loaded STL: " << m_triangles.size() << " triangles" << std::endl;
    
    return true;
}

bool Mesh::load(const std::string& filename) {
    std::string ext = filename.substr(filename.find_last_of('.') + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    
    if (ext == "obj") {
        return loadOBJ(filename);
    } else if (ext == "stl") {
        return loadSTL(filename);
    } else {
        std::cerr << "Unsupported mesh format: " << ext << std::endl;
        return false;
    }
}

void Mesh::computeTriangleNormals() {
    for (auto& tri : m_triangles) {
        // Edge vectors
        float e1[3] = {
            tri.v1[0] - tri.v0[0],
            tri.v1[1] - tri.v0[1],
            tri.v1[2] - tri.v0[2]
        };
        float e2[3] = {
            tri.v2[0] - tri.v0[0],
            tri.v2[1] - tri.v0[1],
            tri.v2[2] - tri.v0[2]
        };
        
        // Cross product
        float nx = e1[1] * e2[2] - e1[2] * e2[1];
        float ny = e1[2] * e2[0] - e1[0] * e2[2];
        float nz = e1[0] * e2[1] - e1[1] * e2[0];
        
        // Normalize
        float len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (len > 1e-8f) {
            tri.normal[0] = nx / len;
            tri.normal[1] = ny / len;
            tri.normal[2] = nz / len;
        } else {
            tri.normal[0] = 0;
            tri.normal[1] = 1;
            tri.normal[2] = 0;
        }
    }
}

void Mesh::translate(float x, float y, float z) {
    for (auto& tri : m_triangles) {
        tri.v0[0] += x; tri.v0[1] += y; tri.v0[2] += z;
        tri.v1[0] += x; tri.v1[1] += y; tri.v1[2] += z;
        tri.v2[0] += x; tri.v2[1] += y; tri.v2[2] += z;
    }
}

void Mesh::scale(float sx, float sy, float sz) {
    for (auto& tri : m_triangles) {
        tri.v0[0] *= sx; tri.v0[1] *= sy; tri.v0[2] *= sz;
        tri.v1[0] *= sx; tri.v1[1] *= sy; tri.v1[2] *= sz;
        tri.v2[0] *= sx; tri.v2[1] *= sy; tri.v2[2] *= sz;
    }
    computeTriangleNormals();
}

void Mesh::rotateX(float angle) {
    float c = std::cos(angle);
    float s = std::sin(angle);
    
    for (auto& tri : m_triangles) {
        float y, z;
        
        y = tri.v0[1]; z = tri.v0[2];
        tri.v0[1] = c * y - s * z;
        tri.v0[2] = s * y + c * z;
        
        y = tri.v1[1]; z = tri.v1[2];
        tri.v1[1] = c * y - s * z;
        tri.v1[2] = s * y + c * z;
        
        y = tri.v2[1]; z = tri.v2[2];
        tri.v2[1] = c * y - s * z;
        tri.v2[2] = s * y + c * z;
    }
    computeTriangleNormals();
}

void Mesh::rotateY(float angle) {
    float c = std::cos(angle);
    float s = std::sin(angle);
    
    for (auto& tri : m_triangles) {
        float x, z;
        
        x = tri.v0[0]; z = tri.v0[2];
        tri.v0[0] = c * x + s * z;
        tri.v0[2] = -s * x + c * z;
        
        x = tri.v1[0]; z = tri.v1[2];
        tri.v1[0] = c * x + s * z;
        tri.v1[2] = -s * x + c * z;
        
        x = tri.v2[0]; z = tri.v2[2];
        tri.v2[0] = c * x + s * z;
        tri.v2[2] = -s * x + c * z;
    }
    computeTriangleNormals();
}

void Mesh::rotateZ(float angle) {
    float c = std::cos(angle);
    float s = std::sin(angle);
    
    for (auto& tri : m_triangles) {
        float x, y;
        
        x = tri.v0[0]; y = tri.v0[1];
        tri.v0[0] = c * x - s * y;
        tri.v0[1] = s * x + c * y;
        
        x = tri.v1[0]; y = tri.v1[1];
        tri.v1[0] = c * x - s * y;
        tri.v1[1] = s * x + c * y;
        
        x = tri.v2[0]; y = tri.v2[1];
        tri.v2[0] = c * x - s * y;
        tri.v2[1] = s * x + c * y;
    }
    computeTriangleNormals();
}

void Mesh::getBoundingBox(float& minX, float& maxX, 
                          float& minY, float& maxY, 
                          float& minZ, float& maxZ) const {
    minX = minY = minZ = std::numeric_limits<float>::max();
    maxX = maxY = maxZ = std::numeric_limits<float>::lowest();
    
    for (const auto& tri : m_triangles) {
        for (const float* v : {tri.v0, tri.v1, tri.v2}) {
            minX = std::min(minX, v[0]);
            maxX = std::max(maxX, v[0]);
            minY = std::min(minY, v[1]);
            maxY = std::max(maxY, v[1]);
            minZ = std::min(minZ, v[2]);
            maxZ = std::max(maxZ, v[2]);
        }
    }
}

void Mesh::centerAtOrigin() {
    float minX, maxX, minY, maxY, minZ, maxZ;
    getBoundingBox(minX, maxX, minY, maxY, minZ, maxZ);
    
    float cx = (minX + maxX) / 2.0f;
    float cy = (minY + maxY) / 2.0f;
    float cz = (minZ + maxZ) / 2.0f;
    
    translate(-cx, -cy, -cz);
}

void Mesh::fitToBounds(float minX, float maxX, float minY, float maxY, float minZ, float maxZ) {
    float meshMinX, meshMaxX, meshMinY, meshMaxY, meshMinZ, meshMaxZ;
    getBoundingBox(meshMinX, meshMaxX, meshMinY, meshMaxY, meshMinZ, meshMaxZ);
    
    float meshSizeX = meshMaxX - meshMinX;
    float meshSizeY = meshMaxY - meshMinY;
    float meshSizeZ = meshMaxZ - meshMinZ;
    
    float targetSizeX = maxX - minX;
    float targetSizeY = maxY - minY;
    float targetSizeZ = maxZ - minZ;
    
    // Find uniform scale to fit within bounds
    float scaleX = (meshSizeX > 0) ? targetSizeX / meshSizeX : 1.0f;
    float scaleY = (meshSizeY > 0) ? targetSizeY / meshSizeY : 1.0f;
    float scaleZ = (meshSizeZ > 0) ? targetSizeZ / meshSizeZ : 1.0f;
    float uniformScale = std::min({scaleX, scaleY, scaleZ});
    
    // Center, scale, then translate to target center
    centerAtOrigin();
    scale(uniformScale);
    
    float targetCenterX = (minX + maxX) / 2.0f;
    float targetCenterY = (minY + maxY) / 2.0f;
    float targetCenterZ = (minZ + maxZ) / 2.0f;
    translate(targetCenterX, targetCenterY, targetCenterZ);
}

// Create a sphere using icosahedron subdivision
Mesh Mesh::createSphere(float radius, int subdivisions) {
    Mesh mesh;
    
    // Golden ratio
    const float phi = (1.0f + std::sqrt(5.0f)) / 2.0f;
    const float scale = radius / std::sqrt(1.0f + phi * phi);
    
    // Icosahedron vertices
    std::vector<std::array<float, 3>> vertices = {
        {-1,  phi, 0}, { 1,  phi, 0}, {-1, -phi, 0}, { 1, -phi, 0},
        { 0, -1,  phi}, { 0,  1,  phi}, { 0, -1, -phi}, { 0,  1, -phi},
        { phi, 0, -1}, { phi, 0,  1}, {-phi, 0, -1}, {-phi, 0,  1}
    };
    
    // Normalize and scale
    for (auto& v : vertices) {
        float len = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        v[0] = v[0] / len * radius;
        v[1] = v[1] / len * radius;
        v[2] = v[2] / len * radius;
    }
    
    // Icosahedron faces
    std::vector<std::array<int, 3>> faces = {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
    };
    
    // Subdivide
    for (int s = 0; s < subdivisions; s++) {
        std::vector<std::array<int, 3>> newFaces;
        std::map<std::pair<int, int>, int> midpointCache;
        
        auto getMidpoint = [&](int i1, int i2) -> int {
            auto key = std::make_pair(std::min(i1, i2), std::max(i1, i2));
            auto it = midpointCache.find(key);
            if (it != midpointCache.end()) {
                return it->second;
            }
            
            // Create midpoint
            std::array<float, 3> mid = {
                (vertices[i1][0] + vertices[i2][0]) / 2.0f,
                (vertices[i1][1] + vertices[i2][1]) / 2.0f,
                (vertices[i1][2] + vertices[i2][2]) / 2.0f
            };
            
            // Normalize to sphere surface
            float len = std::sqrt(mid[0]*mid[0] + mid[1]*mid[1] + mid[2]*mid[2]);
            mid[0] = mid[0] / len * radius;
            mid[1] = mid[1] / len * radius;
            mid[2] = mid[2] / len * radius;
            
            int idx = vertices.size();
            vertices.push_back(mid);
            midpointCache[key] = idx;
            return idx;
        };
        
        for (const auto& face : faces) {
            int a = getMidpoint(face[0], face[1]);
            int b = getMidpoint(face[1], face[2]);
            int c = getMidpoint(face[2], face[0]);
            
            newFaces.push_back({face[0], a, c});
            newFaces.push_back({face[1], b, a});
            newFaces.push_back({face[2], c, b});
            newFaces.push_back({a, b, c});
        }
        
        faces = newFaces;
    }
    
    // Convert to triangles
    mesh.m_triangles.resize(faces.size());
    for (size_t i = 0; i < faces.size(); i++) {
        const auto& v0 = vertices[faces[i][0]];
        const auto& v1 = vertices[faces[i][1]];
        const auto& v2 = vertices[faces[i][2]];
        
        mesh.m_triangles[i].v0[0] = v0[0];
        mesh.m_triangles[i].v0[1] = v0[1];
        mesh.m_triangles[i].v0[2] = v0[2];
        
        mesh.m_triangles[i].v1[0] = v1[0];
        mesh.m_triangles[i].v1[1] = v1[1];
        mesh.m_triangles[i].v1[2] = v1[2];
        
        mesh.m_triangles[i].v2[0] = v2[0];
        mesh.m_triangles[i].v2[1] = v2[1];
        mesh.m_triangles[i].v2[2] = v2[2];
    }
    
    mesh.computeTriangleNormals();
    return mesh;
}

Mesh Mesh::createBox(float width, float height, float depth) {
    Mesh mesh;
    
    float hw = width / 2.0f;
    float hh = height / 2.0f;
    float hd = depth / 2.0f;
    
    // 8 vertices of the box
    std::array<float, 3> v[8] = {
        {-hw, -hh, -hd}, { hw, -hh, -hd}, { hw,  hh, -hd}, {-hw,  hh, -hd},
        {-hw, -hh,  hd}, { hw, -hh,  hd}, { hw,  hh,  hd}, {-hw,  hh,  hd}
    };
    
    // 12 triangles (2 per face)
    int faces[12][3] = {
        // Front
        {4, 5, 6}, {4, 6, 7},
        // Back
        {1, 0, 3}, {1, 3, 2},
        // Top
        {7, 6, 2}, {7, 2, 3},
        // Bottom
        {0, 1, 5}, {0, 5, 4},
        // Right
        {5, 1, 2}, {5, 2, 6},
        // Left
        {0, 4, 7}, {0, 7, 3}
    };
    
    mesh.m_triangles.resize(12);
    for (int i = 0; i < 12; i++) {
        mesh.m_triangles[i].v0[0] = v[faces[i][0]][0];
        mesh.m_triangles[i].v0[1] = v[faces[i][0]][1];
        mesh.m_triangles[i].v0[2] = v[faces[i][0]][2];
        
        mesh.m_triangles[i].v1[0] = v[faces[i][1]][0];
        mesh.m_triangles[i].v1[1] = v[faces[i][1]][1];
        mesh.m_triangles[i].v1[2] = v[faces[i][1]][2];
        
        mesh.m_triangles[i].v2[0] = v[faces[i][2]][0];
        mesh.m_triangles[i].v2[1] = v[faces[i][2]][1];
        mesh.m_triangles[i].v2[2] = v[faces[i][2]][2];
    }
    
    mesh.computeTriangleNormals();
    return mesh;
}

Mesh Mesh::createCylinder(float radius, float height, int segments) {
    Mesh mesh;
    
    float halfHeight = height / 2.0f;
    
    // Create vertices for top and bottom circles
    std::vector<std::array<float, 3>> topVerts, bottomVerts;
    
    for (int i = 0; i < segments; i++) {
        float angle = 2.0f * M_PI * i / segments;
        float x = radius * std::cos(angle);
        float z = radius * std::sin(angle);
        
        topVerts.push_back({x, halfHeight, z});
        bottomVerts.push_back({x, -halfHeight, z});
    }
    
    // Side triangles
    for (int i = 0; i < segments; i++) {
        int next = (i + 1) % segments;
        
        Triangle t1, t2;
        
        // First triangle
        t1.v0[0] = bottomVerts[i][0];
        t1.v0[1] = bottomVerts[i][1];
        t1.v0[2] = bottomVerts[i][2];
        
        t1.v1[0] = topVerts[i][0];
        t1.v1[1] = topVerts[i][1];
        t1.v1[2] = topVerts[i][2];
        
        t1.v2[0] = topVerts[next][0];
        t1.v2[1] = topVerts[next][1];
        t1.v2[2] = topVerts[next][2];
        
        mesh.m_triangles.push_back(t1);
        
        // Second triangle
        t2.v0[0] = bottomVerts[i][0];
        t2.v0[1] = bottomVerts[i][1];
        t2.v0[2] = bottomVerts[i][2];
        
        t2.v1[0] = topVerts[next][0];
        t2.v1[1] = topVerts[next][1];
        t2.v1[2] = topVerts[next][2];
        
        t2.v2[0] = bottomVerts[next][0];
        t2.v2[1] = bottomVerts[next][1];
        t2.v2[2] = bottomVerts[next][2];
        
        mesh.m_triangles.push_back(t2);
    }
    
    // Top cap (fan triangulation)
    for (int i = 1; i < segments - 1; i++) {
        Triangle t;
        t.v0[0] = topVerts[0][0];
        t.v0[1] = topVerts[0][1];
        t.v0[2] = topVerts[0][2];
        
        t.v1[0] = topVerts[i][0];
        t.v1[1] = topVerts[i][1];
        t.v1[2] = topVerts[i][2];
        
        t.v2[0] = topVerts[i + 1][0];
        t.v2[1] = topVerts[i + 1][1];
        t.v2[2] = topVerts[i + 1][2];
        
        mesh.m_triangles.push_back(t);
    }
    
    // Bottom cap (fan triangulation, reversed winding)
    for (int i = 1; i < segments - 1; i++) {
        Triangle t;
        t.v0[0] = bottomVerts[0][0];
        t.v0[1] = bottomVerts[0][1];
        t.v0[2] = bottomVerts[0][2];
        
        t.v1[0] = bottomVerts[i + 1][0];
        t.v1[1] = bottomVerts[i + 1][1];
        t.v1[2] = bottomVerts[i + 1][2];
        
        t.v2[0] = bottomVerts[i][0];
        t.v2[1] = bottomVerts[i][1];
        t.v2[2] = bottomVerts[i][2];
        
        mesh.m_triangles.push_back(t);
    }
    
    mesh.computeTriangleNormals();
    return mesh;
}
