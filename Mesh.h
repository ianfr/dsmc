//
// Mesh.h
// 3D Mesh loader for DSMC
//

#ifndef DSMC_MESH_H
#define DSMC_MESH_H

#include <vector>
#include <string>
#include <array>

// Triangle structure for mesh data (must match MetalCompute.h)
struct Triangle {
    float v0[3];
    float v1[3];
    float v2[3];
    float normal[3];
};

class Mesh {
public:
    Mesh() = default;
    ~Mesh() = default;
    
    // Load mesh from file (OBJ format supported)
    bool loadOBJ(const std::string& filename);
    
    // Load mesh from STL file
    bool loadSTL(const std::string& filename);
    
    // Auto-detect format and load
    bool load(const std::string& filename);
    
    // Transform mesh
    void translate(float x, float y, float z);
    void scale(float sx, float sy, float sz);
    void scale(float s) { scale(s, s, s); }
    void rotateX(float angleRadians);
    void rotateY(float angleRadians);
    void rotateZ(float angleRadians);
    
    // Center mesh at origin
    void centerAtOrigin();
    
    // Fit mesh within a bounding box
    void fitToBounds(float minX, float maxX, float minY, float maxY, float minZ, float maxZ);
    
    // Get triangles for GPU
    const std::vector<Triangle>& getTriangles() const { return m_triangles; }
    std::vector<Triangle>& getTriangles() { return m_triangles; }
    
    // Mesh info
    size_t getTriangleCount() const { return m_triangles.size(); }
    size_t getVertexCount() const { return m_vertices.size(); }
    
    // Bounding box
    void getBoundingBox(float& minX, float& maxX, 
                        float& minY, float& maxY, 
                        float& minZ, float& maxZ) const;
    
    // Create primitive shapes for testing
    static Mesh createSphere(float radius, int subdivisions = 3);
    static Mesh createBox(float width, float height, float depth);
    static Mesh createCylinder(float radius, float height, int segments = 32);
    
private:
    std::vector<std::array<float, 3>> m_vertices;
    std::vector<std::array<float, 3>> m_normals;
    std::vector<Triangle> m_triangles;
    
    void computeTriangleNormals();
    void recomputeTriangles();
};

#endif // DSMC_MESH_H
