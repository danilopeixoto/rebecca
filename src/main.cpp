// Copyright (c) 2019, Danilo Peixoto. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <aurora/Math.h>
#include <aurora/Utility.h>
#include <aurora/Vector.h>
#include <aurora/Color.h>
#include <aurora/Matrix.h>
#include <aurora/TriangleMesh.h>

#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <vector>

AURORA_NAMESPACE_USING

struct Ray {
    Vector3 origin;
    Vector3 direction;
    
    Ray() {}
    Ray(const Vector3 & origin, const Vector3 & direction) {
        this->origin = origin;
        this->direction = direction;
    }
    
    Vector3 point(float distance) const {
        return origin + direction * distance;
    }
};

struct Intersection {
    bool hit;
    float distance;
    size_t index;
    
    Intersection() {
        hit = false;
        distance = AURORA_INFINITY;
        index = std::numeric_limits<size_t>::max();
    }
    Intersection(bool hit, float distance, size_t index) {
        this->hit = hit;
        this->distance = distance;
        this->index = index;
    }
};

struct ShaderGlobals {
    Vector3 point;
    Vector3 normal;
    Vector2 textureCoordinate;
    Vector2 uv;
    Vector3 tangentU;
    Vector3 tangentV;
    Vector3 viewDirection;
    Vector3 lightDirection;
    Vector3 lightPoint;
    Vector3 lightNormal;
    
    ShaderGlobals() {}
    ShaderGlobals(
            const Vector3 & point,
            const Vector3 & normal, const Vector2 & textureCoordinate,
            const Vector2 & uv, const Vector3 & tangentU, const Vector3 & tangentV,
            const Vector3 & viewDirection, const Vector3 & lightDirection,
            const Vector3 & lightPoint, const Vector3 & lightNormal) {
        this->point = point;
        this->normal = normal;
        this->textureCoordinate = textureCoordinate;
        this->uv = uv;
        this->tangentU = tangentU;
        this->tangentV = tangentV;
        this->viewDirection = viewDirection;
        this->lightDirection = lightDirection;
        this->lightPoint = lightPoint;
        this->lightNormal = lightNormal;
    }
};

enum BSDFType {
    Light = 0,
    Diffuse,
    Specular,
    None
};

struct BSDF {
    BSDFType type;
    Color3 color;
    
    BSDF() {}
    BSDF(BSDFType type, const Color3 & color) {
        this->type = type;
        this->color = color;
    }
};

struct Shape {
    BSDF * bsdf;
    
    Shape() {}
    Shape(BSDF * bsdf) {
        this->bsdf = bsdf;
    }
    
    virtual bool intersects(
            const Ray & ray, Intersection & intersection) const = 0;
    virtual void calculateShaderGlobals(
            const Ray & ray, const Intersection & intersection,
            ShaderGlobals & shaderGlobals) const = 0;
    virtual float surfaceArea() const = 0;
};

struct Sphere : Shape {
    Vector3 position;
    float radius;
    
    Sphere() : Shape() {}
    Sphere(const Vector3 & position, float radius, BSDF * bsdf) : Shape(bsdf) {
        this->position = position;
        this->radius = radius;
    }
    
    virtual bool intersects(
            const Ray & ray, Intersection & intersection) const {
        Vector3 l = position - ray.origin;
        float t = l.dot(ray.direction);
        
        if (t < 0.0f)
            return false;
            
        float d2 = l.dot(l) - t * t;
        float r2 = radius * radius;
        
        if (d2 > r2)
            return false;
        
        float dt = std::sqrt(r2 - d2);
        
        float t0 = t - dt;
        float t1 = t + dt;
        
        if (t0 > t1)
            std::swap(t0, t1);
        
        if (t0 < 0.0f) {
            t0 = t1;
            
            if (t0 < 0.0f)
                return false;
        }
        
        intersection.hit = true;
        intersection.distance = t0;
        
        return true;
    }
    virtual void calculateShaderGlobals(
            const Ray & ray, const Intersection & intersection,
            ShaderGlobals & shaderGlobals) const {
        shaderGlobals.point = ray.point(intersection.distance);
        shaderGlobals.normal = (shaderGlobals.point - position).normalize();
        
        float theta = std::atan2(shaderGlobals.normal.x, shaderGlobals.normal.z);
        float phi = std::acos(shaderGlobals.normal.y);
        
        shaderGlobals.uv.x = theta * AURORA_INV_PI * 0.5f;
        shaderGlobals.uv.y = phi * AURORA_INV_PI;
        
        shaderGlobals.textureCoordinate.x = shaderGlobals.uv.x;
        shaderGlobals.textureCoordinate.y = shaderGlobals.uv.y;
        
        shaderGlobals.tangentU.x = std::cos(theta);
        shaderGlobals.tangentU.y = 0.0f;
        shaderGlobals.tangentU.z = -std::sin(theta);
        
        shaderGlobals.tangentV.x = -std::sin(theta) * std::cos(phi);
        shaderGlobals.tangentV.y = std::sin(phi);
        shaderGlobals.tangentV.z = -std::cos(theta) * std::cos(phi);
        
        shaderGlobals.viewDirection = -ray.direction;
    }
    virtual float surfaceArea() const {
        return 4.0f * AURORA_PI * radius * radius;
    }
};

struct Vertex {
    Vector3 position;
    Vector3 normal;
    Vector2 textureCoordinate;
    
    Vertex() {}
    Vertex(
            const Vector3 & position,
            const Vector3 & normal,
            const Vector2 & textureCoordinate) {
        this->position = position;
        this->normal = normal;
        this->textureCoordinate = textureCoordinate;
    }
};

struct Triangle : Shape {
    Vertex vertices[3];
    float area;
    
    Triangle() : Shape() {}
    Triangle(
            const Vertex & vertex0,
            const Vertex & vertex1,
            const Vertex & vertex2,
            BSDF * bsdf) : Shape(bsdf) {
        this->vertices[0] = vertex0;
        this->vertices[1] = vertex1;
        this->vertices[2] = vertex2;
        
        updateSurfaceArea();
    }
    
    virtual bool intersects(
            const Ray & ray, Intersection & intersection) const {
        const Vector3 & v0 = vertices[0].position;
        const Vector3 & v1 = vertices[1].position;
        const Vector3 & v2 = vertices[2].position;
        
        Vector3 u = v1 - v0;
        Vector3 v = v2 - v0;
        
        Vector3 p = ray.direction.cross(v);
        float d = u.dot(p);
        
        if (std::abs(d) < AURORA_EPSILON)
            return false;
        
        Vector3 t = ray.origin - v0;
        float inverseD = 1.0f / d;
        
        float alpha = t.dot(p) * inverseD;
        
        if (alpha < 0.0f || alpha > 1.0f)
            return false;
            
        Vector3 q = t.cross(u);
        
        float beta = ray.direction.dot(q) * inverseD;
        
        if (beta < 0.0f || alpha + beta > 1.0f)
            return false;
        
        float t0 = v.dot(q) * inverseD;
        
        if (t0 < AURORA_EPSILON)
            return false;
        
        intersection.hit = true;
        intersection.distance = t0;
        
        return true;
    }
    virtual void calculateShaderGlobals(
            const Ray & ray, const Intersection & intersection,
            ShaderGlobals & shaderGlobals) const {
        const Vector3 & p0 = vertices[0].position;
        const Vector3 & p1 = vertices[1].position;
        const Vector3 & p2 = vertices[2].position;
        
        const Vector3 & n0 = vertices[0].normal;
        const Vector3 & n1 = vertices[1].normal;
        const Vector3 & n2 = vertices[2].normal;
        
        const Vector2 & t0 = vertices[0].textureCoordinate;
        const Vector2 & t1 = vertices[1].textureCoordinate;
        const Vector2 & t2 = vertices[2].textureCoordinate;
        
        Vector3 b = barycentric(shaderGlobals.point, p0, p1, p2);
        
        shaderGlobals.point = ray.point(intersection.distance);
        
        shaderGlobals.normal = (n0 * b.x + n1 * b.y + n2 * b.z).normalize();
        shaderGlobals.textureCoordinate = t0 * b.x + t1 * b.y + t2 * b.z;
        
        shaderGlobals.uv = Vector2(b.x, b.y);
        
        calculateTangents(
            shaderGlobals.normal,
            shaderGlobals.tangentU,
            shaderGlobals.tangentV);
        
        shaderGlobals.viewDirection = -ray.direction;
    }
    virtual float surfaceArea() const {
        return area;
    }
    
    Triangle & updateSurfaceArea() {
        const Vector3 & v0 = vertices[0].position;
        const Vector3 & v1 = vertices[1].position;
        const Vector3 & v2 = vertices[2].position;
        
        Vector3 normal = calculateVectorArea(v0, v1, v2);
        area = normal.length() * 0.5f;
        
        return *this;
    }
};

struct Scene {
    std::vector<Shape *> shapes;
    
    Scene() {}
    Scene(const std::vector<Shape *> & shapes) {
        this->shapes = shapes;
    }
    
    Scene & fromMesh(const TriangleMesh & triangleMesh, BSDF * bsdf = nullptr) {
        size_t triangleCount = triangleMesh.getTriangleCount();
        
        bool hasNormals = triangleMesh.hasNormals();
        bool hasTextureCoordinates = triangleMesh.hasTextureCoordinates();
        
        shapes.reserve(shapes.size() + triangleCount);
        
        for (size_t i = 0; i < triangleCount; i++) {
            size_t indices[3];
            
            triangleMesh.getVertexIndices(i, indices[0], indices[1], indices[2]);
            
            Vertex vertex0;
            Vertex vertex1;
            Vertex vertex2;
            
            vertex0.position = triangleMesh.getVertex(indices[0]);
            vertex1.position = triangleMesh.getVertex(indices[1]);
            vertex2.position = triangleMesh.getVertex(indices[2]);
            
            if (hasNormals) {
                triangleMesh.getNormalIndices(i, indices[0], indices[1], indices[2]);
                
                vertex0.normal = triangleMesh.getNormal(indices[0]);
                vertex1.normal = triangleMesh.getNormal(indices[1]);
                vertex2.normal = triangleMesh.getNormal(indices[2]);
            }
            else {
                vertex0.normal = calculateVectorArea(
                    vertex0.position, vertex1.position, vertex2.position).normalize();
                
                vertex1.normal = vertex0.normal;
                vertex2.normal = vertex0.normal;
            }
            
            if (hasTextureCoordinates) {
                triangleMesh.getTextureIndices(i, indices[0], indices[1], indices[2]);
                
                vertex0.textureCoordinate = triangleMesh.getTextureCoordinates(indices[0]);
                vertex1.textureCoordinate = triangleMesh.getTextureCoordinates(indices[1]);
                vertex2.textureCoordinate = triangleMesh.getTextureCoordinates(indices[2]);
            }
            else {
                vertex0.textureCoordinate = Vector2(0.0f, 0.0f);
                vertex1.textureCoordinate = Vector2(1.0f, 0.0f);
                vertex2.textureCoordinate = Vector2(0.0f, 1.0f);
            }
            
            shapes.push_back(new Triangle(vertex0, vertex1, vertex2, bsdf));
        }
        
        return *this;
    }
    
    bool intersects(const Ray & ray, Intersection & intersection) const {
        for (size_t i = 0; i < shapes.size(); i++) {
            Shape * shape = shapes[i];
            
            Intersection temp;
            shape->intersects(ray, temp);
            
            if (temp.hit && temp.distance < intersection.distance) {
                intersection.hit = temp.hit;
                intersection.distance = temp.distance;
                intersection.index = i;
            }
        }
        
        return intersection.hit;
    }
};

struct Film {
    float width;
    float height;
    
    Film() {}
    Film(float width, float height) {
        this->width = width;
        this->height = height;
    }
    
    float aspectRatio() const {
        return width / height;
    }
};

struct Camera {
    float fieldOfView;
    float focalLength;
    float aperture;
    Film * film;
    Matrix4 worldMatrix;
    
    Camera() {}
    Camera(
            float fieldOfView, float focalLength, float aperture,
            Film * film, const Matrix4 & worldMatrix) {
        this->fieldOfView = fieldOfView;
        this->focalLength = focalLength;
        this->aperture = aperture;
        this->film = film;
        this->worldMatrix = worldMatrix;
    }
    
    void lookAt(
            const Vector3 & position,
            const Vector3 & target, const Vector3 & up) {
        Vector3 w = (position - target).normalize();
        Vector3 u = up.cross(w).normalize();
        Vector3 v = w.cross(u);
        
        worldMatrix[0][0] = u.x;
        worldMatrix[0][1] = u.y;
        worldMatrix[0][2] = u.z;
        worldMatrix[0][3] = 0.0f;
        
        worldMatrix[1][0] = v.x;
        worldMatrix[1][1] = v.y;
        worldMatrix[1][2] = v.z;
        worldMatrix[1][3] = 0.0f;
        
        worldMatrix[2][0] = w.x;
        worldMatrix[2][1] = w.y;
        worldMatrix[2][2] = w.z;
        worldMatrix[2][3] = 0.0f;
        
        worldMatrix[3][0] = position.x;
        worldMatrix[3][1] = position.y;
        worldMatrix[3][2] = position.z;
        worldMatrix[3][3] = 1.0f;
    }
    void generateRay(
            float x, float y, const Vector2 & sample,
            Ray & ray) const {
        float scale = focalLength * std::tan(fieldOfView * 0.5f);
        
        Vector3 position;
        
        if (aperture > 0.0f) {
            float radius = 0.5f * focalLength / aperture;
            Vector2 diskSample = radius * concentricSampleDisk(sample);
            
            position.x = diskSample.x;
            position.y = diskSample.y;
            position.z = 0.0f;
        }
        
        position *= worldMatrix;
        
        Vector3 pixel;
        
        pixel.x = (2.0f * (x + sample.x + 0.5f) / film->width - 1.0f) * scale * film->aspectRatio();
        pixel.y = (1.0f - 2.0f * (y + sample.y + 0.5f) / film->height) * scale;
        pixel.z = -focalLength;
        
        pixel *= worldMatrix;
        
        Vector3 direction = (pixel - position).normalize();
        
        ray.origin = position;
        ray.direction = direction;
    }
};

int main(int argc, char ** argv) {
    BSDF diffuse(BSDFType::Diffuse, Color3(1.0f, 1.0f, 1.0f));
    
    Sphere sphere(Vector3(0.0f, 0.0f, 0.0f), 1.0f, &diffuse);
    
    Vertex vertex0(Vector3(-0.5f, -0.5f, 0.0f), Vector3(0.0f, 0.0f, 1.0f), Vector2(0.0f, 0.0f));
    Vertex vertex1(Vector3(0.5f, -0.5f, 0.0f), Vector3(0.0f, 0.0f, 1.0f), Vector2(1.0f, 0.0f));
    Vertex vertex2(Vector3(0.0f, 1.0f, 0.0f), Vector3(0.0f, 0.0f, 1.0f), Vector2(0.0f, 1.0f));
    
    Triangle triangle(vertex0, vertex1, vertex2, &diffuse);
    
    std::vector<Shape *> shapes;
    
    shapes.push_back(&sphere);
    shapes.push_back(&triangle);
    
    Scene scene(shapes);
    
    Film film(800, 600);
    
    Matrix4 worldMatrix(
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f);
    
    Camera camera(radians(20.0f), 1.0f, 0.0f, &film, worldMatrix);
    camera.lookAt(Vector3(0.0f, 0.0f, 35.0f), Vector3(0.0f, 0.0f, 0.0f), Vector3(0.0f, 1.0f, 0.0f));
    
    Ray ray;
    camera.generateRay(400.0f, 300.0f, Vector2(), ray);
    
    std::cout << "Ray origin: " << ray.origin << std::endl;
    std::cout << "Ray direction: " << ray.direction << std::endl;
    
    Intersection intersection;
    scene.intersects(ray, intersection);
    
    std::cout << "Hit: " << intersection.hit << std::endl;
    std::cout << "Distance: " << intersection.distance << std::endl;
    std::cout << "Index: " << intersection.index << std::endl;
    
    if (intersection.hit) {
        Shape * shape = scene.shapes[intersection.index];
        
        ShaderGlobals shaderGlobals;
        shape->calculateShaderGlobals(ray, intersection, shaderGlobals);
        
        std::cout << "Point: " << shaderGlobals.point << std::endl;
        std::cout << "Normal: " << shaderGlobals.normal << std::endl;
        std::cout << "Texture coordinate: " << shaderGlobals.textureCoordinate << std::endl;
        std::cout << "UV: " << shaderGlobals.uv << std::endl;
        std::cout << "Tangent U: " << shaderGlobals.tangentU << std::endl;
        std::cout << "Tangent V: " << shaderGlobals.tangentV << std::endl;
        std::cout << "View direction: " << shaderGlobals.viewDirection << std::endl;
        std::cout << "Light direction: " << shaderGlobals.lightDirection << std::endl;
        std::cout << "Light point: " << shaderGlobals.lightPoint << std::endl;
        std::cout << "Light normal: " << shaderGlobals.lightNormal << std::endl;
    }
    
    return 0;
}
