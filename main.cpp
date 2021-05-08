#include "pr.hpp"
#include <iostream>
#include <memory>
#include <algorithm>

#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/variantSets.h>
#include <pxr/usd/usd/editContext.h>
#include <pxr/usd/usd/modelApi.h>
#include <pxr/usd/sdf/types.h>
#include <pxr/usd/usdGeom/xform.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/sphere.h>
#include <pxr/usd/usdGeom/xformCommonAPI.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/base/gf/matrix4f.h>
//void partitionsize(int& L, int& R, int partitionSize)
//{
//    int nodeCount = partitionSize / 2;
//    int H = floor(log2(nodeCount));
//    int s = pow(2, H - 1) - 1;
//    int S = pow(2, H) - 1;
//    int O = std::max(0, (nodeCount - 1) - s - S);
//    R = 2 * (s + O);
//    L = 2 * (nodeCount - 1) - R;
//}
uint32_t log2i(uint32_t x)
{
    uint32_t v = 0;
    while (x >>= 1) 
        ++v;
    return v;
}

uint32_t partitionsL( uint32_t n )
{
    uint32_t h = log2i((n + 1) / 2);
    uint32_t baseleafcnt = 1 << h;
    uint32_t basenodecnt = baseleafcnt * 2 - 1;
    uint32_t L = (baseleafcnt - 1) /* base left node count */ + std::min((n - basenodecnt) /* remind */, baseleafcnt);
    return L;
}

// 0 based heap traversal
inline int heapParent(int i) { return (i - 1) / 2; }
inline int heapChildL(int i) { return 2 * i + 1; }
inline int heapChildR(int i) { return 2 * (i + 1); }


struct Poly
{
    int indices[4];
};
std::vector<glm::vec3> g_vertices;
std::vector<Poly> g_polygons;
std::vector<int> g_polygonIndices;

// is inplace building possible?
std::vector<int> g_polygonIndicesBvh;

// debug
std::set<int> g_stored;


inline float minss(float a, float b)
{
    return (a < b) ? a : b;
}
inline float maxss(float a, float b)
{
    return (b < a) ? a : b;
}

void buildBvh( int node, int begIndex, int endIndex, int depth, const std::vector<glm::vec3>& vs, const std::vector<Poly>& polygons, std::vector<int>& polygonIndices, std::vector<int>& polygonIndicesOut )
{
    int axis = depth % 3;
    int nElement = endIndex - begIndex;

    if( nElement == 0 )
    {
        return;
    }

    assert(node * 2 < polygonIndicesOut.size());
    assert(node * 2 + 1 < polygonIndicesOut.size());

    if( nElement == 2 )
    {
        polygonIndicesOut[node * 2] = polygonIndices[begIndex];
        polygonIndicesOut[node * 2 + 1] = polygonIndices[begIndex + 1];
        return;
    }

    assert( 2 < nElement );

    float minValue = +FLT_MAX;
    int minIndex = -1;
    for (int i = begIndex; i < endIndex; ++i)
    {
        const Poly& p = polygons[polygonIndices[i]];
        for (int j = 0; j < 4; ++j)
        {
            float value = vs[p.indices[j]][axis];
            if( value < minValue )
            {
                minValue = value;
                minIndex = i;
            }
        }
    }
    std::swap( polygonIndices[begIndex], polygonIndices[minIndex] );

    float maxValue = -FLT_MAX;
    int maxIndex = -1;
    for (int i = begIndex + 1; i < endIndex; ++i)
    {
        const Poly& p = polygons[polygonIndices[i]];
        for (int j = 0; j < 4; ++j)
        {
            float value = vs[p.indices[j]][axis];
            if( maxValue < value ) 
            {
                maxValue = value;
                maxIndex = i;
            }
        }
    }
    std::swap( polygonIndices[begIndex + 1], polygonIndices[maxIndex] );

    polygonIndicesOut[node * 2] = polygonIndices[begIndex];
    polygonIndicesOut[node * 2 + 1] = polygonIndices[begIndex + 1];

    // validation
    for (int i = begIndex + 2; i < endIndex; ++i)
    {
        const Poly& p = polygons[polygonIndices[i]];
        for (int j = 0; j < 4; ++j)
        {
            float value = vs[p.indices[j]][axis];
            assert(minValue <= value);
            assert(value <= maxValue);
        }
    }

    uint32_t L = partitionsL( nElement / 2 ) * 2;
    std::nth_element( polygonIndices.begin() + begIndex + 2, polygonIndices.begin() + begIndex + 2 + L, polygonIndices.begin() + endIndex, [&]( int a, int b ) {
        float aValue = 0.0f;
        float bValue = 0.0f;
        for( int j = 0; j < 4; ++j )
        {
            int ai = polygons[a].indices[j];
            int bi = polygons[b].indices[j];
            aValue += vs[ai][axis];
            bValue += vs[bi][axis];
        }
        return aValue < bValue;
    });

    buildBvh(heapChildL(node), begIndex + 2,     begIndex + 2 + L, depth + 1, vs, polygons, polygonIndices, polygonIndicesOut);
    buildBvh(heapChildR(node), begIndex + 2 + L, endIndex,         depth + 1, vs, polygons, polygonIndices, polygonIndicesOut);
}

struct Hit
{
    // a hit is fired when tmin < t < tmax
    float tmin;
    float tmax;

    // geometric normal in world space
    float ng[3];
};


inline bool intersect_ray_triangle(glm::vec3 O, glm::vec3 D, glm::vec3 tri_v0, glm::vec3 tri_v1, glm::vec3 tri_v2, float* tfar, glm::vec2* uv, glm::vec3* ng)
{
    glm::vec3 v0 = tri_v0 - O;
    glm::vec3 v1 = tri_v1 - O;
    glm::vec3 v2 = tri_v2 - O;

    glm::vec3 e0 = v0 - v2;
    glm::vec3 e1 = v1 - v0;
    glm::vec3 e2 = v2 - v1;

    float U = glm::dot(glm::cross(e0, v0 + v2), D);
    float V = glm::dot(glm::cross(e1, v1 + v0), D);
    float W = glm::dot(glm::cross(e2, v2 + v1), D);

    bool hit =
        (U <= 0.0f && V <= 0.0f && W <= 0.0f) ||
        (0.0f <= U && 0.0f <= V && 0.0f <= W);

    if (hit == false)
    {
        return false;
    }

    glm::vec3 Ng = glm::cross(e2, e0);
    float denom = glm::dot(Ng, D);

    if (fabs(denom) < FLT_EPSILON)
    {
        return false;
    }

    float t = glm::dot(v0, Ng) / denom;
    float UVW = U + V + W;
    *tfar = t;
    *uv = glm::vec2(U, V) / UVW;
    *ng = Ng;

    return true;
}

void intersect_poly( Hit* h, const glm::vec3& ro, const glm::vec3& rd, const std::vector<glm::vec3>& vs, const Poly& poly )
{
    uint32_t v0i = poly.indices[0];
    uint32_t v1i = poly.indices[1];
    uint32_t v2i = poly.indices[2];
    uint32_t v3i = poly.indices[3];
    glm::vec3 v0 = vs[v0i];
    glm::vec3 v1 = vs[v1i];
    glm::vec3 v2 = vs[v2i];
    glm::vec3 v3 = vs[v3i];

    float t;
    glm::vec2 uv;
    glm::vec3 Ng;

    /*
    v0, v1, v3,
    v2, v3, v1
    */
    if (intersect_ray_triangle(ro, rd, v0, v1, v3, &t, &uv, &Ng))
    {
        if (h->tmin < t && t < h->tmax)
        {
            h->tmax = t;

            h->ng[0] = Ng.x;
            h->ng[1] = Ng.y;
            h->ng[2] = Ng.z;
        }
    }
    if (v2i != v3i && intersect_ray_triangle(ro, rd, v2, v3, v1, &t, &uv, &Ng))
    {
        if (h->tmin < t && t < h->tmax)
        {
            h->tmax = t;

            h->ng[0] = Ng.x;
            h->ng[1] = Ng.y;
            h->ng[2] = Ng.z;
        }
    }
}

void traverseBvh( Hit *hit, const glm::vec3& ro, const glm::vec3& rd, const glm::vec3& one_over_rd, int node, int depth, float tmin, float tmax, const std::vector<glm::vec3>& vs, const std::vector<Poly>& polygons, std::vector<int>& polygonIndices )
{
    int polygonSrc0 = node * 2;
    int polygonSrc1 = node * 2 + 1;
    if( polygonIndices.size() <= polygonSrc1)
    {
        return;
    }

    int polygon0 = polygonIndices[polygonSrc0];
    int polygon1 = polygonIndices[polygonSrc1];

    // intersect slab
    int axis = depth % 3;
    float minS = +FLT_MAX;
    float maxS = -FLT_MAX;
    for( int i = 0; i < 4; ++i )
    {
        int ai = polygons[polygon0].indices[i];
        int bi = polygons[polygon1].indices[i];
        minS = minss(minS, vs[ai][axis]);
        minS = minss(minS, vs[bi][axis]);
        maxS = maxss(maxS, vs[ai][axis]);
        maxS = maxss(maxS, vs[bi][axis]);
    }
    float t0 = (minS - ro[axis]) * one_over_rd[axis];
    float t1 = (maxS - ro[axis]) * one_over_rd[axis];

    // slab (region_min, region_max)
    float region_min = minss(t0, t1);
    float region_max = maxss(t0, t1);
    
    // update
    tmin = maxss(tmin, region_min);
    tmax = minss(tmax, region_max);

    if( tmax < tmin || hit->tmax < tmin /* slab is backside of hit */ )
    {
        return;
    }

    intersect_poly( hit, ro, rd, vs, polygons[polygon0] );
    intersect_poly( hit, ro, rd, vs, polygons[polygon1] );

    int lChild = heapChildL(node);
    int rChild = heapChildR(node);
    if( ro[axis] < 0 )
    {
        std::swap( lChild, rChild );
    }

    traverseBvh(hit, ro, rd, one_over_rd, lChild, depth + 1, tmin, tmax, vs, polygons, polygonIndices);
    traverseBvh(hit, ro, rd, one_over_rd, rChild, depth + 1, tmin, tmax, vs, polygons, polygonIndices);
}

int main() {
    using namespace pr;
    //for( int i = 0; i < 40; ++i )
    //{
    //    uint32_t n = i;
    //    //uint32_t h = log2i(n);
    //    //uint32_t base = 1 << h;
    //    //uint32_t L = std::min( base / 2 + (n - base) /* remind */, base );

    //    uint32_t h = log2i((n + 1) / 2);
    //    uint32_t baseleafcnt = 1 << h;
    //    uint32_t basenodecnt = baseleafcnt * 2 - 1;
    //    uint32_t L = ( baseleafcnt - 1 ) /* base left node count */ + std::min( (n - basenodecnt) /* remind */, baseleafcnt );
    //    printf("%d, baseleafcnt=%d, %d-%d\n", n, baseleafcnt, L, n - L - 1);
    //    // printf("%d, base %d %d - %d\n", n, base, L, n - L);
    //    printf("");
    //}

    Config config;
    config.ScreenWidth = 960;
    config.ScreenHeight = 960;
    config.SwapInterval = 1;
    Initialize(config);

    Camera3D camera;
    camera.origin = { 0, 0, 4 };
    camera.lookat = { 0, 0, 0 };
    camera.zUp = false;

    pr::SetDataDir(pr::ExecutableDir());
    auto pStage = pxr::UsdStage::Open(pr::GetDataPath("data/teapot.usda"));

    for (pxr::UsdPrim p : pStage->Traverse())
    {
        if (p.IsA<pxr::UsdGeomMesh>())
        {
            pxr::UsdGeomMesh mesh(p);
            pxr::VtArray<pxr::GfVec3f> points;
            pxr::VtArray<int> vertexCounts;
            pxr::VtArray<int> indices;
            mesh.GetPointsAttr().Get(&points);
            mesh.GetFaceVertexCountsAttr().Get(&vertexCounts);
            mesh.GetFaceVertexIndicesAttr().Get(&indices);

            g_vertices.resize(points.size());

            for (int i = 0; i < points.size(); ++i)
            {
                pxr::GfVec3f point = points[i];
                g_vertices[i] = glm::vec3( point[0], point[1], point[2] );
            }

            g_polygons.resize(vertexCounts.size());
            g_polygonIndices.resize(vertexCounts.size());
            g_polygonIndicesBvh.resize(vertexCounts.size());

            int indexStart = 0;
            for (int i = 0; i < vertexCounts.size(); ++i)
            {
                int vn = vertexCounts[i];
                int nLoopV = std::min(vn, 4);

                Poly p;
                for (int i = 0; i < nLoopV; ++i)
                {
                    p.indices[i] = indices[indexStart + i];
                }
                if( vn == 3 )
                {
                    p.indices[3] = p.indices[2];
                }
                g_polygons[i] = p;
                g_polygonIndices[i] = i;
                indexStart += vn;
            }
        }
    }

    assert(g_polygons.size() % 2 == 0);
    

    // TODO odd polygons
    buildBvh( 0, 0, g_polygonIndices.size(), 0, g_vertices, g_polygons, g_polygonIndices, g_polygonIndicesBvh );


    pr::ITexture* texture = CreateTexture();

    bool drawwire = false;

    SetDepthTest(true);

    double e = GetElapsedTime();

    while (pr::NextFrame() == false) {
        if (IsImGuiUsingMouse() == false) {
            UpdateCameraBlenderLike(&camera);
        }

        BeginCamera(camera);

        int _stride = 2;
        Image2DRGBA8 image;
        image.allocate(GetScreenWidth() / _stride, GetScreenHeight() / _stride);

        CameraRayGenerator rayGenerator(GetCurrentViewMatrix(), GetCurrentProjMatrix(), image.width(), image.height());

        pr::Stopwatch sw;

        // SerialFor ParallelFor
        ParallelFor(image.height(), [&](int j) {
            for (int i = 0; i < image.width(); ++i)
            {
                glm::vec3 ro, rd;
                rayGenerator.shoot(&ro, &rd, i, j, 0.5f, 0.5f);
         
                glm::vec3 one_over_rd = glm::vec3(1.0f) / rd;

                Hit h = {};
                h.tmin = 0.0f;
                h.tmax = FLT_MAX;
                traverseBvh( &h, ro, rd, one_over_rd, 0, 0, 0.0f, FLT_MAX, g_vertices, g_polygons, g_polygonIndicesBvh );

                glm::vec3 ng = glm::vec3(-1, -1, -1);
                if( h.tmax != FLT_MAX )
                {
                    ng = glm::normalize(glm::vec3(h.ng[0], h.ng[1], h.ng[2]));
                }
                glm::vec3 ncolor = (ng + glm::vec3(1.0f)) * 0.5f;
                image(i, j) = { 255 * ncolor.r, 255 * ncolor.g, 255 * ncolor.b, 255 };
            }
         });

        //image.saveAsPng("hair.png");
        //break;

        texture->upload(image);
        ClearBackground(texture);

        PushGraphicState();

        DrawGrid(GridAxis::XZ, 1.0f, 10, { 128, 128, 128 });
        DrawXYZAxis(1.0f);

        if(drawwire)
        for (pxr::UsdPrim p : pStage->Traverse())
        {
            if (p.IsA<pxr::UsdGeomMesh>())
            {
                pxr::UsdGeomMesh mesh(p);
                pxr::VtArray<pxr::GfVec3f> points;
                pxr::VtArray<int> vertexCounts;
                pxr::VtArray<int> indices;
                mesh.GetPointsAttr().Get(&points);
                mesh.GetFaceVertexCountsAttr().Get(&vertexCounts);
                mesh.GetFaceVertexIndicesAttr().Get(&indices);

                //pxr::GfMatrix4f matrix = pxr::GfMatrix4f(mesh.ComputeLocalToWorldTransform(pxr::UsdTimeCode::Default()));
                //glm::mat4 m;
                //memcpy(glm::value_ptr(m), matrix.data(), sizeof(glm::mat4));
                //pr::SetObjectTransform(m);

                for (int i = 0; i < points.size(); ++i)
                {
                    pxr::GfVec3f point = points[i];
                    DrawPoint({ point[0], point[1], point[2] }, { 255, 0, 0 }, 2);
                }

                int indexStart = 0;
                for (int i = 0; i < vertexCounts.size(); ++i)
                {
                    int vn = vertexCounts[i];
                    int nLoopV = std::min(vn, 4);
                    for (int i = 0; i < nLoopV; ++i)
                    {
                        pxr::GfVec3f p0 = points[indices[indexStart + i]];
                        pxr::GfVec3f p1 = points[indices[indexStart + (i + 1) % nLoopV]];
                        DrawLine({ p0[0], p0[1], p0[2] }, { p1[0], p1[1], p1[2] }, { 255, 255, 255 }, 1);
                    }
                    indexStart += vn;
                }

                // pr::SetObjectIdentify();
            }
        }

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());
        ImGui::Checkbox("drawwire", &drawwire);

        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
