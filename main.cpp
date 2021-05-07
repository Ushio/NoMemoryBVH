#include "pr.hpp"
#include <iostream>
#include <memory>

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

int main() {
    using namespace pr;

    Config config;
    config.ScreenWidth = 1920;
    config.ScreenHeight = 1080;
    config.SwapInterval = 1;
    Initialize(config);

    Camera3D camera;
    camera.origin = { 4, 4, 4 };
    camera.lookat = { 0, 0, 0 };
    camera.zUp = false;

    pr::SetDataDir(pr::ExecutableDir());
    auto pStage = pxr::UsdStage::Open(pr::GetDataPath("data/teapots.usda"));

    double e = GetElapsedTime();

    while (pr::NextFrame() == false) {
        if (IsImGuiUsingMouse() == false) {
            UpdateCameraBlenderLike(&camera);
        }

        ClearBackground(0.1f, 0.1f, 0.1f, 1);

        BeginCamera(camera);

        PushGraphicState();

        DrawGrid(GridAxis::XZ, 1.0f, 10, { 128, 128, 128 });
        DrawXYZAxis(1.0f);

        for (pxr::UsdPrim p : pStage->Traverse())
        {
            if (p.IsA<pxr::UsdGeomMesh>())
            {
                pxr::UsdGeomMesh mesh(p);
                pxr::UsdAttribute pointAttrib = mesh.GetPointsAttr();
                pxr::VtArray<pxr::GfVec3f> points;
                pointAttrib.Get(&points);
                pxr::GfMatrix4f matrix = pxr::GfMatrix4f(mesh.ComputeLocalToWorldTransform(pxr::UsdTimeCode::Default()));
                glm::mat4 m;
                memcpy(glm::value_ptr(m), matrix.data(), sizeof(glm::mat4));
                pr::SetObjectTransform(m);
                for (int i = 0; i < points.size(); ++i)
                {
                    pxr::GfVec3f point = points[i];
                    DrawPoint({ point[0], point[1], point[2] }, { 255, 255, 255 }, 2);
                }
                pr::SetObjectIdentify();
            }
        }

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());

        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
