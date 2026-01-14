#include "../include/GLProgram.h"
#include "MyNurbs.h"
#include <string>
#include <random>
#include <iomanip>
#include <numeric>
#include <fstream>

#define WINDOW_WIDTH 1600
#define WINDOW_HEIGHT 1200

// declare static members for use in callback functions
int GLProgram::windowWidth = WINDOW_WIDTH;
int GLProgram::windowHeight = WINDOW_HEIGHT;
Camera GLProgram::camera;
bool GLProgram::mousePressed = false;
double GLProgram::prevMouseX, GLProgram::prevMouseY;
glm::mat4 GLProgram::modelMatrix = glm::mat4(1.0f);



// 辅助函数：生成随机点 (基于正弦波加噪声，保证曲线比较平滑好看)
std::vector<glm::dvec3> GenerateRandomPoints(int count, double scale, double noiseLevel) {
    std::vector<glm::dvec3> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> noise(-noiseLevel, noiseLevel);

    for (int i = 0; i < count; ++i) {
        double t = (double)i / (count - 1) * 10.0; // x范围 0 到 10
        double x = t * scale;
        // y = sin(x) + random
        double y = std::sin(t) * scale * 2.0 + noise(gen);
        // z 稍微波动一下
        double z = noise(gen) * 2.0;

        points.push_back(glm::dvec3(x, y, 5.0 + z));
    }
    return points;
}

// 辅助函数：平移曲线控制点 
void TranslateCurve(MyNurbs::RationalCurve* curve, glm::dvec3 offset) {
    if (!curve) return;
    for (auto& p : curve->control_points) {
        p += offset;
    }
}

// 测试函数：生成、评估、排版
void GenerateAndAnalyzeCurves(std::vector<MyNurbs::RationalCurve>& outCurves) {
    using namespace MyNurbs;

    // Part 1: 统计评估 (500+ Sets) & 导出 CSV
    int numTests = 500;
    int pointsCount = 18; // 题目要求至少15个点
    int integrationSteps = 100; // 优化计算速度

    // 打开 CSV 文件准备写入
    std::ofstream csvFile("energy_data.csv");
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open energy_data.csv for writing!" << std::endl;
        return;
    }

    // 写入 CSV 表头
    // 格式：ID, Uniform_Stretch, Uniform_Bending, Chord_Stretch, Chord_Bending, ...
    csvFile << "TestID,"
        << "Uniform_Stretch,Uniform_Bending,"
        << "Chord_Stretch,Chord_Bending,"
        << "Centripetal_Stretch,Centripetal_Bending,"
        << "Universal_Stretch,Universal_Bending\n";

    std::cout << "Running assessment on " << numTests << " random datasets..." << std::endl;
    std::cout << "Exporting data to 'energy_data.csv'..." << std::endl;
    std::cout << "Progress: 0/" << numTests << std::flush;

    // 累加器（用于最后算平均值）
    double totalStretch[4] = { 0 };
    double totalBending[4] = { 0 };

    RationalCurve::InterpolateMethod methods[4] = {
        RationalCurve::Uniform,
        RationalCurve::ChordLength,
        RationalCurve::Centripetal,
        RationalCurve::Universal
    };
    std::string names[4] = { "Uniform", "Chord Length", "Centripetal", "Universal" };

    for (int i = 0; i < numTests; i++) {
        // 进度条
        if (i % 10 == 0) {
            std::cout << "\rProgress: " << i << "/" << numTests << "   " << std::flush;
        }

        // 1. 生成一组随机点
        auto pts = GenerateRandomPoints(pointsCount, 1.0, 0.5);

        // 2. 写入当前测试ID
        csvFile << i << ",";

        // 3. 对该组点应用 4 种方法
        for (int m = 0; m < 4; m++) {
            RationalCurve* curve = nullptr;
            double sE = 0.0;
            double bE = 0.0;

            // 尝试插值
            if (RationalCurve::Interpolate(methods[m], 3, pts, curve)) {
                // 计算能量
                sE = curve->stretchEnergy(integrationSteps);
                bE = curve->bendingEnergy(integrationSteps);
                delete curve;
            }

            // 累加到总和（用于控制台平均值）
            totalStretch[m] += sE;
            totalBending[m] += bE;

            // 写入 CSV (当前方法的 Stretch 和 Bending)
            csvFile << sE << "," << bE;

            // 如果不是最后一个方法，加逗号
            if (m < 3) csvFile << ",";
        }
        // 换行，开始下一组数据
        csvFile << "\n";
    }

    csvFile.close(); // 关闭文件
    std::cout << "\rProgress: " << numTests << "/" << numTests << " (Done!)" << std::endl;
    std::cout << "Data saved to 'energy_data.csv'." << std::endl;

    // 打印控制台汇总报表
    std::cout << "\n=============================================================\n";
    std::cout << "  Performance Assessment (Average over " << numTests << " runs)\n";
    std::cout << "=============================================================\n";
    std::cout << std::left << std::setw(20) << "Method"
        << std::setw(20) << "Avg Stretch Energy"
        << std::setw(20) << "Avg Bending Energy" << std::endl;
    std::cout << "-------------------------------------------------------------\n";

    for (int m = 0; m < 4; m++) {
        std::cout << std::left << std::setw(20) << names[m]
            << std::setw(20) << (totalStretch[m] / numTests)
            << std::setw(20) << (totalBending[m] / numTests) << std::endl;
    }
    std::cout << "=============================================================\n\n";

    // Part 2: 生成可视化曲线 (1 Set, Spatially Separated)
    std::cout << "Generating visualization curves..." << std::endl;
    auto vizPoints = GenerateRandomPoints(20, 1.5, 0.2);

    for (int m = 0; m < 4; m++) {
        RationalCurve* curve = nullptr;
        RationalCurve::Interpolate(methods[m], 3, vizPoints, curve);
        if (curve) {
            double yOffset = m * 6.0;
            TranslateCurve(curve, glm::dvec3(0, yOffset, 0));
            outCurves.push_back(*curve);
            delete curve;
        }
    }
    std::cout << "Visualization curves generated." << std::endl;
}

int main() {
    GLProgram program;

    tinynurbs::RationalCurve<double> crv; // Planar curve using float32
    crv.control_points = { glm::vec3(0, 0, 5), // std::vector of 3D points
                          glm::vec3(0, 3, 4),
                          glm::vec3(0, 4, 3),
                          glm::vec3(0, 5, 6),
                          glm::vec3(0, 7, 5),
                          glm::vec3(0, 11, 3),
                          glm::vec3(0, 13, 4),
                          glm::vec3(0, 15, 5),

    };
    crv.knots = { 0, 0, 0, 0, 1.0 / 5, 2.0 / 5, 3.0 / 5, 4.0 / 5, 1, 1, 1, 1 }; // std::vector of floats

    crv.degree = 3;

    crv.weights = { 0.2, 0.3, 0.4 , 1, 0.5, 0.6, 0.7, 0.8 };

    //MyNurbs::RationalCurve* mycrv = nullptr;
    //MyNurbs::RationalCurve::Create(crv.degree, crv.control_points, crv.knots, crv.weights, mycrv);
    //for (int i = 0; i < 20; i++) {
    //    auto point = mycrv->Evaluate(i * 0.05);
    //    std::cout << "u=" << i * 0.05 << "  Point: " << point.x << ", " << point.y << ", " << point.z << ", " << std::endl;
    //}

    //// 求直到2阶导数 (点, 速度, 加速度)
    //int derivOrder = 2;
    //std::vector<glm::dvec3> derivs = mycrv->Derivatives(0.5, derivOrder);

    //std::cout << "Position: " << derivs[0].x << ", " << derivs[0].y << ", " << derivs[0].z << std::endl;
    //std::cout << "Velocity: " << derivs[1].x << ", " << derivs[1].y << ", " << derivs[1].z << std::endl;
    //std::cout << "Accel   : " << derivs[2].x << ", " << derivs[2].y << ", " << derivs[2].z << std::endl;


    tinynurbs::RationalCurve<double> crv1; // Planar curve using float32
    crv1.control_points = {glm::fvec3(0,0,5),glm::fvec3(0,3,5),glm::fvec3(0,4,5),glm::fvec3(0,5,5),glm::fvec3(0,7,5),glm::fvec3(0,11,5),glm::fvec3(0,13,5),glm::fvec3(0,15,5)};
    crv1.knots = { 0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1 };
    crv1.degree = 3;
    crv1.weights = { 1, 1, 1 , 1, 1, 1, 1, 1 };

    tinynurbs::RationalCurve<double> crv2; // Planar curve using float32
    crv2.control_points = { glm::fvec3(2,0,5),glm::fvec3(2,3,4),glm::fvec3(2,4,3),glm::fvec3(2,5,4),glm::fvec3(2,7,6),glm::fvec3(2,11,5),glm::fvec3(2,13,3),glm::fvec3(2,15,5) };
    crv2.knots = { 0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1 };
    crv2.degree = 3;
    crv2.weights = { 1, 1, 1 , 1, 1, 1, 1, 1 };

    tinynurbs::RationalCurve<double> crv3; // Planar curve using float32
    crv3.control_points = { glm::fvec3(3,0,5),glm::fvec3(3,3,3),glm::fvec3(3,4,1),glm::fvec3(3,5,4),glm::fvec3(3,7,7),glm::fvec3(3,11,9),glm::fvec3(3,13,7),glm::fvec3(3,15,5) };
    crv3.knots = { 0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1 };
    crv3.degree = 3;
    crv3.weights = { 1, 1, 1 , 1, 1, 1, 1, 1 };

    tinynurbs::RationalCurve<double> crv4; // Planar curve using float32
    crv4.control_points = { glm::fvec3(5,0,5),glm::fvec3(5,3,4),glm::fvec3(5,4,3),glm::fvec3(5,5,3),glm::fvec3(5,7,2),glm::fvec3(5,11,3),glm::fvec3(5,13,4),glm::fvec3(5,15,5) };
    crv4.knots = { 0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1 };
    crv4.degree = 3;
    crv4.weights = { 1, 1, 1 , 1, 1, 1, 1, 1 };

    tinynurbs::RationalCurve<double> crv5; // Planar curve using float32
    crv5.control_points = { glm::fvec3(8,0,5),glm::fvec3(8,3,5),glm::fvec3(8,4,5),glm::fvec3(8,5,5),glm::fvec3(8,7,5),glm::fvec3(8,11,5),glm::fvec3(8,13,5),glm::fvec3(8,15,5) };
    crv5.knots = { 0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1 };
    crv5.degree = 3;
    crv5.weights = { 1, 1, 1 , 1, 1, 1, 1, 1 };


   // vector<tinynurbs::RationalCurve<double>> curves;
   //// curves.push_back(crv);
   // curves.push_back(crv1);
   // curves.push_back(crv2);
   // curves.push_back(crv3);
   // curves.push_back(crv4);
   // curves.push_back(crv5);



    MyNurbs::RationalCurve* mycrv1 = nullptr;
    MyNurbs::RationalCurve::Create(crv1.degree, crv1.control_points, crv1.knots, crv1.weights, mycrv1);

    MyNurbs::RationalCurve* mycrv2 = nullptr;
    MyNurbs::RationalCurve::Create(crv2.degree, crv2.control_points, crv2.knots, crv2.weights, mycrv2);

    MyNurbs::RationalCurve* mycrv3 = nullptr;
    MyNurbs::RationalCurve::Create(crv3.degree, crv3.control_points, crv3.knots, crv3.weights, mycrv3);

    MyNurbs::RationalCurve* mycrv4 = nullptr;
    MyNurbs::RationalCurve::Create(crv4.degree, crv4.control_points, crv4.knots, crv4.weights, mycrv4);

    MyNurbs::RationalCurve* mycrv5 = nullptr;
    MyNurbs::RationalCurve::Create(crv5.degree, crv5.control_points, crv5.knots, crv5.weights, mycrv5);

    vector<glm::dvec3> testPoints = {
        glm::fvec3(0,0,5),
        glm::fvec3(0,3,4),
        glm::fvec3(0,4,3),
        glm::fvec3(0,4.5,6),
        glm::fvec3(0,5,2),
        glm::fvec3(0,13,6),
        glm::fvec3(0,14,1),
        glm::fvec3(0,14.5,5)
    };

    MyNurbs::RationalCurve* testInterpolateCurve[4] = {nullptr, nullptr, nullptr, nullptr};

    MyNurbs::RationalCurve::Interpolate(
        MyNurbs::RationalCurve::InterpolateMethod::Uniform,
        3,
        testPoints,
        testInterpolateCurve[0]);

    for (int i = 0; i < testPoints.size(); i++) {
        testPoints[i].x += 2;
    }
    MyNurbs::RationalCurve::Interpolate(
        MyNurbs::RationalCurve::InterpolateMethod::ChordLength,
        3,
        testPoints,
        testInterpolateCurve[1]);

    for (int i = 0; i < testPoints.size(); i++) {
        testPoints[i].x += 1;
    }
    MyNurbs::RationalCurve::Interpolate(
        MyNurbs::RationalCurve::InterpolateMethod::Centripetal,
        3,
        testPoints,
        testInterpolateCurve[2]);


    for (int i = 0; i < testPoints.size(); i++) {
        testPoints[i].x += 1;
    }
    MyNurbs::RationalCurve::Interpolate(
        MyNurbs::RationalCurve::InterpolateMethod::Universal,
        3,
        testPoints,
        testInterpolateCurve[3]);
        
    vector<MyNurbs::RationalCurve> MyCurves;
    // curves.push_back(crv);

    //MyCurves.push_back(*mycrv1);
    //MyCurves.push_back(*mycrv2);
    //MyCurves.push_back(*mycrv3);
    //MyCurves.push_back(*mycrv4);
    //MyCurves.push_back(*mycrv5);


    MyCurves.push_back(*testInterpolateCurve[0]);
    MyCurves.push_back(*testInterpolateCurve[1]);
    MyCurves.push_back(*testInterpolateCurve[2]);
    MyCurves.push_back(*testInterpolateCurve[3]);
    string name[4] = {"Uniform", "Chord Length", "Centripetal", "Universal"};
    double stretchEnergy[4];
    double bendingEnergy[4];

    for (int i = 0; i < 4; i++) {
        stretchEnergy[i] = testInterpolateCurve[i]->stretchEnergy();
        bendingEnergy[i] = testInterpolateCurve[i]->bendingEnergy();
        cout << name[i] << " stretch energy " << stretchEnergy[i] << "  bending energy " << bendingEnergy[i] << endl;
    }



    tinynurbs::RationalSurface<double> srf;

    srf.degree_u = 3;
    srf.degree_v = 3;
    srf.knots_u = { 0, 0, 0, 0, 1, 1, 1, 1 };
    srf.knots_v = { 0, 0, 0, 0, 1, 1, 1, 1 };

    // 2D array of control points using tinynurbs::array2<T> container
    // Example from geometrictools.com/Documentation/NURBSCircleSphere.pdf
    srf.control_points = { 4, 4,
                          {glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1),
                           glm::vec3(2, 0, 1), glm::vec3(2, 4, 1),  glm::vec3(-2, 4, 1),  glm::vec3(-2, 0, 1),
                           glm::vec3(2, 0, -1), glm::vec3(2, 4, -1), glm::vec3(-2, 4, -1), glm::vec3(-2, 0, -1),
                           glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1)
                          }
    };
    srf.weights = { 4, 4,
                   {1,       1.f / 3.f, 1.f / 3.f, 1,
                   1.f / 3.f, 1.f / 9.f, 1.f / 9.f, 1.f / 3.f,
                   1.f / 3.f, 1.f / 9.f, 1.f / 9.f, 1.f / 3.f,
                   1,       1.f / 3.f, 1.f / 3.f, 1
                   }
    };


    tinynurbs::RationalSurface<double> srf1;

    srf1.degree_u = 3;
    srf1.degree_v = 3;
    srf1.knots_u = { 0, 0, 0, 0, 1, 1, 1, 1 };
    srf1.knots_v = { 0, 0, 0, 0, 1, 1, 1, 1 };

    // 2D array of control points using tinynurbs::array2<T> container
    // Example from geometrictools.com/Documentation/NURBSCircleSphere.pdf
    srf1.control_points = { 4, 4,
                          {glm::vec3(-0.75, -0.75, -0.5), glm::vec3(-0.75, -0.25, -0.75), glm::vec3(-0.75, 0.25, 0.0), glm::vec3(-0.75, 0.75, -0.5),
                           glm::vec3(-0.25, -0.75, 0.0), glm::vec3(-0.25, -0.25, 0.5),  glm::vec3(-0.25, 0.25, -0.5),  glm::vec3(-0.25, 0.75, -1.0),
                           glm::vec3(0.25, -0.75, 0.0), glm::vec3(0.25, -0.25, 0.5), glm::vec3(0.25, 0.25, -0.5), glm::vec3(0.25, 0.75, -1.0),
                           glm::vec3(0.75, -0.75, -0.5), glm::vec3(0.75, -0.25, -0.75), glm::vec3(0.75, 0.25, 0.0), glm::vec3(0.75, 0.75, -0.5)
                          }
    };
    srf1.weights = { 4, 4,
                   {1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1
                   }
    };



    tinynurbs::RationalSurface<double> srf3;

    srf3.degree_u = 3;
    srf3.degree_v = 3;
    srf3.knots_u = { 0, 0, 0, 0, 0.2, 0.4 ,0.6 ,0.8 ,1 ,1 ,1 ,1 };
    srf3.knots_v = { 0, 0, 0, 0, 0.5, 1, 1, 1, 1 };

    srf3.control_points =
    { 8, 5,
        {
        glm::fvec3(0,0,5),glm::fvec3(2,0,5),glm::fvec3(3,0,5),glm::fvec3(5,0,5),glm::fvec3(8,0,5),
        glm::fvec3(0,3,5),glm::fvec3(2,3,4),glm::fvec3(3,3,3),glm::fvec3(5,3,4),glm::fvec3(8,3,5),
        glm::fvec3(0,4,5),glm::fvec3(2,4,3),glm::fvec3(3,4,1),glm::fvec3(5,4,3),glm::fvec3(8,4,5),
        glm::fvec3(0,5,5),glm::fvec3(2,5,4),glm::fvec3(3,5,4),glm::fvec3(5,5,3),glm::fvec3(8,5,5),
        glm::fvec3(0,7,5),glm::fvec3(2,7,6),glm::fvec3(3,7,7),glm::fvec3(5,7,2),glm::fvec3(8,7,5),
        glm::fvec3(0,11,5),glm::fvec3(2,11,5),glm::fvec3(3,11,9),glm::fvec3(5,11,3),glm::fvec3(8,11,5),
        glm::fvec3(0,13,5),glm::fvec3(2,13,3),glm::fvec3(3,13,7),glm::fvec3(5,13,4),glm::fvec3(8,13,5),
        glm::fvec3(0,15,5),glm::fvec3(2,15,5),glm::fvec3(3,15,5),glm::fvec3(5,15,5),glm::fvec3(8,15,5)
        }
    };

    srf3.weights =
    { 8,5,
        {
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1
        }
    };



    MyNurbs::RationalSurface* mysrf3 = nullptr;
    MyNurbs::RationalSurface::Create(srf3.degree_u, srf3.degree_v,
        {
        glm::fvec3(0,0,5),glm::fvec3(2,0,5),glm::fvec3(3,0,5),glm::fvec3(5,0,5),glm::fvec3(8,0,5),
        glm::fvec3(0,3,5),glm::fvec3(2,3,4),glm::fvec3(3,3,3),glm::fvec3(5,3,4),glm::fvec3(8,3,5),
        glm::fvec3(0,4,5),glm::fvec3(2,4,3),glm::fvec3(3,4,1),glm::fvec3(5,4,3),glm::fvec3(8,4,5),
        glm::fvec3(0,5,5),glm::fvec3(2,5,4),glm::fvec3(3,5,4),glm::fvec3(5,5,3),glm::fvec3(8,5,5),
        glm::fvec3(0,7,5),glm::fvec3(2,7,6),glm::fvec3(3,7,7),glm::fvec3(5,7,2),glm::fvec3(8,7,5),
        glm::fvec3(0,11,5),glm::fvec3(2,11,5),glm::fvec3(3,11,9),glm::fvec3(5,11,3),glm::fvec3(8,11,5),
        glm::fvec3(0,13,5),glm::fvec3(2,13,3),glm::fvec3(3,13,7),glm::fvec3(5,13,4),glm::fvec3(8,13,5),
        glm::fvec3(0,15,5),glm::fvec3(2,15,5),glm::fvec3(3,15,5),glm::fvec3(5,15,5),glm::fvec3(8,15,5)
        },
        srf3.knots_u,
        srf3.knots_v,
        {
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1
        },
        mysrf3);
    



    if (!tinynurbs::surfaceIsValid(srf3)) {

        std::cout << "test";
        // check if degree, knots and control points are configured as per
        // #knots == #control points + degree + 1
    }

    vector<tinynurbs::RationalSurface<double>> surfaces;
    surfaces.push_back(srf3);


    vector<MyNurbs::RationalSurface> MySurfaces;
    MySurfaces.push_back(*mysrf3);



    std::vector<MyNurbs::RationalCurve> testCurves;
    GenerateAndAnalyzeCurves(testCurves);


    //program.init(curves, surfaces);
    program.init(testCurves, MySurfaces);
    program.setClearColor(0.05f, 0.18f, 0.25f, 1.0f);
    program.run(testCurves, MySurfaces);
    program.cleanup();
    return 0;
}
