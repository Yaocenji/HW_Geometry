#include "../include/GLProgram.h"
#include <string>

#define WINDOW_WIDTH 1600
#define WINDOW_HEIGHT 1200

// declare static members for use in callback functions
int GLProgram::windowWidth = WINDOW_WIDTH;
int GLProgram::windowHeight = WINDOW_HEIGHT;
Camera GLProgram::camera;
bool GLProgram::mousePressed = false;
double GLProgram::prevMouseX, GLProgram::prevMouseY;
glm::mat4 GLProgram::modelMatrix = glm::mat4(1.0f);

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


    vector<tinynurbs::RationalCurve<double>> curves;
   // curves.push_back(crv);
    curves.push_back(crv1);
    curves.push_back(crv2);
    curves.push_back(crv3);
    curves.push_back(crv4);
    curves.push_back(crv5);


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



    if (!tinynurbs::surfaceIsValid(srf3)) {

        std::cout << "test";
        // check if degree, knots and control points are configured as per
        // #knots == #control points + degree + 1
    }

    vector<tinynurbs::RationalSurface<double>> surfaces;
    surfaces.push_back(srf3);


    program.init(curves, surfaces);
    program.setClearColor(0.05f, 0.18f, 0.25f, 1.0f);
    program.run(curves, surfaces);
    program.cleanup();
    return 0;
}