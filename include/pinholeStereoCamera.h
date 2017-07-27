/*****************************************************************************
**   PL-SLAM: stereo visual SLAM with points and line segment features  	**
******************************************************************************
**																			**
**	Copyright(c) 2017, Ruben Gomez-Ojeda, University of Malaga              **
**	Copyright(c) 2017, MAPIR group, University of Malaga					**
**																			**
**  This program is free software: you can redistribute it and/or modify	**
**  it under the terms of the GNU General Public License (version 3) as		**
**	published by the Free Software Foundation.								**
**																			**
**  This program is distributed in the hope that it will be useful, but		**
**	WITHOUT ANY WARRANTY; without even the implied warranty of				**
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the			**
**  GNU General Public License for more details.							**
**																			**
**  You should have received a copy of the GNU General Public License		**
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.	**
**																			**
*****************************************************************************/

#pragma once

using namespace std;

#include <opencv/cv.h>
#include <opencv2/calib3d.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <line_descriptor_custom.hpp>
using namespace cv;
using namespace line_descriptor;

#include <eigen3/Eigen/Core>
using namespace Eigen;

// Pinhole model for a Stereo Camera in an ideal configuration (horizontal)
class PinholeStereoCamera
{

private:
    const int           width, height;
    double              fx, fy;
    double              cx, cy;
    const double        b;
    Matrix3d            K;
    bool                dist;
    Matrix<double,5,1>  d;
    Mat                 Kl, Kr, Dl, Dr, Rl, Rr, Pl, Pr, R, t, Q;
    Mat                 undistmap1l, undistmap2l, undistmap1r, undistmap2r;

public:

    PinholeStereoCamera( int width_, int height_, double fx_, double fy_, double cx_, double cy_, double b_,
                         double d0 = 0.0, double d1 = 0.0, double d2 = 0.0, double d3 = 0.0, double d4 = 0.0);
    PinholeStereoCamera( int width_, int height_, double fx_, double fy_, double cx_, double cy_, double b_, Mat Rl, Mat Rr,
                         double d0 = 0.0, double d1 = 0.0, double d2 = 0.0, double d3 = 0.0, double d4 = 0.0);

    //PinholeStereoCamera(int width_, int height_, double b_, Mat Kl_, Mat Kr_, Mat Rl_, Mat Rr_, Mat Dl_, Mat Dr_, bool equi );
    PinholeStereoCamera(int width_, int height_, double b_, Mat Kl_, Mat Kr_, Mat R_, Mat t_, Mat Dl_, Mat Dr_, bool equi );

    ~PinholeStereoCamera();

    // Image rectification
    void rectifyImage( const Mat& img_src, Mat& img_rec);
    void rectifyImagesLR( const Mat& img_src_l, Mat& img_rec_l, const Mat& img_src_r, Mat& img_rec_r );

    // Proyection and Back-projection
    Vector3d backProjection_unit(const double &u, const double &v, const double &disp, double &depth);
    Vector3d backProjection(const double &u, const double &v, const double &disp);
    Vector2d projection(Vector3d P);
    Vector3d projectionNH(Vector3d P);
    Vector2d nonHomogeneous( Vector3d x);

    // Getters
    inline const int getWidth()             const { return width; };
    inline const int getHeight()            const { return height; };
    inline const Matrix3d&    getK()        const { return K; };
    inline const double       getB()        const { return b; };
    inline const Matrix<double,5,1> getD()  const { return d; };
    inline const double getFx()             const { return fx; };
    inline const double getFy()             const { return fy; };
    inline const double getCx()             const { return cx; };
    inline const double getCy()             const { return cy; };

};

