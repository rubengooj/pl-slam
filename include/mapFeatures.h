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
#include <iostream>
#include <vector>
#include <list>
#include <opencv/cv.h>
#include <eigen3/Eigen/Core>

using namespace cv;
using namespace std;
using namespace Eigen;


typedef Matrix<int,5,1>    Vector5i;
typedef Matrix<double,6,1> Vector6d;
typedef Matrix<double,6,6> Matrix6d;

namespace PLSLAM{

class MapPoint
{

public:

    MapPoint(){};
    MapPoint(int idx_, Vector3d point3D_, Mat desc_, int kf_obs_, Vector2d obs_, Vector3d dir_, double sigma2_ = 1.f);
    ~MapPoint(){};

    void addMapPointObservation(Mat desc_, int kf_obs_, Vector2d obs_, Vector3d dir_,  double sigma2_ = 1.f);
    void updateAverageDescDir();

    int            idx;

    bool           inlier;
    bool           local;
    Vector3d       point3D;
    Vector3d       med_obs_dir;
    Mat            med_desc;

    vector<Mat>      desc_list;       // list with the descriptor of each observation
    vector<Vector2d> obs_list;        // list with the coordinates of each observation
    vector<Vector3d> dir_list;        // list with the direction unit vector of each observation
    vector<int>      kf_obs_list;     // list with KF index where the feature has been observed
    vector<double>   sigma_list;      // list with the sigma scale of each observation

};

class MapLine
{

public:

    MapLine(){};
    MapLine(int idx_, Vector6d line3D_, Mat desc_, int kf_obs_, Vector3d obs_, Vector3d dir_, Vector4d pts_,  double sigma2_ = 1.f);
    ~MapLine(){};

    void addMapLineObservation(Mat desc_, int kf_obs_, Vector3d obs_, Vector3d dir_, Vector4d pts_,  double sigma2_ = 1.f);
    void updateAverageDescDir();

    int            idx;

    bool           inlier;
    bool           local;
    Vector6d       line3D;            // 3D endpoints of the line segment
    Vector3d       med_obs_dir;
    Mat            med_desc;

    vector<Mat>      desc_list;       // list with the descriptor of each observation
    vector<Vector3d> obs_list;        // list with the coordinates of each observation ( 2D line equation, normalized by sqrt(lx2+ly2) )
    vector<Vector4d> pts_list;        // list with the coordinates of each endpoint (four coordinates)
    vector<Vector3d> dir_list;        // list with the direction unit vector of each observation (middle point)
    vector<int>      kf_obs_list;     // list with KF index where the feature has been observed
    vector<double>   sigma_list;      // list with the sigma scale of each observation

};

}
