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
#include <stereoFrame.h>
#include <stereoFeatures.h>

typedef Matrix<double,6,6> Matrix6d;
typedef Matrix<double,6,1> Vector6d;

class StereoFrame;

namespace StVO{

class StereoFrameHandler
{

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    StereoFrameHandler( PinholeStereoCamera* cam_ );
    ~StereoFrameHandler();

    void initialize( const Mat img_l_, const Mat img_r_, const int idx_);
    void updateFrame();
    void insertStereoPair(const Mat img_l_, const Mat img_r_, const int idx_);
    void f2fTracking();
    void optimizePose();
    void optimizePose(Matrix4d DT_ini);
    double lineSegmentOverlap( Vector2d spl_obs, Vector2d epl_obs, Vector2d spl_proj, Vector2d epl_proj  );

    // adaptative fast
    int orb_fast_th;

    // slam-specific functions
    bool needNewKF();
    void currFrameIsKF();

    list<PointFeature*> matched_pt;
    list<LineFeature*>  matched_ls;

    StereoFrame* prev_frame;
    StereoFrame* curr_frame;
    PinholeStereoCamera* cam;

    int  n_inliers, n_inliers_pt, n_inliers_ls;

    // slam-specific variables
    bool     prev_f_iskf;
    double   entropy_first_prevKF;
    Matrix4d T_prevKF;
    Matrix6d cov_prevKF_currF;

private:

    void removeOutliers( Matrix4d DT );
    void gaussNewtonOptimization(Matrix4d &DT, Matrix6d &DT_cov, double &err_, int max_iters);
    void optimizeFunctions(Matrix4d DT, Matrix6d &H, Vector6d &g, double &e);

};

}
