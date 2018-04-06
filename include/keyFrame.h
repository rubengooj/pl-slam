/*****************************************************************************
**      Stereo VO and SLAM by combining point and line segment features     **
******************************************************************************
**                                                                          **
**  Copyright(c) 2016-2018, Ruben Gomez-Ojeda, University of Malaga         **
**  Copyright(c) 2016-2018, David Zuñiga-Noël, University of Malaga         **
**  Copyright(c) 2016-2018, MAPIR group, University of Malaga               **
**                                                                          **
**  This program is free software: you can redistribute it and/or modify    **
**  it under the terms of the GNU General Public License (version 3) as     **
**  published by the Free Software Foundation.                              **
**                                                                          **
**  This program is distributed in the hope that it will be useful, but     **
**  WITHOUT ANY WARRANTY; without even the implied warranty of              **
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            **
**  GNU General Public License for more details.                            **
**                                                                          **
**  You should have received a copy of the GNU General Public License       **
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.   **
**                                                                          **
*****************************************************************************/

#pragma once
#include <vector>
#include <eigen3/Eigen/Core>
#include <opencv/cv.h>

#include <DBoW2/TemplatedVocabulary.h>
#include <DBoW2/FORB.h>
#include <DBoW2/BowVector.h>
#include <DBoW2/FClass.h>
#include <DBoW2/FeatureVector.h>
#include <DBoW2/ScoringObject.h>

#include <auxiliar.h>
#include <stereoFeatures.h>
#include <stereoFrame.h>

using namespace cv;
using namespace std;
using namespace Eigen;
using namespace StVO;

typedef Matrix<double,6,1> Vector6d;
typedef Matrix<double,6,6> Matrix6d;

namespace PLSLAM{

class KeyFrame
{

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    KeyFrame() { }
    KeyFrame( const StVO::StereoFrame* sf );
    KeyFrame( const StVO::StereoFrame* sf, int kf_idx_ );
    ~KeyFrame();

    Mat plotKeyFrame();

    bool     local;

    int       f_idx;
    string    img_name;

    int      kf_idx;
    Matrix4d T_kf_w;
    Vector6d x_kf_w;
    Matrix6d xcov_kf_w;

    DBoW2::BowVector descDBoW_P, descDBoW_L;

    StVO::StereoFrame* stereo_frame;



};

}
