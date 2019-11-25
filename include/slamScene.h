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

#include <iomanip>
using namespace std;

#include <mrpt/opengl.h>
#include <mrpt/opengl/CPointCloudColoured.h>
#include <mrpt/gui.h>
#include <mrpt/config/CConfigFile.h>
#include <mrpt/config/CConfigFileBase.h>
using namespace mrpt;
using namespace mrpt::gui;
using namespace mrpt::poses;
using namespace mrpt::config;
using namespace mrpt::math;
using namespace mrpt::opengl;

#include <opencv/cv.h>
using namespace cv;

#include <eigen3/Eigen/Core>
using namespace Eigen;

#include <mapHandler.h>

namespace PLSLAM{

class slamScene{

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    slamScene();
    slamScene(string configFile);
    ~slamScene();

    void initializeScene(Matrix4d x_0);
    void initViewports(int W, int H);
    bool updateScene();
    bool updateSceneVO( Matrix4d T_last_kf );
    bool updateScene(const MapHandler* map);
    bool updateSceneSafe(const MapHandler* map);
    void updateSceneGraphs(const MapHandler* map);

    void setText(int frame_, float time_, int nPoints_, int nPointsH_, int nLines_, int nLinesH_);
    void setText(int frame_, float time_, int nPoints_, int nLines_ );
    void setPose(Matrix4d x_);
    void setGT(Matrix4d xgt_);
    void setComparison(Matrix4d xcomp_);
    void setImage(const Mat &image_);
    void setImage(const string &image_);
    void setLegend();
    void setHelp();
    void setPoints(CMatrixFloat pData_);
    void setLines(CMatrixFloat lData_);
    void setStereoCalibration(Matrix3d K_, float b_);
    void setKF(Matrix4d Tfw);

    bool waitUntilClose();
    bool isOpen();
    bool getYPR(float &yaw, float &pitch, float &roll);
    bool getPose(Matrix4d &T);

//private:

    CMatrixDouble getPoseFormat(Matrix4d T);
    CMatrixDouble33 getCovFormat(MatrixXd cov_);
    CPose3D getPoseXYZ(VectorXd x);

    CDisplayWindow3D*           win;
    COpenGLScene::Ptr             theScene;
    COpenGLViewport::Ptr          image, legend, help;
    opengl::CSetOfObjects::Ptr    bbObj, bbObj1, srefObj, srefObj1, gtObj, srefObjGT, elliObjL, elliObjP;
    opengl::CEllipsoid::Ptr       elliObj;
    opengl::CFrustum::Ptr         frustObj, frustObj1;
    opengl::CAxis::Ptr            axesObj;

    opengl::CSetOfLines::Ptr      lineObj, lineObj_local, kfsLinesObj, voLinesObj;
    opengl::CPointCloud::Ptr      pointObj, pointObj_local;
    opengl::CSetOfObjects::Ptr    kfsObj;

    float           sbb, saxis, srad, sref, sline, sfreq, szoom, selli, selev, sazim, sfrust, slinef;
    CVectorFixedDouble<6> v_aux, v_aux_, v_aux1, v_aux1_, v_auxgt, gt_aux_, v_auxgt_;
    CPose3D         pose, pose_0, pose_gt, pose_ini, ellPose, pose1,  change, frustumL_, frustumR_;
    Matrix4d        x_ini;
    mrptKeyModifier kmods;
    int             key;
    CMatrixDouble33 cov3D;
    bool            hasText, hasCov, hasGT, hasImg, hasLines, hasPoints, hasFrustum, hasComparison, hasLegend, hasHelp, hasAxes, hasTraj, isKitti;

    Matrix4d        x, xgt, xcomp;
    MatrixXd        cov, W;
    unsigned int    frame, nPoints, nPointsH, nLines, nLinesH;
    float           time;
    string          img, img_legend, img_help;
    CMatrixFloat    lData, pData;
    mrpt::img::CImage img_mrpt_legend, img_mrpt_image, img_mrpt_help;

    float           b, sigmaP, sigmaL, f, cx, cy, bsigmaL, bsigmaP;

};

}
