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
#include <list>
#include <map>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/CholmodSupport>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/Jacobi>
#include <Eigen/src/Core/MatrixBase.h>

#include <g2o/types/slam3d/vertex_se3.h>
#include <g2o/types/slam3d/edge_se3.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/solver.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/solvers/structure_only/structure_only_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/types/sba/types_six_dof_expmap.h>

#include "../3rdparty/DBoW2/DBoW2/TemplatedVocabulary.h"
#include "../3rdparty/DBoW2/DBoW2/FORB.h"
#include "../3rdparty/DBoW2/DBoW2/BowVector.h"
#include "../3rdparty/DBoW2/DBoW2/FClass.h"
#include "../3rdparty/DBoW2/DBoW2/FeatureVector.h"
#include "../3rdparty/DBoW2/DBoW2/ScoringObject.h"

#ifdef HAS_MRPT
#include <mrpt/utils/CTicTac.h>
#endif
#include <config.h>
#include <stereoFrame.h>
#include <stereoFrameHandler.h>
#include <keyFrame.h>
#include <mapFeatures.h>

using namespace std;
using namespace Eigen;
using namespace DBoW2;

typedef Matrix<int,  6,1> Vector6i;
typedef Matrix<float,6,1> Vector6f;
typedef Matrix<float,6,6> Matrix6f;
typedef Matrix<float,7,1> Vector7f;
typedef DBoW2::TemplatedVocabulary<DBoW2::FORB::TDescriptor, DBoW2::FORB> Vocabulary;

namespace PLSLAM
{

class MapHandler
{

public:

    MapHandler(PinholeStereoCamera* cam_);
    ~MapHandler(){};

    void initialize(KeyFrame* kf0);
    void finishSLAM();
    void addKeyFrame(KeyFrame *curr_kf);

    void addKeyFrame_multiThread(KeyFrame *curr_kf);
    void loopClosureThread();
    void localMappingThread();

    void lookForCommonMatches(KeyFrame *kf0, KeyFrame *&kf1);
    void expandGraphs();
    void formLocalMap();
    void formLocalMap_old();
    void fuseLandmarksFromLocalMap();
    void removeBadMapLandmarks();
    void removeRedundantKFs();
    void loopClosure();
    bool lookForLoopCandidates(int kf_idx_curr, int &kf_idx_prev);
    void insertKFBowVectorP(KeyFrame *&kf);
    void insertKFBowVectorL(KeyFrame *&kf);
    void insertKFBowVectorPL(KeyFrame *&kf);
    bool isLoopClosure(KeyFrame* kf0, KeyFrame* kf1, Vector6d &pose_inc,
                       vector<Vector4i> &lc_pt_idx, vector<Vector4i> &lc_ls_idx,
                       vector<PointFeature*> &lc_points, vector<LineFeature*>  &lc_lines);
    bool computeRelativePoseGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                Vector6d &pose_inc );
    bool computeRelativePoseRobustGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                Vector6d &pose_inc );
    bool loopClosureOptimizationEssGraphG2O();
    bool loopClosureOptimizationCovGraphG2O();
    void loopClosureFuseLandmarks();

    void localBundleAdjustment();
    void levMarquardtOptimizationLBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  );

    void globalBundleAdjustment();
    void levMarquardtOptimizationGBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  );



    PinholeStereoCamera* cam;

    vector<KeyFrame*> map_keyframes;
    vector<MapPoint*> map_points;
    vector<MapLine*>  map_lines;

    map<int,vector<int>> map_points_kf_idx; // base KF list from which the LM is observed
    map<int,vector<int>> map_lines_kf_idx;

    vector< vector<unsigned int> > full_graph;

    vector< vector<float> > conf_matrix;
    Vocabulary              dbow_voc_p, dbow_voc_l;

    unsigned int max_pt_idx, max_ls_idx, max_kf_idx ;

    // experiment variables
    Vector7f time;
    mrpt::utils::CTicTac clock;

    // lba variables
    vector<int> lba_kfs;

    // lc variables
    enum LCStatus{
        LC_IDLE,
        LC_ACTIVE,
        LC_READY
    };
    LCStatus lc_status;
    vector< Vector3i > lc_idxs,  lc_idx_list;
    vector< Vector6d > lc_poses, lc_pose_list;
    vector< vector<Vector4i> > lc_pt_idxs;
    vector< vector<Vector4i> > lc_ls_idxs;
};

}
