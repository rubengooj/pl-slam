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
#include <mutex>
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

#include <DBoW2/TemplatedVocabulary.h>
#include <DBoW2/FORB.h>
#include <DBoW2/BowVector.h>
#include <DBoW2/FClass.h>
#include <DBoW2/FeatureVector.h>
#include <DBoW2/ScoringObject.h>

#include <slamConfig.h>
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

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    MapHandler(PinholeStereoCamera* cam_);
    ~MapHandler() { }

    void initialize(KeyFrame* kf0);
    void finishSLAM();
    void addKeyFrame(KeyFrame *curr_kf);

    void addKeyFrame_multiThread(KeyFrame *curr_kf, KeyFrame *prev_kf);
    void handlerThread();

    void startThreads();
    void killThreads();

    void loopClosureThread();
    void localMappingThread();

    int matchKF2KFPoints(KeyFrame *prev_kf, KeyFrame *curr_kf);
    int matchMap2KFPoints();

    int matchKF2KFLines(KeyFrame *prev_kf, KeyFrame *curr_kf);
    int matchMap2KFLines();

    void lookForCommonMatches(KeyFrame *kf0, KeyFrame *&kf1);

    void expandGraphs();
    void formLocalMap();
    void formLocalMap( KeyFrame * kf );
    void formLocalMap_old();
    void removeBadMapLandmarks();
    void removeRedundantKFs();
    void loopClosure();
    bool lookForLoopCandidates(int kf_idx_curr, int &kf_idx_prev);
    void insertKFBowVectorP(KeyFrame *kf);
    void insertKFBowVectorL(KeyFrame *kf);
    void insertKFBowVectorPL(KeyFrame *kf);
    bool isLoopClosure(const KeyFrame* kf0, const KeyFrame* kf1, Vector6d &pose_inc,
                       vector<Vector4i> &lc_pt_idx, vector<Vector4i> &lc_ls_idx,
                       vector<PointFeature*> &lc_points, vector<LineFeature*>  &lc_lines);
    bool computeRelativePoseGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                Vector6d &pose_inc ) const;
    bool computeRelativePoseRobustGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                Vector6d &pose_inc ) const;
    bool loopClosureOptimizationEssGraphG2O();
    bool loopClosureOptimizationCovGraphG2O();
    void loopClosureFuseLandmarks();

    int localBundleAdjustment();
    int levMarquardtOptimizationLBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  );

    void globalBundleAdjustment();
    void levMarquardtOptimizationGBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  );

    PinholeStereoCamera* cam;

    vector<KeyFrame*> map_keyframes;
    vector<MapPoint*> map_points;
    vector<MapLine*>  map_lines;

    list<PointFeature*> matched_pt;
    list<LineFeature*>  matched_ls;

    map<int,vector<int>> map_points_kf_idx; // base KF list from which the LM is observed
    map<int,vector<int>> map_lines_kf_idx;

    vector< vector<unsigned int> > full_graph;

    vector< vector<float> > conf_matrix;
    Vocabulary              dbow_voc_p, dbow_voc_l;

    unsigned int max_pt_idx, max_ls_idx, max_kf_idx ;

    KeyFrame *prev_kf, *curr_kf;
    Matrix4d Twf, DT;

    // experiment variables
    Vector7f time;

    // VO status
    mutex m_insert_kf;
    enum VOStatus{
        VO_PROCESSING,
        VO_INSERTING_KF
    };
    VOStatus vo_status;

    // status of the LBA thread
    vector<int> lba_kfs;
    enum LBAState{
        LBA_IDLE,
        LBA_ACTIVE,
        LBA_READY,
        LBA_TERMINATED
    };
    LBAState lba_thread_status;

    // Local Mapping
    std::mutex lba_mutex;
    std::condition_variable lba_start, lba_join;

    vector< Vector3i > lc_idxs,  lc_idx_list;
    vector< Vector6d > lc_poses, lc_pose_list;
    vector< vector<Vector4i> > lc_pt_idxs;
    vector< vector<Vector4i> > lc_ls_idxs;

    std::mutex lc_mutex;
    std::condition_variable lc_start, lc_join;

    enum LCState{
        LC_IDLE,
        LC_ACTIVE,
        LC_READY,
        LC_TERMINATED
    };
    LCState lc_state, lc_thread_status;



    // KF queue
    std::list<pair<KeyFrame*,KeyFrame*>> kf_queue;  // list of curr_kf_mt and prev_kf_mt
    std::mutex kf_queue_mutex;
    std::condition_variable new_kf;
    KeyFrame* curr_kf_mt;
    KeyFrame* prev_kf_mt;

    std::mutex cout_mutex;
    void print_msg(const std::string &msg);

private:

    bool threads_started;

};

}
