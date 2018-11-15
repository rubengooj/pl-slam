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

#include "slamConfig.h"

//STL
#include <iostream>

//Boost
#include <boost/filesystem.hpp>

//YAML
#include <yaml-cpp/yaml.h>

using namespace std;

SlamConfig::SlamConfig()
    : Config()
{
    // Note: do not set default values for StVO Config here! (they can be overriten by the default class constructor...)

    // SLAM parameters
    // -----------------------------------------------------------------------------------------------------
    fast_matching         = false;      // allow for the fast matching (window-based) of the map features
    has_refinement        = false;      // refine the pose between keyframes (disabled as it is also performed by the LBA)
    mutithread_slam       = true;       // if true the system runs with both the VO, LBA and LC in parallel threads

    // lm numbers and errors
    min_lm_obs            = 5;          // min number of observations for a landmark to be considered as inlier
    max_common_fts_kf     = 0.9;        // max number of common features for a keyframe to be considered redundant (disabled)

    max_kf_epip_p         = 1.0;        // max epip distance for points in LBA
    max_kf_epip_l         = 1.0;        // max epip distance for line segments in LBA

    max_point_point_error = 0.1;        // maximum distance to fuse landmarks from local map (points)
    max_point_line_error  = 0.1;        // maximum distance to fuse landmarks from local map (lines)
    max_dir_line_error    = 0.1;        // maximum angle in line direction to fuse landmarks from local map
    kf_inlier_ratio       = 30.0;       // ratio of inliers to fuse local landmarks

    // graphs parameters
    min_lm_ess_graph      = 150;        // minimum number of landmarks for connectivity in Essential graph
    min_lm_cov_graph      = 75;         // minimum number of landmarks for connectivity in Covisibility graph
    min_kf_local_map      = 3;          // min number of landmarks for the local mapping

    // LBA
    lambda_lba_lm         = 0.00001;    // (if auto, this is the initial tau)
    lambda_lba_k          = 10.0;       // lambda_k for LM method in LBA
    max_iters_lba         = 15;         // maximum number of iterations

    // Loop closure
    vocabulary_p          = "/home/ruben/code/pl-slam-dev/vocabulary/mapir_orb.yml";
    vocabulary_l          = "/home/ruben/code/pl-slam-dev/vocabulary/mapir_lsd.yml";

    lc_res                = 1.0;        // maximum residue in relative pose estimation
    lc_unc                = 0.01;       // maximum uncertainty value in relative pose estimation
    lc_inl                = 0.3;        // minimum inlier ratio in relative pose estimation
    lc_trs                = 1.5;        // maximum translation in relative pose estimation
    lc_rot                = 35.0;       // maximum rotation in relative pose estimation

    max_iters_pgo         = 100;        // maximum number of iterations of the PGO
    lc_kf_dist            = 50;         // minimum number of KFs from prev LC
    lc_kf_max_dist        = 50;         // max distance from last LC KF
    lc_nkf_closest        = 4;          // number of KFs closest to the match to consider it as positive
    lc_inlier_ratio       = 30.0;       // inlier ratio to consider or not a loop closure

    min_pt_matches        = 10;         // min number of point observations
    min_ls_matches        = 6;          // min number of line segment observations
}

SlamConfig::~SlamConfig(){}

SlamConfig& SlamConfig::getInstance()
{
    static SlamConfig instance; // Instantiated on first use and guaranteed to be destroyed
    return instance;
}

template<typename T>
inline T loadSafe(const YAML::Node &config, std::string param, T default_value = T()) {

    if (YAML::Node parameter = config[param])
        return parameter.as<T>();
    else
        return default_value;
}

void SlamConfig::loadFromFile( const string &config_file )
{
    if (!boost::filesystem::exists(config_file) || !boost::filesystem::is_regular(config_file)) {
        cout << "[SlamConfig->loadFromFile] Invalid config file, keeping default params..." << endl;
        return;
    }

    // StVO-PL options
    // -----------------------------------------------------------------------------------------------------
    Config::loadFromFile(config_file);

    // PL-SLAM options
    // -----------------------------------------------------------------------------------------------------
    YAML::Node config = YAML::LoadFile(config_file);

    SlamConfig::minLMObs() = loadSafe(config, "min_lm_obs", SlamConfig::minLMObs());
    SlamConfig::maxCommonFtsKF() = loadSafe(config, "max_common_fts_kf", SlamConfig::maxCommonFtsKF());

    SlamConfig::maxKFEpipP() = loadSafe(config, "max_kf_epip_p", SlamConfig::maxKFEpipP());
    SlamConfig::maxKFEpipL() = loadSafe(config, "max_kf_epip_l", SlamConfig::maxKFEpipL());

    SlamConfig::maxLM3DErr() = loadSafe(config, "max_lm_3d_err", SlamConfig::maxLM3DErr());
    SlamConfig::maxLMDirErr() = loadSafe(config, "max_lm_dir_err", SlamConfig::maxLMDirErr());
    SlamConfig::maxPointPointError() = loadSafe(config, "max_point_point_error", SlamConfig::maxPointPointError());
    SlamConfig::maxPointLineError() = loadSafe(config, "max_point_line_error", SlamConfig::maxPointLineError());
    SlamConfig::maxDirLineError() = loadSafe(config, "max_dir_line_error", SlamConfig::maxDirLineError());

    SlamConfig::minLMEssGraph() = loadSafe(config, "min_lm_ess_graph", SlamConfig::minLMEssGraph());
    SlamConfig::minLMCovGraph() = loadSafe(config, "min_lm_cov_graph", SlamConfig::minLMCovGraph());
    SlamConfig::minKFLocalMap() = loadSafe(config, "min_kf_local_map", SlamConfig::minKFLocalMap());

    SlamConfig::lambdaLbaLM() = loadSafe(config, "lambda_lba_lm", SlamConfig::lambdaLbaLM());
    SlamConfig::lambdaLbaK() = loadSafe(config, "lambda_lba_k", SlamConfig::lambdaLbaK());
    SlamConfig::maxItersLba() = loadSafe(config, "max_iters_lba", SlamConfig::maxItersLba());

    SlamConfig::dbowVocP() = loadSafe(config, "vocabulary_p", SlamConfig::dbowVocP());
    SlamConfig::dbowVocL() = loadSafe(config, "vocabulary_l", SlamConfig::dbowVocL());

    SlamConfig::lcMat() = loadSafe(config, "lc_mat", SlamConfig::lcMat());
    SlamConfig::lcRes() = loadSafe(config, "lc_res", SlamConfig::lcRes());
    SlamConfig::lcUnc() = loadSafe(config, "lc_unc", SlamConfig::lcUnc());
    SlamConfig::lcInl() = loadSafe(config, "lc_inl", SlamConfig::lcInl());
    SlamConfig::lcTrs() = loadSafe(config, "lc_trs", SlamConfig::lcTrs());
    SlamConfig::lcRot() = loadSafe(config, "lc_rot", SlamConfig::lcRot());

    SlamConfig::maxItersPGO() = loadSafe(config, "max_iters_pgo", SlamConfig::maxItersPGO());
    SlamConfig::lcKFDist() = loadSafe(config, "lc_kf_dist", SlamConfig::lcKFDist());
    SlamConfig::lcKFMaxDist() = loadSafe(config, "lc_kf_max_dist", SlamConfig::lcKFMaxDist());
    SlamConfig::lcNKFClosest() = loadSafe(config, "lc_nkf_closest", SlamConfig::lcNKFClosest());
    SlamConfig::lcInlierRatio() = loadSafe(config, "lc_inlier_ratio", SlamConfig::lcInlierRatio());

    SlamConfig::minPointMatches() = loadSafe(config, "min_pt_matches", SlamConfig::minPointMatches());
    SlamConfig::minLineMatches() = loadSafe(config, "min_ls_matches", SlamConfig::minLineMatches());
    SlamConfig::kfInlierRatio() = loadSafe(config, "kf_inlier_ratio", SlamConfig::kfInlierRatio());

    SlamConfig::fastMatching() = loadSafe(config, "fast_matching", SlamConfig::fastMatching());
    SlamConfig::hasRefinement() = loadSafe(config, "has_refinement", SlamConfig::hasRefinement());
    SlamConfig::multithreadSLAM() = loadSafe(config, "mutithread_slam", SlamConfig::multithreadSLAM());
}
