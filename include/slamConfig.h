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

#include <string>
#include <config.h>

class SlamConfig : public Config
{

public:

    SlamConfig();
    ~SlamConfig();

    static void loadFromFile( const std::string &config_file );

    static SlamConfig& getInstance();

    // SLAM parameters
    static double&  maxKFEpipP()        { return getInstance().max_kf_epip_p; }
    static double&  maxKFEpipL()        { return getInstance().max_kf_epip_l; }
    static int&     minLMEssGraph()     { return getInstance().min_lm_ess_graph; }
    static int&     minLMCovGraph()     { return getInstance().min_lm_cov_graph; }
    static int&     minKFLocalMap()     { return getInstance().min_kf_local_map; }
    static double&  maxLM3DErr()        { return getInstance().max_lm_3d_err; }
    static double&  maxLMDirErr()       { return getInstance().max_lm_dir_err; }
    static double&  lambdaLbaLM()       { return getInstance().lambda_lba_lm; }
    static double&  lambdaLbaK()        { return getInstance().lambda_lba_k; }
    static int&     maxItersLba()       { return getInstance().max_iters_lba; }
    static int&     minLMObs()          { return getInstance().min_lm_obs; }
    static double&  maxCommonFtsKF()    { return getInstance().max_common_fts_kf; }
    static double&  maxDirLineError()   { return getInstance().max_dir_line_error; }
    static double&  maxPointLineError() { return getInstance().max_point_line_error; }
    static double&  maxPointPointError(){ return getInstance().max_point_point_error; }
    static std::string&  dbowVocP()     { return getInstance().vocabulary_p; }
    static std::string&  dbowVocL()     { return getInstance().vocabulary_l; }
    static double&  lcRes()             { return getInstance().lc_res; }
    static double&  lcUnc()             { return getInstance().lc_unc; }
    static double&  lcInl()             { return getInstance().lc_inl; }
    static double&  lcTrs()             { return getInstance().lc_trs; }
    static double&  lcRot()             { return getInstance().lc_rot; }
    static double&  lcMat()             { return getInstance().lc_mat; }
    static int&     maxItersPGO()       { return getInstance().max_iters_pgo; }
    static int&     lcKFDist()          { return getInstance().lc_kf_dist; }
    static int&     lcKFMaxDist()       { return getInstance().lc_kf_max_dist; }
    static int&     lcNKFClosest()      { return getInstance().lc_nkf_closest; }
    static double&  lcInlierRatio()     { return getInstance().lc_inlier_ratio; }
    static int&     minPointMatches()   { return getInstance().min_pt_matches; }
    static int&     minLineMatches()    { return getInstance().min_ls_matches; }
    static double&  kfInlierRatio()     { return getInstance().kf_inlier_ratio; }
    static bool&    fastMatching()      { return getInstance().fast_matching; }
    static bool&    hasRefinement()     { return getInstance().has_refinement; }
    static bool&    multithreadSLAM()   { return getInstance().mutithread_slam; }

    // SLAM parameters
    int    max_kf_num_frames;
    double min_kf_t_dist;
    double min_kf_r_dist;
    int    min_kf_n_feats;
    double max_kf_epip_p;
    double max_kf_epip_l;
    int    min_lm_cov_graph;
    int    min_lm_ess_graph;
    int    min_kf_local_map;
    double max_lm_3d_err;
    double max_lm_dir_err;
    double max_dir_line_error;
    double max_point_line_error;
    double max_point_point_error;
    double lambda_lba_lm;
    double lambda_lba_k;
    int    max_iters_lba;
    int    min_lm_obs;
    double max_common_fts_kf;
    std::string vocabulary_p, vocabulary_l;
    double lc_res;
    double lc_unc;
    double lc_inl;
    double lc_trs;
    double lc_rot;
    double lc_mat;
    int    max_iters_pgo;
    int    lc_kf_dist;
    int    lc_kf_max_dist;
    int    lc_nkf_closest;
    double lc_inlier_ratio;
    int    min_pt_matches;
    int    min_ls_matches;
    double kf_inlier_ratio;
    bool   fast_matching;
    bool   has_refinement;
    bool   mutithread_slam;

};

