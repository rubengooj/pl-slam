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

#include <mapHandler.h>

namespace PLSLAM
{

MapHandler::MapHandler(PinholeStereoCamera* cam_)
{
    cam = cam_;
    // load vocabulary
    if( Config::hasPoints() )
        dbow_voc_p.load( Config::dbowVocP() );
    if( Config::hasLines() )
        dbow_voc_l.load( Config::dbowVocL() );

}

void MapHandler::initialize( KeyFrame *kf0 )
{

    // reset information from the map
    map_keyframes.clear();
    map_points.clear();
    map_points_kf_idx.clear();
    map_lines.clear();
    map_lines_kf_idx.clear();
    full_graph.clear();
    conf_matrix.clear();
    lc_idx_list.clear();
    lc_pose_list.clear();
    max_pt_idx = 0;
    max_ls_idx = 0;
    max_kf_idx = 0;

    // initialize graphs
    full_graph.resize(1);
    full_graph[0].resize(1);
    full_graph[0][0] = 0;
    conf_matrix.resize(1);
    conf_matrix[0].resize(1);
    conf_matrix[0][0] = 1.0;

    // reset indices
    for( vector<PointFeature*>::iterator it = kf0->stereo_frame->stereo_pt.begin(); it != kf0->stereo_frame->stereo_pt.end(); it++)
        (*it)->idx = -1;
    for( vector<LineFeature*>::iterator  it = kf0->stereo_frame->stereo_ls.begin(); it != kf0->stereo_frame->stereo_ls.end(); it++)
        (*it)->idx = -1;

    // initialize DBoW descriptor vector and LC status
    vector<Mat> curr_desc;
    if( Config::hasPoints() )
    {
        curr_desc.reserve( kf0->stereo_frame->pdesc_l.rows );
        for ( int i = 0; i < kf0->stereo_frame->pdesc_l.rows; i++ )
            curr_desc.push_back( kf0->stereo_frame->pdesc_l.row(i) );
        dbow_voc_p.transform( curr_desc, kf0->descDBoW_P );
        curr_desc.clear();
    }
    if( Config::hasLines() )
    {
        curr_desc.reserve( kf0->stereo_frame->ldesc_l.rows );
        for ( int i = 0; i < kf0->stereo_frame->ldesc_l.rows; i++ )
            curr_desc.push_back( kf0->stereo_frame->ldesc_l.row(i) );
        dbow_voc_l.transform( curr_desc, kf0->descDBoW_L );
    }
    lc_status = LC_IDLE;

    // insert keyframe and add to map of indexes
    vector<int> aux_vec;
    map_keyframes.push_back( kf0 );
    map_points_kf_idx.insert( std::pair<int,vector<int>>(kf0->kf_idx,aux_vec) );
    map_lines_kf_idx.insert(  std::pair<int,vector<int>>(kf0->kf_idx,aux_vec) );

}

void MapHandler::finishSLAM()
{

    // close loop if there is an active LC
    if( lc_status == LC_ACTIVE || lc_status == LC_READY )
    {
        loopClosureOptimizationCovGraphG2O();
        // if succesfully closed
        lc_status = LC_IDLE;
        lc_idxs.clear();
        lc_poses.clear();
        lc_pt_idxs.clear();
        lc_ls_idxs.clear();
    }

}

void MapHandler::addKeyFrame( KeyFrame *curr_kf )
{


    // expand graphs
    expandGraphs();

    // select previous keyframe
    KeyFrame* prev_kf;
    prev_kf = map_keyframes.back();
    max_kf_idx++;
    curr_kf->kf_idx = max_kf_idx;
    curr_kf->local  = true;

    // update pose of current keyframe wrt previous one (in case of LC)
    Matrix4d T_curr_w = prev_kf->T_kf_w * curr_kf->T_kf_w;
    curr_kf->T_kf_w = T_curr_w;
    curr_kf->x_kf_w = logmap_se3(T_curr_w);

    // reset indices
    for( vector<PointFeature*>::iterator it = curr_kf->stereo_frame->stereo_pt.begin(); it != curr_kf->stereo_frame->stereo_pt.end(); it++)
        (*it)->idx = -1;
    for( vector<LineFeature*>::iterator  it = curr_kf->stereo_frame->stereo_ls.begin(); it != curr_kf->stereo_frame->stereo_ls.end(); it++)
        (*it)->idx = -1;

    // look for common matches and update the full graph
    lookForCommonMatches( prev_kf, curr_kf );
    if( Config::hasPoints() && Config::hasLines() )
        insertKFBowVectorPL(curr_kf);
    else if( Config::hasPoints() && !Config::hasLines() )
        insertKFBowVectorP(curr_kf);
    else if( !Config::hasPoints() && Config::hasLines() )
        insertKFBowVectorL(curr_kf);

    // insert keyframe and add to map of indexes
    clock.Tic();
    vector<int> aux_vec;
    map_keyframes.push_back( curr_kf );
    map_points_kf_idx.insert( std::pair<int,vector<int>>(curr_kf->kf_idx,aux_vec) );
    map_lines_kf_idx.insert(  std::pair<int,vector<int>>(curr_kf->kf_idx,aux_vec) );

    // form local map
    formLocalMap();

    // fuse similar 3D landmarks from local map
    //fuseLandmarksFromLocalMap();

    // perform local bundle adjustment
    localBundleAdjustment();

    // Recent map LMs culling (implement filters for line segments, which seems to be unaccurate)
    removeBadMapLandmarks();

    // Check for loop closures
    loopClosure();

}

void MapHandler::lookForCommonMatches( KeyFrame* kf0, KeyFrame* &kf1 )
{

    // grab frame number
    int kf0_idx = kf0->kf_idx;
    int kf1_idx = kf1->kf_idx;

    // estimates pose increment
    Matrix4d DT = inverse_se3( kf1->T_kf_w ) * kf0->T_kf_w;

    // Estimates Twf
    Matrix4d Twf = inverse_se3( kf1->T_kf_w );

    bool has_refinement = false;
    list<PointFeature*> matched_pt;
    list<LineFeature*>  matched_ls;

    // ---------------------------------------------------
    // find matches between prev_keyframe and curr_frame
    // ---------------------------------------------------
    // points f2f tracking
    int common_pt = 0;
    if( Config::hasPoints() && !(kf1->stereo_frame->stereo_pt.size()==0) && !(kf0->stereo_frame->stereo_pt.size()==0)  )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        Mat pdesc_l1, pdesc_l2;
        vector<vector<DMatch>> pmatches_12, pmatches_21;
        // 12 and 21 matches
        pdesc_l1 = kf0->stereo_frame->pdesc_l;
        pdesc_l2 = kf1->stereo_frame->pdesc_l;
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StVO::StereoFrame::matchPointFeatures,
                                      kf0->stereo_frame, bfm, pdesc_l1, pdesc_l2, ref(pmatches_12) );
                auto match_r = async( launch::async, &StVO::StereoFrame::matchPointFeatures,
                                      kf0->stereo_frame, bfm, pdesc_l2, pdesc_l1, ref(pmatches_21) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( pdesc_l1, pdesc_l2, pmatches_12, 2);
                bfm->knnMatch( pdesc_l2, pdesc_l1, pmatches_21, 2);
            }
        }
        else
            bfm->knnMatch( pdesc_l1, pdesc_l2, pmatches_12, 2);
        // sort matches by the distance between the best and second best matches
        double nn12_dist_th = Config::minRatio12P();
        // resort according to the queryIdx
        sort( pmatches_12.begin(), pmatches_12.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( pmatches_21.begin(), pmatches_21.end(), sort_descriptor_by_queryIdx() );
        // bucle around pmatches
        for( int i = 0; i < pmatches_12.size(); i++ )
        {
            // check if they are mutual best matches
            int lr_qdx = pmatches_12[i][0].queryIdx;
            int lr_tdx = pmatches_12[i][0].trainIdx;
            int rl_tdx;
            if( Config::bestLRMatches() )
                rl_tdx = pmatches_21[lr_tdx][0].trainIdx;
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_nn = pmatches_12[i][0].distance;
            double dist_12 = pmatches_12[i][0].distance / pmatches_12[i][1].distance;
            // check the f2f max disparity condition
            if( lr_qdx == rl_tdx  && dist_12 > nn12_dist_th )
            {
                Vector3d P_ = DT.block(0,0,3,3) * kf0->stereo_frame->stereo_pt[lr_qdx]->P + DT.col(3).head(3);
                Vector2d pl_proj = cam->projection( P_ );
                double   error = ( pl_proj - kf1->stereo_frame->stereo_pt[lr_tdx]->pl ).norm() * sqrt(kf0->stereo_frame->stereo_pt[lr_qdx]->sigma2);
                if( error < sqrt(7.815) )
                //if( error < Config::maxKFEpipP() )
                {
                    common_pt++;
                    // new 3D landmark
                    if( kf0->stereo_frame->stereo_pt[lr_qdx]->idx == -1 )
                    {
                        // assign indices
                        kf0->stereo_frame->stereo_pt[lr_qdx]->idx = max_pt_idx;
                        kf1->stereo_frame->stereo_pt[lr_tdx]->idx = max_pt_idx;
                        // create new 3D landmark with the observation from previous KF
                        Matrix4d Tfw = ( kf0->T_kf_w );
                        Vector3d P3d = Tfw.block(0,0,3,3) * kf0->stereo_frame->stereo_pt[lr_qdx]->P + Tfw.col(3).head(3);
                        //Vector3d dir = kf0->stereo_frame->stereo_pt[lr_qdx]->P / kf0->stereo_frame->stereo_pt[lr_qdx]->P.norm();
                        Vector3d dir = P3d / P3d.norm();
                        MapPoint* map_point = new MapPoint(max_pt_idx,P3d,kf0->stereo_frame->pdesc_l.row(lr_qdx),kf0->kf_idx,kf0->stereo_frame->stereo_pt[lr_qdx]->pl,dir);
                        // add new 3D landmark to kf_idx where it was first observed
                        map_points_kf_idx.at( kf0->kf_idx ).push_back( max_pt_idx );
                        // add observation of the 3D landmark from current KF
                        P3d = kf1->T_kf_w.block(0,0,3,3) * kf1->stereo_frame->stereo_pt[lr_tdx]->P + kf1->T_kf_w.col(3).head(3);
                        //dir = kf1->stereo_frame->stereo_pt[lr_tdx]->P / kf1->stereo_frame->stereo_pt[lr_tdx]->P.norm();
                        dir = P3d / P3d.norm();
                        map_point->addMapPointObservation(kf1->stereo_frame->pdesc_l.row(lr_tdx),kf1->kf_idx,kf1->stereo_frame->stereo_pt[lr_tdx]->pl,dir);
                        // add 3D landmark to map
                        map_points.push_back(map_point);
                        max_pt_idx++;
                        // update full graph (new feature)
                        full_graph[kf1_idx][kf0_idx]++;
                        full_graph[kf0_idx][kf1_idx]++;
                        // if has refine pose:
                        if( has_refinement )
                        {
                            PointFeature* pt = kf0->stereo_frame->stereo_pt[lr_qdx];
                            pt->pl_obs = kf1->stereo_frame->stereo_pt[lr_tdx]->pl;
                            pt->inlier = true;
                            matched_pt.push_back(pt);
                        }
                    }
                    // 3D landmark exists: copy idx && add observation to map landmark
                    else
                    {
                        // copy idx
                        int lm_idx = kf0->stereo_frame->stereo_pt[lr_qdx]->idx;
                        if( map_points[lm_idx] != NULL )
                        {
                            kf1->stereo_frame->stereo_pt[lr_tdx]->idx = lm_idx;
                            // add observation of the 3D landmark from current KF
                            Vector3d p3d = kf1->T_kf_w.block(0,0,3,3) * kf1->stereo_frame->stereo_pt[lr_tdx]->P + kf1->T_kf_w.col(3).head(3);
                            //Vector3d dir = kf1->stereo_frame->stereo_pt[lr_tdx]->P / kf1->stereo_frame->stereo_pt[lr_tdx]->P.norm();
                            Vector3d dir = p3d / p3d.norm();
                            map_points[lm_idx]->addMapPointObservation(kf1->stereo_frame->pdesc_l.row(lr_tdx),kf1->kf_idx,kf1->stereo_frame->stereo_pt[lr_tdx]->pl,dir);
                            // update full graph (previously observed feature)
                            for( vector<int>::iterator obs = map_points[lm_idx]->kf_obs_list.begin(); obs != map_points[lm_idx]->kf_obs_list.end(); obs++ )
                            {
                                if( (*obs) != kf1_idx )
                                {
                                    full_graph[kf1_idx][(*obs)]++;
                                    full_graph[(*obs)][kf1_idx]++;
                                }
                            }
                            // if has refine pose:
                            if( has_refinement )
                            {
                                PointFeature* pt = kf0->stereo_frame->stereo_pt[lr_qdx];
                                pt->pl_obs = kf1->stereo_frame->stereo_pt[lr_tdx]->pl;
                                pt->inlier = true;
                                matched_pt.push_back(pt);
                            }
                        }
                    }
                }
            }
        }
    }

    // line segments f2f tracking
    int common_ls = 0;
    if( Config::hasLines() && !(kf1->stereo_frame->stereo_ls.size()==0) && !(kf0->stereo_frame->stereo_ls.size()==0)  )
    {
        Ptr<BinaryDescriptorMatcher> bdm = BinaryDescriptorMatcher::createBinaryDescriptorMatcher();
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        Mat ldesc_l1, ldesc_l2;
        vector<vector<DMatch>> lmatches_12, lmatches_21;
        // 12 and 21 matches
        ldesc_l1 = kf0->stereo_frame->ldesc_l;
        ldesc_l2 = kf1->stereo_frame->ldesc_l;
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StVO::StereoFrame::matchLineFeatures, kf0->stereo_frame, bfm, ldesc_l1, ldesc_l2, ref(lmatches_12) );
                auto match_r = async( launch::async, &StVO::StereoFrame::matchLineFeatures, kf0->stereo_frame, bfm, ldesc_l2, ldesc_l1, ref(lmatches_21) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( ldesc_l1, ldesc_l2, lmatches_12, 2);
                bfm->knnMatch( ldesc_l2, ldesc_l1, lmatches_21, 2);
            }
        }
        else
            bfm->knnMatch( ldesc_l1, ldesc_l2, lmatches_12, 2);
        // sort matches by the distance between the best and second best matches
        double nn_dist_th, nn12_dist_th;
        kf0->stereo_frame->lineDescriptorMAD(lmatches_12,nn_dist_th, nn12_dist_th);
        nn12_dist_th  = nn12_dist_th * Config::descThL();
        // resort according to the queryIdx
        sort( lmatches_12.begin(), lmatches_12.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( lmatches_21.begin(), lmatches_21.end(), sort_descriptor_by_queryIdx() );
        // bucle around lmatches
        for( int i = 0; i < lmatches_12.size(); i++ )
        {
            // check if they are mutual best matches
            int lr_qdx = lmatches_12[i][0].queryIdx;
            int lr_tdx = lmatches_12[i][0].trainIdx;
            int rl_tdx;
            if( Config::bestLRMatches() )
                rl_tdx = lmatches_21[lr_tdx][0].trainIdx;
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_12 = lmatches_12[i][1].distance - lmatches_12[i][0].distance;            
            // this shouldn't happen
            if( kf0->stereo_frame->stereo_ls[lr_qdx] == NULL || kf1->stereo_frame->stereo_ls[lr_tdx] == NULL )
                continue;
            if( lr_qdx == rl_tdx  && dist_12 > nn12_dist_th )
            {
                Vector3d sP_ = DT.block(0,0,3,3) * kf0->stereo_frame->stereo_ls[lr_qdx]->sP + DT.col(3).head(3);
                Vector2d spl_proj = cam->projection( sP_ );
                Vector3d eP_ = DT.block(0,0,3,3) * kf0->stereo_frame->stereo_ls[lr_qdx]->eP + DT.col(3).head(3);
                Vector2d epl_proj = cam->projection( eP_ );
                Vector3d l_obs = kf0->stereo_frame->stereo_ls[lr_qdx]->le;
                // projection error
                Vector2d err_ls;
                err_ls(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
                err_ls(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
                if( err_ls.norm() *sqrt(kf0->stereo_frame->stereo_ls[lr_qdx]->sigma2) < sqrt(7.815) )
                {
                    common_ls++;
                    // new 3D landmark
                    if( kf0->stereo_frame->stereo_ls[lr_qdx]->idx == -1 )
                    {
                        // assign indices
                        kf0->stereo_frame->stereo_ls[lr_qdx]->idx = max_ls_idx;
                        kf1->stereo_frame->stereo_ls[lr_tdx]->idx = max_ls_idx;
                        // create new 3D landmark with the observation from previous KF
                        Matrix4d Tfw = ( kf0->T_kf_w );
                        Vector3d sP3d = Tfw.block(0,0,3,3) * kf0->stereo_frame->stereo_ls[lr_qdx]->sP + Tfw.col(3).head(3);
                        Vector3d eP3d = Tfw.block(0,0,3,3) * kf0->stereo_frame->stereo_ls[lr_qdx]->eP + Tfw.col(3).head(3);
                        Vector6d L3d;
                        L3d << sP3d, eP3d;
                        Vector3d mP3d = 0.5*(sP3d+eP3d);
                        mP3d = mP3d / mP3d.norm();
                        Vector4d pts;
                        pts << kf0->stereo_frame->stereo_ls[lr_qdx]->spl, kf0->stereo_frame->stereo_ls[lr_qdx]->epl;
                        MapLine* map_line = new MapLine(max_ls_idx,L3d,kf0->stereo_frame->ldesc_l.row(lr_qdx),kf0->kf_idx,kf0->stereo_frame->stereo_ls[lr_qdx]->le,mP3d,pts);
                        mP3d = 0.5*( kf1->stereo_frame->stereo_ls[lr_tdx]->sP + kf1->stereo_frame->stereo_ls[lr_tdx]->eP );
                        mP3d = kf1->T_kf_w.block(0,0,3,3) * mP3d + kf1->T_kf_w.col(3).head(3);
                        mP3d = mP3d / mP3d.norm();
                        pts << kf1->stereo_frame->stereo_ls[lr_tdx]->spl, kf1->stereo_frame->stereo_ls[lr_tdx]->epl;
                        map_line->addMapLineObservation(kf1->stereo_frame->ldesc_l.row(lr_tdx),kf1->kf_idx,kf1->stereo_frame->stereo_ls[lr_tdx]->le,mP3d,pts);
                        // add new 3D landmark to kf_idx where it was first observed
                        map_lines_kf_idx.at( kf0->kf_idx ).push_back( max_ls_idx );
                        // add 3D landmark to map
                        map_lines.push_back(map_line);
                        max_ls_idx++;
                        // update full graph (new feature)
                        full_graph[kf1_idx][kf0_idx]++;
                        full_graph[kf0_idx][kf1_idx]++;
                        // if has refine pose:
                        if( has_refinement )
                        {
                            LineFeature* ls = kf0->stereo_frame->stereo_ls[lr_qdx];
                            ls->sdisp_obs   = kf1->stereo_frame->stereo_ls[lr_tdx]->sdisp;
                            ls->edisp_obs   = kf1->stereo_frame->stereo_ls[lr_tdx]->edisp;
                            ls->spl_obs     = kf1->stereo_frame->stereo_ls[lr_tdx]->spl;
                            ls->epl_obs     = kf1->stereo_frame->stereo_ls[lr_tdx]->epl;
                            ls->le_obs      = kf1->stereo_frame->stereo_ls[lr_tdx]->le;
                            ls->inlier      = true;
                            matched_ls.push_back( ls );
                        }
                    }
                    // 3D landmark exists: copy idx && add observation to map landmark
                    else
                    {
                        // copy idx
                        int lm_idx = kf0->stereo_frame->stereo_ls[lr_qdx]->idx;
                        if( map_lines[lm_idx] != NULL )
                        {
                            kf1->stereo_frame->stereo_ls[lr_tdx]->idx = lm_idx;
                            // add observation of the 3D landmark from current KF
                            Vector3d mP3d = 0.5*(kf1->stereo_frame->stereo_ls[lr_tdx]->sP+kf1->stereo_frame->stereo_ls[lr_tdx]->eP);
                            mP3d = kf1->T_kf_w.block(0,0,3,3) * mP3d + kf1->T_kf_w.col(3).head(3);
                            mP3d = mP3d / mP3d.norm();
                            Vector4d pts;
                            pts << kf1->stereo_frame->stereo_ls[lr_tdx]->spl, kf1->stereo_frame->stereo_ls[lr_tdx]->epl;
                            map_lines[lm_idx]->addMapLineObservation(kf1->stereo_frame->ldesc_l.row(lr_tdx),kf1->kf_idx,kf1->stereo_frame->stereo_ls[lr_tdx]->le,mP3d,pts);
                            // update full graph (previously observed feature)
                            for( vector<int>::iterator obs = map_lines[lm_idx]->kf_obs_list.begin(); obs != map_lines[lm_idx]->kf_obs_list.end(); obs++ )
                            {
                                if( (*obs) != kf1_idx )
                                {
                                    full_graph[kf1_idx][(*obs)]++;
                                    full_graph[(*obs)][kf1_idx]++;
                                }
                            }
                            // if has refine pose:
                            if( has_refinement )
                            {
                                LineFeature* ls = kf0->stereo_frame->stereo_ls[lr_qdx];
                                ls->sdisp_obs   = kf1->stereo_frame->stereo_ls[lr_tdx]->sdisp;
                                ls->edisp_obs   = kf1->stereo_frame->stereo_ls[lr_tdx]->edisp;
                                ls->spl_obs     = kf1->stereo_frame->stereo_ls[lr_tdx]->spl;
                                ls->epl_obs     = kf1->stereo_frame->stereo_ls[lr_tdx]->epl;
                                ls->le_obs      = kf1->stereo_frame->stereo_ls[lr_tdx]->le;
                                ls->inlier      = true;
                                matched_ls.push_back( ls );
                            }
                        }
                    }
                }
            }
        }
    }

    // ---------------------------------------------------
    // refine pose between kf0 and kf1
    // ---------------------------------------------------
    if( has_refinement )
    {
        StVO::StereoFrameHandler* stf = new StereoFrameHandler( cam );
        stf->matched_pt = matched_pt;
        stf->matched_ls = matched_ls;
        stf->prev_frame = kf0->stereo_frame;
        stf->curr_frame = kf1->stereo_frame;
        stf->optimizePose( DT );
        Matrix4d DT_         = stf->curr_frame->DT;
        Vector6d Dt_cov_eig_ = stf->curr_frame->DT_cov_eig;
        if( Dt_cov_eig_(0) <= 0.01 )
            kf1->T_kf_w = kf0->T_kf_w * inverse_se3( DT_ );
    }

    // ---------------------------------------------------
    // find point matches between local map and curr_frame
    // ---------------------------------------------------
    // select local map
    vector<MapPoint*> map_local_points;
    Mat map_lpt_desc;
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++ )
    {
        if((*pt_it)!=NULL)
        {
            // if it is local and not found in the current KF
            if( (*pt_it)->local && (*pt_it)->kf_obs_list.back() != kf1_idx )
            {
                // if the LM is projected inside the current image
                Vector3d Pf  = Twf.block(0,0,3,3) * (*pt_it)->point3D + Twf.col(3).head(3);
                Vector2d pf   = cam->projection( Pf );
                if( pf(0) > 0 && pf(0) < cam->getWidth() && pf(1) > 0 && pf(1) < cam->getHeight() && Pf(2) > 0.0 )
                {
                    // add the point and its representative descriptor
                    map_local_points.push_back( (*pt_it) );
                    map_lpt_desc.push_back( (*pt_it)->med_desc.row(0) );
                }
            }
        }
    }

    // select unmatched points
    vector<PointFeature*> unmatched_points;
    Mat unmatched_pt_desc;
    int it = 0;
    for( vector<PointFeature*>::iterator pt_it = kf1->stereo_frame->stereo_pt.begin(); pt_it != kf1->stereo_frame->stereo_pt.end(); pt_it++, it++ )
    {
        if( (*pt_it)->idx == -1 )
        {
            unmatched_points.push_back( (*pt_it) );
            unmatched_pt_desc.push_back( kf1->stereo_frame->pdesc_l.row(it) );
        }
    }

    // track points from local map
    if( Config::hasPoints() && !(unmatched_points.size()==0) && !(map_local_points.size()==0)  )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        vector<vector<DMatch>> pmatches_12, pmatches_21;
        // 12 and 21 matches
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StVO::StereoFrame::matchPointFeatures, kf0->stereo_frame, bfm, map_lpt_desc, unmatched_pt_desc, ref(pmatches_12) );
                auto match_r = async( launch::async, &StVO::StereoFrame::matchPointFeatures, kf0->stereo_frame, bfm, unmatched_pt_desc, map_lpt_desc, ref(pmatches_21) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( map_lpt_desc, unmatched_pt_desc, pmatches_12, 2);
                bfm->knnMatch( unmatched_pt_desc, map_lpt_desc, pmatches_21, 2);
            }
        }
        else
            bfm->knnMatch( map_lpt_desc, unmatched_pt_desc, pmatches_12, 2);
        // sort matches by the distance between the best and second best matches
        double nn12_dist_th = Config::minRatio12P();
        // resort according to the queryIdx
        sort( pmatches_12.begin(), pmatches_12.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( pmatches_21.begin(), pmatches_21.end(), sort_descriptor_by_queryIdx() );
        // bucle around pmatches
        for( int i = 0; i < pmatches_12.size(); i++ )
        {
            // check if they are mutual best matches
            int lr_qdx = pmatches_12[i][0].queryIdx;
            int lr_tdx = pmatches_12[i][0].trainIdx;
            int rl_tdx;
            if( Config::bestLRMatches() )
                rl_tdx = pmatches_21[lr_tdx][0].trainIdx;
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_12 = pmatches_12[i][0].distance / pmatches_12[i][1].distance;
            double nn12_dist_th = Config::minRatio12P();
            // check the f2f 12 distance condition
            if( lr_qdx == rl_tdx  && dist_12 > nn12_dist_th )
            {
                Vector3d Pf_map = Twf.block(0,0,3,3) * map_local_points[lr_qdx]->point3D + Twf.col(3).head(3);
                Vector3d Pf_kf  = unmatched_points[lr_tdx]->P;
                // check that the point is projected in front of the camera
                if( Pf_map(2) > 0.0 )
                {
                    // check that the viewing direction condition is also satisfied
                    Vector3d dir_map = Twf.block(0,0,3,3) * map_local_points[lr_qdx]->med_obs_dir + Twf.col(3).head(3);
                    Vector3d dir_kf  = Pf_kf  / Pf_kf.norm();
                    // check that the epipolar constraint is satisfied
                    Vector2d pf_map = cam->projection( Pf_map );
                    Vector2d pf_kf  = unmatched_points[lr_tdx]->pl;
                    double error_epip = ( pf_map - pf_kf ).norm();
                    if( error_epip < Config::maxKFEpipP() )
                    {
                        common_pt++;
                        // copy idx
                        int lm_idx = map_local_points[lr_qdx]->idx;
                        unmatched_points[lr_tdx]->idx = lm_idx;
                        // add observation of the 3D LM from current KF
                        dir_kf = kf1->T_kf_w.block(0,0,3,3) * dir_kf + kf1->T_kf_w.col(3).head(3);
                        map_points[lm_idx]->addMapPointObservation( unmatched_pt_desc.row(lr_tdx), kf1_idx, unmatched_points[lr_tdx]->pl, dir_kf );
                        // update full graph (previously observed feature)
                        for( vector<int>::iterator obs = map_points[lm_idx]->kf_obs_list.begin(); obs != map_points[lm_idx]->kf_obs_list.end(); obs++ )
                        {
                            if( (*obs) != kf1_idx )
                            {
                                full_graph[kf1_idx][(*obs)]++;
                                full_graph[(*obs)][kf1_idx]++;
                            }
                        }
                    }
                }
            }
        }
    }

    // ---------------------------------------------------
    // find line matches between local map and curr_frame
    // ---------------------------------------------------
    // select local map
    vector<MapLine*> map_local_lines;
    Mat map_lls_desc;
    for( vector<MapLine*>::iterator ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++ )
    {
        if((*ls_it)!=NULL)
        {
            // if it is local and not found in the current KF
            if( (*ls_it)->local && (*ls_it)->kf_obs_list.back() != kf1_idx )
            {
                // if the LM is projected inside the current image
                Vector3d sPf  = Twf.block(0,0,3,3) * (*ls_it)->line3D.head(3) + Twf.col(3).head(3);
                Vector2d spf   = cam->projection( sPf );
                Vector3d ePf  = Twf.block(0,0,3,3) * (*ls_it)->line3D.tail(3) + Twf.col(3).head(3);
                Vector2d epf   = cam->projection( ePf );
                if( spf(0) > 0 && spf(0) < cam->getWidth() && spf(1) > 0 && spf(1) < cam->getHeight() && sPf(2) > 0.0 )
                {
                    if( epf(0) > 0 && epf(0) < cam->getWidth() && epf(1) > 0 && epf(1) < cam->getHeight() && ePf(2) > 0.0 )
                    {
                        // add the point and its representative descriptor
                        map_local_lines.push_back( (*ls_it) );
                        map_lls_desc.push_back( (*ls_it)->med_desc.row(0) );
                    }
                }
            }
        }
    }

    // select unmatched line segments
    vector<LineFeature*> unmatched_lines;
    Mat unmatched_ls_desc;
    it = 0;
    for( vector<LineFeature*>::iterator ls_it = kf1->stereo_frame->stereo_ls.begin(); ls_it != kf1->stereo_frame->stereo_ls.end(); ls_it++, it++ )
    {
        if( (*ls_it)->idx == -1 )
        {
            unmatched_lines.push_back( (*ls_it) );
            unmatched_ls_desc.push_back( kf1->stereo_frame->ldesc_l.row(it) );
        }
    }

    // track line segments from local map
    if( Config::hasLines() && !(unmatched_lines.size()==0) && !(map_local_lines.size()==0)  )
    {
        Ptr<BinaryDescriptorMatcher> bdm = BinaryDescriptorMatcher::createBinaryDescriptorMatcher();
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        vector<vector<DMatch>> lmatches_12, lmatches_21;
        // 12 and 21 matches
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StVO::StereoFrame::matchLineFeatures, kf0->stereo_frame, bfm, map_lls_desc, unmatched_ls_desc, ref(lmatches_12) );
                auto match_r = async( launch::async, &StVO::StereoFrame::matchLineFeatures, kf0->stereo_frame, bfm, unmatched_ls_desc, map_lls_desc, ref(lmatches_21) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( map_lls_desc, unmatched_ls_desc, lmatches_12, 2);
                bfm->knnMatch( unmatched_ls_desc, map_lls_desc, lmatches_21, 2);
            }
        }
        else
            bfm->knnMatch( map_lls_desc, unmatched_ls_desc, lmatches_12, 2);

        // sort matches by the distance between the best and second best matches
        double nn_dist_th, nn12_dist_th;
        kf0->stereo_frame->lineDescriptorMAD(lmatches_12,nn_dist_th, nn12_dist_th);
        nn12_dist_th  = nn12_dist_th * Config::descThL();
        // resort according to the queryIdx
        sort( lmatches_12.begin(), lmatches_12.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( lmatches_21.begin(), lmatches_21.end(), sort_descriptor_by_queryIdx() );
        // bucle around lmatches
        for( int i = 0; i < lmatches_12.size(); i++ )
        {
            // check if they are mutual best matches
            int lr_qdx = lmatches_12[i][0].queryIdx;
            int lr_tdx = lmatches_12[i][0].trainIdx;
            if( lr_qdx != -1 && lr_tdx != -1 )
            {
                int rl_tdx;
                if( Config::bestLRMatches() )
                    rl_tdx = lmatches_21[lr_tdx][0].trainIdx;
                else
                    rl_tdx = lr_qdx;
                // check if they are mutual best matches and the minimum distance
                double dist_12 = lmatches_12[i][1].distance - lmatches_12[i][0].distance;
                if( lr_qdx == rl_tdx && dist_12 > nn12_dist_th )
                {
                    Vector3d sP_ = Twf.block(0,0,3,3) * map_local_lines[lr_qdx]->line3D.head(3) + Twf.col(3).head(3);
                    Vector2d spl_proj = cam->projection( sP_ );
                    Vector3d eP_ = Twf.block(0,0,3,3) * map_local_lines[lr_qdx]->line3D.tail(3) + Twf.col(3).head(3);
                    Vector2d epl_proj = cam->projection( eP_ );
                    Vector3d l_obs = unmatched_lines[lr_tdx]->le;
                    if( sP_(2) > 0.0 && eP_(2) > 0.0 )
                    {
                        // check that the features have been observed from a similar point of view
                        Vector3d mP_ = 0.5 * ( sP_ + eP_ );
                        Vector3d dir_map = cam->projectionNH( mP_ );
                        dir_map = Twf.block(0,0,3,3) * map_local_lines[lr_qdx]->med_obs_dir + Twf.col(3).head(3);
                        mP_  = 0.5 * ( unmatched_lines[lr_tdx]->sP + unmatched_lines[lr_tdx]->eP  );
                        Vector3d dir_kf;
                        dir_kf = mP_ / mP_.norm();
                        // check the epipolar constraint
                        Vector2d err_ls;
                        err_ls(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
                        err_ls(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
                        if( err_ls(0) < Config::maxKFEpipL() && err_ls(1) < Config::maxKFEpipL() )
                        {
                            common_ls++;
                            // copy idx
                            int lm_idx = map_local_lines[lr_qdx]->idx;
                            unmatched_lines[lr_tdx]->idx = lm_idx;
                            // add observation of the 3D landmark from current KF
                            Vector3d mP3d = 0.5*(unmatched_lines[lr_tdx]->sP+unmatched_lines[lr_tdx]->eP);
                            mP3d = kf1->T_kf_w.block(0,0,3,3) * mP3d + kf1->T_kf_w.col(3).head(3);
                            mP3d = mP3d / mP3d.norm();
                            Vector4d pts;
                            pts << unmatched_lines[lr_tdx]->spl, unmatched_lines[lr_tdx]->epl;
                            map_lines[lm_idx]->addMapLineObservation( kf1->stereo_frame->ldesc_l.row(lr_tdx), kf1_idx, unmatched_lines[lr_tdx]->le, mP3d, pts );
                            // update full graph (previously observed feature)
                            for( vector<int>::iterator obs = map_lines[lm_idx]->kf_obs_list.begin(); obs != map_lines[lm_idx]->kf_obs_list.end(); obs++ )
                            {
                                if( (*obs) != kf1_idx )
                                {
                                    full_graph[kf1_idx][(*obs)]++;
                                    full_graph[(*obs)][kf1_idx]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}

void MapHandler::expandGraphs()
{
    int g_size = full_graph.size() + 1;
    // full graph
    full_graph.resize( g_size );
    for(unsigned int i = 0; i < g_size; i++ )
    {
        full_graph[i].resize( g_size );
    }
    // confusion matrix
    conf_matrix.resize( g_size );
    for(unsigned int i = 0; i < g_size; i++ )
        conf_matrix[i].resize( g_size );
}

void MapHandler::formLocalMap()
{

    // reset local KFs & LMs
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++ )
    {
        if( (*kf_it) != NULL )
            (*kf_it)->local = false;
    }
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++ )
    {
        if( (*pt_it) != NULL )
            (*pt_it)->local = false;
    }
    for( vector<MapLine*>::iterator  ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++ )
    {
        if( (*ls_it) != NULL )
            (*ls_it)->local = false;
    }

    // set first KF and their associated LMs as local
    map_keyframes.back()->local = true;
    for( vector<PointFeature*>::iterator pt_it = map_keyframes.back()->stereo_frame->stereo_pt.begin(); pt_it != map_keyframes.back()->stereo_frame->stereo_pt.end(); pt_it++ )
    {
        if( (*pt_it) != NULL )
        {
            int lm_idx = (*pt_it)->idx;
            if( lm_idx != -1 && map_points[lm_idx] != NULL )
                map_points[lm_idx]->local = true;
        }
    }
    for( vector<LineFeature*>::iterator ls_it = map_keyframes.back()->stereo_frame->stereo_ls.begin(); ls_it != map_keyframes.back()->stereo_frame->stereo_ls.end(); ls_it++ )
    {
        if( (*ls_it) != NULL )
        {
            int lm_idx = (*ls_it)->idx;
            if( lm_idx != -1 && map_lines[lm_idx] != NULL )
                map_lines[lm_idx]->local = true;
        }
    }

    // loop over covisibility graph / full graph if we want to find more points
    int g_size = full_graph.size()-1;
    for( int i = 0; i < g_size; i++ )
    {
        if( full_graph[g_size][i] >= Config::minLMCovGraph() || abs(g_size-i) <= Config::minKFLocalMap() )
        {
            map_keyframes[i]->local = true;
            // loop over the landmarks seen by KF{i}
            for( vector<PointFeature*>::iterator pt_it = map_keyframes[i]->stereo_frame->stereo_pt.begin(); pt_it != map_keyframes[i]->stereo_frame->stereo_pt.end(); pt_it++ )
            {
                int lm_idx = (*pt_it)->idx;
                if( lm_idx != -1 && map_points[lm_idx] != NULL )
                    map_points[lm_idx]->local = true;
            }
            for( vector<LineFeature*>::iterator ls_it = map_keyframes[i]->stereo_frame->stereo_ls.begin(); ls_it != map_keyframes[i]->stereo_frame->stereo_ls.end(); ls_it++ )
            {
                int lm_idx = (*ls_it)->idx;
                if( lm_idx != -1 && map_lines[lm_idx] != NULL )
                    map_lines[lm_idx]->local = true;
            }
        }
    }

}

void MapHandler::fuseLandmarksFromLocalMap()
{

    // ---------------------------------------------------
    // find point matches between local map and local map
    // ---------------------------------------------------
    // select local map
    vector<MapPoint*> map_local_points;
    vector<int> map_local_points_idx;
    Mat map_pps_desc;
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++ )
    {
        if((*pt_it)!=NULL)
        {
            // if it is local and not found in the current KF
            if( (*pt_it)->local )
            {
                Matrix4d Twf = map_keyframes[ (*pt_it)->kf_obs_list[0] ]->T_kf_w;
                // if the LM is projected inside the current image
                Vector3d Pf  = Twf.block(0,0,3,3) * (*pt_it)->point3D + Twf.col(3).head(3);
                Vector2d pf   = cam->projection( Pf );

                if( pf(0) > 0 && pf(0) < cam->getWidth() && pf(1) > 0 && pf(1) < cam->getHeight() && Pf(2) > 0.0 )
                {
                    // add the point and its representative descriptor
                    map_local_points.push_back( (*pt_it) );
                    map_pps_desc.push_back( (*pt_it)->med_desc.row(0) );
                    map_local_points_idx.push_back( (*pt_it)->idx );
                }
            }
        }
    }

    // track line segments from local map
    vector<Vector2i> map_local_points_to_fuse;
    if( Config::hasPoints() && map_local_points.size() > 1  )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        vector<vector<DMatch>> lmatches_12;
        // 12 and 21 matches
        bfm->knnMatch( map_pps_desc, map_pps_desc, lmatches_12, 2);
        // bucle around lmatches
        for( int i = 0; i < lmatches_12.size(); i++ )
        {
            // check if they are mutual best matches
            int lr_qdx = lmatches_12[i][1].queryIdx;
            int lr_tdx = lmatches_12[i][1].trainIdx;
            if( lr_qdx != -1 && lr_tdx != -1 && lr_qdx != lr_tdx && map_local_points[lr_qdx] != NULL && map_local_points[lr_tdx] != NULL  )
            {
                // compare the directions of the 3D lines
                Vector3d P1 = map_local_points[lr_qdx]->point3D;
                Vector3d P2 = map_local_points[lr_tdx]->point3D;
                double Pd = (P2-P1).norm();
                if( Pd < Config::maxPointPointError() )
                {
                    // add LM indices to list
                    Vector2i idx;
                    idx(0) = lr_qdx;
                    idx(1) = lr_tdx;
                    map_local_points_to_fuse.push_back(idx);
                }
            }
        }
    }

    // fuse landmark list
    int k = 0;
    for( auto it = map_local_points_to_fuse.begin(); it != map_local_points_to_fuse.end(); it++, k++ )
    {
        int idx_1 = (*it)(0);
        int idx_2 = (*it)(1);
        if( map_local_points[idx_1]!=NULL && map_local_points[idx_2]!=NULL && idx_1 != idx_2 )
        {
            MapPoint* P1 = map_local_points[idx_1];
            MapPoint* P2 = map_local_points[idx_2];
            // update graphs
            for( auto ldx1 = P1->kf_obs_list.begin(); ldx1 != P1->kf_obs_list.end(); ldx1++ )
            {
                for( auto ldx2 = P2->kf_obs_list.begin(); ldx2 != P2->kf_obs_list.end(); ldx2++ )
                {
                    full_graph[(*ldx1)][(*ldx2)]++;
                    full_graph[(*ldx2)][(*ldx1)]++;
                }
            }
            // add observations
            for( int i = 0; i < P2->kf_obs_list.size(); i++ )
                P1->addMapPointObservation( P2->desc_list[i],
                                            P2->kf_obs_list[i],
                                            P2->obs_list[i],
                                            P2->dir_list[i] );
            // remove landmark
            map_points[ map_local_points_idx[k] ] = NULL;
        }
    }

    // ---------------------------------------------------
    // find line matches between local map and local map
    // ---------------------------------------------------
    // select local map
    vector<MapLine*> map_local_lines;
    vector<int> map_local_lines_idx;
    Mat map_lls_desc;
    for( vector<MapLine*>::iterator ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++ )
    {
        if((*ls_it)!=NULL)
        {
            // if it is local and not found in the current KF
            if( (*ls_it)->local )
            {
                Matrix4d Twf = map_keyframes[ (*ls_it)->kf_obs_list[0] ]->T_kf_w;
                // if the LM is projected inside the current image
                Vector3d sPf  = Twf.block(0,0,3,3) * (*ls_it)->line3D.head(3) + Twf.col(3).head(3);
                Vector2d spf   = cam->projection( sPf );
                Vector3d ePf  = Twf.block(0,0,3,3) * (*ls_it)->line3D.tail(3) + Twf.col(3).head(3);
                Vector2d epf   = cam->projection( ePf );
                if( spf(0) > 0 && spf(0) < cam->getWidth() && spf(1) > 0 && spf(1) < cam->getHeight() && sPf(2) > 0.0 )
                {
                    if( epf(0) > 0 && epf(0) < cam->getWidth() && epf(1) > 0 && epf(1) < cam->getHeight() && ePf(2) > 0.0 )
                    {
                        // add the point and its representative descriptor
                        map_local_lines.push_back( (*ls_it) );
                        map_lls_desc.push_back( (*ls_it)->med_desc.row(0) );
                        map_local_lines_idx.push_back( (*ls_it)->idx );
                    }
                }
            }
        }
    }

    // track line segments from local map
    vector<Vector2i> map_local_lines_to_fuse;
    if( Config::hasLines() && map_local_lines.size() > 1  )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        vector<vector<DMatch>> lmatches_12;
        // 12 and 21 matches
        bfm->knnMatch( map_lls_desc, map_lls_desc, lmatches_12, 2);
        // bucle around lmatches
        for( int i = 0; i < lmatches_12.size(); i++ )
        {
            // check if they are mutual best matches
            int lr_qdx = lmatches_12[i][1].queryIdx;
            int lr_tdx = lmatches_12[i][1].trainIdx;
            if( lr_qdx != -1 && lr_tdx != -1 && lr_qdx != lr_tdx && map_local_lines[lr_qdx]!=NULL && map_local_lines[lr_tdx]!=NULL )
            {
                /*// compare the directions of the 3D lines
                Vector6d L1 = map_local_lines[lr_qdx]->line3D;
                Vector6d L2 = map_local_lines[lr_tdx]->line3D;
                Vector3d d1 = (L1.tail(3)-L1.head(3)) / (L1.tail(3)-L1.head(3)).norm() ;
                Vector3d d2 = (L2.tail(3)-L2.head(3)) / (L2.tail(3)-L2.head(3)).norm() ;
                double Ld = (d1-d2).norm();
                // compute the distance of the endpoints to the first 3D line
                vector<double> lambda1, lambda2;
                for( int j = 0; j < 3; j++ )
                {
                    lambda1.push_back( (L2(j)-L1(j)) / d1(j)  );
                    lambda1.push_back( (L2(j+3)-L1(j)) / d1(j)  );
                    lambda2.push_back( (L1(j)-L2(j)) / d2(j)  );
                    lambda2.push_back( (L1(j+3)-L2(j)) / d2(j)  );
                }
                double std1 = vector_stdv(lambda1);
                double std2 = vector_stdv(lambda2);
                if( Ld < Config::maxDirLineError() && std1 < Config::maxPointLineError() && std2 < Config::maxPointLineError() )*/

                // compare the distance from the points to the line
                Vector6d L1 = map_local_lines[lr_qdx]->line3D;
                Vector6d L2 = map_local_lines[lr_tdx]->line3D;
                Vector3d d1 = (L1.tail(3)-L1.head(3)) / (L1.tail(3)-L1.head(3)).norm() ;
                Vector3d d2 = (L2.tail(3)-L2.head(3)) / (L2.tail(3)-L2.head(3)).norm() ;
                double Ld = (d1-d2).norm();
                double dP, dQ;
                if( d1.norm() > d2.norm() ) //dist from the points in the 2nd to the 1st line
                {
                    Vector3d x1 = L1.head(3);
                    Vector3d  p = L2.head(3);
                    Vector3d  q = L2.tail(3);
                    dP = ( (p-x1) - (p-x1).dot(d1) * d1 ).norm();
                    dQ = ( (q-x1) - (q-x1).dot(d1) * d1 ).norm();
                }
                else
                {
                    Vector3d x2 = L2.head(3);
                    Vector3d  p = L1.head(3);
                    Vector3d  q = L1.tail(3);
                    dP = ( (p-x2) - (p-x2).dot(d2) * d2 ).norm();
                    dQ = ( (q-x2) - (q-x2).dot(d2) * d2 ).norm();
                }
                if( Ld < Config::maxDirLineError() && dP < Config::maxPointLineError() && dQ < Config::maxPointLineError() )
                {
                    // add LM indices to list
                    Vector2i idx;
                    idx(0) = lr_qdx;
                    idx(1) = lr_tdx;
                    map_local_lines_to_fuse.push_back(idx);
                }
            }
        }
    }

    // fuse landmark list
    k = 0;
    for( auto it = map_local_lines_to_fuse.begin(); it != map_local_lines_to_fuse.end(); it++, k++ )
    {
        int idx_1 = (*it)(0);
        int idx_2 = (*it)(1);
        if( map_local_lines[idx_1]!=NULL && map_local_lines[idx_2]!=NULL && idx_1 != idx_2 )
        {
            // look for the biggest line
            double d1 = ( map_local_lines[idx_1]->line3D.head(3)-map_local_lines[idx_1]->line3D.tail(3) ).norm();
            double d2 = ( map_local_lines[idx_2]->line3D.head(3)-map_local_lines[idx_2]->line3D.tail(3) ).norm();
            MapLine* L1;
            MapLine* L2;
            if( d1 > d2 )
            {
                L1 = map_local_lines[idx_1];
                L2 = map_local_lines[idx_2];
            }
            else
            {
                L1 = map_local_lines[idx_2];
                L2 = map_local_lines[idx_1];
            }
            // update graphs
            for( auto ldx1 = L1->kf_obs_list.begin(); ldx1 != L1->kf_obs_list.end(); ldx1++ )
            {
                for( auto ldx2 = L2->kf_obs_list.begin(); ldx2 != L2->kf_obs_list.end(); ldx2++ )
                {
                    full_graph[(*ldx1)][(*ldx2)]++;
                    full_graph[(*ldx2)][(*ldx1)]++;
                }
            }
            // add observations
            for( int i = 0; i < L2->kf_obs_list.size(); i++ )
                L1->addMapLineObservation( L2->desc_list[i],
                                           L2->kf_obs_list[i],
                                           L2->obs_list[i],
                                           L2->dir_list[i],
                                           L2->pts_list[i] );
            // remove landmark
            map_lines[ map_local_lines_idx[k] ] = NULL;
        }
    }

}

// -----------------------------------------------------------------------------------------------------------------------------
// Local Bundle Adjustment functions
// -----------------------------------------------------------------------------------------------------------------------------

void MapHandler::localBundleAdjustment()
{

    vector<double> X_aux;

    // create list of local keyframes
    vector<int> kf_list;
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++)
    {
        if( (*kf_it)!= NULL )
        {
            if( (*kf_it)->local && (*kf_it)->kf_idx != 0 )
            {
                Vector6d pose_aux = (*kf_it)->x_kf_w;
                for(int i = 0; i < 6; i++)
                    X_aux.push_back( pose_aux(i) );
                kf_list.push_back( (*kf_it)->kf_idx );
            }
        }
    }

    // create list of local point landmarks
    vector<Vector6i> pt_obs_list;
    vector<int> pt_list;
    int lm_local_idx = 0;
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++)
    {
        if( (*pt_it)!= NULL )
        {
            if( (*pt_it)->local )
            {
                Vector3d point_aux = (*pt_it)->point3D;
                for(int i = 0; i < 3; i++)
                    X_aux.push_back( point_aux(i) );
                // gather all observations
                for( int i = 0; i < (*pt_it)->obs_list.size(); i++)
                {
                    Vector6i obs_aux;
                    obs_aux(0) = (*pt_it)->idx; // LM idx
                    obs_aux(1) = lm_local_idx;  // LM local idx
                    obs_aux(2) = i;             // LM obs idx
                    int kf_obs_list_ = (*pt_it)->kf_obs_list[i];
                    obs_aux(3) = kf_obs_list_;  // KF idx
                    obs_aux(4) = -1;            // KF local idx (-1 if not local)
                    obs_aux(5) = 1;             // 1 if the observation is an inlier
                    for( int j = 0; j < kf_list.size(); j++ )
                    {
                        if( kf_list[j] == kf_obs_list_ )
                        {
                            obs_aux(4) = j;
                            break;
                        }
                    }
                    pt_obs_list.push_back( obs_aux );
                }
                lm_local_idx++;
                // pt_list
                pt_list.push_back( (*pt_it)->idx );
            }
        }
    }

    // create list of local line segment landmarks
    vector<Vector6i> ls_obs_list;
    vector<int> ls_list;
    lm_local_idx = 0;
    for( vector<MapLine*>::iterator ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++)
    {
        if( (*ls_it)!= NULL )
        {
            if( (*ls_it)->local )
            {
                Vector6d line_aux = (*ls_it)->line3D;
                for(int i = 0; i < 6; i++)
                    X_aux.push_back( line_aux(i) );
                // gather all observations
                for( int i = 0; i < (*ls_it)->obs_list.size(); i++)
                {
                    Vector6i obs_aux;
                    obs_aux(0) = (*ls_it)->idx; // LM idx
                    obs_aux(1) = lm_local_idx;  // LM local idx
                    obs_aux(2) = i;             // LM obs idx
                    int kf_obs_list_ = (*ls_it)->kf_obs_list[i];
                    obs_aux(3) = kf_obs_list_;  // KF idx
                    obs_aux(4) = -1;            // KF local idx (-1 if not local)
                    obs_aux(5) = 1;             // 1 if the observation is an inlier
                    for( int j = 0; j < kf_list.size(); j++ )
                    {
                        if( kf_list[j] == kf_obs_list_ )
                        {
                            obs_aux(4) = j;
                            break;
                        }
                    }
                    ls_obs_list.push_back( obs_aux );
                }
                lm_local_idx++;
                // ls_list
                ls_list.push_back( (*ls_it)->idx );
            }
        }
    }

    // Levenberg-Marquardt optimization of the local map
    if( pt_obs_list.size() + ls_obs_list.size() != 0 )
        levMarquardtOptimizationLBA(X_aux,kf_list,pt_list,ls_list,pt_obs_list,ls_obs_list);

}

void MapHandler::levMarquardtOptimizationLBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  )
{

    // create Levenberg-Marquardt variables
    int    Nkf = kf_list.size();
    int      N = X_aux.size();

    VectorXd X = VectorXd::Zero(N), DX = VectorXd::Zero(N);
    VectorXd g = VectorXd::Zero(N);
    MatrixXd H = MatrixXd::Zero(N,N);

    for(int i = 0; i < N; i++)
        X(i) = X_aux[i];
    SparseMatrix<double> H_(N,N);

    // create Levenberg-Marquardt parameters
    double err = 0.0, err_prev = 999999999.9;
    double lambda = Config::lambdaLbaLM(), lambda_k = Config::lambdaLbaK();
    int    max_iters = Config::maxItersLba();

    // estimate H and g to precalculate lambda
    //---------------------------------------------------------------------------------------------
    // point observations
    int Npt = 0, Npt_obs = 0;
    if( pt_obs_list.size() != 0 )
        Npt = pt_obs_list.back()(1)+1;
    for( vector<Vector6i>::iterator pt_it = pt_obs_list.begin(); pt_it != pt_obs_list.end(); pt_it++ )
    {
        int lm_idx_map = (*pt_it)(0);
        int lm_idx_loc = (*pt_it)(1);
        int lm_idx_obs = (*pt_it)(2);
        int kf_idx_map = (*pt_it)(3);
        int kf_idx_loc = (*pt_it)(4);
        if( map_points[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
        {
            // grab 3D LM (Xwj)
            Vector3d Xwj   = map_points[lm_idx_map]->point3D;
            // grab 6DoF KF (Tiw)
            Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
            // projection error
            Tiw = inverse_se3( Tiw );
            Vector3d Xwi   = Tiw.block(0,0,3,3) * Xwj + Tiw.block(0,3,3,1);
            Vector2d p_prj = cam->projection( Xwi );
            Vector2d p_obs = map_points[lm_idx_map]->obs_list[lm_idx_obs];
            Vector2d p_err    = p_obs - p_prj;
            double p_err_norm = p_err.norm();
            // estimate useful variables
            double gx   = Xwi(0);
            double gy   = Xwi(1);
            double gz   = Xwi(2);
            double gz2  = gz*gz;
            gz2         = 1.0 / std::max(Config::homogTh(),gz2);
            double fx   = cam->getFx();
            double fy   = cam->getFy();
            double dx   = p_err(0);
            double dy   = p_err(1);
            double fxdx = fx*dx;
            double fydy = fy*dy;
            // estimate Jacobian wrt KF pose
            Vector6d Jij_Tiw = Vector6d::Zero();
            Jij_Tiw << + gz2 * fxdx * gz,
                       + gz2 * fydy * gz,
                       - gz2 * ( fxdx*gx + fydy*gy ),
                       - gz2 * ( fxdx*gx*gy + fydy*gy*gy + fydy*gz*gz ),
                       + gz2 * ( fxdx*gx*gx + fxdx*gz*gz + fydy*gx*gy ),
                       + gz2 * ( fydy*gx*gz - fxdx*gy*gz );
            Jij_Tiw = Jij_Tiw / std::max(Config::homogTh(),p_err_norm);
            // estimate Jacobian wrt LM
            Vector3d Jij_Xwj = Vector3d::Zero();
            Jij_Xwj << + gz2 * fxdx * gz,
                       + gz2 * fydy * gz,
                       - gz2 * ( fxdx*gx + fydy*gy );
            Jij_Xwj = Jij_Xwj.transpose() * Tiw.block(0,0,3,3) / std::max(Config::homogTh(),p_err_norm);
            // if employing robust cost function
            double w = 1.0;
            double s2 = map_points[lm_idx_map]->sigma_list[lm_idx_obs];
            w = 1.0 / ( 1.0 + p_err_norm * p_err_norm * s2 );
            // update hessian, gradient, and error
            VectorXd g_aux = VectorXd::Zero(N);
            MatrixXd Haux  = MatrixXd::Zero(3,6);
            int idx = 6 * kf_idx_loc;
            int jdx = 6*Nkf + 3*lm_idx_loc;
            if( kf_idx_loc == -1 )
            {
                g.block(jdx,0,3,1) += Jij_Xwj * p_err_norm * w ;
                err += p_err_norm * p_err_norm * w ;
                H.block(jdx,jdx,3,3) += Jij_Xwj * Jij_Xwj.transpose() * w ;
            }
            else
            {
                g.block(idx,0,6,1) += Jij_Tiw * p_err_norm * w;
                g.block(jdx,0,3,1) += Jij_Xwj * p_err_norm * w;
                err += p_err_norm * p_err_norm * w;
                Haux = Jij_Xwj * Jij_Tiw.transpose() * w ;
                H.block(idx,idx,6,6) += Jij_Tiw * Jij_Tiw.transpose() * w;
                H.block(jdx,idx,3,6) += Haux;
                H.block(idx,jdx,6,3) += Haux.transpose();
                H.block(jdx,jdx,3,3) += Jij_Xwj * Jij_Xwj.transpose() * w;
            }
        }
    }
    // line segment observations
    int Nls = 0, Nls_obs = 0;
    if( ls_obs_list.size() != 0 )
        Nls = ls_obs_list.back()(1)+1;
    for( vector<Vector6i>::iterator ls_it = ls_obs_list.begin(); ls_it != ls_obs_list.end(); ls_it++ )
    {
        int lm_idx_map = (*ls_it)(0);
        int lm_idx_loc = (*ls_it)(1);
        int lm_idx_obs = (*ls_it)(2);
        int kf_idx_map = (*ls_it)(3);
        int kf_idx_loc = (*ls_it)(4);
        if( map_lines[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
        {
            // grab 3D LM (Pwj and Qwj)
            Vector3d Pwj   = map_lines[lm_idx_map]->line3D.head(3);
            Vector3d Qwj   = map_lines[lm_idx_map]->line3D.tail(3);
            // grab 6DoF KF (Tiw)
            Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
            // projection error
            Tiw = inverse_se3( Tiw );
            Vector3d Pwi   = Tiw.block(0,0,3,3) * Pwj + Tiw.block(0,3,3,1);
            Vector3d Qwi   = Tiw.block(0,0,3,3) * Qwj + Tiw.block(0,3,3,1);
            Vector2d p_prj = cam->projection( Pwi );
            Vector2d q_prj = cam->projection( Qwi );
            Vector3d l_obs = map_lines[lm_idx_map]->obs_list[lm_idx_obs];
            Vector2d l_err;
            l_err(0) = l_obs(0) * p_prj(0) + l_obs(1) * p_prj(1) + l_obs(2);
            l_err(1) = l_obs(0) * q_prj(0) + l_obs(1) * q_prj(1) + l_obs(2);
            double l_err_norm = l_err.norm();
            // start point
            double gx   = Pwi(0);
            double gy   = Pwi(1);
            double gz   = Pwi(2);
            double gz2  = gz*gz;
            gz2         = 1.0 / std::max(Config::homogTh(),gz2);
            double fx   = cam->getFx();
            double fy   = cam->getFy();
            double lx   = l_err(0);
            double ly   = l_err(1);
            double fxlx = fx*lx;
            double fyly = fy*ly;
            // - jac. wrt. KF pose
            Vector6d Jij_Piw = Vector6d::Zero();
            Jij_Piw << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy ),
                       - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                       + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                       + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
            // - jac. wrt. LM
            Vector3d Jij_Pwj = Vector3d::Zero();
            Jij_Pwj << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy );
            Jij_Pwj = Jij_Pwj.transpose() * Tiw.block(0,0,3,3) * l_err(0) / std::max(Config::homogTh(),l_err_norm);
            // end point
            gx   = Qwi(0);
            gy   = Qwi(1);
            gz   = Qwi(2);
            gz2  = gz*gz;
            gz2         = 1.0 / std::max(Config::homogTh(),gz2);
            // - jac. wrt. KF pose
            Vector6d Jij_Qiw = Vector6d::Zero();
            Jij_Qiw << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy ),
                       - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                       + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                       + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
            // - jac. wrt. LM
            Vector3d Jij_Qwj = Vector3d::Zero();
            Jij_Qwj << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy );
            Jij_Qwj = Jij_Qwj.transpose() * Tiw.block(0,0,3,3) * l_err(1) / std::max(Config::homogTh(),l_err_norm);
            // estimate Jacobian wrt KF pose
            Vector6d Jij_Tiw = Vector6d::Zero();
            Jij_Tiw = ( Jij_Piw * l_err(0) + Jij_Qiw * l_err(1) ) / std::max(Config::homogTh(),l_err_norm);
            // estimate Jacobian wrt LM
            Vector6d Jij_Lwj = Vector6d::Zero();
            Jij_Lwj.head(3) = Jij_Pwj;
            Jij_Lwj.tail(3) = Jij_Qwj;
            // if employing robust cost function
            double s2 = map_lines[lm_idx_map]->sigma_list[lm_idx_obs];
            double w = 1.0 / ( 1.0 + l_err_norm * l_err_norm * s2 );
            // update hessian, gradient, and error
            VectorXd g_aux = VectorXd::Zero(N);
            MatrixXd Haux  = MatrixXd::Zero(3,6);
            int idx = 6 * kf_idx_loc;
            int jdx = 6*Nkf + 3*Npt + 6*lm_idx_loc;
            if( kf_idx_loc == -1 )
            {
                g.block(jdx,0,6,1) += Jij_Lwj * l_err_norm * w;
                err += l_err_norm * l_err_norm * w;
                H.block(jdx,jdx,6,6) += Jij_Lwj * Jij_Lwj.transpose() * w;
            }
            else
            {
                g.block(idx,0,6,1) += Jij_Tiw * l_err_norm * w;
                g.block(jdx,0,6,1) += Jij_Lwj * l_err_norm * w;
                err += l_err_norm * l_err_norm * w;
                Haux = Jij_Lwj * Jij_Tiw.transpose() * w;
                H.block(idx,idx,6,6) += Jij_Tiw * Jij_Tiw.transpose() * w;
                H.block(jdx,idx,6,6) += Haux;
                H.block(idx,jdx,6,6) += Haux.transpose();
                H.block(jdx,jdx,6,6) += Jij_Lwj * Jij_Lwj.transpose() * w;
            }
        }
    }
    err /= (Npt_obs+Nls_obs);
    // initial guess of lambda
    int Hmax = 0.0;
    for( int i = 0; i < N; i++)
    {
        if( H(i,i) > Hmax || H(i,i) < -Hmax )
            Hmax = fabs( H(i,i) );
    }
    lambda *= Hmax;
    // solve the first iteration
    for(int i = 0; i < N; i++)
        H(i,i) += lambda * H(i,i) ;
    H_ = H.sparseView();
    SimplicialLDLT< SparseMatrix<double> > solver1(H_);
    DX = solver1.solve( g );
    // update KFs
    for( int i = 0; i < Nkf; i++)
    {
        Matrix4d Tprev = expmap_se3( X.block(6*i,0,6,1) );
        Matrix4d Tcurr = Tprev * inverse_se3( expmap_se3( DX.block(6*i,0,6,1) ) );
        X.block(6*i,0,6,1) = logmap_se3( Tcurr );
    }
    // update LMs
    for( int i = 6*Nkf; i < N; i++)
        X(i) += DX(i);
    // update error
    err_prev = err;

    // LM iterations
    //---------------------------------------------------------------------------------------------
    int iters;
    for( iters = 1; iters < max_iters; iters++)
    {
        // estimate hessian and gradient (reset)
        DX = VectorXd::Zero(N);
        g  = VectorXd::Zero(N);
        H  = MatrixXd::Zero(N,N);
        err = 0.0;        
        // - point observations
        for( vector<Vector6i>::iterator pt_it = pt_obs_list.begin(); pt_it != pt_obs_list.end(); pt_it++ )
        {
            int lm_idx_map = (*pt_it)(0);
            int lm_idx_loc = (*pt_it)(1);
            int lm_idx_obs = (*pt_it)(2);
            int kf_idx_map = (*pt_it)(3);
            int kf_idx_loc = (*pt_it)(4);
            if( map_points[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
            {
                // grab 3D LM (Xwj)
                Vector3d Xwj = X.block(6*Nkf+3*lm_idx_loc,0,3,1);
                // grab 6DoF KF (Tiw)
                Matrix4d Tiw;
                if( kf_idx_loc != -1 )
                    Tiw = expmap_se3( X.block( 6*kf_idx_loc,0,6,1 ) );
                else
                    Tiw = map_keyframes[kf_idx_map]->T_kf_w;
                // projection error
                Tiw = inverse_se3( Tiw );
                Vector3d Xwi   = Tiw.block(0,0,3,3) * Xwj + Tiw.block(0,3,3,1);
                Vector2d p_prj = cam->projection( Xwi );
                Vector2d p_obs = map_points[lm_idx_map]->obs_list[lm_idx_obs];
                Vector2d p_err    = p_obs - p_prj;
                double p_err_norm = p_err.norm();
                // useful variables
                double gx   = Xwi(0);
                double gy   = Xwi(1);
                double gz   = Xwi(2);
                double gz2  = gz*gz;
                gz2         = 1.0 / std::max(Config::homogTh(),gz2);
                double fx   = cam->getFx();
                double fy   = cam->getFy();
                double dx   = p_err(0);
                double dy   = p_err(1);
                double fxdx = fx*dx;
                double fydy = fy*dy;
                // estimate Jacobian wrt KF pose
                Vector6d Jij_Tiw = Vector6d::Zero();
                Jij_Tiw << + gz2 * fxdx * gz,
                           + gz2 * fydy * gz,
                           - gz2 * ( fxdx*gx + fydy*gy ),
                           - gz2 * ( fxdx*gx*gy + fydy*gy*gy + fydy*gz*gz ),
                           + gz2 * ( fxdx*gx*gx + fxdx*gz*gz + fydy*gx*gy ),
                           + gz2 * ( fydy*gx*gz - fxdx*gy*gz );
                Jij_Tiw = Jij_Tiw / std::max(Config::homogTh(),p_err_norm);
                // estimate Jacobian wrt LM
                Vector3d Jij_Xwj = Vector3d::Zero();
                Jij_Xwj << + gz2 * fxdx * gz,
                           + gz2 * fydy * gz,
                           - gz2 * ( fxdx*gx + fydy*gy );
                Jij_Xwj = Jij_Xwj.transpose() * Tiw.block(0,0,3,3) / std::max(Config::homogTh(),p_err_norm);
                // if employing robust cost function
                double s2 = map_points[lm_idx_map]->sigma_list[lm_idx_obs];
                double w = 1.0 / ( 1.0 + p_err_norm * p_err_norm * s2 );
                // update hessian, gradient, and error
                VectorXd g_aux = VectorXd::Zero(N);
                MatrixXd Haux  = MatrixXd::Zero(3,6);
                int idx = 6 * kf_idx_loc;
                int jdx = 6*Nkf + 3*lm_idx_loc;
                if( kf_idx_loc == -1 )
                {
                    g.block(jdx,0,3,1) += Jij_Xwj * p_err_norm * w;
                    err += p_err_norm * p_err_norm * w;
                    H.block(jdx,jdx,3,3) += Jij_Xwj * Jij_Xwj.transpose() * w;
                }
                else
                {
                    g.block(idx,0,6,1) += Jij_Tiw * p_err_norm * w;
                    g.block(jdx,0,3,1) += Jij_Xwj * p_err_norm * w;
                    err += p_err_norm * p_err_norm * w;
                    Haux = Jij_Xwj * Jij_Tiw.transpose() * w ;
                    H.block(idx,idx,6,6) += Jij_Tiw * Jij_Tiw.transpose() * w ;
                    H.block(jdx,idx,3,6) += Haux;
                    H.block(idx,jdx,6,3) += Haux.transpose();
                    H.block(jdx,jdx,3,3) += Jij_Xwj * Jij_Xwj.transpose() * w ;
                }
            }
        }
        // - line segment observations
        for( vector<Vector6i>::iterator ls_it = ls_obs_list.begin(); ls_it != ls_obs_list.end(); ls_it++ )
        {
            int lm_idx_map = (*ls_it)(0);
            int lm_idx_loc = (*ls_it)(1);
            int lm_idx_obs = (*ls_it)(2);
            int kf_idx_map = (*ls_it)(3);
            int kf_idx_loc = (*ls_it)(4);
            if( map_lines[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
            {
                // grab 3D LM (Pwj and Qwj)
                Vector3d Pwj   = map_lines[lm_idx_map]->line3D.head(3);
                Vector3d Qwj   = map_lines[lm_idx_map]->line3D.tail(3);
                // grab 6DoF KF (Tiw)
                Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
                // projection error
                Tiw = inverse_se3( Tiw );
                Vector3d Pwi   = Tiw.block(0,0,3,3) * Pwj + Tiw.block(0,3,3,1);
                Vector3d Qwi   = Tiw.block(0,0,3,3) * Qwj + Tiw.block(0,3,3,1);
                Vector2d p_prj = cam->projection( Pwi );
                Vector2d q_prj = cam->projection( Qwi );
                Vector3d l_obs = map_lines[lm_idx_map]->obs_list[lm_idx_obs];
                Vector2d l_err;
                l_err(0) = l_obs(0) * p_prj(0) + l_obs(1) * p_prj(1) + l_obs(2);
                l_err(1) = l_obs(0) * q_prj(0) + l_obs(1) * q_prj(1) + l_obs(2);
                double l_err_norm = l_err.norm();
                // start point
                double gx   = Pwi(0);
                double gy   = Pwi(1);
                double gz   = Pwi(2);
                double gz2  = gz*gz;
                gz2         = 1.0 / std::max(0.0000001,gz2);
                double fx   = cam->getFx();
                double fy   = cam->getFy();
                double lx   = l_err(0);
                double ly   = l_err(1);
                double fxlx = fx*lx;
                double fyly = fy*ly;
                // - jac. wrt. KF pose
                Vector6d Jij_Piw = Vector6d::Zero();
                Jij_Piw << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy ),
                           - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                           + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                           + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
                // - jac. wrt. LM
                Vector3d Jij_Pwj = Vector3d::Zero();
                Jij_Pwj << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy );
                Jij_Pwj = Jij_Pwj.transpose() * Tiw.block(0,0,3,3) * l_err(0) / std::max(0.0000001,l_err_norm);
                // end point
                gx   = Qwi(0);
                gy   = Qwi(1);
                gz   = Qwi(2);
                gz2  = gz*gz;
                gz2         = 1.0 / std::max(0.0000001,gz2);
                // - jac. wrt. KF pose
                Vector6d Jij_Qiw = Vector6d::Zero();
                Jij_Qiw << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy ),
                           - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                           + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                           + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
                // - jac. wrt. LM
                Vector3d Jij_Qwj = Vector3d::Zero();
                Jij_Qwj << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy );
                Jij_Qwj = Jij_Qwj.transpose() * Tiw.block(0,0,3,3) * l_err(1) / std::max(0.0000001,l_err_norm);
                // estimate Jacobian wrt KF pose
                Vector6d Jij_Tiw = Vector6d::Zero();
                Jij_Tiw = ( Jij_Piw * l_err(0) + Jij_Qiw * l_err(1) ) / std::max(0.0000001,l_err_norm);
                // estimate Jacobian wrt LM
                Vector6d Jij_Lwj = Vector6d::Zero();
                Jij_Lwj.head(3) = Jij_Pwj;
                Jij_Lwj.tail(3) = Jij_Qwj;
                // if employing robust cost function
                double s2 = map_lines[lm_idx_map]->sigma_list[lm_idx_obs];
                double w = 1.0 / ( 1.0 + l_err_norm * l_err_norm * s2 );
                // update hessian, gradient, and error
                VectorXd g_aux = VectorXd::Zero(N);
                MatrixXd Haux  = MatrixXd::Zero(3,6);
                int idx = 6 * kf_idx_loc;
                int jdx = 6*Nkf + 3*Npt + 6*lm_idx_loc;
                if( kf_idx_loc == -1 )
                {
                    g.block(jdx,0,6,1) += Jij_Lwj * l_err_norm * w;
                    err += l_err_norm * l_err_norm * w;
                    H.block(jdx,jdx,6,6) += Jij_Lwj * Jij_Lwj.transpose() * w;
                }
                else
                {
                    g.block(idx,0,6,1) += Jij_Tiw * l_err_norm * w;
                    g.block(jdx,0,6,1) += Jij_Lwj * l_err_norm * w;
                    err += l_err_norm * l_err_norm * w;
                    Haux = Jij_Lwj * Jij_Tiw.transpose() * w;
                    H.block(idx,idx,6,6) += Jij_Tiw * Jij_Tiw.transpose() * w;
                    H.block(jdx,idx,6,6) += Haux;
                    H.block(idx,jdx,6,6) += Haux.transpose();
                    H.block(jdx,jdx,6,6) += Jij_Lwj * Jij_Lwj.transpose() * w;
                }
            }
        }
        err /= (Npt+Nls);
        // if the difference is very small stop
        if( abs(err-err_prev) < numeric_limits<double>::epsilon() || err < numeric_limits<double>::epsilon() )
            break;
        // add lambda to hessian
        for(int i = 0; i < N; i++)
            H(i,i) += lambda * H(i,i) ;
        // solve iteration
        H_ = H.sparseView();
        SimplicialLDLT< SparseMatrix<double> > solver1(H_);
        DX = solver1.solve( g );
        // update lambda
        if( err > err_prev ){
            lambda /= lambda_k;
        }
        else
        {
            lambda *= lambda_k;
            // update KFs
            for( int i = 0; i < Nkf; i++)
            {
                Matrix4d Tprev = expmap_se3( X.block(6*i,0,6,1) );
                Matrix4d Tcurr = Tprev * inverse_se3( expmap_se3( DX.block(6*i,0,6,1) ) );
                X.block(6*i,0,6,1) = logmap_se3( Tcurr );
            }
            // update LMs
            for( int i = 6*Nkf; i < N; i++)
                X(i) += DX(i);
        }
        // if the parameter change is small stop
        if( DX.norm() < numeric_limits<double>::epsilon() )
            break;
        // update previous values
        err_prev = err;

    }
    cout << endl << "LBA Iterations:    " << iters << endl << endl;

    // Update KFs and LMs
    //---------------------------------------------------------------------------------------------
    // update KFs
    for( int i = 0; i < Nkf; i++)
    {
        Matrix4d Test = expmap_se3( X.block( 6*i,0,6,1 ) );
        map_keyframes[ kf_list[i] ]->T_kf_w = Test;
    }
    // update point LMs
    for( int i = 0; i < Npt; i++)
    {
        map_points[ pt_list[i] ]->point3D(0) = X(6*Nkf+3*i);
        map_points[ pt_list[i] ]->point3D(1) = X(6*Nkf+3*i+1);
        map_points[ pt_list[i] ]->point3D(2) = X(6*Nkf+3*i+2);
    }
    // update line segment LMs
    for( int i = 0; i < Nls; i++)
    {
        map_lines[ ls_list[i] ]->line3D(0) = X(6*Nkf+3*Npt+6*i);
        map_lines[ ls_list[i] ]->line3D(1) = X(6*Nkf+3*Npt+6*i+1);
        map_lines[ ls_list[i] ]->line3D(2) = X(6*Nkf+3*Npt+6*i+2);
        map_lines[ ls_list[i] ]->line3D(3) = X(6*Nkf+3*Npt+6*i+3);
        map_lines[ ls_list[i] ]->line3D(4) = X(6*Nkf+3*Npt+6*i+4);
        map_lines[ ls_list[i] ]->line3D(5) = X(6*Nkf+3*Npt+6*i+5);
    }

    // Remove bad observations
    //---------------------------------------------------------------------------------------------
    for( vector<Vector6i>::iterator pt_it = pt_obs_list.end(); pt_it != pt_obs_list.begin(); pt_it-- )
    {
        if( (*pt_it)(5) == -1 )
        {
            int lm_idx_map = (*pt_it)(0);
            int lm_idx_obs = (*pt_it)(2);
            if( map_points[lm_idx_map] != NULL )
            {
                int kf_obs = map_points[lm_idx_map]->kf_obs_list[lm_idx_obs];
                // remove observations from map_points
                if( map_points[lm_idx_map]->obs_list.size() > 1 )
                {
                    // if it is the first observation, update it from map_points_kf_idx
                    if( lm_idx_obs == 0 )
                    {
                        // delete observation from map_points_kf_idx
                        for( auto it = map_points_kf_idx.at(kf_obs).begin(); it != map_points_kf_idx.at(kf_obs).end(); it++)
                        {
                            if( (*it) == lm_idx_map )
                            {
                                int new_kf_base = map_points[(*it)]->kf_obs_list[1];
                                map_points_kf_idx.at(new_kf_base).push_back( (*it) );
                                break;
                            }
                        }
                    }
                    // remove observations from map points
                    map_points[lm_idx_map]->desc_list.erase( map_points[lm_idx_map]->desc_list.begin() + lm_idx_obs );
                    map_points[lm_idx_map]->obs_list.erase( map_points[lm_idx_map]->obs_list.begin() + lm_idx_obs );
                    map_points[lm_idx_map]->dir_list.erase( map_points[lm_idx_map]->dir_list.begin() + lm_idx_obs );
                    map_points[lm_idx_map]->kf_obs_list.erase( map_points[lm_idx_map]->kf_obs_list.begin() + lm_idx_obs );
                    // remove idx from KeyFrame stereo points
                    for(vector<PointFeature*>::iterator st_pt = map_keyframes[kf_obs]->stereo_frame->stereo_pt.begin();
                        st_pt != map_keyframes[kf_obs]->stereo_frame->stereo_pt.end(); st_pt++ )
                    {
                        if( (*st_pt)->idx == lm_idx_map )
                        {
                            (*st_pt)->idx = -1;
                            st_pt = map_keyframes[kf_obs]->stereo_frame->stereo_pt.end()-1;
                        }
                    }
                    // update main descriptor and direction
                    map_points[lm_idx_map]->updateAverageDescDir();
                    // update graphs
                    for( int i = 0; i < map_points[lm_idx_map]->kf_obs_list.size(); i++ )
                    {
                        int idx = map_points[lm_idx_map]->kf_obs_list[i];
                        if( kf_obs != idx )
                        {
                            full_graph[kf_obs][idx]--;
                            full_graph[idx][kf_obs]--;
                        }
                    }
                }
                else
                    map_points[lm_idx_map]->inlier = false;
            }
        }
    }

    for( vector<Vector6i>::iterator ls_it = ls_obs_list.end(); ls_it != ls_obs_list.begin(); ls_it-- )
    {
        if( (*ls_it)(5) == -1 )
        {
            int lm_idx_map = (*ls_it)(0);
            int lm_idx_obs = (*ls_it)(2);
            if( map_lines[lm_idx_map] != NULL )
            {
                int kf_obs = map_lines[lm_idx_map]->kf_obs_list[lm_idx_obs];
                // remove observations from map_points
                if( map_lines[lm_idx_map]->obs_list.size() > 1 )
                {
                    // if it is the first observation, update it from map_points_kf_idx
                    if( lm_idx_obs == 0 )
                    {
                        // delete observation from map_points_kf_idx
                        for( auto it = map_lines_kf_idx.at(kf_obs).begin(); it != map_lines_kf_idx.at(kf_obs).end(); it++)
                        {
                            if( (*it) == lm_idx_map )
                            {
                                int new_kf_base = map_lines[(*it)]->kf_obs_list[1];
                                map_lines_kf_idx.at(new_kf_base).push_back( (*it) );
                                break;
                            }
                        }
                    }
                    // remove observations
                    map_lines[lm_idx_map]->desc_list.erase( map_lines[lm_idx_map]->desc_list.begin() + lm_idx_obs );
                    map_lines[lm_idx_map]->obs_list.erase( map_lines[lm_idx_map]->obs_list.begin() + lm_idx_obs );
                    map_lines[lm_idx_map]->pts_list.erase( map_lines[lm_idx_map]->pts_list.begin() + lm_idx_obs );
                    map_lines[lm_idx_map]->dir_list.erase( map_lines[lm_idx_map]->dir_list.begin() + lm_idx_obs );
                    map_lines[lm_idx_map]->kf_obs_list.erase( map_lines[lm_idx_map]->kf_obs_list.begin() + lm_idx_obs );
                    // remove idx from KeyFrame stereo lines
                    for(vector<LineFeature*>::iterator st_ls = map_keyframes[kf_obs]->stereo_frame->stereo_ls.begin();
                        st_ls != map_keyframes[kf_obs]->stereo_frame->stereo_ls.end(); st_ls++ )
                    {
                        if( (*st_ls)->idx == lm_idx_map )
                        {
                            (*st_ls)->idx = -1;
                            st_ls = map_keyframes[kf_obs]->stereo_frame->stereo_ls.end()-1;
                        }
                    }
                    // update main descriptor and direction
                    map_lines[lm_idx_map]->updateAverageDescDir();
                    // update graphs
                    for( int i = 0; i < map_lines[lm_idx_map]->kf_obs_list.size(); i++ )
                    {
                        int idx = map_lines[lm_idx_map]->kf_obs_list[i];
                        if( kf_obs != idx )
                        {
                            full_graph[kf_obs][idx]--;
                            full_graph[idx][kf_obs]--;
                        }
                    }
                }
                else
                    map_lines[lm_idx_map]->inlier = false;
            }

        }
    }

}

// -----------------------------------------------------------------------------------------------------------------------------
// Global Bundle Adjustment functions
// -----------------------------------------------------------------------------------------------------------------------------

void MapHandler::globalBundleAdjustment()
{

    vector<double> X_aux;

    // create list of keyframes
    vector<int> kf_list;
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++)
    {
        if( (*kf_it)!= NULL )
        {
            if( (*kf_it)->kf_idx != 0 )
            {
                Vector6d pose_aux = (*kf_it)->x_kf_w;
                for(int i = 0; i < 6; i++)
                    X_aux.push_back( pose_aux(i) );
                kf_list.push_back( (*kf_it)->kf_idx );
            }
        }
    }

    // create list of point landmarks
    vector<Vector6i> pt_obs_list;
    vector<int> pt_list;
    int lm_local_idx = 0;
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++)
    {
        if( (*pt_it)!= NULL )
        {
            Vector3d point_aux = (*pt_it)->point3D;
            for(int i = 0; i < 3; i++)
                X_aux.push_back( point_aux(i) );
            // gather all observations
            for( int i = 0; i < (*pt_it)->obs_list.size(); i++)
            {
                Vector6i obs_aux;
                obs_aux(0) = (*pt_it)->idx; // LM idx
                obs_aux(1) = lm_local_idx;  // LM local idx
                obs_aux(2) = i;             // LM obs idx
                int kf_obs_list_ = (*pt_it)->kf_obs_list[i];
                obs_aux(3) = kf_obs_list_;  // KF idx
                obs_aux(4) = -1;            // KF local idx (-1 if not local)
                obs_aux(5) = 1;             // 1 if the observation is an inlier
                for( int j = 0; j < kf_list.size(); j++ )
                {
                    if( kf_list[j] == kf_obs_list_ )
                    {
                        obs_aux(4) = j;
                        break;
                    }
                }
                pt_obs_list.push_back( obs_aux );
            }
            lm_local_idx++;
            // pt_list
            pt_list.push_back( (*pt_it)->idx );
        }
    }

    // create list of line segment landmarks
    vector<Vector6i> ls_obs_list;
    vector<int> ls_list;
    lm_local_idx = 0;
    for( vector<MapLine*>::iterator ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++)
    {
        if( (*ls_it)!= NULL )
        {
            Vector6d line_aux = (*ls_it)->line3D;
            for(int i = 0; i < 6; i++)
                X_aux.push_back( line_aux(i) );
            // gather all observations
            for( int i = 0; i < (*ls_it)->obs_list.size(); i++)
            {
                Vector6i obs_aux;
                obs_aux(0) = (*ls_it)->idx; // LM idx
                obs_aux(1) = lm_local_idx;  // LM local idx
                obs_aux(2) = i;             // LM obs idx
                int kf_obs_list_ = (*ls_it)->kf_obs_list[i];
                obs_aux(3) = kf_obs_list_;  // KF idx
                obs_aux(4) = -1;            // KF local idx (-1 if not local)
                obs_aux(5) = 1;             // 1 if the observation is an inlier
                for( int j = 0; j < kf_list.size(); j++ )
                {
                    if( kf_list[j] == kf_obs_list_ )
                    {
                        obs_aux(4) = j;
                        break;
                    }
                }
                ls_obs_list.push_back( obs_aux );
            }
            lm_local_idx++;
            // ls_list
            ls_list.push_back( (*ls_it)->idx );
        }
    }

    // Levenberg-Marquardt optimization
    levMarquardtOptimizationGBA(X_aux,kf_list,pt_list,ls_list,pt_obs_list,ls_obs_list);

    // -------------------------------------------------------------------------------------------------------------------

    // Recent map LMs culling (implement filters for line segments, which seems to be unaccurate)
    // Recent KFs culling
}

void MapHandler::levMarquardtOptimizationGBA( vector<double> X_aux, vector<int> kf_list, vector<int> pt_list, vector<int> ls_list, vector<Vector6i> pt_obs_list, vector<Vector6i> ls_obs_list  )
{

    // create Levenberg-Marquardt variables
    int    Nkf = kf_list.size();
    int      N = X_aux.size();
    VectorXd DX, X(N), gdense(N);
    SparseVector<double> g(N);
    SparseMatrix<double> H(N,N);
    for(int i = 0; i < N; i++)
        X.coeffRef(i) = X_aux[i];
    H.reserve( VectorXi::Constant(N,5000) );    //[TODO: change the allocation]

    // create Levenberg-Marquardt parameters
    double err, err_prev = 999999999.9;
    double lambda = Config::lambdaLbaLM(), lambda_k = Config::lambdaLbaK(); // TODO: change variable name in case they are equal with EBA
    int    max_iters = Config::maxItersLba();

    // estimate H and g to precalculate lambda
    //---------------------------------------------------------------------------------------------
    // point observations
    int Npt = 0, Npt_obs = 0;
    if( pt_obs_list.size() != 0 )
        Npt = pt_obs_list.back()(1)+1;
    for( vector<Vector6i>::iterator pt_it = pt_obs_list.begin(); pt_it != pt_obs_list.end(); pt_it++ )
    {
        int lm_idx_map = (*pt_it)(0);
        int lm_idx_loc = (*pt_it)(1);
        int lm_idx_obs = (*pt_it)(2);
        int kf_idx_map = (*pt_it)(3);
        int kf_idx_loc = (*pt_it)(4);
        if( map_points[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
        {
            // grab 3D LM (Xwj)
            Vector3d Xwj   = map_points[lm_idx_map]->point3D;
            // grab 6DoF KF (Tiw)
            Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
            // projection error
            Tiw = inverse_se3( Tiw );
            Vector3d Xwi   = Tiw.block(0,0,3,3) * Xwj + Tiw.block(0,3,3,1);
            Vector2d p_prj = cam->projection( Xwi );
            Vector2d p_obs = map_points[lm_idx_map]->obs_list[lm_idx_obs];
            Vector2d p_err    = p_obs - p_prj;
            double p_err_norm = p_err.norm();
            double gx   = Xwi(0);
            double gy   = Xwi(1);
            double gz   = Xwi(2);
            double gz2  = gz*gz;
            gz2         = 1.0 / std::max(Config::homogTh(),gz2);
            double fx   = cam->getFx();
            double fy   = cam->getFy();
            double dx   = p_err(0);
            double dy   = p_err(1);
            double fxdx = fx*dx;
            double fydy = fy*dy;
            // estimate Jacobian wrt KF pose
            Vector6d Jij_Tiw = Vector6d::Zero();
            Jij_Tiw << + gz2 * fxdx * gz,
                       + gz2 * fydy * gz,
                       - gz2 * ( fxdx*gx + fydy*gy ),
                       - gz2 * ( fxdx*gx*gy + fydy*gy*gy + fydy*gz*gz ),
                       + gz2 * ( fxdx*gx*gx + fxdx*gz*gz + fydy*gx*gy ),
                       + gz2 * ( fydy*gx*gz - fxdx*gy*gz );
            Jij_Tiw = Jij_Tiw / std::max(Config::homogTh(),p_err_norm);
            // estimate Jacobian wrt LM
            Vector3d Jij_Xwj = Vector3d::Zero();
            Jij_Xwj << + gz2 * fxdx * gz,
                       + gz2 * fydy * gz,
                       - gz2 * ( fxdx*gx + fydy*gy );
            Jij_Xwj = Jij_Xwj.transpose() * Tiw.block(0,0,3,3) / std::max(Config::homogTh(),p_err_norm);
            // if employing robust cost function
            double s2 = map_points[lm_idx_map]->sigma_list[lm_idx_obs];
            double w = 1.0 / ( 1.0 + p_err_norm * p_err_norm * s2 );
            // update hessian, gradient, and error
            VectorXd g_aux = VectorXd::Zero(N);
            int idx = 6 * kf_idx_loc;
            int jdx = 6*Nkf + 3*lm_idx_loc;
            if( kf_idx_loc == -1 )
            {
                err += p_err_norm * p_err_norm * w;
                Vector3d gi;
                Matrix3d Hjj;
                Hjj = Jij_Xwj * Jij_Xwj.transpose() * w;
                gi  = Jij_Xwj * p_err_norm * w;
                for(int i = 0; i < 3; i++)
                {
                    g.coeffRef(jdx+i) += gi(i);
                    for(int j = 0; j < 3; j++)
                        H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                }
            }
            else
            {
                err += p_err_norm * p_err_norm * w;
                Vector3d gj;
                Vector6d gi;
                Matrix6d Hii;
                gi = Jij_Tiw * p_err_norm * w;
                gj = Jij_Xwj * p_err_norm * w;
                Hii = Jij_Tiw * Jij_Tiw.transpose() * w;
                for(int i = 0; i < 6; i++)
                {
                    g.coeffRef(i+idx) += gi(i);
                    for(int j = 0; j < 6; j++)
                        H.coeffRef(i+idx,j+idx) += Hii(i,j);
                }
                Matrix3d Hjj;
                Hjj = Jij_Xwj * Jij_Xwj.transpose() * w;
                for(int i = 0; i < 3; i++)
                {
                    g.coeffRef(i+jdx) += gj(i);
                    for(int j = 0; j < 3; j++)
                        H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                }
                MatrixXd Hij  = MatrixXd::Zero(3,6);
                Hij = Jij_Xwj * Jij_Tiw.transpose() * w;
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 6; j++)
                    {
                        H.coeffRef(i+jdx,j+idx) += Hij(i,j);
                        H.coeffRef(j+idx,i+jdx) += Hij(i,j);
                    }
                }
            }
        }
    }
    // line segment observations
    int Nls = 0, Nls_obs = 0;
    if( ls_obs_list.size() != 0 )
        Nls = ls_obs_list.back()(1)+1;
    for( vector<Vector6i>::iterator ls_it = ls_obs_list.begin(); ls_it != ls_obs_list.end(); ls_it++ )
    {
        int lm_idx_map = (*ls_it)(0);
        int lm_idx_loc = (*ls_it)(1);
        int lm_idx_obs = (*ls_it)(2);
        int kf_idx_map = (*ls_it)(3);
        int kf_idx_loc = (*ls_it)(4);
        if( map_lines[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
        {
            // grab 3D LM (Pwj and Qwj)
            Vector3d Pwj   = map_lines[lm_idx_map]->line3D.head(3);
            Vector3d Qwj   = map_lines[lm_idx_map]->line3D.tail(3);
            // grab 6DoF KF (Tiw)
            Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
            // projection error
            Tiw = inverse_se3( Tiw );
            Vector3d Pwi   = Tiw.block(0,0,3,3) * Pwj + Tiw.block(0,3,3,1);
            Vector3d Qwi   = Tiw.block(0,0,3,3) * Qwj + Tiw.block(0,3,3,1);
            Vector2d p_prj = cam->projection( Pwi );
            Vector2d q_prj = cam->projection( Qwi );
            Vector3d l_obs = map_lines[lm_idx_map]->obs_list[lm_idx_obs];
            Vector2d l_err;
            l_err(0) = l_obs(0) * p_prj(0) + l_obs(1) * p_prj(1) + l_obs(2);
            l_err(1) = l_obs(0) * q_prj(0) + l_obs(1) * q_prj(1) + l_obs(2);
            double l_err_norm = l_err.norm();
            // start point
            double gx   = Pwi(0);
            double gy   = Pwi(1);
            double gz   = Pwi(2);
            double gz2  = gz*gz;
            gz2         = 1.0 / std::max(Config::homogTh(),gz2);
            double fx   = cam->getFx();
            double fy   = cam->getFy();
            double lx   = l_err(0);
            double ly   = l_err(1);
            double fxlx = fx*lx;
            double fyly = fy*ly;
            // - jac. wrt. KF pose
            Vector6d Jij_Piw = Vector6d::Zero();
            Jij_Piw << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy ),
                       - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                       + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                       + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
            // - jac. wrt. LM
            Vector3d Jij_Pwj = Vector3d::Zero();
            Jij_Pwj << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy );
            Jij_Pwj = Jij_Pwj.transpose() * Tiw.block(0,0,3,3) * l_err(0) / std::max(Config::homogTh(),l_err_norm);
            // end point
            gx   = Qwi(0);
            gy   = Qwi(1);
            gz   = Qwi(2);
            gz2  = gz*gz;
            gz2  = 1.0 / std::max(Config::homogTh(),gz2);
            // - jac. wrt. KF pose
            Vector6d Jij_Qiw = Vector6d::Zero();
            Jij_Qiw << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy ),
                       - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                       + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                       + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
            // - jac. wrt. LM
            Vector3d Jij_Qwj = Vector3d::Zero();
            Jij_Qwj << + gz2 * fxlx * gz,
                       + gz2 * fyly * gz,
                       - gz2 * ( fxlx*gx + fyly*gy );
            Jij_Qwj = Jij_Qwj.transpose() * Tiw.block(0,0,3,3) * l_err(1) / std::max(Config::homogTh(),l_err_norm);
            // estimate Jacobian wrt KF pose
            Vector6d Jij_Tiw = Vector6d::Zero();
            Jij_Tiw = ( Jij_Piw * l_err(0) + Jij_Qiw * l_err(1) ) / std::max(Config::homogTh(),l_err_norm);
            // estimate Jacobian wrt LM
            Vector6d Jij_Lwj = Vector6d::Zero();
            Jij_Lwj.head(3) = Jij_Pwj;
            Jij_Lwj.tail(3) = Jij_Qwj;
            // if employing robust cost function
            double s2 = map_lines[lm_idx_map]->sigma_list[lm_idx_obs];
            double w = 1.0 / ( 1.0 + l_err_norm * l_err_norm * s2 );
            // update hessian, gradient, and error
            VectorXd gi = VectorXd::Zero(N), gj = VectorXd::Zero(N);
            int idx = 6 * kf_idx_loc;
            int jdx = 6*Nkf + 3*Npt + 6*lm_idx_loc;
            if( kf_idx_loc == -1 )
            {
                gj = Jij_Lwj * l_err_norm * w;
                err += l_err_norm * l_err_norm * w;
                Matrix6d Hjj;
                Hjj = Jij_Lwj * Jij_Lwj.transpose() * w;
                for(int i = 0; i < 6; i++)
                {
                    g.coeffRef(jdx+i) += gj(i);
                    for(int j = 0; j < 6; j++)
                        H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                }
            }
            else
            {
                gi = Jij_Tiw * l_err_norm * w;
                gj = Jij_Lwj * l_err_norm * w;
                err += l_err_norm * l_err_norm * w;
                Matrix6d Hii, Hjj;
                MatrixXd Hij   = MatrixXd::Zero(3,6);
                Hii = Jij_Tiw * Jij_Tiw.transpose() * w;
                Hjj = Jij_Lwj * Jij_Lwj.transpose() * w;
                Hij = Jij_Lwj * Jij_Tiw.transpose() * w;
                for(int i = 0; i < 6; i++)
                {
                    g.coeffRef(i+idx) += gi(i);
                    g.coeffRef(i+jdx) += gj(i);
                    for(int j = 0; j < 6; j++)
                    {
                        H.coeffRef(idx+i,idx+j) += Hii(i,j);
                        H.coeffRef(jdx+i,jdx+j) += Hjj(i,j);
                        H.coeffRef(idx+i,jdx+j) += Hij(i,j);
                        H.coeffRef(jdx+j,idx+i) += Hij(i,j);
                    }
                }
            }
        }
    }
    err /= (Npt_obs+Nls_obs);
    // initial guess of lambda
    int Hmax = 0.0;
    for( int i = 0; i < N; i++)
    {
        if( H.coeffRef(i,i) > Hmax || H.coeffRef(i,i) < -Hmax )
            Hmax = fabs( H.coeffRef(i,i) );
    }
    lambda *= Hmax;
    // solve the first iteration
    H.makeCompressed();
    for(int i = 0; i < N; i++)
    {
        H.coeffRef(i,i) += lambda * H.coeffRef(i,i);
        gdense(i) = g.coeffRef(i);
    }
    SimplicialLDLT< SparseMatrix<double> > solver1(H);
    DX = solver1.solve( gdense );
    // update KFs
    for( int i = 0; i < Nkf; i++)
    {
        Matrix4d Tprev = expmap_se3( X.block(6*i,0,6,1) );
        Matrix4d Tcurr = Tprev * inverse_se3( expmap_se3( DX.block(6*i,0,6,1) ) );       
        X.block(6*i,0,6,1) = logmap_se3( Tcurr );
    }
    // update LMs
    for( int i = 6*Nkf; i < N; i++)
        X(i) += DX(i);
    // update error
    err_prev = err;

    // LM iterations
    //---------------------------------------------------------------------------------------------
    int iters;
    for( iters = 1; iters < Config::maxItersLba(); iters++ )
    {
        // estimate hessian and gradient (reset)
        DX.setZero();
        g.setZero();
        H.setZero();
        H.reserve( VectorXi::Constant(N,5000) );
        err = 0.0;
        // - point observations
        for( vector<Vector6i>::iterator pt_it = pt_obs_list.begin(); pt_it != pt_obs_list.end(); pt_it++ )
        {
            int lm_idx_map = (*pt_it)(0);
            int lm_idx_loc = (*pt_it)(1);
            int lm_idx_obs = (*pt_it)(2);
            int kf_idx_map = (*pt_it)(3);
            int kf_idx_loc = (*pt_it)(4);
            if( map_points[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
            {
                // grab 3D LM (Xwj)
                Vector3d Xwj = X.block(6*Nkf+3*lm_idx_loc,0,3,1);
                // grab 6DoF KF (Tiw)
                Matrix4d Tiw;
                if( kf_idx_loc != -1 )
                    Tiw = expmap_se3( X.block( 6*kf_idx_loc,0,6,1 ) );
                else
                    Tiw = map_keyframes[kf_idx_map]->T_kf_w;
                // projection error
                Tiw = inverse_se3( Tiw );
                Vector3d Xwi   = Tiw.block(0,0,3,3) * Xwj + Tiw.block(0,3,3,1);
                Vector2d p_prj = cam->projection( Xwi );
                Vector2d p_obs = map_points[lm_idx_map]->obs_list[lm_idx_obs];
                Vector2d p_err    = p_obs - p_prj;
                double p_err_norm = p_err.norm();
                // useful variables
                double gx   = Xwi(0);
                double gy   = Xwi(1);
                double gz   = Xwi(2);
                double gz2  = gz*gz;
                gz2         = 1.0 / std::max(Config::homogTh(),gz2);
                double fx   = cam->getFx();
                double fy   = cam->getFy();
                double dx   = p_err(0);
                double dy   = p_err(1);
                double fxdx = fx*dx;
                double fydy = fy*dy;
                // estimate Jacobian wrt KF pose
                Vector6d Jij_Tiw = Vector6d::Zero();
                Jij_Tiw << + gz2 * fxdx * gz,
                           + gz2 * fydy * gz,
                           - gz2 * ( fxdx*gx + fydy*gy ),
                           - gz2 * ( fxdx*gx*gy + fydy*gy*gy + fydy*gz*gz ),
                           + gz2 * ( fxdx*gx*gx + fxdx*gz*gz + fydy*gx*gy ),
                           + gz2 * ( fydy*gx*gz - fxdx*gy*gz );
                Jij_Tiw = Jij_Tiw / std::max(Config::homogTh(),p_err_norm);
                // estimate Jacobian wrt LM
                Vector3d Jij_Xwj = Vector3d::Zero();
                Jij_Xwj << + gz2 * fxdx * gz,
                           + gz2 * fydy * gz,
                           - gz2 * ( fxdx*gx + fydy*gy );
                Jij_Xwj = Jij_Xwj.transpose() * Tiw.block(0,0,3,3) / std::max(Config::homogTh(),p_err_norm);
                // if employing robust cost function
                double s2 = map_points[lm_idx_map]->sigma_list[lm_idx_obs];
                double w = 1.0 / ( 1.0 + p_err_norm * p_err_norm * s2 );
                // update hessian, gradient, and error
                VectorXd g_aux = VectorXd::Zero(N);
                int idx = 6 * kf_idx_loc;
                int jdx = 6*Nkf + 3*lm_idx_loc;
                if( kf_idx_loc == -1 )
                {
                    err += p_err_norm * p_err_norm * w;
                    Vector3d gi;
                    Matrix3d Hjj;
                    Hjj = Jij_Xwj * Jij_Xwj.transpose() * w;
                    gi  = Jij_Xwj * p_err_norm * w;
                    for(int i = 0; i < 3; i++)
                    {
                        g.coeffRef(jdx+i) += gi(i);
                        for(int j = 0; j < 3; j++)
                            H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                    }
                }
                else
                {
                    err += p_err_norm * p_err_norm * w;
                    Vector3d gj;
                    Vector6d gi;
                    Matrix6d Hii;
                    gi = Jij_Tiw * p_err_norm * w;
                    gj = Jij_Xwj * p_err_norm * w;
                    Hii = Jij_Tiw * Jij_Tiw.transpose() * w;
                    for(int i = 0; i < 6; i++)
                    {
                        g.coeffRef(i+idx) += gi(i);
                        for(int j = 0; j < 6; j++)
                            H.coeffRef(i+idx,j+idx) += Hii(i,j);
                    }
                    Matrix3d Hjj;
                    Hjj = Jij_Xwj * Jij_Xwj.transpose() * w;
                    for(int i = 0; i < 3; i++)
                    {
                        g.coeffRef(i+jdx) += gj(i);
                        for(int j = 0; j < 3; j++)
                            H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                    }
                    MatrixXd Hij  = MatrixXd::Zero(3,6);
                    Hij = Jij_Xwj * Jij_Tiw.transpose() * w;
                    for(int i = 0; i < 3; i++)
                    {
                        for(int j = 0; j < 6; j++)
                        {
                            H.coeffRef(i+jdx,j+idx) += Hij(i,j);
                            H.coeffRef(j+idx,i+jdx) += Hij(i,j);
                        }
                    }
                }
            }
        }
        // - line segment observations
        for( vector<Vector6i>::iterator ls_it = ls_obs_list.begin(); ls_it != ls_obs_list.end(); ls_it++ )
        {
            int lm_idx_map = (*ls_it)(0);
            int lm_idx_loc = (*ls_it)(1);
            int lm_idx_obs = (*ls_it)(2);
            int kf_idx_map = (*ls_it)(3);
            int kf_idx_loc = (*ls_it)(4);
            if( map_lines[lm_idx_map] != NULL && map_keyframes[kf_idx_map] != NULL)
            {
                // grab 3D LM (Pwj and Qwj)
                Vector3d Pwj   = map_lines[lm_idx_map]->line3D.head(3);
                Vector3d Qwj   = map_lines[lm_idx_map]->line3D.tail(3);
                // grab 6DoF KF (Tiw)
                Matrix4d Tiw   = map_keyframes[kf_idx_map]->T_kf_w;
                // projection error
                Tiw = inverse_se3( Tiw );
                Vector3d Pwi   = Tiw.block(0,0,3,3) * Pwj + Tiw.block(0,3,3,1);
                Vector3d Qwi   = Tiw.block(0,0,3,3) * Qwj + Tiw.block(0,3,3,1);
                Vector2d p_prj = cam->projection( Pwi );
                Vector2d q_prj = cam->projection( Qwi );
                Vector3d l_obs = map_lines[lm_idx_map]->obs_list[lm_idx_obs];
                Vector2d l_err;
                l_err(0) = l_obs(0) * p_prj(0) + l_obs(1) * p_prj(1) + l_obs(2);
                l_err(1) = l_obs(0) * q_prj(0) + l_obs(1) * q_prj(1) + l_obs(2);
                double l_err_norm = l_err.norm();
                // start point
                double gx   = Pwi(0);
                double gy   = Pwi(1);
                double gz   = Pwi(2);
                double gz2  = gz*gz;
                gz2         = 1.0 / std::max(Config::homogTh(),gz2);
                double fx   = cam->getFx();
                double fy   = cam->getFy();
                double lx   = l_err(0);
                double ly   = l_err(1);
                double fxlx = fx*lx;
                double fyly = fy*ly;
                // - jac. wrt. KF pose
                Vector6d Jij_Piw = Vector6d::Zero();
                Jij_Piw << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy ),
                           - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                           + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                           + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
                // - jac. wrt. LM
                Vector3d Jij_Pwj = Vector3d::Zero();
                Jij_Pwj << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy );
                Jij_Pwj = Jij_Pwj.transpose() * Tiw.block(0,0,3,3) * l_err(0) / std::max(Config::homogTh(),l_err_norm);
                // end point
                gx   = Qwi(0);
                gy   = Qwi(1);
                gz   = Qwi(2);
                gz2  = gz*gz;
                gz2  = 1.0 / std::max(Config::homogTh(),gz2);
                // - jac. wrt. KF pose
                Vector6d Jij_Qiw = Vector6d::Zero();
                Jij_Qiw << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy ),
                           - gz2 * ( fxlx*gx*gy + fyly*gy*gy + fyly*gz*gz ),
                           + gz2 * ( fxlx*gx*gx + fxlx*gz*gz + fyly*gx*gy ),
                           + gz2 * ( fyly*gx*gz - fxlx*gy*gz );
                // - jac. wrt. LM
                Vector3d Jij_Qwj = Vector3d::Zero();
                Jij_Qwj << + gz2 * fxlx * gz,
                           + gz2 * fyly * gz,
                           - gz2 * ( fxlx*gx + fyly*gy );
                Jij_Qwj = Jij_Qwj.transpose() * Tiw.block(0,0,3,3) * l_err(1) / std::max(Config::homogTh(),l_err_norm);
                // estimate Jacobian wrt KF pose
                Vector6d Jij_Tiw = Vector6d::Zero();
                Jij_Tiw = ( Jij_Piw * l_err(0) + Jij_Qiw * l_err(1) ) / std::max(Config::homogTh(),l_err_norm);
                // estimate Jacobian wrt LM
                Vector6d Jij_Lwj = Vector6d::Zero();
                Jij_Lwj.head(3) = Jij_Pwj;
                Jij_Lwj.tail(3) = Jij_Qwj;
                // if employing robust cost function
                double s2 = map_lines[lm_idx_map]->sigma_list[lm_idx_obs];
                double w = 1.0 / ( 1.0 + l_err_norm * l_err_norm * s2 );
                // update hessian, gradient, and error
                VectorXd gi = VectorXd::Zero(N), gj = VectorXd::Zero(N);
                int idx = 6 * kf_idx_loc;
                int jdx = 6*Nkf + 3*Npt + 6*lm_idx_loc;
                if( kf_idx_loc == -1 )
                {
                    gj = Jij_Lwj * l_err_norm * w;
                    err += l_err_norm * l_err_norm * w;
                    Matrix6d Hjj;
                    Hjj = Jij_Lwj * Jij_Lwj.transpose() * w;
                    for(int i = 0; i < 6; i++)
                    {
                        g.coeffRef(jdx+i) += gj(i);
                        for(int j = 0; j < 6; j++)
                            H.coeffRef(i+jdx,j+jdx) += Hjj(i,j);
                    }
                }
                else
                {
                    gi = Jij_Tiw * l_err_norm * w;
                    gj = Jij_Lwj * l_err_norm * w;
                    err += l_err_norm * l_err_norm * w;
                    Matrix6d Hii, Hjj;
                    MatrixXd Hij   = MatrixXd::Zero(3,6);
                    Hii = Jij_Tiw * Jij_Tiw.transpose() * w;
                    Hjj = Jij_Lwj * Jij_Lwj.transpose() * w;
                    Hij = Jij_Lwj * Jij_Tiw.transpose() * w;
                    for(int i = 0; i < 6; i++)
                    {
                        g.coeffRef(i+idx) += gi(i);
                        g.coeffRef(i+jdx) += gj(i);
                        for(int j = 0; j < 6; j++)
                        {
                            H.coeffRef(idx+i,idx+j) += Hii(i,j);
                            H.coeffRef(jdx+i,jdx+j) += Hjj(i,j);
                            H.coeffRef(idx+i,jdx+j) += Hij(i,j);
                            H.coeffRef(jdx+j,idx+i) += Hij(j,i);
                        }
                    }
                }
            }
        }
        err /= (Npt_obs+Nls_obs);
        // if the difference is very small stop
        if( abs(err-err_prev) < numeric_limits<double>::epsilon() || err < numeric_limits<double>::epsilon() )
            break;
        // add lambda to diagonal
        for(int i = 0; i < N; i++)
            H.coeffRef(i,i) += lambda * H.coeffRef(i,i) ;
        // solve iteration
        H.makeCompressed();
        for( int i = 0; i < N; i++)
            gdense(i) = g.coeffRef(i);
        SimplicialLDLT< SparseMatrix<double> > solver1(H);
        DX = solver1.solve( gdense );
        // update lambda
        if( err > err_prev ){
            lambda /= lambda_k;
        }
        else
        {
            lambda *= lambda_k;
            // update KFs
            for( int i = 0; i < Nkf; i++)
            {
                Matrix4d Tprev = expmap_se3( X.block(6*i,0,6,1) );
                Matrix4d Tcurr = Tprev * inverse_se3( expmap_se3( DX.block(6*i,0,6,1) ) );
                X.block(6*i,0,6,1) = logmap_se3( Tcurr );
            }
            // update LMs
            for( int i = 6*Nkf; i < N; i++)
                X(i) += DX(i);
        }
        // if the parameter change is small stop (TODO: change with two parameters, one for R and another one for t)
        if( DX.norm() < numeric_limits<double>::epsilon() )
            break;
        // update previous values
        err_prev = err;
    }
    cout << endl << endl << "GBA iterations: " << iters << endl;

    // Update KFs and LMs
    //---------------------------------------------------------------------------------------------
    // update KFs
    for( int i = 0; i < Nkf; i++)
    {
        Matrix4d Test = expmap_se3( X.block( 6*i,0,6,1 ) );
        map_keyframes[ kf_list[i] ]->T_kf_w = Test;
    }
    // update point LMs
    for( int i = 0; i < Npt; i++)
    {
        map_points[ pt_list[i] ]->point3D(0) = X(6*Nkf+3*i);
        map_points[ pt_list[i] ]->point3D(1) = X(6*Nkf+3*i+1);
        map_points[ pt_list[i] ]->point3D(2) = X(6*Nkf+3*i+2);
    }
    // update line segment LMs
    for( int i = 0; i < Nls; i++)
    {
        map_lines[ ls_list[i] ]->line3D(0) = X(6*Nkf+3*Npt+6*i);
        map_lines[ ls_list[i] ]->line3D(1) = X(6*Nkf+3*Npt+6*i+1);
        map_lines[ ls_list[i] ]->line3D(2) = X(6*Nkf+3*Npt+6*i+2);
        map_lines[ ls_list[i] ]->line3D(3) = X(6*Nkf+3*Npt+6*i+3);
        map_lines[ ls_list[i] ]->line3D(4) = X(6*Nkf+3*Npt+6*i+4);
        map_lines[ ls_list[i] ]->line3D(5) = X(6*Nkf+3*Npt+6*i+5);
    }

}

// -----------------------------------------------------------------------------------------------------------------------------
// Culling functions
// -----------------------------------------------------------------------------------------------------------------------------

void MapHandler::removeBadMapLandmarks()
{

    // point features
    for( vector<MapPoint*>::iterator pt_it = map_points.begin(); pt_it != map_points.end(); pt_it++)
    {
        if( (*pt_it)!=NULL )
        {
            if( (*pt_it)->local == false && max_kf_idx - (*pt_it)->kf_obs_list[0] > 10 )
            {
                if( (*pt_it)->inlier == false || (*pt_it)->obs_list.size() < Config::minLMObs() )
                {
                    int kf_obs = (*pt_it)->kf_obs_list[0];
                    int lm_idx = (*pt_it)->idx;
                    // remove idx from KeyFrame stereo points
                    for(vector<PointFeature*>::iterator st_pt = map_keyframes[kf_obs]->stereo_frame->stereo_pt.begin();
                        st_pt != map_keyframes[kf_obs]->stereo_frame->stereo_pt.end(); st_pt++ )
                    {
                        if( (*st_pt)->idx == lm_idx )
                        {
                            (*st_pt)->idx = -1;
                            st_pt = map_keyframes[kf_obs]->stereo_frame->stereo_pt.end()-1;
                        }
                    }
                    // remove from map_points_kf_idx
                    int iter = 0;
                    for( auto it = map_points_kf_idx.at(kf_obs).begin(); it != map_points_kf_idx.at(kf_obs).end(); it++, iter++)
                    {
                        if( (*it) == lm_idx )
                        {
                            map_points_kf_idx.at(kf_obs).erase( map_points_kf_idx.at(kf_obs).begin() + iter );
                            break;
                        }
                    }
                    // remove LM
                    (*pt_it) = NULL;
                }
            }
        }
    }

    // line features
    for( vector<MapLine*>::iterator ls_it = map_lines.begin(); ls_it != map_lines.end(); ls_it++)
    {
        if( (*ls_it)!=NULL )
        {
            if( (*ls_it)->local == false && max_kf_idx-(*ls_it)->kf_obs_list[0] > 10 )
            {
                if( (*ls_it)->inlier == false || (*ls_it)->obs_list.size() < Config::minLMObs() )
                {
                    int kf_obs = (*ls_it)->kf_obs_list[0];
                    int lm_idx = (*ls_it)->idx;
                    // remove idx from KeyFrame stereo points
                    for(vector<LineFeature*>::iterator st_ls = map_keyframes[kf_obs]->stereo_frame->stereo_ls.begin();
                        st_ls != map_keyframes[kf_obs]->stereo_frame->stereo_ls.end(); st_ls++ )
                    {
                        if( (*st_ls)->idx == lm_idx )
                        {
                            (*st_ls)->idx = -1;
                            st_ls = map_keyframes[kf_obs]->stereo_frame->stereo_ls.end()-1;
                        }
                    }
                    // remove from map_points_kf_idx
                    int iter = 0;
                    for( auto it = map_lines_kf_idx.at(kf_obs).begin(); it != map_lines_kf_idx.at(kf_obs).end(); it++, iter++)
                    {
                        if( (*it) == lm_idx )
                        {
                            map_lines_kf_idx.at(kf_obs).erase( map_lines_kf_idx.at(kf_obs).begin() + iter );
                            break;
                        }
                    }

                    // remove LM
                    (*ls_it) = NULL;
                }
            }
        }
    }

}

void MapHandler::removeRedundantKFs()
{

    // select which KFs are going to remove
    vector<int> kf_idxs;
    for( vector<KeyFrame*>::iterator kf_it = map_keyframes.begin(); kf_it != map_keyframes.end(); kf_it++)
    {
        if( (*kf_it)!= NULL )
        {
            int kf_idx = (*kf_it)->kf_idx;
            if( !(*kf_it)->local && kf_idx > 1 && kf_idx < max_kf_idx )
            {
                // estimate number of landmarks observed by this KF
                int n_feats = 0;
                for( vector<PointFeature*>::iterator pt_it = (*kf_it)->stereo_frame->stereo_pt.begin();
                     pt_it != (*kf_it)->stereo_frame->stereo_pt.end(); pt_it++ )
                {
                    if( (*pt_it)->idx != -1 )
                        n_feats++;
                }
                for( vector<LineFeature*>::iterator ls_it = (*kf_it)->stereo_frame->stereo_ls.begin();
                     ls_it != (*kf_it)->stereo_frame->stereo_ls.end(); ls_it++)
                {
                    if( (*ls_it)->idx != -1 )
                        n_feats++;
                }
                int max_n_feats = int( Config::maxCommonFtsKF() * double(n_feats) );
                // check if the KF is redundant
                int n_graph = full_graph.size();
                for(int i = 0; i < n_graph-1, i != kf_idx; i++)
                {
                    if( map_keyframes[i] != NULL )
                    {
                        if( full_graph[kf_idx][i] > n_feats )
                            kf_idxs.push_back(kf_idx);
                    }
                }
            }
        }
    }
    if( kf_idxs.size() != 0 )
    {
        cout << endl << "KFs to be erased: ";
        for(int i = 0; i < kf_idxs.size(); i++)
            cout << kf_idxs[i] << " " ;
    }
    cout << endl ;

    // eliminate KFs, LMs observed only by this KFs, and all observations from this KF
    for( int i = 0; i < kf_idxs.size(); i++)
    {
        int kf_idx = kf_idxs[i];
        if( map_keyframes[kf_idx] != NULL )
        {
            // delete observation from map_points_kf_idx
            if( !(map_points_kf_idx.find(kf_idx) == map_points_kf_idx.end()) )
            {
                for( auto it = map_points_kf_idx.at(kf_idx).begin(); it != map_points_kf_idx.at(kf_idx).end(); it++)
                {
                    if( map_points[(*it)]!= NULL &&  map_points[(*it)]->kf_obs_list.size() <= 1 )
                    {
                        bool found = false;
                        for( int k = 1; k < map_points[(*it)]->kf_obs_list.size(); k++)
                        {
                            map_points_kf_idx.find( map_points[(*it)]->kf_obs_list[k] );
                            int new_kf_base = map_points[(*it)]->kf_obs_list[k];    // if the second dont exist...
                            if( !(map_points_kf_idx.find(new_kf_base)==map_points_kf_idx.end()) )
                            {
                                map_points_kf_idx.at(new_kf_base).push_back( (*it) );
                                found = true;
                            }
                        }
                        if( !found )
                            map_points[(*it)]->inlier = false;
                    }
                }
                map_points_kf_idx.erase(kf_idx);
            }

            // delete observation from map_lines_kf_idx
            for( auto it = map_lines_kf_idx.at(kf_idx).begin(); it != map_lines_kf_idx.at(kf_idx).end(); it++)
            {
                if( map_points[(*it)]!= NULL )
                {
                    int new_kf_base = map_lines[(*it)]->kf_obs_list[1];
                    map_lines_kf_idx.at(new_kf_base).push_back( (*it) );
                }
            }
            map_lines_kf_idx.erase(kf_idx);

            // iterate over point features
            for( vector<PointFeature*>::iterator pt_it = map_keyframes[kf_idx]->stereo_frame->stereo_pt.begin();
                 pt_it != map_keyframes[kf_idx]->stereo_frame->stereo_pt.end(); pt_it++)
            {
                int pt_idx = (*pt_it)->idx;
                if( pt_idx != -1 )
                {
                    if( map_points[pt_idx] != NULL )
                    {
                        for( int j = 0; j < map_points[pt_idx]->obs_list.size(); j++)
                        {
                            if( map_points[pt_idx]->kf_obs_list[j] == kf_idx )
                            {
                                // delete observation
                                map_points[pt_idx]->desc_list.erase( map_points[pt_idx]->desc_list.begin() + j );
                                map_points[pt_idx]->obs_list.erase( map_points[pt_idx]->obs_list.begin() + j );
                                map_points[pt_idx]->dir_list.erase( map_points[pt_idx]->dir_list.begin() + j );
                                map_points[pt_idx]->kf_obs_list.erase( map_points[pt_idx]->kf_obs_list.begin() + j );
                                // update main descriptor and direction
                                map_points[pt_idx]->updateAverageDescDir();
                            }
                        }
                    }
                }
            }

            // iterate over line segment features
            for( vector<LineFeature*>::iterator ls_it = map_keyframes[kf_idx]->stereo_frame->stereo_ls.begin();
                 ls_it != map_keyframes[kf_idx]->stereo_frame->stereo_ls.end(); ls_it++)
            {
                int ls_idx = (*ls_it)->idx;
                if( ls_idx != -1 )
                {
                    if( map_lines[ls_idx] != NULL )
                    {
                        for( int j = 0; j < map_lines[ls_idx]->obs_list.size(); j++)
                        {
                            if( map_lines[ls_idx]->kf_obs_list[j] == kf_idx )
                            {
                                // delete observation
                                map_lines[ls_idx]->desc_list.erase( map_lines[ls_idx]->desc_list.begin() + j );
                                map_lines[ls_idx]->obs_list.erase( map_lines[ls_idx]->obs_list.begin() + j );
                                map_lines[ls_idx]->dir_list.erase( map_lines[ls_idx]->dir_list.begin() + j );
                                map_lines[ls_idx]->kf_obs_list.erase( map_lines[ls_idx]->kf_obs_list.begin() + j );
                                map_lines[ls_idx]->pts_list.erase( map_lines[ls_idx]->pts_list.begin() + j );
                                // update main descriptor and direction
                                map_lines[ls_idx]->updateAverageDescDir();
                            }
                        }
                    }
                }
            }

            // update full graph
            for( int k = 0; k < full_graph.size()-1; k++ )
            {
                full_graph[kf_idx][k] = 0;
                full_graph[k][kf_idx] = 0;
            }

            // erase KF
            map_keyframes[kf_idx] = NULL;
        }
    }

    if( kf_idxs.size() != 0 )
        cout << endl << "Erased, now update graph" << endl;

    // update graphs
    if( kf_idxs.size() != 0 )
        cout << endl << "Success erasing KFs" << endl;

}

// -----------------------------------------------------------------------------------------------------------------------------
// Loop Closure functions
// -----------------------------------------------------------------------------------------------------------------------------

void MapHandler::loopClosure()
{

    // look for loop closure candidates
    int kf_prev_idx, kf_curr_idx;
    kf_curr_idx = max_kf_idx;
    clock.Tic();
    bool is_lc_candidate = lookForLoopCandidates(kf_curr_idx,kf_prev_idx);
    time(4) = 1000 * clock.Tac(); //ms

    // compute relative transformation if it is LC candidate
    if( is_lc_candidate )
    {
        vector<Vector4i> lc_pt_idx, lc_ls_idx;
        vector<PointFeature*> lc_points;
        vector<LineFeature*>  lc_lines;
        Vector6d pose_inc;
        clock.Tic();
        bool isLC = isLoopClosure( map_keyframes[kf_prev_idx], map_keyframes[kf_curr_idx], pose_inc, lc_pt_idx, lc_ls_idx, lc_points, lc_lines );
        time(5) = 1000 * clock.Tac(); //ms
        // if it is loop closure, add information and update status
        if( isLC )
        {
            lc_pt_idxs.push_back( lc_pt_idx );
            lc_ls_idxs.push_back( lc_ls_idx );
            lc_poses.push_back( pose_inc );
            lc_pose_list.push_back( pose_inc );
            Vector3i lc_idx;
            lc_idx(0) = kf_prev_idx;
            lc_idx(1) = kf_curr_idx;
            lc_idx(2) = 1;
            lc_idxs.push_back( lc_idx );
            lc_idx_list.push_back(lc_idx);
            if( lc_status == LC_IDLE )
            {
                cout << endl << "Lopp closure is now active." << endl;
                lc_status = LC_ACTIVE;
            }
        }
        else
        {
            if( lc_status == LC_ACTIVE )
                lc_status = LC_READY;
        }
        cout << "isLC      = " << isLC << endl;
        cout << "LC_status = " << lc_status << endl;
    }
    else
    {
        if( lc_status == LC_ACTIVE )
            lc_status = LC_READY;
    }

    // LC computation
    if( lc_status == LC_READY )  // add condition indicating that the LC has finished (i.e. the car pass by an already visited street)
    {
        clock.Tic();
        loopClosureOptimizationCovGraphG2O();
        time(6) = 1000 * clock.Tac(); //ms
        lc_status = LC_IDLE;
    }

}

void MapHandler::insertKFBowVectorP( KeyFrame* &kf  )
{

    // transform Mat to vector<Mat>
    vector<Mat> curr_desc;
    curr_desc.reserve( kf->stereo_frame->pdesc_l.rows );
    for ( int i = 0; i < kf->stereo_frame->pdesc_l.rows; i++ )
    {
        curr_desc.push_back( kf->stereo_frame->pdesc_l.row(i) );
    }
    // transform to DBoW2::BowVector
    BowVector kf_bv;
    dbow_voc_p.transform( curr_desc, kf_bv );
    kf->descDBoW_P = kf_bv;
    // fill the confusion matrix for the new KF
    int idx = kf->kf_idx;
    int max_kf = map_keyframes.size();
    for( int i = 0; i < idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            double score = dbow_voc_p.score( kf_bv, map_keyframes[i]->descDBoW_P );
            conf_matrix[idx][i] = score;
            conf_matrix[i][idx] = score;
        }
    }
    conf_matrix[idx][idx] = dbow_voc_p.score( kf_bv, kf_bv );

}

void MapHandler::insertKFBowVectorL( KeyFrame* &kf  )
{

    // transform Mat to vector<Mat>
    vector<Mat> curr_desc;
    curr_desc.reserve( kf->stereo_frame->ldesc_l.rows );
    for ( int i = 0; i < kf->stereo_frame->ldesc_l.rows; i++ )
    {
        curr_desc.push_back( kf->stereo_frame->ldesc_l.row(i) );
    }
    // transform to DBoW2::BowVector
    BowVector kf_bv;
    dbow_voc_l.transform( curr_desc, kf_bv );
    kf->descDBoW_L = kf_bv;
    // fill the confusion matrix for the new KF
    int idx = kf->kf_idx;
    int max_kf = map_keyframes.size();
    for( int i = 0; i < idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            double score = dbow_voc_l.score( kf_bv, map_keyframes[i]->descDBoW_L );
            conf_matrix[idx][i] = score;
            conf_matrix[i][idx] = score;
        }
    }
    conf_matrix[idx][idx] = dbow_voc_l.score( kf_bv, kf_bv );

}

void MapHandler::insertKFBowVectorPL( KeyFrame* &kf  )
{

    // Point Features
    // --------------------------------------------------------------
    // transform Mat to vector<Mat>
    vector<Mat> curr_desc;
    curr_desc.reserve( kf->stereo_frame->pdesc_l.rows );
    for ( int i = 0; i < kf->stereo_frame->pdesc_l.rows; i++ )
        curr_desc.push_back( kf->stereo_frame->pdesc_l.row(i) );
    // transform to DBoW2::BowVector
    BowVector kf_bv_p;
    dbow_voc_p.transform( curr_desc, kf_bv_p );
    kf->descDBoW_P = kf_bv_p;
    // estimate dispersion of point features
    vector<double> pt_x, pt_y;
    for( auto it = kf->stereo_frame->stereo_pt.begin(); it != kf->stereo_frame->stereo_pt.end(); it++ )
    {
        pt_x.push_back( (*it)->pl(0) );
        pt_y.push_back( (*it)->pl(1) );
    }
    double std_pt = vector_stdv( pt_x ) + vector_stdv( pt_y );
    int    n_pt   = pt_x.size();

    // Line Segment Features
    // --------------------------------------------------------------
    // transform Mat to vector<Mat>
    curr_desc.clear();
    curr_desc.reserve( kf->stereo_frame->ldesc_l.rows );
    for ( int i = 0; i < kf->stereo_frame->ldesc_l.rows; i++ )
        curr_desc.push_back( kf->stereo_frame->ldesc_l.row(i) );
    // transform to DBoW2::BowVector
    BowVector kf_bv_l;
    dbow_voc_l.transform( curr_desc, kf_bv_l );
    kf->descDBoW_L = kf_bv_l;
    // estimate dispersion of point features
    vector<double> ls_x, ls_y;
    for( auto it = kf->stereo_frame->stereo_ls.begin(); it != kf->stereo_frame->stereo_ls.end(); it++ )
    {
        Vector2d mp;
        mp << ((*it)->spl + (*it)->epl)*0.5;
        ls_x.push_back( mp(0) );
        ls_y.push_back( mp(1) );
    }
    double std_ls = vector_stdv( ls_x ) + vector_stdv( ls_y );
    double std_pl = std_ls + std_pt;
    int    n_ls  = ls_x.size();
    int    n_pl  = n_pt + n_ls;

    // Combined Approach
    // --------------------------------------------------------------
    // fill the confusion matrix for the new KF
    double score, score_p, score_l;
    int idx = kf->kf_idx;
    int max_kf = map_keyframes.size();
    for( int i = 0; i < idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            score_p = dbow_voc_p.score( kf_bv_p, map_keyframes[i]->descDBoW_P );
            score_l = dbow_voc_l.score( kf_bv_l, map_keyframes[i]->descDBoW_L );
            score = 0.0;
            score += ( score_p * n_pt   + score_l * n_ls   ) / n_pl;    // strategy#1
            score += ( score_p * std_pt + score_l * std_ls ) / std_pl;  // strategy#2
            conf_matrix[idx][i] = score;
            conf_matrix[i][idx] = score;
        }
    }
    score_p = dbow_voc_p.score( kf_bv_p, kf_bv_p );
    score_l = dbow_voc_l.score( kf_bv_l, kf_bv_l );
    score = 0.0;
    score += ( score_p * n_pt   + score_l * n_ls   ) / n_pl;    // strategy#1
    score += ( score_p * std_pt + score_l * std_ls ) / std_pl;  // strategy#2
    conf_matrix[idx][idx] = score;

}

bool MapHandler::lookForLoopCandidates( int kf_curr_idx, int &kf_prev_idx )
{

    // variables to be included in config ( besides K which is used in loopClosure() )
    bool is_lc_candidate = false;
    kf_prev_idx = -1;

    // find the best matches
    vector<Vector2d> max_confmat;
    for(int i = 0; i < kf_curr_idx - Config::lcKFDist(); i++)
    {
        if( map_keyframes[i]!=NULL )
        {
            Vector2d aux;
            aux(0) = i;
            aux(1) = conf_matrix[i][kf_curr_idx];
            max_confmat.push_back( aux );
        }
    }
    sort( max_confmat.begin(), max_confmat.end(), sort_confmat_by_score() );

    // find the minimum score in the covisibility graph (and/or 3 previous keyframes)
    double lc_min_score = 1.0;
    for( int i = 0; i < kf_curr_idx; i++ )
    {
        if( full_graph[i][kf_curr_idx] >= Config::minLMCovGraph() || kf_curr_idx - i <= Config::minKFLocalMap()+3 )
        {
            double score_i = conf_matrix[i][kf_curr_idx];
            if( score_i < lc_min_score && score_i > 0.001 )
                lc_min_score = score_i;
        }
    }

    // if there are enough matches..
    if( max_confmat.size() > Config::lcKFMaxDist() )
    {

        // the best match must has an score above lc_dbow_score_max
        int idx_max = int( max_confmat[0](0) );
        int Nkf_closest = 0;
        if( max_confmat[0](1) >= lc_min_score )
        {
            // there must be at least lc_nkf_closest KFs conected to the LC candidate with a score above lc_dbow_score_min
            for( int i = 1; i < max_confmat.size(); i++ )
            {
                int idx = int( max_confmat[i](0) );
                // frame closest && connected by the cov_graph && score > lc_dbow_score_min
                if( abs(idx-idx_max) <= Config::lcKFMaxDist() &&
                    //full_graph[idx][idx_max] >= Config::minLMCovGraph() &&
                    max_confmat[i](1) >= lc_min_score * 0.8 )
                    Nkf_closest++;
            }
        }

        // update in case of being loop closure candidate
        if( Nkf_closest >= Config::lcNKFClosest() )
        {
            is_lc_candidate = true;
            kf_prev_idx     = idx_max;
        }


        // ****************************************************************** //
        cout << endl << "lc_min_score: " << lc_min_score        ;
        cout << endl << "Nkf_closest:  " << Nkf_closest         ;
        cout << endl << "idx_max:  "     << idx_max  << endl;
        for( int i = 0; i < Config::lcKFMaxDist(); i++ )
            cout <<  max_confmat[i](0) << "\t" <<  max_confmat[i](1) << endl;
        // ****************************************************************** //

    }

    return is_lc_candidate;

}

bool MapHandler::isLoopClosure( KeyFrame* kf0, KeyFrame* kf1, Vector6d &pose_inc,
                                vector<Vector4i> &lc_pt_idx, vector<Vector4i> &lc_ls_idx,
                                vector<PointFeature*> &lc_points, vector<LineFeature*>  &lc_lines)
{

    // grab frame number
    int kf0_idx = kf0->kf_idx;
    int kf1_idx = kf1->kf_idx;

    // number of stereo matches
    int n_pt_0 = kf0->stereo_frame->stereo_pt.size();
    int n_pt_1 = kf1->stereo_frame->stereo_pt.size();
    int n_ls_0 = kf0->stereo_frame->stereo_ls.size();
    int n_ls_1 = kf1->stereo_frame->stereo_ls.size();

    // initial pose increment
    Matrix4d DT = Matrix4d::Identity();

    // structures to optimize camera drift
    lc_points.clear();
    lc_lines.clear();

    // find matches between both KFs
    // ---------------------------------------------------
    // points f2f tracking
    int common_pt = 0;
    if( Config::hasPoints() && !(kf1->stereo_frame->stereo_pt.size()==0) && !(kf0->stereo_frame->stereo_pt.size()==0)  )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        Mat pdesc_l1, pdesc_l2;
        vector<vector<DMatch>> pmatches_12, pmatches_21;
        // 12 and 21 matches
        pdesc_l1 = kf0->stereo_frame->pdesc_l;
        pdesc_l2 = kf1->stereo_frame->pdesc_l;
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StVO::StereoFrame::matchPointFeatures, kf0->stereo_frame, bfm, pdesc_l1, pdesc_l2, ref(pmatches_12) );
                auto match_r = async( launch::async, &StVO::StereoFrame::matchPointFeatures, kf0->stereo_frame, bfm, pdesc_l2, pdesc_l1, ref(pmatches_21) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( pdesc_l1, pdesc_l2, pmatches_12, 2);
                bfm->knnMatch( pdesc_l2, pdesc_l1, pmatches_21, 2);
            }
        }
        else
            bfm->knnMatch( pdesc_l1, pdesc_l2, pmatches_12, 2);
        // sort matches by the distance between the best and second best matches
        double nn12_dist_th = Config::minRatio12P();
        // resort according to the queryIdx
        sort( pmatches_12.begin(), pmatches_12.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( pmatches_21.begin(), pmatches_21.end(), sort_descriptor_by_queryIdx() );
        // bucle around pmatches
        for( int i = 0; i < pmatches_12.size(); i++ )
        {
            // check if they are mutual best matches
            int lr_qdx = pmatches_12[i][0].queryIdx;
            int lr_tdx = pmatches_12[i][0].trainIdx;
            int rl_tdx;
            if( Config::bestLRMatches() )
                rl_tdx = pmatches_21[lr_tdx][0].trainIdx;
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_nn = pmatches_12[i][0].distance;
            double dist_12 = pmatches_12[i][0].distance / pmatches_12[i][1].distance;
            // check the f2f max disparity condition
            if( lr_qdx == rl_tdx  )
            {
                common_pt++;
                // save data for optimization
                Vector3d P       = kf0->stereo_frame->stereo_pt[lr_qdx]->P;
                Vector2d pl_obs  = kf1->stereo_frame->stereo_pt[lr_tdx]->pl;
                PointFeature* pt = new PointFeature( P, pl_obs );
                lc_points.push_back(pt);
                // save indices for fusing LMs
                Vector4i idx;
                idx(0) = kf0->stereo_frame->stereo_pt[lr_qdx]->idx;
                idx(1) = lr_qdx;
                idx(2) = kf1->stereo_frame->stereo_pt[lr_tdx]->idx;
                idx(3) = lr_tdx;
                lc_pt_idx.push_back(idx);
            }
        }
    }

    // line segments f2f tracking
    int common_ls = 0;
    if( Config::hasLines() && !(kf1->stereo_frame->stereo_ls.size()==0) && !(kf0->stereo_frame->stereo_ls.size()==0)  )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        Mat ldesc_l1, ldesc_l2;
        vector<vector<DMatch>> lmatches_12, lmatches_21;
        // 12 and 21 matches
        ldesc_l1 = kf0->stereo_frame->ldesc_l;
        ldesc_l2 = kf1->stereo_frame->ldesc_l;
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StVO::StereoFrame::matchLineFeatures, kf0->stereo_frame, bfm, ldesc_l1, ldesc_l2, ref(lmatches_12) );
                auto match_r = async( launch::async, &StVO::StereoFrame::matchLineFeatures, kf0->stereo_frame, bfm, ldesc_l2, ldesc_l1, ref(lmatches_21) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( ldesc_l1, ldesc_l2, lmatches_12, 2);
                bfm->knnMatch( ldesc_l2, ldesc_l1, lmatches_21, 2);
            }
        }
        else
            bfm->knnMatch( ldesc_l1, ldesc_l2, lmatches_12, 2);

        // sort matches by the distance between the best and second best matches
        double nn_dist_th, nn12_dist_th;
        kf0->stereo_frame->lineDescriptorMAD(lmatches_12,nn_dist_th, nn12_dist_th);
        nn12_dist_th  = nn12_dist_th * Config::descThL();
        // resort according to the queryIdx
        sort( lmatches_12.begin(), lmatches_12.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( lmatches_21.begin(), lmatches_21.end(), sort_descriptor_by_queryIdx() );
        // bucle around lmatches
        for( int i = 0; i < lmatches_12.size(); i++ )
        {
            // check if they are mutual best matches
            int lr_qdx = lmatches_12[i][0].queryIdx;
            int lr_tdx = lmatches_12[i][0].trainIdx;
            int rl_tdx;
            if( Config::bestLRMatches() )
                rl_tdx = lmatches_21[lr_tdx][0].trainIdx;
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_12 = lmatches_12[i][1].distance - lmatches_12[i][0].distance;
            // f2f angle diff and flow
            double a1 = kf0->stereo_frame->stereo_ls[lr_qdx]->angle;
            double a2 = kf1->stereo_frame->stereo_ls[lr_tdx]->angle;
            Vector2d x1 = (kf0->stereo_frame->stereo_ls[lr_qdx]->spl + kf0->stereo_frame->stereo_ls[lr_qdx]->epl);
            Vector2d x2 = (kf1->stereo_frame->stereo_ls[lr_tdx]->spl + kf1->stereo_frame->stereo_ls[lr_tdx]->epl);
            if( lr_qdx == rl_tdx  )
            {
                common_ls++;
                // save data for optimization
                Vector3d sP     = kf0->stereo_frame->stereo_ls[lr_qdx]->sP;
                Vector3d eP     = kf0->stereo_frame->stereo_ls[lr_qdx]->eP;
                Vector3d le_obs = kf1->stereo_frame->stereo_ls[lr_tdx]->le;
                LineFeature* ls = new LineFeature( sP, eP, le_obs, kf1->stereo_frame->stereo_ls[lr_tdx]->spl, kf1->stereo_frame->stereo_ls[lr_tdx]->epl );
                lc_lines.push_back(ls);
                // save indices for fusing LMs
                Vector4i idx;
                idx(0) = kf0->stereo_frame->stereo_ls[lr_qdx]->idx;
                idx(1) = lr_qdx;
                idx(2) = kf1->stereo_frame->stereo_ls[lr_tdx]->idx;
                idx(3) = lr_tdx;
                lc_ls_idx.push_back(idx);
            }
        }

    }

    cout << endl << "[Checking Loop Closure]:" << endl;
    cout << "%Found points: " << 100.0 * common_pt / n_pt_0 << " / " << 100.0 * common_pt / n_pt_1 << endl;
    cout << "%Found lines:  " << 100.0 * common_ls / n_ls_0 << " / " << 100.0 * common_ls / n_ls_1 << endl;
    cout << "Common points (without map): " << common_pt << " \t Common lines (without map): " << common_ls << endl ;

    // estimate relative pose between both KFs
    // ---------------------------------------------------
    double inl_ratio_pt = max( 100.0 * common_pt / n_pt_0, 100.0 * common_pt / n_pt_1 );
    double inl_ratio_ls = max( 100.0 * common_ls / n_ls_0, 100.0 * common_ls / n_ls_1 );
    bool inl_ratio_condition = false;

    // TODO: ignore this condition?
    inl_ratio_condition = true;
    if( Config::hasPoints() && Config::hasLines() )
    {
        if( inl_ratio_pt > Config::lcInlierRatio() && inl_ratio_ls >Config::lcInlierRatio() )
            inl_ratio_condition = true;
    }
    else if( Config::hasPoints() && !Config::hasLines() )
    {
        if( inl_ratio_pt > Config::lcInlierRatio() )
            inl_ratio_condition = true;
    }
    else if( !Config::hasPoints() && Config::hasLines() )
    {
        if( inl_ratio_ls > Config::lcInlierRatio() )
            inl_ratio_condition = true;
    }
    if( inl_ratio_condition )
        return computeRelativePoseGN(lc_points,lc_lines,lc_pt_idx,lc_ls_idx,pose_inc);
    else
        return false;

}

bool MapHandler::computeRelativePoseGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                        vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                        Vector6d &pose_inc)
{

    // create GN variables
    Vector6d x_inc = Vector6d::Zero(), x_prev = Vector6d::Zero();
    Matrix4d T_inc = Matrix4d::Identity(), T_prev = Matrix4d::Zero();
    Matrix6d H_l, H_p, H;
    Vector6d g_l, g_p, g;
    H = Matrix6d::Zero();
    g = Vector6d::Zero();
    H_l = H; H_p = H;
    g_l = g; g_p = g;
    double   e_l = 0.0, e_p = 0.0, e = 0.0, S_l, S_p;

    // create GN parameters
    double err, err_prev = 999999999.9;
    int    max_iters_first = Config::maxIters(), max_iters = Config::maxItersRef();

    // GN iterations
    //---------------------------------------------------------------------------------------------
    for( int iters = 0; iters < max_iters_first; iters++)
    {
        H = Matrix6d::Zero(); H_l = H; H_p = H;
        g = Vector6d::Zero(); g_l = g; g_p = g;
        e = 0.0; e_p = 0.0; e_l = 0.0;
        // point observations
        int N_p = 0;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
        {
            if( (*pt_it)->inlier )
            {
                Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
                Vector2d pl_proj = cam->projection( P_ );
                // projection error
                Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
                double err_i_norm = err_i.norm();
                double gx   = P_(0);
                double gy   = P_(1);
                double gz   = P_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                double dx   = err_i(0);
                double dy   = err_i(1);
                // jacobian
                Vector6d J_aux;
                J_aux << + fgz2 * dx * gz,
                         + fgz2 * dy * gz,
                         - fgz2 * ( gx*dx + gy*dy ),
                         - fgz2 * ( gx*gy*dx + gy*gy*dy + gz*gz*dy ),
                         + fgz2 * ( gx*gx*dx + gz*gz*dx + gx*gy*dy ),
                         + fgz2 * ( gx*gz*dy - gy*gz*dx );
                J_aux = J_aux / std::max(Config::homogTh(),err_i_norm);
                // if employing robust cost function
                double s2 = (*pt_it)->sigma2;
                double w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
                // update hessian, gradient, and error
                H_p += J_aux * J_aux.transpose() * w;
                g_p += J_aux * err_i_norm * w;
                e_p += err_i_norm * err_i_norm * w;
                N_p++;
            }
        }
        // line segment features
        int N_l = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
        {
            if( (*ls_it)->inlier )
            {
                Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
                Vector2d spl_proj = cam->projection( sP_ );
                Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
                Vector2d epl_proj = cam->projection( eP_ );
                Vector3d l_obs = (*ls_it)->le_obs;
                // projection error
                Vector2d err_i;
                err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
                err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
                double err_i_norm = err_i.norm();
                // start point
                double gx   = sP_(0);
                double gy   = sP_(1);
                double gz   = sP_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                double ds   = err_i(0);
                double de   = err_i(1);
                double lx   = l_obs(0);
                double ly   = l_obs(1);
                Vector6d Js_aux;
                Js_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // end point
                gx   = eP_(0);
                gy   = eP_(1);
                gz   = eP_(2);
                gz2  = gz*gz;
                fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                Vector6d Je_aux, J_aux;
                Je_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // jacobian
                J_aux = ( Js_aux * ds + Je_aux * de ) / std::max(Config::homogTh(),err_i_norm);
                // if employing robust cost function
                double s2 = (*ls_it)->sigma2;
                double w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
                // update hessian, gradient, and error
                H_l += J_aux * J_aux.transpose() * w;
                g_l += J_aux * err_i_norm * w;
                e_l += err_i_norm * err_i_norm * w;
                N_l++;
            }
        }
        // sum H, g and err from both points and lines
        H = H_p + H_l;
        g = g_p + g_l;
        e = e_p + e_l;
        // normalize error
        e /= (N_l+N_p);
        // if the difference is very small stop
        if( abs(e-err_prev) < numeric_limits<double>::epsilon() || e < numeric_limits<double>::epsilon() )
            break;
        // solve
        LDLT<MatrixXd> solver(H);
        x_inc = solver.solve( g );
        T_inc = T_inc * inverse_se3( expmap_se3(x_inc) );
        // if the parameter change is small stop (TODO: change with two parameters, one for R and another one for t)
        if( x_inc.norm() < numeric_limits<double>::epsilon() )
            break;
        // update error
        err_prev = e;
    }
    x_inc = logmap_se3(T_inc);

    // Remove outliers
    //---------------------------------------------------------------------------------------------
    // point observations
    for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
    {
        if( (*pt_it)->inlier )
        {
            Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
            Vector2d pl_proj = cam->projection( P_ );
            // projection error
            Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
            double s2 = (*pt_it)->sigma2;
            if( err_i.norm() * sqrt(s2) > sqrt(7.815) )
                (*pt_it)->inlier = false;
        }
    }
    // line segments observations
    for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
    {
        if((*ls_it)->inlier)
        {
            Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
            Vector2d spl_proj = cam->projection( sP_ );
            Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
            Vector2d epl_proj = cam->projection( eP_ );
            Vector3d l_obs = (*ls_it)->le_obs;
            // projection error
            Vector2d err_i;
            err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
            err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
            double s2 = sqrt((*ls_it)->sigma2);
            if( err_i.norm() * s2 > sqrt(7.815) )
                (*ls_it)->inlier = false;
        }
    }

    // Check whether it is Loop Closure or not
    //---------------------------------------------------------------------------------------------
    // Residue value
    bool lc_res = ( e < Config::lcRes() );

    // Uncertainty value
    Matrix6d DT_cov;
    DT_cov = H.inverse();
    SelfAdjointEigenSolver<Matrix6d> eigensolver(DT_cov);
    Vector6d DT_cov_eig = eigensolver.eigenvalues();
    bool lc_unc = ( DT_cov_eig(5) < Config::lcUnc() );

    // Ratio of outliers
    int N = lc_points.size() + lc_lines.size(), N_inl = 0;
    for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
    {
        if( (*pt_it)->inlier )
            N_inl++;
    }
    for( vector<LineFeature*>::iterator  ls_it = lc_lines.begin();  ls_it != lc_lines.end();  ls_it++ )
    {
        if( (*ls_it)->inlier )
            N_inl++;
    }
    double ratio_inliers = double(N_inl) / double(N);
    bool lc_inl = ( ratio_inliers > Config::lcInl() );
    //lc_inl = true;

    // Geometry condition
    double t = x_inc.head(3).norm();
    double r = x_inc.tail(3).norm() * 180.f / CV_PI;
    bool lc_trs = ( t < Config::lcTrs() );
    bool lc_rot = ( r < Config::lcRot() );

    // Verbose
    cout << endl << "Residue:       (" << lc_res << ")  " << e << " px" ;
    cout << endl << "Uncertainty:   (" << lc_unc << ")  " << DT_cov_eig(5) ;
    cout << endl << "Inliers ratio: (" << lc_inl << ")  " << ratio_inliers << "\t" << N_inl << " " << N ;
    cout << endl << "Translation:   (" << lc_trs << ")  " << t ;
    cout << endl << "Rotation:      (" << lc_rot << ")  " << r << endl << endl;

    // Decision
    if( lc_res && lc_unc && lc_inl && lc_trs && lc_rot )
    {
        // erase outliers from      lc_pt_idx & lc_ls_idx
        int iter = 0;
        vector<Vector4i> lc_pt_idx_, lc_ls_idx_;
        vector<PointFeature*> lc_points_;
        vector<LineFeature*>  lc_lines_;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++, iter++ )
        {
            if( (*pt_it)!= NULL )
            {
                if( (*pt_it)->inlier )
                {
                    lc_points_.push_back( (*pt_it) );
                    lc_pt_idx_.push_back( lc_pt_idx[iter] );
                }
            }
        }
        iter = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++, iter++ )
        {
            if( (*ls_it)!= NULL )
            {
                if( (*ls_it)->inlier )
                {
                    lc_lines_.push_back( (*ls_it) );
                    lc_ls_idx_.push_back( lc_ls_idx[iter] );
                }
            }
        }
        lc_pt_idx.clear();
        lc_pt_idx = lc_pt_idx_;
        lc_ls_idx.clear();
        lc_ls_idx = lc_ls_idx_;
        lc_points.clear();
        lc_points = lc_points_;
        lc_lines.clear();
        lc_points = lc_points_;
        // assign pose increment
        pose_inc = logmap_se3( inverse_se3( T_inc ) );
        return true;
    }
    else
        return false;

}

bool MapHandler::computeRelativePoseRobustGN( vector<PointFeature*> &lc_points, vector<LineFeature*> &lc_lines,
                                        vector<Vector4i>      &lc_pt_idx, vector<Vector4i>     &lc_ls_idx,
                                        Vector6d &pose_inc)
{

    // create GN variables
    Vector6d x_inc = Vector6d::Zero(), x_prev = Vector6d::Zero();
    Matrix4d T_inc = Matrix4d::Identity(), T_prev = Matrix4d::Zero();
    Matrix6d H_l, H_p, H;
    Vector6d g_l, g_p, g;
    H = Matrix6d::Zero();
    g = Vector6d::Zero();
    H_l = H; H_p = H;
    g_l = g; g_p = g;
    double   e_l = 0.0, e_p = 0.0, e = 0.0, S_l, S_p;

    // create GN parameters
    double err, err_prev = 999999999.9;
    int    max_iters_first = Config::maxIters(), max_iters = Config::maxItersRef();

    // GN iterations
    //---------------------------------------------------------------------------------------------
    for( int iters = 0; iters < max_iters_first; iters++)
    {
        H = Matrix6d::Zero(); H_l = H; H_p = H;
        g = Vector6d::Zero(); g_l = g; g_p = g;
        e = 0.0; e_p = 0.0; e_l = 0.0;
        // point observations
        int N_p = 0;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
        {
            if( (*pt_it)->inlier )
            {
                Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
                Vector2d pl_proj = cam->projection( P_ );
                // projection error
                Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
                double err_i_norm = err_i.norm();
                double gx   = P_(0);
                double gy   = P_(1);
                double gz   = P_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                double dx   = err_i(0);
                double dy   = err_i(1);
                // jacobian
                Vector6d J_aux;
                J_aux << + fgz2 * dx * gz,
                         + fgz2 * dy * gz,
                         - fgz2 * ( gx*dx + gy*dy ),
                         - fgz2 * ( gx*gy*dx + gy*gy*dy + gz*gz*dy ),
                         + fgz2 * ( gx*gx*dx + gz*gz*dx + gx*gy*dy ),
                         + fgz2 * ( gx*gz*dy - gy*gz*dx );
                J_aux = J_aux / std::max(Config::homogTh(),err_i_norm);
                // if employing robust cost function
                double s2 = (*pt_it)->sigma2;
                double w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
                // update hessian, gradient, and error
                H_p += J_aux * J_aux.transpose() * w;
                g_p += J_aux * err_i_norm * w;
                e_p += err_i_norm * err_i_norm * w;
                N_p++;
            }
        }
        // line segment features
        int N_l = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
        {
            if( (*ls_it)->inlier )
            {
                Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
                Vector2d spl_proj = cam->projection( sP_ );
                Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
                Vector2d epl_proj = cam->projection( eP_ );
                Vector3d l_obs = (*ls_it)->le_obs;
                // projection error
                Vector2d err_i;
                err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
                err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
                double err_i_norm = err_i.norm();
                // start point
                double gx   = sP_(0);
                double gy   = sP_(1);
                double gz   = sP_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                double ds   = err_i(0);
                double de   = err_i(1);
                double lx   = l_obs(0);
                double ly   = l_obs(1);
                Vector6d Js_aux;
                Js_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // end point
                gx   = eP_(0);
                gy   = eP_(1);
                gz   = eP_(2);
                gz2  = gz*gz;
                fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                Vector6d Je_aux, J_aux;
                Je_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // jacobian
                J_aux = ( Js_aux * ds + Je_aux * de ) / std::max(Config::homogTh(),err_i_norm);
                // if employing robust cost function
                double s2 = (*ls_it)->sigma2;
                double w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
                // update hessian, gradient, and error
                H_l += J_aux * J_aux.transpose() * w;
                g_l += J_aux * err_i_norm * w;
                e_l += err_i_norm * err_i_norm * w;
                N_l++;
            }
        }
        // sum H, g and err from both points and lines
        H = H_p + H_l;
        g = g_p + g_l;
        e = e_p + e_l;
        // normalize error
        e /= (N_l+N_p);

        // if the difference is very small stop
        if( abs(e-err_prev) < numeric_limits<double>::epsilon() || e < numeric_limits<double>::epsilon() )
            break;

        // solve
        LDLT<MatrixXd> solver(H);
        x_inc = solver.solve( g );
        T_inc = T_inc * inverse_se3( expmap_se3(x_inc) );

        // if the parameter change is small stop (TODO: change with two parameters, one for R and another one for t)
        if( x_inc.norm() < numeric_limits<double>::epsilon() )
            break;

        // update error
        err_prev = e;

    }

    // Remove outliers
    //---------------------------------------------------------------------------------------------
    // point observations
    for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
    {
        if( (*pt_it)->inlier )
        {
            Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
            Vector2d pl_proj = cam->projection( P_ );
            // projection error
            Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
            double s2 = (*pt_it)->sigma2;
            if( err_i.norm() * sqrt(s2) > sqrt(7.815) )
                (*pt_it)->inlier = false;
        }
    }
    // line segments observations
    for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
    {
        if((*ls_it)->inlier)
        {
            Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
            Vector2d spl_proj = cam->projection( sP_ );
            Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
            Vector2d epl_proj = cam->projection( eP_ );
            Vector3d l_obs = (*ls_it)->le_obs;
            // projection error
            Vector2d err_i;
            err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
            err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
            double s2 = sqrt((*ls_it)->sigma2);
            if( err_i.norm() * s2 > sqrt(7.815) )
                (*ls_it)->inlier = false;
        }
    }

    // GN refinement
    //---------------------------------------------------------------------------------------------
    for( int iters = 0; iters < max_iters; iters++)
    {
        H = Matrix6d::Zero(); H_l = H; H_p = H;
        g = Vector6d::Zero(); g_l = g; g_p = g;
        e = 0.0; e_p = 0.0; e_l = 0.0;
        // point observations
        int N_p = 0;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
        {
            if( (*pt_it)->inlier )
            {
                Vector3d P_ = T_inc.block(0,0,3,3) * (*pt_it)->P + T_inc.col(3).head(3);
                Vector2d pl_proj = cam->projection( P_ );
                // projection error
                Vector2d err_i    = pl_proj - (*pt_it)->pl_obs;
                double err_i_norm = err_i.norm();
                double gx   = P_(0);
                double gy   = P_(1);
                double gz   = P_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                double dx   = err_i(0);
                double dy   = err_i(1);
                // jacobian
                Vector6d J_aux;
                J_aux << + fgz2 * dx * gz,
                         + fgz2 * dy * gz,
                         - fgz2 * ( gx*dx + gy*dy ),
                         - fgz2 * ( gx*gy*dx + gy*gy*dy + gz*gz*dy ),
                         + fgz2 * ( gx*gx*dx + gz*gz*dx + gx*gy*dy ),
                         + fgz2 * ( gx*gz*dy - gy*gz*dx );
                J_aux = J_aux / std::max(Config::homogTh(),err_i_norm);
                // if employing robust cost function
                double s2 = (*pt_it)->sigma2;
                double w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
                // update hessian, gradient, and error
                H_p += J_aux * J_aux.transpose() * w;
                g_p += J_aux * err_i_norm * w;
                e_p += err_i_norm * err_i_norm * w;
                N_p++;
            }
        }
        // line segment features
        int N_l = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++ )
        {
            if( (*ls_it)->inlier )
            {
                Vector3d sP_ = T_inc.block(0,0,3,3) * (*ls_it)->sP + T_inc.col(3).head(3);
                Vector2d spl_proj = cam->projection( sP_ );
                Vector3d eP_ = T_inc.block(0,0,3,3) * (*ls_it)->eP + T_inc.col(3).head(3);
                Vector2d epl_proj = cam->projection( eP_ );
                Vector3d l_obs = (*ls_it)->le_obs;
                // projection error
                Vector2d err_i;
                err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
                err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
                double err_i_norm = err_i.norm();
                // start point
                double gx   = sP_(0);
                double gy   = sP_(1);
                double gz   = sP_(2);
                double gz2  = gz*gz;
                double fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                double ds   = err_i(0);
                double de   = err_i(1);
                double lx   = l_obs(0);
                double ly   = l_obs(1);
                Vector6d Js_aux;
                Js_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // end point
                gx   = eP_(0);
                gy   = eP_(1);
                gz   = eP_(2);
                gz2  = gz*gz;
                fgz2 = cam->getFx() / std::max(Config::homogTh(),gz2);
                Vector6d Je_aux, J_aux;
                Je_aux << + fgz2 * lx * gz,
                          + fgz2 * ly * gz,
                          - fgz2 * ( gx*lx + gy*ly ),
                          - fgz2 * ( gx*gy*lx + gy*gy*ly + gz*gz*ly ),
                          + fgz2 * ( gx*gx*lx + gz*gz*lx + gx*gy*ly ),
                          + fgz2 * ( gx*gz*ly - gy*gz*lx );
                // jacobian
                J_aux = ( Js_aux * ds + Je_aux * de ) / std::max(Config::homogTh(),err_i_norm);
                // if employing robust cost function
                double s2 = (*ls_it)->sigma2;
                double w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
                // update hessian, gradient, and error
                H_l += J_aux * J_aux.transpose() * w;
                g_l += J_aux * err_i_norm * w;
                e_l += err_i_norm * err_i_norm * w;
                N_l++;
            }
        }
        // sum H, g and err from both points and lines
        H = H_p + H_l;
        g = g_p + g_l;
        e = e_p + e_l;
        // normalize error
        e /= (N_l+N_p);

        // if the difference is very small stop
        if( abs(e-err_prev) < numeric_limits<double>::epsilon() || e < numeric_limits<double>::epsilon() )
            break;

        // solve
        LDLT<MatrixXd> solver(H);
        x_inc = solver.solve( g );
        T_inc = T_inc * inverse_se3( expmap_se3(x_inc) );

        // if the parameter change is small stop (TODO: change with two parameters, one for R and another one for t)
        if( x_inc.norm() < numeric_limits<double>::epsilon() )
            break;

        // update error
        err_prev = e;
    }

    x_inc = logmap_se3(T_inc);

    // Check whether it is Loop Closure or not
    //---------------------------------------------------------------------------------------------
    // Residue value
    bool lc_res = ( e < Config::lcRes() );

    // Uncertainty value
    Matrix6d DT_cov;
    DT_cov = H.inverse();
    SelfAdjointEigenSolver<Matrix6d> eigensolver(DT_cov);
    Vector6d DT_cov_eig = eigensolver.eigenvalues();
    bool lc_unc = ( DT_cov_eig(5) < Config::lcUnc() );

    // Ratio of outliers
    int N = lc_points.size() + lc_lines.size(), N_inl = 0;
    for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++ )
    {
        if( (*pt_it)->inlier )
            N_inl++;
    }
    for( vector<LineFeature*>::iterator  ls_it = lc_lines.begin();  ls_it != lc_lines.end();  ls_it++ )
    {
        if( (*ls_it)->inlier )
            N_inl++;
    }
    double ratio_inliers = double(N_inl) / double(N);
    bool lc_inl = ( ratio_inliers > Config::lcInl() );

    lc_inl = true;

    // Geometry condition
    double t = x_inc.head(3).norm();
    double r = x_inc.tail(3).norm() * 180.f / CV_PI;
    bool lc_trs = ( t < Config::lcTrs() );
    bool lc_rot = ( r < Config::lcRot() );

    // Verbose
    cout << endl << "Residue:       (" << lc_res << ")  " << e << " px" ;
    cout << endl << "Uncertainty:   (" << lc_unc << ")  " << DT_cov_eig(5) ;
    cout << endl << "Inliers ratio: (" << lc_inl << ")  " << ratio_inliers << "\t" << N_inl << " " << N ;
    cout << endl << "Translation:   (" << lc_trs << ")  " << t ;
    cout << endl << "Rotation:      (" << lc_rot << ")  " << r << endl << endl;

    cout << endl << T_inc << endl;

    // Decision
    if( lc_res && lc_unc && lc_inl && lc_trs && lc_rot )
    {
        // erase outliers from      lc_pt_idx & lc_ls_idx
        int iter = 0;
        vector<Vector4i> lc_pt_idx_, lc_ls_idx_;
        vector<PointFeature*> lc_points_;
        vector<LineFeature*>  lc_lines_;
        for( vector<PointFeature*>::iterator pt_it = lc_points.begin(); pt_it != lc_points.end(); pt_it++, iter++ )
        {
            if( (*pt_it)!= NULL )
            {
                if( (*pt_it)->inlier )
                {
                    lc_points_.push_back( (*pt_it) );
                    lc_pt_idx_.push_back( lc_pt_idx[iter] );
                }
            }
        }
        iter = 0;
        for( vector<LineFeature*>::iterator ls_it = lc_lines.begin(); ls_it != lc_lines.end(); ls_it++, iter++ )
        {
            if( (*ls_it)!= NULL )
            {
                if( (*ls_it)->inlier )
                {
                    lc_lines_.push_back( (*ls_it) );
                    lc_ls_idx_.push_back( lc_ls_idx[iter] );
                }
            }
        }
        lc_pt_idx.clear();
        lc_pt_idx = lc_pt_idx_;
        lc_ls_idx.clear();
        lc_ls_idx = lc_ls_idx_;
        lc_points.clear();
        lc_points = lc_points_;
        lc_lines.clear();
        lc_points = lc_points_;
        // assign pose increment
        pose_inc = logmap_se3( inverse_se3( expmap_se3(x_inc) ) );
        return true;
    }
    else
        return false;

}

bool MapHandler::loopClosureOptimizationEssGraphG2O()
{

    // define G2O variables
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(true);
    g2o::BlockSolver_6_3::LinearSolverType* linearSolver;
    linearSolver = new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>();
    g2o::BlockSolver_6_3* solver_ptr = new g2o::BlockSolver_6_3(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    solver->setUserLambdaInit(1e-10);
    optimizer.setAlgorithm(solver);

    // select min and max KF indices
    int kf_prev_idx =  2 * max_kf_idx;
    int kf_curr_idx = -1;
    for( auto it = lc_idxs.begin(); it != lc_idxs.end(); it++)
    {
        if( (*it)(0) < kf_prev_idx )
            kf_prev_idx = (*it)(0);
        if( (*it)(1) > kf_curr_idx )
            kf_curr_idx = (*it)(1);
    }

    // grab the KFs included in the optimization
    vector<int> kf_list;
    for( int i = kf_prev_idx; i <= kf_curr_idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            // check if it is a LC vertex
            bool is_lc_i = false;
            bool is_lc_j = false;
            int id = 0;
            for( auto it = lc_idxs.begin(); it != lc_idxs.end(); it++, id++ )
            {
                if( (*it)(0) == i )
                {
                    is_lc_i = true;
                    break;
                }
                if( (*it)(1) == i )
                {
                    is_lc_j = true;
                    break;
                }
            }
            kf_list.push_back(i);
            // create SE3 vertex
            g2o::VertexSE3* v_se3 = new g2o::VertexSE3();
            v_se3->setId(i);
            v_se3->setMarginalized( false );
            if( is_lc_j )
            {
                // update pose of LC vertex
                v_se3->setFixed(true);
                v_se3->setEstimate( g2o::SE3Quat::exp( reverse_se3(logmap_se3( (expmap_se3(lc_pose_list[id])) * map_keyframes[lc_idxs[id](0)]->T_kf_w )) ) );
            }
            else
            {
                v_se3->setEstimate( g2o::SE3Quat::exp( reverse_se3(map_keyframes[i]->x_kf_w) ) );
                if( is_lc_i || i == 0 )
                    v_se3->setFixed(true);
            }
            optimizer.addVertex( v_se3 );
        }
    }

    // introduce edges
    for( int i = kf_prev_idx; i <= kf_curr_idx; i++ )
    {
        for( int j = i+1; j <= kf_curr_idx; j++ )
        {
            if( map_keyframes[i] != NULL && map_keyframes[j] != NULL &&
                ( full_graph[i][j] >= Config::minLMEssGraph() || abs(i-j) == 1  ) )
            {
                // kf2kf constraint
                Matrix4d T_ji_constraint = inverse_se3( map_keyframes[i]->T_kf_w ) * map_keyframes[j]->T_kf_w;
                // add edge
                g2o::EdgeSE3* e_se3 = new g2o::EdgeSE3();
                e_se3->setVertex( 0, optimizer.vertex(i) );
                e_se3->setVertex( 1, optimizer.vertex(j) );
                Vector6d x;
                x = reverse_se3(logmap_se3(T_ji_constraint) );
                e_se3->setMeasurement( g2o::SE3Quat::exp(x) );
                //e_se3->information() = map_keyframes[j]->xcov_kf_w;
                e_se3->setInformation( Matrix6d::Identity() );
                optimizer.addEdge( e_se3 );
            }
        }
    }

    // introduce loop closure edges
    int id = 0;
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++, id++ )
    {
        // add edge
        g2o::EdgeSE3* e_se3 = new g2o::EdgeSE3();
        e_se3->setVertex( 0, optimizer.vertex((*it)(0)) );
        e_se3->setVertex( 1, optimizer.vertex((*it)(1)) );
        Vector6d x;
        x = reverse_se3( lc_pose_list[id] );
        e_se3->setMeasurement( g2o::SE3Quat::exp(x) );
        //e_se3->information() = kf_Tij_cov_constraint[id];
        e_se3->information() = Matrix6d::Identity();
        optimizer.addEdge( e_se3 );
    }

    // optimize graph
    optimizer.initializeOptimization();
    optimizer.computeInitialGuess();
    optimizer.computeActiveErrors();
    optimizer.optimize(100);

    // recover pose and update map
    Matrix4d Tkfw_corr;
    for( auto kf_it = kf_list.begin(); kf_it != kf_list.end(); kf_it++)
    {
        g2o::VertexSE3* v_se3 = static_cast<g2o::VertexSE3*>(optimizer.vertex( (*kf_it) ));
        g2o::SE3Quat Tiw_corr =  v_se3->estimateAsSE3Quat();
        Vector6d x;
        Matrix4d Tkfw, Tkfw_prev;
        x = reverse_se3(Tiw_corr.log());
        Tkfw = expmap_se3( x );
        Tkfw_prev = map_keyframes[ (*kf_it) ]->T_kf_w;
        map_keyframes[ (*kf_it) ]->T_kf_w = Tkfw;
        map_keyframes[ (*kf_it) ]->x_kf_w = logmap_se3(Tkfw);
        // update map
        Tkfw_corr = Tkfw * inverse_se3( Tkfw_prev );
        for( auto it = map_points_kf_idx.at((*kf_it)).begin(); it != map_points_kf_idx.at((*kf_it)).end(); it++ )
        {
           if( map_points[(*it)] != NULL )
           {
               // update 3D position
               Vector3d point3D = map_points[(*it)]->point3D;
               map_points[(*it)]->point3D = Tkfw_corr.block(0,0,3,3) * point3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_points[(*it)]->med_obs_dir;
               map_points[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // [CHECK] update direction of each observation
               for( auto dir_it = map_points[(*it)]->dir_list.begin(); dir_it != map_points[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
        for( auto it = map_lines_kf_idx.at((*kf_it)).begin(); it != map_lines_kf_idx.at((*kf_it)).end(); it++ )
        {
           if( map_lines[(*it)] != NULL )
           {
               // update 3D position
               Vector3d sP3D = map_lines[(*it)]->line3D.head(3);
               Vector3d eP3D = map_lines[(*it)]->line3D.tail(3);
               map_lines[(*it)]->line3D.head(3) = Tkfw_corr.block(0,0,3,3) * sP3D + Tkfw_corr.block(0,3,3,1);
               map_lines[(*it)]->line3D.tail(3) = Tkfw_corr.block(0,0,3,3) * eP3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_lines[(*it)]->med_obs_dir;
               map_lines[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // [CHECK] update direction of each observation
               for( auto dir_it = map_lines[(*it)]->dir_list.begin(); dir_it != map_lines[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
    }

    // update pose and map of the rest of frames
    for( int i = kf_curr_idx + 1; i < map_keyframes.size(); i++ )
    {
        // update pose
        map_keyframes[i]->T_kf_w = Tkfw_corr * map_keyframes[i]->T_kf_w;
        map_keyframes[i]->x_kf_w = logmap_se3(map_keyframes[i]->T_kf_w);
        // update landmarks
        for( auto it = map_points_kf_idx.at(i).begin(); it != map_points_kf_idx.at(i).end(); it++ )
        {
           if( map_points[(*it)] != NULL )
           {
               // update 3D position
               Vector3d point3D = map_points[(*it)]->point3D;
               map_points[(*it)]->point3D = Tkfw_corr.block(0,0,3,3) * point3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_points[(*it)]->med_obs_dir;
               map_points[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // [CHECK] update direction of each observation
               for( auto dir_it = map_points[(*it)]->dir_list.begin(); dir_it != map_points[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
        for( auto it = map_lines_kf_idx.at(i).begin(); it != map_lines_kf_idx.at(i).end(); it++ )
        {
           if( map_lines[(*it)] != NULL )
           {
               // update 3D position
               Vector3d sP3D = map_lines[(*it)]->line3D.head(3);
               Vector3d eP3D = map_lines[(*it)]->line3D.tail(3);
               map_lines[(*it)]->line3D.head(3) = Tkfw_corr.block(0,0,3,3) * sP3D + Tkfw_corr.block(0,3,3,1);
               map_lines[(*it)]->line3D.tail(3) = Tkfw_corr.block(0,0,3,3) * eP3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_lines[(*it)]->med_obs_dir;
               map_lines[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // [CHECK] update direction of each observation
               for( auto dir_it = map_lines[(*it)]->dir_list.begin(); dir_it != map_lines[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
    }

    // mark as optimized the lc_idx_list edges
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++)
        (*it)(2) = 0;

    // fuse local map from both sides of the loop and update graphs
    loopClosureFuseLandmarks();

    return true;

}

bool MapHandler::loopClosureOptimizationCovGraphG2O()
{

    // define G2O variables
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(true);
    g2o::BlockSolver_6_3::LinearSolverType* linearSolver;
    linearSolver = new g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>();
    g2o::BlockSolver_6_3* solver_ptr = new g2o::BlockSolver_6_3(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    solver->setUserLambdaInit(1e-10);
    optimizer.setAlgorithm(solver);

    // select min and max KF indices
    int kf_prev_idx =  2 * max_kf_idx;
    int kf_curr_idx = -1;
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++)
    {
        if( (*it)(0) < kf_prev_idx )
            kf_prev_idx = (*it)(0);
        if( (*it)(1) > kf_curr_idx )
            kf_curr_idx = (*it)(1);
    }
    kf_prev_idx = 0;

    // grab the KFs included in the optimization
    vector<int> kf_list;
    for( int i = kf_prev_idx; i <= kf_curr_idx; i++)
    {
        if( map_keyframes[i] != NULL )
        {
            // check if it is a LC vertex
            bool is_lc_i = false;
            bool is_lc_j = false;
            int id = 0;
            for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++, id++ )
            {
                if( (*it)(0) == i )
                {
                    is_lc_i = true;
                    break;
                }
                if( (*it)(1) == i )
                {
                    is_lc_j = true;
                    break;
                }
            }
            kf_list.push_back(i);
            // create SE3 vertex
            g2o::VertexSE3* v_se3 = new g2o::VertexSE3();
            v_se3->setId(i);
            v_se3->setMarginalized( false );
            if( is_lc_j )
            {
                // update pose of LC vertex
                v_se3->setFixed(false);
                v_se3->setEstimate( g2o::SE3Quat::exp( reverse_se3(logmap_se3( (expmap_se3(lc_pose_list[id])) * map_keyframes[lc_idx_list[id](0)]->T_kf_w )) ) );
            }
            else
            {
                v_se3->setEstimate( g2o::SE3Quat::exp( reverse_se3(map_keyframes[i]->x_kf_w) ) );
                //if( is_lc_i || i == 0 )
                if( i == 0 )
                    v_se3->setFixed(true);
            }
            optimizer.addVertex( v_se3 );
        }
    }

    // introduce edges
    for( int i = kf_prev_idx; i <= kf_curr_idx; i++ )
    {
        for( int j = i+1; j <= kf_curr_idx; j++ )
        {
            if( map_keyframes[i] != NULL && map_keyframes[j] != NULL &&
                ( full_graph[i][j] >= Config::minLMEssGraph() || full_graph[i][j] >= Config::minLMCovGraph() || abs(i-j) == 1  ) )
            {
                // kf2kf constraint
                Matrix4d T_ji_constraint = inverse_se3( map_keyframes[i]->T_kf_w ) * map_keyframes[j]->T_kf_w;
                // add edge
                g2o::EdgeSE3* e_se3 = new g2o::EdgeSE3();
                e_se3->setVertex( 0, optimizer.vertex(i) );
                e_se3->setVertex( 1, optimizer.vertex(j) );
                Vector6d x;
                x = reverse_se3(logmap_se3(T_ji_constraint) );
                e_se3->setMeasurement( g2o::SE3Quat::exp(x) );
                //e_se3->information() = map_keyframes[j]->xcov_kf_w;
                e_se3->setInformation( Matrix6d::Identity() );
                optimizer.addEdge( e_se3 );
            }
        }
    }

    // introduce loop closure edges
    int id = 0;
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++, id++ )
    {
        // add edge
        g2o::EdgeSE3* e_se3 = new g2o::EdgeSE3();
        e_se3->setVertex( 0, optimizer.vertex((*it)(0)) );
        e_se3->setVertex( 1, optimizer.vertex((*it)(1)) );
        Vector6d x;
        x = reverse_se3( lc_pose_list[id] );
        e_se3->setMeasurement( g2o::SE3Quat::exp(x) );
        //e_se3->information() = kf_Tij_cov_constraint[id];
        e_se3->information() = Matrix6d::Identity();
        optimizer.addEdge( e_se3 );
    }

    // optimize graph
    optimizer.initializeOptimization();
    optimizer.computeInitialGuess();
    optimizer.computeActiveErrors();
    optimizer.optimize(100);

    // recover pose and update map
    Matrix4d Tkfw_corr;
    for( auto kf_it = kf_list.begin(); kf_it != kf_list.end(); kf_it++)
    {
        g2o::VertexSE3* v_se3 = static_cast<g2o::VertexSE3*>(optimizer.vertex( (*kf_it) ));
        g2o::SE3Quat Tiw_corr =  v_se3->estimateAsSE3Quat();
        Vector6d x;
        Matrix4d Tkfw, Tkfw_prev;
        x = reverse_se3(Tiw_corr.log());
        Tkfw = expmap_se3( x );
        Tkfw_prev = map_keyframes[ (*kf_it) ]->T_kf_w;
        map_keyframes[ (*kf_it) ]->T_kf_w = Tkfw;
        map_keyframes[ (*kf_it) ]->x_kf_w = logmap_se3(Tkfw);
        // update map
        Tkfw_corr = Tkfw * inverse_se3( Tkfw_prev );
        for( auto it = map_points_kf_idx.at((*kf_it)).begin(); it != map_points_kf_idx.at((*kf_it)).end(); it++ )
        {
           if( map_points[(*it)] != NULL )
           {
               // update 3D position
               Vector3d point3D = map_points[(*it)]->point3D;
               map_points[(*it)]->point3D = Tkfw_corr.block(0,0,3,3) * point3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_points[(*it)]->med_obs_dir;
               map_points[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // [CHECK] update direction of each observation
               for( auto dir_it = map_points[(*it)]->dir_list.begin(); dir_it != map_points[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
        for( auto it = map_lines_kf_idx.at((*kf_it)).begin(); it != map_lines_kf_idx.at((*kf_it)).end(); it++ )
        {
           if( map_lines[(*it)] != NULL )
           {
               // update 3D position
               Vector3d sP3D = map_lines[(*it)]->line3D.head(3);
               Vector3d eP3D = map_lines[(*it)]->line3D.tail(3);
               map_lines[(*it)]->line3D.head(3) = Tkfw_corr.block(0,0,3,3) * sP3D + Tkfw_corr.block(0,3,3,1);
               map_lines[(*it)]->line3D.tail(3) = Tkfw_corr.block(0,0,3,3) * eP3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_lines[(*it)]->med_obs_dir;
               map_lines[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // [CHECK] update direction of each observation
               for( auto dir_it = map_lines[(*it)]->dir_list.begin(); dir_it != map_lines[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
    }

    // update pose and map of the rest of frames
    for( int i = kf_curr_idx + 1; i < map_keyframes.size(); i++ )
    {
        // update pose
        map_keyframes[i]->T_kf_w = Tkfw_corr * map_keyframes[i]->T_kf_w;
        map_keyframes[i]->x_kf_w = logmap_se3(map_keyframes[i]->T_kf_w);
        // update landmarks
        for( auto it = map_points_kf_idx.at(i).begin(); it != map_points_kf_idx.at(i).end(); it++ )
        {
           if( map_points[(*it)] != NULL )
           {
               // update 3D position
               Vector3d point3D = map_points[(*it)]->point3D;
               map_points[(*it)]->point3D = Tkfw_corr.block(0,0,3,3) * point3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_points[(*it)]->med_obs_dir;
               map_points[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // [CHECK] update direction of each observation
               for( auto dir_it = map_points[(*it)]->dir_list.begin(); dir_it != map_points[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
        for( auto it = map_lines_kf_idx.at(i).begin(); it != map_lines_kf_idx.at(i).end(); it++ )
        {
           if( map_lines[(*it)] != NULL )
           {
               // update 3D position
               Vector3d sP3D = map_lines[(*it)]->line3D.head(3);
               Vector3d eP3D = map_lines[(*it)]->line3D.tail(3);
               map_lines[(*it)]->line3D.head(3) = Tkfw_corr.block(0,0,3,3) * sP3D + Tkfw_corr.block(0,3,3,1);
               map_lines[(*it)]->line3D.tail(3) = Tkfw_corr.block(0,0,3,3) * eP3D + Tkfw_corr.block(0,3,3,1);
               // update direction of observation
               Vector3d obs_dir = map_lines[(*it)]->med_obs_dir;
               map_lines[(*it)]->med_obs_dir = Tkfw_corr.block(0,0,3,3) * obs_dir + Tkfw_corr.block(0,3,3,1);
               // [CHECK] update direction of each observation
               for( auto dir_it = map_lines[(*it)]->dir_list.begin(); dir_it != map_lines[(*it)]->dir_list.end(); dir_it++)
               {
                   Vector3d dir_list_ = (*dir_it);
                   (*dir_it) = Tkfw_corr.block(0,0,3,3) * dir_list_ + Tkfw_corr.block(0,3,3,1);
               }
           }
        }
    }

    // mark as optimized the lc_idx_list edges
    for( auto it = lc_idx_list.begin(); it != lc_idx_list.end(); it++)
        (*it)(2) = 0;

    // fuse local map from both sides of the loop and update graphs
    loopClosureFuseLandmarks();

    return true;

}

void MapHandler::loopClosureFuseLandmarks()
{

    // point matches
    int lc_idx = 0;
    for( auto idx_it = lc_pt_idxs.begin(); idx_it != lc_pt_idxs.end(); idx_it++, lc_idx++ )
    {
        if( lc_idx_list[lc_idx](2) == 1 )   // if not already optimized
        {
            int kf_prev_idx = lc_idx_list[lc_idx](0);
            int kf_curr_idx = lc_idx_list[lc_idx](1);
            /*if( map_keyframes[kf_prev_idx] == NULL || map_keyframes[kf_curr_idx] == NULL )
                continue;*/
            for( auto lm_it = (*idx_it).begin(); lm_it != (*idx_it).end(); lm_it++ )
            {
                // grab indices
                int lm_idx0  = (*lm_it)(0);
                int lm_ldx0  = (*lm_it)(1); // lr_qdx
                int lm_idx1  = (*lm_it)(2);
                int lm_ldx1  = (*lm_it)(3); // lr_tdx
                // if the LM exists just once, add observation
                if( lm_idx0 == -1 && lm_idx1 != -1 )
                {
                    if( map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0] != NULL && map_points[lm_idx1] != NULL )
                    {
                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->idx = lm_idx1;
                        Vector3d dir  = map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->P / map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->P.norm();
                        map_points[lm_idx1]->addMapPointObservation( map_keyframes[kf_prev_idx]->stereo_frame->pdesc_l.row(lm_ldx0), map_keyframes[kf_prev_idx]->kf_idx, map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->pl, dir );
                        // increase full graph for each KF that has already observed this LM
                        for( auto kf_it = map_points[lm_idx1]->kf_obs_list.begin(); kf_it != map_points[lm_idx1]->kf_obs_list.end(); kf_it++)
                        {
                            full_graph[(*kf_it)][kf_curr_idx]++;
                            full_graph[kf_curr_idx][(*kf_it)]++;
                        }
                    }
                }
                if( lm_idx0 != -1 && lm_idx1 == -1 )
                {
                    if( map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1] != NULL && map_points[lm_idx0] != NULL )
                    {
                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->idx = lm_idx0;
                        Vector3d dir  = map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->P / map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->P.norm();
                        map_points[lm_idx0]->addMapPointObservation( map_keyframes[kf_curr_idx]->stereo_frame->pdesc_l.row(lm_ldx1),
                                                                     map_keyframes[kf_curr_idx]->kf_idx,
                                                                     map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->pl, dir );
                        // increase full graph for each KF that has already observed this LM
                        for( auto kf_it = map_points[lm_idx0]->kf_obs_list.begin(); kf_it != map_points[lm_idx0]->kf_obs_list.end(); kf_it++)
                        {
                            full_graph[(*kf_it)][kf_prev_idx]++;
                            full_graph[kf_prev_idx][(*kf_it)]++;
                        }
                    }
                }
                // if not, create LM and add both observations
                if( lm_idx0 == -1 && lm_idx1 == -1 )
                {
                    if( map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0] != NULL && map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1] != NULL )
                    {
                        // assign indices
                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->idx = max_pt_idx;
                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->idx = max_pt_idx;
                        // create new 3D landmark with the observation from previous KF
                        Matrix4d Tfw = ( map_keyframes[kf_prev_idx]->T_kf_w );
                        Vector3d P3d = Tfw.block(0,0,3,3) * map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->P + Tfw.col(3).head(3);
                        //Vector3d dir = map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->P / map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->P.norm();
                        Vector3d dir = P3d / P3d.norm();
                        MapPoint* map_point = new MapPoint(max_pt_idx,P3d,map_keyframes[kf_prev_idx]->stereo_frame->pdesc_l.row(lm_ldx0),map_keyframes[kf_prev_idx]->kf_idx,map_keyframes[kf_prev_idx]->stereo_frame->stereo_pt[lm_ldx0]->pl,dir);
                        // add new 3D landmark to kf_idx where it was first observed
                        map_points_kf_idx.at( kf_prev_idx ).push_back( max_pt_idx );
                        // add observation of the 3D landmark from current KF
                        P3d = map_keyframes[kf_curr_idx]->T_kf_w.block(0,0,3,3) *  map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->P + map_keyframes[kf_curr_idx]->T_kf_w.col(3).head(3);
                        //dir = map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->P / map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->P.norm();
                        dir = P3d / P3d.norm();
                        map_point->addMapPointObservation(map_keyframes[kf_curr_idx]->stereo_frame->pdesc_l.row(lm_ldx1),map_keyframes[kf_curr_idx]->kf_idx,map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->pl,dir);
                        // add 3D landmark to map
                        map_points.push_back(map_point);
                        // update full graph (new feature)
                        max_pt_idx++;
                        full_graph[kf_prev_idx][kf_curr_idx]++;
                        full_graph[kf_curr_idx][kf_prev_idx]++;
                    }
                }
                // if the LM observed is different in each KF, then fuse them and erase the old one
                if( lm_idx0 != -1 && lm_idx1 != -1 )
                {
                    if( map_points[lm_idx0] != NULL && map_points[lm_idx1] != NULL
                        && map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1] != NULL )
                    {
                        int Nobs_lm_prev = map_points[lm_idx0]->kf_obs_list.size();
                        // fuse LMs while updating the full graph
                        int iter = 0;
                        for( auto it = map_points[lm_idx1]->desc_list.begin(); it != map_points[lm_idx1]->desc_list.end(); it++, iter++)
                        {
                            // concatenate desc, obs, dir, and kf_obs lists
                            map_points[lm_idx0]->desc_list.push_back(   (*it) );
                            map_points[lm_idx0]->obs_list.push_back(    map_points[lm_idx1]->obs_list[iter]    );
                            map_points[lm_idx0]->dir_list.push_back(    map_points[lm_idx1]->dir_list[iter]    );
                            map_points[lm_idx0]->kf_obs_list.push_back( map_points[lm_idx1]->kf_obs_list[iter] );
                            // update full graph
                            int jdx = map_points[lm_idx1]->kf_obs_list[iter];
                            for( int i = 0; i < Nobs_lm_prev; i++ )
                            {
                                int idx = map_points[lm_idx0]->kf_obs_list[i];
                                full_graph[idx][jdx]++;
                                full_graph[jdx][idx]++;
                            }
                            // update average descriptor and direction of observation
                            map_points[lm_idx0]->updateAverageDescDir();
                            // change idx in stereo_pt
                            map_keyframes[kf_curr_idx]->stereo_frame->stereo_pt[lm_ldx1]->idx = lm_idx0;
                        }
                        // remove from map_points_kf_idx
                        iter = 0;
                        int kf_lm_obs = map_points[lm_idx1]->kf_obs_list[0];
                        for( auto it = map_points_kf_idx.at(kf_lm_obs).begin(); it != map_points_kf_idx.at(kf_lm_obs).end(); it++, iter++)
                        {
                            if( (*it) == lm_idx1 )
                            {
                                map_points_kf_idx.at(kf_lm_obs).erase( map_points_kf_idx.at(kf_lm_obs).begin() + iter );
                                break;
                            }
                        }
                        // erase old landmark
                        map_points[lm_idx1] = NULL;
                    }
                }
            }
        }
    }

    // line segment matches
    lc_idx = 0;
    for( auto idx_it = lc_ls_idxs.begin(); idx_it != lc_ls_idxs.end(); idx_it++, lc_idx++ )
    {
        if( lc_idx_list[lc_idx](2) == 1 )   // if not already optimized
        {
            int kf_prev_idx = lc_idx_list[lc_idx](0);
            int kf_curr_idx = lc_idx_list[lc_idx](1);
            /*if( map_keyframes[kf_prev_idx] == NULL || map_keyframes[kf_curr_idx] == NULL )
                continue;*/
            for( auto lm_it = (*idx_it).begin(); lm_it != (*idx_it).end(); lm_it++ )
            {
                // grab indices
                int lm_idx0  = (*lm_it)(0);
                int lm_ldx0  = (*lm_it)(1); // lr_qdx
                int lm_idx1  = (*lm_it)(2);
                int lm_ldx1  = (*lm_it)(3); // lr_tdx
                cout << endl << (*lm_it).transpose() ;
                // if the LM exists just once, add observation
                if( lm_idx0 == -1 && lm_idx1 != -1 )
                {

                    cout << "\t" << map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls.size();

                    if( map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0] != NULL && map_lines[lm_idx1] != NULL )
                    {
                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->idx = lm_idx1;
                        Vector3d dir  = ( map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->sP + map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->eP )
                                      / ( map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->sP + map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->eP ).norm();
                        Vector4d pts;
                        pts.head(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->spl_obs;
                        pts.tail(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->epl_obs;
                        map_lines[lm_idx1]->addMapLineObservation( map_keyframes[kf_prev_idx]->stereo_frame->ldesc_l.row(lm_ldx0),
                                                                   map_keyframes[kf_prev_idx]->kf_idx,
                                                                   map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->le,
                                                                   dir, pts );
                        // increase full graph for each KF that has already observed this LM
                        for( auto kf_it = map_lines[lm_idx1]->kf_obs_list.begin(); kf_it != map_lines[lm_idx1]->kf_obs_list.end(); kf_it++)
                        {
                            full_graph[(*kf_it)][kf_curr_idx]++;
                            full_graph[kf_curr_idx][(*kf_it)]++;
                        }
                    }
                }
                if( lm_idx0 != -1 && lm_idx1 == -1 )
                {
                    if( map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1] != NULL && map_lines[lm_idx0] != NULL )
                    {

                        cout << "\t" << map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls.size();

                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->idx = lm_idx0;
                        Vector3d dir  = (map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->sP + map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->eP)
                                / (map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->sP+map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->eP).norm();
                        Vector4d pts;
                        pts.head(2) = map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->spl_obs;
                        pts.tail(2) = map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->epl_obs;
                        map_lines[lm_idx0]->addMapLineObservation( map_keyframes[kf_curr_idx]->stereo_frame->ldesc_l.row(lm_ldx1),
                                                                   map_keyframes[kf_curr_idx]->kf_idx,
                                                                   map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->le, dir, pts );
                        // increase full graph for each KF that has already observed this LM
                        for( auto kf_it = map_lines[lm_idx0]->kf_obs_list.begin(); kf_it != map_lines[lm_idx0]->kf_obs_list.end(); kf_it++)
                        {
                            full_graph[(*kf_it)][kf_prev_idx]++;
                            full_graph[kf_prev_idx][(*kf_it)]++;
                        }
                    }
                }
                // if not, create LM and add both observations
                if( lm_idx0 == -1 && lm_idx1 == -1 )
                {
                    if( map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0] != NULL && map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1] != NULL )
                    {
                        // assign indices
                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->idx = max_ls_idx;
                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->idx = max_ls_idx;
                        // create new 3D landmark with the observation from previous KF
                        Matrix4d Tfw  = ( map_keyframes[kf_prev_idx]->T_kf_w );
                        Vector3d sP3d = Tfw.block(0,0,3,3) * map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->sP + Tfw.col(3).head(3);
                        Vector3d eP3d = Tfw.block(0,0,3,3) * map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->eP + Tfw.col(3).head(3);
                        Vector3d mP3d = 0.5 * ( sP3d + eP3d );
                        Vector3d dir = mP3d / mP3d.norm();
                        Vector6d L3D; L3D.head(3) = sP3d; L3D.tail(3) = eP3d;
                        Vector4d pts;
                        pts.head(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->spl;
                        pts.tail(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->epl;
                        MapLine* map_line = new MapLine(max_ls_idx,L3D,map_keyframes[kf_prev_idx]->stereo_frame->ldesc_l.row(lm_ldx0),map_keyframes[kf_prev_idx]->kf_idx,
                                                        map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx0]->le,dir,pts);
                        // add new 3D landmark to kf_idx where it was first observed
                        map_lines_kf_idx.at( kf_prev_idx ).push_back( max_ls_idx );
                        // add observation of the 3D landmark from current KF
                        sP3d = map_keyframes[kf_curr_idx]->T_kf_w.block(0,0,3,3) *  map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->sP + map_keyframes[kf_curr_idx]->T_kf_w.col(3).head(3);
                        eP3d = map_keyframes[kf_curr_idx]->T_kf_w.block(0,0,3,3) *  map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->eP + map_keyframes[kf_curr_idx]->T_kf_w.col(3).head(3);
                        mP3d = 0.5 * ( sP3d + eP3d );
                        dir = mP3d / mP3d.norm();
                        pts.head(2) = map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->spl;
                        pts.tail(2) = map_keyframes[kf_prev_idx]->stereo_frame->stereo_ls[lm_ldx1]->epl;
                        map_line->addMapLineObservation(map_keyframes[kf_curr_idx]->stereo_frame->ldesc_l.row(lm_ldx1),map_keyframes[kf_curr_idx]->kf_idx,
                                                        map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->le,dir,pts);
                        // add 3D landmark to map
                        map_lines.push_back(map_line);
                        // update full graph (new feature)
                        max_ls_idx++;
                        full_graph[kf_prev_idx][kf_curr_idx]++;
                        full_graph[kf_curr_idx][kf_prev_idx]++;
                    }
                }
                // if the LM observed is different in each KF, then fuse them and erase the old one
                if( lm_idx0 != -1 && lm_idx1 != -1 )
                {
                    if( map_lines[lm_idx0] != NULL && map_lines[lm_idx1] != NULL
                        && map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1] != NULL )
                    {
                        int Nobs_lm_prev = map_lines[lm_idx0]->kf_obs_list.size();
                        // fuse LMs while updating the full graph
                        int iter = 0;
                        for( auto it = map_lines[lm_idx1]->desc_list.begin(); it != map_lines[lm_idx1]->desc_list.end(); it++, iter++)
                        {
                            // concatenate desc, obs, dir, pts, and kf_obs lists
                            map_lines[lm_idx0]->desc_list.push_back(   (*it) );
                            map_lines[lm_idx0]->obs_list.push_back(    map_lines[lm_idx1]->obs_list[iter]    );
                            map_lines[lm_idx0]->dir_list.push_back(    map_lines[lm_idx1]->dir_list[iter]    );
                            map_lines[lm_idx0]->pts_list.push_back(    map_lines[lm_idx1]->pts_list[iter]    );
                            map_lines[lm_idx0]->kf_obs_list.push_back( map_lines[lm_idx1]->kf_obs_list[iter] );
                            // update full graph
                            int jdx = map_lines[lm_idx1]->kf_obs_list[iter];
                            for( int i = 0; i < Nobs_lm_prev; i++ )
                            {
                                int idx = map_lines[lm_idx0]->kf_obs_list[i];
                                full_graph[idx][jdx]++;
                                full_graph[jdx][idx]++;
                            }
                            // update average descriptor and direction of observation
                            map_lines[lm_idx0]->updateAverageDescDir();
                            // change idx in stereo_ls
                            map_keyframes[kf_curr_idx]->stereo_frame->stereo_ls[lm_ldx1]->idx = lm_idx0;
                        }
                        // remove from map_points_kf_idx
                        iter = 0;
                        int kf_lm_obs = map_lines[lm_idx1]->kf_obs_list[0];
                        for( auto it = map_lines_kf_idx.at(kf_lm_obs).begin(); it != map_lines_kf_idx.at(kf_lm_obs).end(); it++, iter++)
                        {
                            if( (*it) == lm_idx1 )
                            {
                                map_lines_kf_idx.at(kf_lm_obs).erase( map_lines_kf_idx.at(kf_lm_obs).begin() + iter );
                                break;
                            }
                        }
                        // erase old landmark
                        map_lines[lm_idx1] = NULL;
                    }
                }
            }
        }
    }

}

}

