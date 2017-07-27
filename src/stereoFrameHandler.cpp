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

#include <stereoFrameHandler.h>

namespace StVO{

StereoFrameHandler::StereoFrameHandler( PinholeStereoCamera *cam_ ) : cam(cam_) {}

StereoFrameHandler::~StereoFrameHandler(){}

void StereoFrameHandler::initialize(const Mat img_l_, const Mat img_r_ , const int idx_)
{
    prev_frame = new StereoFrame( img_l_, img_r_, idx_, cam );
    prev_frame->extractInitialStereoFeatures();
    prev_frame->Tfw = Matrix4d::Identity();
    prev_frame->Tfw_cov = Matrix6d::Identity();
    prev_frame->DT  = Matrix4d::Identity();
    // variables for adaptative FAST
    orb_fast_th = Config::orbFastTh();
    // SLAM variables for KF decision
    T_prevKF         = Matrix4d::Identity();
    cov_prevKF_currF = Matrix6d::Zero();
    prev_f_iskf      = true;
}

void StereoFrameHandler::insertStereoPair(const Mat img_l_, const Mat img_r_ , const int idx_)
{
    curr_frame = new StereoFrame( img_l_, img_r_, idx_, cam );
    curr_frame->extractStereoFeatures( orb_fast_th );
    f2fTracking();
}

void StereoFrameHandler::f2fTracking()
{

    // points f2f tracking
    matched_pt.clear();
    if( Config::hasPoints() && !(curr_frame->stereo_pt.size()==0) && !(prev_frame->stereo_pt.size()==0)  )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        Mat pdesc_l1, pdesc_l2;
        vector<vector<DMatch>> pmatches_12, pmatches_21;
        // 12 and 21 matches
        pdesc_l1 = prev_frame->pdesc_l;
        pdesc_l2 = curr_frame->pdesc_l;        
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StereoFrame::matchPointFeatures, prev_frame, bfm, pdesc_l1, pdesc_l2, ref(pmatches_12) );
                auto match_r = async( launch::async, &StereoFrame::matchPointFeatures, prev_frame, bfm, pdesc_l2, pdesc_l1, ref(pmatches_21) );
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
            if( lr_qdx == rl_tdx  && dist_12 > nn12_dist_th )
            {
                PointFeature* point_ = prev_frame->stereo_pt[lr_qdx];
                point_->pl_obs = curr_frame->stereo_pt[lr_tdx]->pl;
                point_->inlier = true;
                matched_pt.push_back( point_ );                
                curr_frame->stereo_pt[lr_tdx]->idx = prev_frame->stereo_pt[lr_qdx]->idx; // prev idx
            }
        }
    }

    // line segments f2f tracking
    matched_ls.clear();
    if( Config::hasLines() && !(curr_frame->stereo_ls.size()==0) && !(prev_frame->stereo_ls.size()==0)  )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );    // cross-check
        Mat ldesc_l1, ldesc_l2;
        vector<vector<DMatch>> lmatches_12, lmatches_21;
        // 12 and 21 matches
        ldesc_l1 = prev_frame->ldesc_l;
        ldesc_l2 = curr_frame->ldesc_l;
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StereoFrame::matchLineFeatures, prev_frame, bfm, ldesc_l1, ldesc_l2, ref(lmatches_12) );
                auto match_r = async( launch::async, &StereoFrame::matchLineFeatures, prev_frame, bfm, ldesc_l2, ldesc_l1, ref(lmatches_21) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( ldesc_l1,ldesc_l2, lmatches_12, 2);
                bfm->knnMatch( ldesc_l2,ldesc_l1, lmatches_21, 2);
            }
        }
        else
            bfm->knnMatch( ldesc_l1,ldesc_l2, lmatches_12, 2);
        // sort matches by the distance between the best and second best matches
        double nn_dist_th, nn12_dist_th;
        curr_frame->lineDescriptorMAD(lmatches_12,nn_dist_th, nn12_dist_th);
        nn12_dist_th  = nn12_dist_th * Config::descThL();
        // resort according to the queryIdx
        sort( lmatches_12.begin(), lmatches_12.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( lmatches_21.begin(), lmatches_21.end(), sort_descriptor_by_queryIdx() );
        // bucle around pmatches
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
            if( lr_qdx == rl_tdx  && dist_12 > nn12_dist_th )
            {
                LineFeature* line_ = prev_frame->stereo_ls[lr_qdx];
                line_->sdisp_obs = curr_frame->stereo_ls[lr_tdx]->sdisp;
                line_->edisp_obs = curr_frame->stereo_ls[lr_tdx]->edisp;
                line_->spl_obs = curr_frame->stereo_ls[lr_tdx]->spl;
                line_->epl_obs = curr_frame->stereo_ls[lr_tdx]->epl;
                line_->le_obs  = curr_frame->stereo_ls[lr_tdx]->le;               
                line_->inlier  = true;
                matched_ls.push_back( line_ );
                curr_frame->stereo_ls[lr_tdx]->idx = prev_frame->stereo_ls[lr_qdx]->idx; // prev idx
            }
        }
    }

    n_inliers_pt = matched_pt.size();
    n_inliers_ls = matched_ls.size();
    n_inliers    = n_inliers_pt + n_inliers_ls;

}

void StereoFrameHandler::updateFrame()
{

    // update FAST threshold for the keypoint detection
    if( Config::adaptativeFAST() )
    {
        int min_fast  = Config::fastMinTh();
        int max_fast  = Config::fastMaxTh();
        int fast_inc  = Config::fastIncTh();
        int feat_th   = Config::fastFeatTh();
        float err_th  = Config::fastErrTh();

        // if bad optimization, -= 2*fast_inc
        if( curr_frame->DT == Matrix4d::Identity() || curr_frame->err_norm > err_th )
            orb_fast_th = std::max( min_fast, orb_fast_th - 2*fast_inc );
        // elif number of features ...
        else if( n_inliers_pt < feat_th )
            orb_fast_th = std::max( min_fast, orb_fast_th - 2*fast_inc );
        else if( n_inliers < feat_th * 2 )
            orb_fast_th = std::max( min_fast, orb_fast_th - fast_inc );
        else if( n_inliers > feat_th * 3 )
            orb_fast_th = std::min( max_fast, orb_fast_th + fast_inc );
        else if( n_inliers > feat_th * 4 )
            orb_fast_th = std::min( max_fast, orb_fast_th + 2*fast_inc );
    }

    // clean and update variables
    matched_pt.clear();
    matched_ls.clear();
    delete prev_frame;
    prev_frame = curr_frame;
    curr_frame = NULL;

}

void StereoFrameHandler::optimizePose()
{

    // definitions
    Matrix6d DT_cov;
    Matrix4d DT, DT_;
    Vector6d DT_cov_eig;
    double   err;

    // set init pose
    DT     = prev_frame->DT;
    DT_cov = prev_frame->DT_cov;

    // solver
    if( n_inliers > Config::minFeatures() )
    {
        // optimize
        DT_ = DT;
        gaussNewtonOptimization(DT_,DT_cov,err,Config::maxIters());
        // remove outliers (implement some logic based on the covariance's eigenvalues and optim error)
        if( is_finite(DT_) )
        {
            removeOutliers(DT_);
            // refine without outliers
            if( n_inliers > Config::minFeatures() )
                gaussNewtonOptimization(DT,DT_cov,err,Config::maxItersRef());
            else
            {
                DT     = Matrix4d::Identity();
                DT_cov = Matrix6d::Zero();
            }
        }
        else
        {
            DT     = Matrix4d::Identity();
            DT_cov = Matrix6d::Zero();
        }
    }
    else
    {
        DT     = Matrix4d::Identity();
        DT_cov = Matrix6d::Zero();
    }

    // set estimated pose
    if( is_finite(DT) )
    {
        curr_frame->DT     = inverse_se3( DT );
        curr_frame->Tfw    = prev_frame->Tfw * curr_frame->DT;
        curr_frame->DT_cov = DT_cov;
        SelfAdjointEigenSolver<Matrix6d> eigensolver(DT_cov);
        curr_frame->DT_cov_eig = eigensolver.eigenvalues();
        curr_frame->Tfw_cov = unccomp_se3( prev_frame->Tfw, prev_frame->Tfw_cov, DT_cov );
        curr_frame->err_norm   = err;
    }
    else
    {
        curr_frame->DT     = Matrix4d::Identity();
        curr_frame->Tfw    = prev_frame->Tfw;
        curr_frame->Tfw_cov= prev_frame->Tfw_cov;
        curr_frame->DT_cov = DT_cov;
        SelfAdjointEigenSolver<Matrix6d> eigensolver(DT_cov);
        curr_frame->DT_cov_eig = eigensolver.eigenvalues();
        curr_frame->err_norm   = -1.0;
    }

}

void StereoFrameHandler::optimizePose(Matrix4d DT_ini)
{

    // definitions
    Matrix6d DT_cov;
    Matrix4d DT, DT_;
    Vector6d DT_cov_eig;
    double   err;

    // set init pose
    DT     = DT_ini;
    DT_cov = prev_frame->DT_cov;

    // solver
    if( n_inliers > Config::minFeatures() )
    {
        // optimize
        DT_ = DT;
        gaussNewtonOptimization(DT_,DT_cov,err,Config::maxIters());
        // remove outliers (implement some logic based on the covariance's eigenvalues and optim error)
        if( is_finite(DT_) )
        {
            removeOutliers(DT_);
            // refine without outliers
            if( n_inliers > Config::minFeatures() )
                gaussNewtonOptimization(DT,DT_cov,err,Config::maxItersRef());
            else
            {
                DT     = Matrix4d::Identity();
                DT_cov = Matrix6d::Zero();
            }
        }
        else
        {
            DT     = Matrix4d::Identity();
            DT_cov = Matrix6d::Zero();
        }
    }
    else
    {
        DT     = Matrix4d::Identity();
        DT_cov = Matrix6d::Zero();
    }

    // set estimated pose
    if( is_finite(DT) )
    {
        curr_frame->DT     = inverse_se3( DT );
        curr_frame->Tfw    = prev_frame->Tfw * curr_frame->DT;
        curr_frame->DT_cov = DT_cov;
        SelfAdjointEigenSolver<Matrix6d> eigensolver(DT_cov);
        curr_frame->DT_cov_eig = eigensolver.eigenvalues();
        curr_frame->Tfw_cov = unccomp_se3( prev_frame->Tfw, prev_frame->Tfw_cov, DT_cov );
        curr_frame->err_norm   = err;
    }
    else
    {
        curr_frame->DT     = Matrix4d::Identity();
        curr_frame->Tfw    = prev_frame->Tfw;
        curr_frame->Tfw_cov= prev_frame->Tfw_cov;
        curr_frame->DT_cov = DT_cov;
        SelfAdjointEigenSolver<Matrix6d> eigensolver(DT_cov);
        curr_frame->DT_cov_eig = eigensolver.eigenvalues();
        curr_frame->err_norm   = -1.0;
    }




}

void StereoFrameHandler::gaussNewtonOptimization(Matrix4d &DT, Matrix6d &DT_cov, double &err_, int max_iters)
{
    Matrix6d H;
    Vector6d g, DT_inc;
    double err, err_prev = 999999999.9;
    for( int iters = 0; iters < max_iters; iters++)
    {
        // estimate hessian and gradient (select)
        optimizeFunctions( DT, H, g, err );
        // if the difference is very small stop
        if( ( abs(err-err_prev) < Config::minErrorChange() ) || ( err < Config::minError()) )
            break;
        // update step
        LDLT<Matrix6d> solver(H);
        DT_inc = solver.solve(g);
        DT  << DT * inverse_se3( expmap_se3(DT_inc) );
        // if the parameter change is small stop (TODO: change with two parameters, one for R and another one for t)
        if( DT_inc.norm() < numeric_limits<double>::epsilon() )
            break;
        // update previous values
        err_prev = err;
    }
    DT_cov = H.inverse();
    err_   = err;
}

void StereoFrameHandler::removeOutliers(Matrix4d DT)
{

    vector<double> res_p, res_l, ove_l;

    // point features
    int iter = 0;
    for( list<PointFeature*>::iterator it = matched_pt.begin(); it!=matched_pt.end(); it++, iter++)
    {
        // projection error
        Vector3d P_ = DT.block(0,0,3,3) * (*it)->P + DT.col(3).head(3);
        Vector2d pl_proj = cam->projection( P_ );
        res_p.push_back( ( pl_proj - (*it)->pl_obs ).norm() * sqrt((*it)->sigma2) );
    }

    // line segment features
    for( list<LineFeature*>::iterator it = matched_ls.begin(); it!=matched_ls.end(); it++, iter++)
    {
        // projection error
        Vector3d sP_ = DT.block(0,0,3,3) * (*it)->sP + DT.col(3).head(3);
        Vector3d eP_ = DT.block(0,0,3,3) * (*it)->eP + DT.col(3).head(3);
        Vector2d spl_proj = cam->projection( sP_ );
        Vector2d epl_proj = cam->projection( eP_ );
        Vector3d l_obs    = (*it)->le_obs;
        Vector2d err_li;
        err_li(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
        err_li(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
        res_l.push_back( err_li.norm() * sqrt((*it)->sigma2) );
    }

    // estimate mad standard deviation
    double inlier_th_p =  Config::inlierK() * vector_stdv_mad( res_p );
    double inlier_th_l =  Config::inlierK() * vector_stdv_mad( res_l );
    //inlier_th_p = sqrt(7.815);
    //inlier_th_l = sqrt(7.815);

    // filter outliers
    iter = 0;
    for( list<PointFeature*>::iterator it = matched_pt.begin(); it!=matched_pt.end(); it++, iter++)
    {
        if( res_p[iter] > inlier_th_p )
        {
            (*it)->inlier = false;
            n_inliers--;
            n_inliers_pt--;
        }
    }
    iter = 0;
    for( list<LineFeature*>::iterator it = matched_ls.begin(); it!=matched_ls.end(); it++, iter++)
    {
        if( res_l[iter] > inlier_th_l )
        {
            (*it)->inlier = false;
            n_inliers--;
            n_inliers_ls--;
        }
    }

}

void StereoFrameHandler::optimizeFunctions(Matrix4d DT, Matrix6d &H, Vector6d &g, double &e )
{

    // define hessians, gradients, and residuals
    Matrix6d H_l, H_p;
    Vector6d g_l, g_p;
    double   e_l = 0.0, e_p = 0.0, S_l, S_p;
    H   = Matrix6d::Zero(); H_l = H; H_p = H;
    g   = Vector6d::Zero(); g_l = g; g_p = g;
    e   = 0.0;

    // point features
    int N_p = 0;
    for( list<PointFeature*>::iterator it = matched_pt.begin(); it!=matched_pt.end(); it++)
    {
        if( (*it)->inlier )
        {
            Vector3d P_ = DT.block(0,0,3,3) * (*it)->P + DT.col(3).head(3);
            Vector2d pl_proj = cam->projection( P_ );
            // projection error
            Vector2d err_i    = pl_proj - (*it)->pl_obs;
            double err_i_norm = err_i.norm();
            // estimate variables for J, H, and g
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
            double w  = 1.0;
            double s2 = (*it)->sigma2;
            w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
            // update hessian, gradient, and error
            H_p += J_aux * J_aux.transpose() * w;
            g_p += J_aux * err_i_norm * w;
            e_p += err_i_norm * err_i_norm * w;
            N_p++;
        }
    }

    // line segment features
    int N_l = 0;
    for( list<LineFeature*>::iterator it = matched_ls.begin(); it!=matched_ls.end(); it++)
    {
        if( (*it)->inlier )
        {
            Vector3d sP_ = DT.block(0,0,3,3) * (*it)->sP + DT.col(3).head(3);
            Vector2d spl_proj = cam->projection( sP_ );
            Vector3d eP_ = DT.block(0,0,3,3) * (*it)->eP + DT.col(3).head(3);
            Vector2d epl_proj = cam->projection( eP_ );
            Vector3d l_obs = (*it)->le_obs;
            // projection error
            Vector2d err_i;
            err_i(0) = l_obs(0) * spl_proj(0) + l_obs(1) * spl_proj(1) + l_obs(2);
            err_i(1) = l_obs(0) * epl_proj(0) + l_obs(1) * epl_proj(1) + l_obs(2);
            double err_i_norm = err_i.norm();
            // estimate variables for J, H, and g
            // -- start point
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
            // -- end point
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
            double w  = 1.0;
            double s2 = (*it)->sigma2;
            w = 1.0 / ( 1.0 + err_i_norm * err_i_norm * s2 );
            // estimating overlap between line segments
            /*bool has_overlap = false;
            double overlap = 1.0;
            if( has_overlap )
                overlap = lineSegmentOverlap( (*it)->spl, (*it)->epl, spl_proj, epl_proj );
            w *= overlap;*/
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

}

double StereoFrameHandler::lineSegmentOverlap( Vector2d spl_obs, Vector2d epl_obs, Vector2d spl_proj, Vector2d epl_proj  )
{

    Vector2d l = spl_obs - epl_obs;
    double lxx  = l(0)*l(0);
    double lyy  = l(1)*l(1);
    double lxy  = l(0)*l(1);
    double lxy2 = 1.f/(lxx+lyy);

    Matrix2d u;
    u << lxx, lxy, lxy, lyy;
    u = u * lxy2;
    Vector2d v;
    Matrix2d vm;
    vm << lxx, -lxy, -lxy, lyy;
    v = lxy2 * vm * spl_obs;

    Vector2d sp;
    sp << u * spl_proj + v;
    Vector2d ep;
    ep << u * epl_proj + v;

    double lnorm  = 1.f / l.norm();
    double seno   = -l(0)*lnorm;
    double coseno = -l(1)*lnorm;

    Matrix2d rot; rot << coseno, -seno, seno, coseno;
    Vector2d sl, el;

    sl     << rot * spl_obs;
    el     << rot * epl_obs;
    sp     << rot * sp;
    ep     << rot * ep;

    double sln    = min(sl(1), el(1));
    double eln    = max(sl(1), el(1));
    double spn    = min(sp(1), ep(1));
    double epn    = max(sp(1), ep(1));

    double length = eln-spn;

    double overlap;
    if ( (epn < sln) || (spn > eln) )
        overlap = 0.f;
    else{
        if ( (epn>eln) && (spn<sln) )
            overlap = eln-sln;
        else
            overlap = min(eln,epn) - max(sln,spn);
    }

    if(length>0.01f)
        overlap = overlap / length;
    else
        overlap = 0.f;

    return overlap;

}

/*  slam functions  */

bool StereoFrameHandler::needNewKF()
{

    // if the previous KF was a KF, update the entropy_first_prevKF value
    if( prev_f_iskf )
    {
        entropy_first_prevKF = 3.0*(1.0+log(2.0*acos(-1))) + 0.5*log( curr_frame->DT_cov.determinant() );
        prev_f_iskf = false;
    }

    // check geometric distances from previous KF
    Matrix4d DT = inverse_se3( curr_frame->Tfw ) * T_prevKF;
    Vector6d dX = logmap_se3( DT );
    double t = dX.head(3).norm();
    double r = dX.tail(3).norm() * 180.f / CV_PI;

    // check cumulated covariance from previous KF
    Matrix6d adjTprevkf = adjoint_se3( T_prevKF );
    Matrix6d covDTinv   = uncTinv_se3( curr_frame->DT, curr_frame->DT_cov );
    cov_prevKF_currF += adjTprevkf * covDTinv * adjTprevkf.transpose();
    double entropy_curr  = 3.0*(1.0+log(2.0*acos(-1))) + 0.5*log( cov_prevKF_currF.determinant() );
    double entropy_ratio = entropy_curr / entropy_first_prevKF;

    // decide if a new KF is needed
    if( entropy_ratio < Config::minEntropyRatio() || std::isnan(entropy_ratio) || std::isinf(entropy_ratio) ||
        ( curr_frame->DT_cov == Matrix6d::Zero() && curr_frame->DT == Matrix4d::Identity() ) )
    {
        cout << endl << "Entropy ratio: " << entropy_ratio   << endl;
        return true;
    }
    else
    {
        cout << endl << "No new KF needed" << endl << endl;
        return false;
    }

}

void StereoFrameHandler::currFrameIsKF()
{

    // restart point indices
    int idx_pt = 0;
    for( vector<PointFeature*>::iterator it = curr_frame->stereo_pt.begin(); it != curr_frame->stereo_pt.end(); it++)
    {
        (*it)->idx = idx_pt;
        idx_pt++;
    }

    // restart line indices
    int idx_ls = 0;
    for( vector<LineFeature*>::iterator it = curr_frame->stereo_ls.begin(); it != curr_frame->stereo_ls.end(); it++)
    {
        (*it)->idx = idx_ls;
        idx_ls++;
    }

    // update KF
    curr_frame->Tfw     = Matrix4d::Identity();
    curr_frame->Tfw_cov = Matrix6d::Identity();

    // update SLAM variables for KF decision
    T_prevKF = curr_frame->Tfw;
    cov_prevKF_currF = Matrix6d::Zero();
    prev_f_iskf = true;

}

}
