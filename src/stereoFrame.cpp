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

#include <stereoFrame.h>

namespace StVO{

StereoFrame::StereoFrame(){}

StereoFrame::StereoFrame(const Mat img_l_, const Mat img_r_ , const int idx_, PinholeStereoCamera *cam_) :
    img_l(img_l_), img_r(img_r_), frame_idx(idx_), cam(cam_) {}

StereoFrame::~StereoFrame(){}

void StereoFrame::extractInitialStereoFeatures( int fast_th )
{

    // Feature detection and description
    vector<KeyPoint> points_l, points_r;
    vector<KeyLine>  lines_l, lines_r;
    double min_line_length_th = Config::minLineLength() * std::min( cam->getWidth(), cam->getHeight() );
    if( Config::lrInParallel() )
    {
        auto detect_l = async(launch::async, &StereoFrame::detectFeatures, this, img_l, ref(points_l), ref(pdesc_l), ref(lines_l), ref(ldesc_l), min_line_length_th, fast_th  );
        auto detect_r = async(launch::async, &StereoFrame::detectFeatures, this, img_r, ref(points_r), ref(pdesc_r), ref(lines_r), ref(ldesc_r), min_line_length_th, fast_th );
        detect_l.wait();
        detect_r.wait();
    }
    else
    {
        detectFeatures(img_l,points_l,pdesc_l,lines_l,ldesc_l,min_line_length_th,fast_th);
        detectFeatures(img_r,points_r,pdesc_r,lines_r,ldesc_r,min_line_length_th,fast_th);
    }

    // Points stereo matching
    if( Config::hasPoints() && !(points_l.size()==0) && !(points_r.size()==0) )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );
        vector<vector<DMatch>> pmatches_lr, pmatches_rl, pmatches_lr_;
        Mat pdesc_l_;
        stereo_pt.clear();
        // LR and RL matches
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StereoFrame::matchPointFeatures, this, bfm, pdesc_l, pdesc_r, ref(pmatches_lr) );
                auto match_r = async( launch::async, &StereoFrame::matchPointFeatures, this, bfm, pdesc_r, pdesc_l, ref(pmatches_rl) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( pdesc_l, pdesc_r, pmatches_lr, 2);
                bfm->knnMatch( pdesc_r, pdesc_l, pmatches_rl, 2);
            }
        }
        else
            bfm->knnMatch( pdesc_l, pdesc_r, pmatches_lr, 2);
        // sort matches by the distance between the best and second best matches
        double nn12_dist_th  = Config::minRatio12P();
        // resort according to the queryIdx
        sort( pmatches_lr.begin(), pmatches_lr.end(), sort_descriptor_by_queryIdx() );
        if(Config::bestLRMatches())
            sort( pmatches_rl.begin(), pmatches_rl.end(), sort_descriptor_by_queryIdx() );
        // bucle around pmatches
        int pt_idx = 0;
        for( int i = 0; i < pmatches_lr.size(); i++ )
        {
            int lr_qdx, lr_tdx, rl_tdx;
            lr_qdx = pmatches_lr[i][0].queryIdx;
            lr_tdx = pmatches_lr[i][0].trainIdx;
            if( Config::bestLRMatches() )
            {
                // check if they are mutual best matches
                rl_tdx = pmatches_rl[lr_tdx][0].trainIdx;
            }
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_12 = pmatches_lr[i][0].distance / pmatches_lr[i][1].distance;
            if( lr_qdx == rl_tdx  && dist_12 > nn12_dist_th )
            {
                // check stereo epipolar constraint
                if( fabsf( points_l[lr_qdx].pt.y-points_r[lr_tdx].pt.y) <= Config::maxDistEpip() )
                {
                    // check minimal disparity
                    double disp_ = points_l[lr_qdx].pt.x - points_r[lr_tdx].pt.x;
                    if( disp_ >= Config::minDisp() ){
                        pdesc_l_.push_back( pdesc_l.row(lr_qdx) );
                        PointFeature* point_;
                        Vector2d pl_; pl_ << points_l[lr_qdx].pt.x, points_l[lr_qdx].pt.y;
                        Vector3d P_;  P_ = cam->backProjection( pl_(0), pl_(1), disp_);
                        stereo_pt.push_back( new PointFeature(pl_,disp_,P_,pt_idx,points_l[lr_qdx].octave) );
                        pt_idx++;
                    }
                }
            }
        }
        pdesc_l_.copyTo(pdesc_l);
    }

    // Line segments stereo matching
    if( Config::hasLines() && !lines_l.empty() && !lines_r.empty() )
    {
        stereo_ls.clear();
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );
        vector<vector<DMatch>> lmatches_lr, lmatches_rl;
        Mat ldesc_l_;
        // LR and RL matches
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StereoFrame::matchLineFeatures, this, bfm, ldesc_l, ldesc_r, ref(lmatches_lr) );
                auto match_r = async( launch::async, &StereoFrame::matchLineFeatures, this, bfm, ldesc_r, ldesc_l, ref(lmatches_rl) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( ldesc_l,ldesc_r, lmatches_lr, 2);
                bfm->knnMatch( ldesc_r,ldesc_l, lmatches_rl, 2);
            }
        }
        else
            bfm->knnMatch( ldesc_l,ldesc_r, lmatches_lr, 2);
        // sort matches by the distance between the best and second best matches
        double nn_dist_th, nn12_dist_th;
        lineDescriptorMAD(lmatches_lr,nn_dist_th, nn12_dist_th);        
        nn12_dist_th  = nn12_dist_th * Config::descThL();
        // bucle around pmatches
        sort( lmatches_lr.begin(), lmatches_lr.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( lmatches_rl.begin(), lmatches_rl.end(), sort_descriptor_by_queryIdx() );
        int n_matches;
        if( Config::bestLRMatches() )
            n_matches = min(lmatches_lr.size(),lmatches_rl.size());
        else
            n_matches = lmatches_lr.size();
        int ls_idx = 0;
        for( int i = 0; i < n_matches; i++ )
        {
            // check if they are mutual best matches ( if bestLRMatches() )
            int lr_qdx = lmatches_lr[i][0].queryIdx;
            int lr_tdx = lmatches_lr[i][0].trainIdx;
            int rl_tdx;
            if( Config::bestLRMatches() )
                rl_tdx = lmatches_rl[lr_tdx][0].trainIdx;
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_12 = lmatches_lr[i][1].distance - lmatches_lr[i][0].distance;
            double length  = lines_r[lr_tdx].lineLength;
            if( lr_qdx == rl_tdx && dist_12 > nn12_dist_th )
            {
                // estimate the disparity of the endpoints
                Vector3d sp_l; sp_l << lines_l[lr_qdx].startPointX, lines_l[lr_qdx].startPointY, 1.0;
                Vector3d ep_l; ep_l << lines_l[lr_qdx].endPointX,   lines_l[lr_qdx].endPointY,   1.0;
                Vector3d le_l; le_l << sp_l.cross(ep_l); le_l = le_l / sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) );
                Vector3d sp_r; sp_r << lines_r[lr_tdx].startPointX, lines_r[lr_tdx].startPointY, 1.0;
                Vector3d ep_r; ep_r << lines_r[lr_tdx].endPointX,   lines_r[lr_tdx].endPointY,   1.0;
                Vector3d le_r; le_r << sp_r.cross(ep_r);                
                double overlap = lineSegmentOverlapStereo( sp_l(1), ep_l(1), sp_r(1), ep_r(1) );
                sp_r << - (le_r(2)+le_r(1)*lines_l[lr_qdx].startPointY )/le_r(0) , lines_l[lr_qdx].startPointY ,  1.0;
                ep_r << - (le_r(2)+le_r(1)*lines_l[lr_qdx].endPointY   )/le_r(0) , lines_l[lr_qdx].endPointY ,    1.0;
                double disp_s = lines_l[lr_qdx].startPointX - sp_r(0);
                double disp_e = lines_l[lr_qdx].endPointX   - ep_r(0);
                // check minimal disparity
                if( disp_s >= Config::minDisp() && disp_e >= Config::minDisp()
                        && fabsf(le_r(0)) > Config::lineHorizTh()
                        && overlap > Config::stereoOverlapTh() )
                {
                    ldesc_l_.push_back( ldesc_l.row(lr_qdx) );
                    Vector3d sP_; sP_ = cam->backProjection( sp_l(0), sp_l(1), disp_s);
                    Vector3d eP_; eP_ = cam->backProjection( ep_l(0), ep_l(1), disp_e);
                    double angle_l = lines_l[lr_qdx].angle;
                    stereo_ls.push_back( new LineFeature(Vector2d(sp_l(0),sp_l(1)),disp_s,sP_,Vector2d(ep_l(0),ep_l(1)),
                                                         disp_e,eP_,le_l,angle_l,ls_idx, lines_l[lr_qdx].octave) );
                    ls_idx++;
                }
            }
        }
        ldesc_l_.copyTo(ldesc_l);
    }

}

void StereoFrame::extractStereoFeatures( int fast_th)
{

    // Feature detection and description
    vector<KeyPoint> points_l, points_r;
    vector<KeyLine>  lines_l, lines_r;
    double min_line_length_th = Config::minLineLength() * std::min( cam->getWidth(), cam->getHeight() ) ;
    if( Config::lrInParallel() )
    {
        auto detect_l = async(launch::async, &StereoFrame::detectFeatures, this, img_l, ref(points_l), ref(pdesc_l), ref(lines_l), ref(ldesc_l), min_line_length_th, fast_th );
        auto detect_r = async(launch::async, &StereoFrame::detectFeatures, this, img_r, ref(points_r), ref(pdesc_r), ref(lines_r), ref(ldesc_r), min_line_length_th, fast_th );
        detect_l.wait();
        detect_r.wait();
    }
    else
    {
        detectFeatures(img_l,points_l,pdesc_l,lines_l,ldesc_l,min_line_length_th,fast_th);
        detectFeatures(img_r,points_r,pdesc_r,lines_r,ldesc_r,min_line_length_th,fast_th);
    }

    // Points stereo matching
    if( Config::hasPoints() && !(points_l.size()==0) && !(points_r.size()==0) )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );
        vector<vector<DMatch>> pmatches_lr, pmatches_rl, pmatches_lr_;
        Mat pdesc_l_;
        stereo_pt.clear();
        // LR and RL matches
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StereoFrame::matchPointFeatures, this, bfm, pdesc_l, pdesc_r, ref(pmatches_lr) );
                auto match_r = async( launch::async, &StereoFrame::matchPointFeatures, this, bfm, pdesc_r, pdesc_l, ref(pmatches_rl) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( pdesc_l, pdesc_r, pmatches_lr, 2);
                bfm->knnMatch( pdesc_r, pdesc_l, pmatches_rl, 2);
            }
        }
        else
            bfm->knnMatch( pdesc_l, pdesc_r, pmatches_lr, 2);
        // sort matches by the distance between the best and second best matches
        double nn12_dist_th  = Config::minRatio12P();
        // resort according to the queryIdx
        sort( pmatches_lr.begin(), pmatches_lr.end(), sort_descriptor_by_queryIdx() );
        if(Config::bestLRMatches())
            sort( pmatches_rl.begin(), pmatches_rl.end(), sort_descriptor_by_queryIdx() );
        // bucle around pmatches
        for( int i = 0; i < pmatches_lr.size(); i++ )
        {
            int lr_qdx, lr_tdx, rl_tdx;
            lr_qdx = pmatches_lr[i][0].queryIdx;
            lr_tdx = pmatches_lr[i][0].trainIdx;
            if( Config::bestLRMatches() )
            {
                // check if they are mutual best matches
                rl_tdx = pmatches_rl[lr_tdx][0].trainIdx;
            }
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_12 = pmatches_lr[i][0].distance / pmatches_lr[i][1].distance;
            if( lr_qdx == rl_tdx  && dist_12 > nn12_dist_th )
            {
                // check stereo epipolar constraint
                if( fabsf( points_l[lr_qdx].pt.y-points_r[lr_tdx].pt.y) <= Config::maxDistEpip() )
                {
                    // check minimal disparity
                    double disp_ = points_l[lr_qdx].pt.x - points_r[lr_tdx].pt.x;
                    if( disp_ >= Config::minDisp() ){
                        PointFeature* point_;
                        Vector2d pl_; pl_ << points_l[lr_qdx].pt.x, points_l[lr_qdx].pt.y;
                        Vector3d P_;  P_ = cam->backProjection( pl_(0), pl_(1), disp_);
                        pdesc_l_.push_back( pdesc_l.row(lr_qdx) );
                        stereo_pt.push_back( new PointFeature(pl_,disp_,P_,-1, points_l[lr_qdx].octave) );
                    }
                }
            }
        }
        pdesc_l_.copyTo(pdesc_l);
    }

    // Line segments stereo matching
    if( Config::hasLines() && !lines_l.empty() && !lines_r.empty() )
    {
        stereo_ls.clear();
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING2, false );
        vector<vector<DMatch>> lmatches_lr, lmatches_rl;
        Mat ldesc_l_;
        // LR and RL matches
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StereoFrame::matchLineFeatures, this, bfm, ldesc_l, ldesc_r, ref(lmatches_lr) );
                auto match_r = async( launch::async, &StereoFrame::matchLineFeatures, this, bfm, ldesc_r, ldesc_l, ref(lmatches_rl) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( ldesc_l,ldesc_r, lmatches_lr, 2);
                bfm->knnMatch( ldesc_r,ldesc_l, lmatches_rl, 2);
            }
        }
        else
            bfm->knnMatch( ldesc_l,ldesc_r, lmatches_lr, 2);
        // sort matches by the distance between the best and second best matches
        double nn_dist_th, nn12_dist_th;
        lineDescriptorMAD(lmatches_lr,nn_dist_th, nn12_dist_th);
        nn12_dist_th  = nn12_dist_th * Config::descThL();
        // bucle around pmatches
        sort( lmatches_lr.begin(), lmatches_lr.end(), sort_descriptor_by_queryIdx() );
        if( Config::bestLRMatches() )
            sort( lmatches_rl.begin(), lmatches_rl.end(), sort_descriptor_by_queryIdx() );
        int n_matches;
        if( Config::bestLRMatches() )
            n_matches = min(lmatches_lr.size(),lmatches_rl.size());
        else
            n_matches = lmatches_lr.size();
        for( int i = 0; i < n_matches; i++ )
        {
            // check if they are mutual best matches ( if bestLRMatches() )
            int lr_qdx = lmatches_lr[i][0].queryIdx;
            int lr_tdx = lmatches_lr[i][0].trainIdx;
            int rl_tdx;
            if( Config::bestLRMatches() )
                rl_tdx = lmatches_rl[lr_tdx][0].trainIdx;
            else
                rl_tdx = lr_qdx;
            // check if they are mutual best matches and the minimum distance
            double dist_12 = lmatches_lr[i][1].distance - lmatches_lr[i][0].distance;
            double length  = lines_r[lr_tdx].lineLength;
            if( lr_qdx == rl_tdx && dist_12 > nn12_dist_th )
            {
                // estimate the disparity of the endpoints
                Vector3d sp_l; sp_l << lines_l[lr_qdx].startPointX, lines_l[lr_qdx].startPointY, 1.0;
                Vector3d ep_l; ep_l << lines_l[lr_qdx].endPointX,   lines_l[lr_qdx].endPointY,   1.0;
                Vector3d le_l; le_l << sp_l.cross(ep_l); le_l = le_l / sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) );
                Vector3d sp_r; sp_r << lines_r[lr_tdx].startPointX, lines_r[lr_tdx].startPointY, 1.0;
                Vector3d ep_r; ep_r << lines_r[lr_tdx].endPointX,   lines_r[lr_tdx].endPointY,   1.0;
                Vector3d le_r; le_r << sp_r.cross(ep_r);
                double overlap = lineSegmentOverlapStereo( sp_l(1), ep_l(1), sp_r(1), ep_r(1) );
                sp_r << - (le_r(2)+le_r(1)*lines_l[lr_qdx].startPointY )/le_r(0) , lines_l[lr_qdx].startPointY ,  1.0;
                ep_r << - (le_r(2)+le_r(1)*lines_l[lr_qdx].endPointY   )/le_r(0) , lines_l[lr_qdx].endPointY ,    1.0;
                double disp_s = lines_l[lr_qdx].startPointX - sp_r(0);
                double disp_e = lines_l[lr_qdx].endPointX   - ep_r(0);
                // check minimal disparity
                if( disp_s >= Config::minDisp() && disp_e >= Config::minDisp()
                        && fabsf(le_r(0)) > Config::lineHorizTh()
                        && overlap > Config::stereoOverlapTh() )
                {
                    Vector3d sP_; sP_ = cam->backProjection( sp_l(0), sp_l(1), disp_s);
                    Vector3d eP_; eP_ = cam->backProjection( ep_l(0), ep_l(1), disp_e);
                    double angle_l = lines_l[lr_qdx].angle;
                    //----------------- DEBUG: 24/05/2016 ----------------------
                    // estimate the uncertainty of the endpoints
                    double cx = cam->getCx();
                    double cy = cam->getCy();
                    double f  = cam->getFx();
                    // - start point
                    double px_hat = sp_l(0) - cx;
                    double py_hat = sp_l(1) - cy;
                    double disp   = disp_s;
                    double disp2  = disp * disp;
                    Matrix3d covS_an;
                    covS_an(0,0) = disp2+2.f*px_hat*px_hat;
                    covS_an(0,1) = 2.f*px_hat*py_hat;
                    covS_an(0,2) = 2.f*f*px_hat;
                    covS_an(1,1) = disp2+2.f*py_hat*py_hat;
                    covS_an(1,2) = 2.f*f*py_hat;
                    covS_an(2,2) = 2.f*f*f;
                    covS_an(1,0) = covS_an(0,1);
                    covS_an(2,0) = covS_an(0,2);
                    covS_an(2,1) = covS_an(1,2);
                    covS_an << covS_an * cam->getB() * cam->getB() / (disp2*disp2);
                    // - end point
                    px_hat = ep_l(0) - cam->getCx();
                    py_hat = ep_l(1) - cam->getCy();
                    disp   = disp_e;
                    disp2  = disp * disp;
                    Matrix3d covE_an;
                    covE_an(0,0) = disp2+2.f*px_hat*px_hat;
                    covE_an(0,1) = 2.f*px_hat*py_hat;
                    covE_an(0,2) = 2.f*f*px_hat;
                    covE_an(1,1) = disp2+2.f*py_hat*py_hat;
                    covE_an(1,2) = 2.f*f*py_hat;
                    covE_an(2,2) = 2.f*f*f;
                    covE_an(1,0) = covE_an(0,1);
                    covE_an(2,0) = covE_an(0,2);
                    covE_an(2,1) = covE_an(1,2);
                    covE_an << covE_an * cam->getB() * cam->getB() / (disp2*disp2);
                    // - estimate eigenvalues
                    Vector3d S_eigen, E_eigen;
                    SelfAdjointEigenSolver<Matrix3d> eigensolver_s(covS_an);
                    S_eigen = eigensolver_s.eigenvalues();
                    SelfAdjointEigenSolver<Matrix3d> eigensolver_e(covE_an);
                    E_eigen = eigensolver_e.eigenvalues();
                    double max_eig = max( S_eigen(2),E_eigen(2) );
                    // - dbg plot
                    if(max_eig < Config::lineCovTh())
                    {
                        ldesc_l_.push_back( ldesc_l.row(lr_qdx) );
                        stereo_ls.push_back( new LineFeature(Vector2d(sp_l(0),sp_l(1)),disp_s,sP_,
                                                             Vector2d(ep_l(0),ep_l(1)),disp_e,eP_,
                                                             le_l,angle_l,-1,lines_l[lr_qdx].octave) );
                    }
                    //----------------------------------------------------------
                }
            }
        }
        ldesc_l_.copyTo(ldesc_l);
    }

}

void StereoFrame::detectFeatures(Mat img, vector<KeyPoint> &points, Mat &pdesc, vector<KeyLine> &lines, Mat &ldesc, double min_line_length, int fast_th)
{

    if( Config::plInParallel() && Config::hasPoints() && Config::hasLines() )
    {
        auto detect_l = async(launch::async, &StereoFrame::detectLineFeatures,  this, img, ref(lines),  ref(ldesc), min_line_length );
        auto detect_p = async(launch::async, &StereoFrame::detectPointFeatures, this, img, ref(points), ref(pdesc), fast_th );
        detect_l.wait();
        detect_p.wait();
    }
    else
    {
        // Detect point features
        if( Config::hasPoints() )
        {
            int fast_th_ = Config::orbFastTh();
            if( fast_th != 0 )
                fast_th_ = fast_th;
            Ptr<ORB> orb = ORB::create( Config::orbNFeatures(), Config::orbScaleFactor(), Config::orbNLevels(),
                                        Config::orbEdgeTh(), 0, Config::orbWtaK(), Config::orbScore(),
                                        Config::orbPatchSize(), fast_th_ );
            orb->detectAndCompute( img, Mat(), points, pdesc, false);
        }
        // Detect line features
        Ptr<BinaryDescriptor>   lbd = BinaryDescriptor::createBinaryDescriptor();
        if( Config::hasLines() )
        {
            Ptr<line_descriptor::LSDDetectorC> lsd = line_descriptor::LSDDetectorC::createLSDDetectorC();
            // lsd parameters
            line_descriptor::LSDDetectorC::LSDOptions opts;
            opts.refine       = Config::lsdRefine();
            opts.scale        = Config::lsdScale();
            opts.sigma_scale  = Config::lsdSigmaScale();
            opts.quant        = Config::lsdQuant();
            opts.ang_th       = Config::lsdAngTh();
            opts.log_eps      = Config::lsdLogEps();
            opts.density_th   = Config::lsdDensityTh();
            opts.n_bins       = Config::lsdNBins();
            opts.min_length   = min_line_length;
            lsd->detect( img, lines, Config::lsdScale(), 1, opts);
            // filter lines
            if( lines.size()>Config::lsdNFeatures() && Config::lsdNFeatures()!=0  )
            {
                // sort lines by their response
                sort( lines.begin(), lines.end(), sort_lines_by_response() );
                //sort( lines.begin(), lines.end(), sort_lines_by_length() );
                lines.resize(Config::lsdNFeatures());
                // reassign index
                for( int i = 0; i < Config::lsdNFeatures(); i++  )
                    lines[i].class_id = i;
            }
            lbd->compute( img, lines, ldesc);
        }
    }

}

void StereoFrame::detectPointFeatures( Mat img, vector<KeyPoint> &points, Mat &pdesc, int fast_th )
{
    // Detect point features
    if( Config::hasPoints() )
    {
        int fast_th_ = Config::orbFastTh();
        if( fast_th != 0 )
            fast_th_ = fast_th;
        Ptr<ORB> orb = ORB::create( Config::orbNFeatures(), Config::orbScaleFactor(), Config::orbNLevels(),
                                    Config::orbEdgeTh(), 0, Config::orbWtaK(), Config::orbScore(),
                                    Config::orbPatchSize(), fast_th_ );
        orb->detectAndCompute( img, Mat(), points, pdesc, false);
    }
}

void StereoFrame::detectLineFeatures( Mat img, vector<KeyLine> &lines, Mat &ldesc, double min_line_length )
{
    // Detect line features
    Ptr<BinaryDescriptor>   lbd = BinaryDescriptor::createBinaryDescriptor();
    if( Config::hasLines() )
    {
        Ptr<line_descriptor::LSDDetectorC> lsd = line_descriptor::LSDDetectorC::createLSDDetectorC();
        // lsd parameters
        line_descriptor::LSDDetectorC::LSDOptions opts;
        opts.refine       = Config::lsdRefine();
        opts.scale        = Config::lsdScale();
        opts.sigma_scale  = Config::lsdSigmaScale();
        opts.quant        = Config::lsdQuant();
        opts.ang_th       = Config::lsdAngTh();
        opts.log_eps      = Config::lsdLogEps();
        opts.density_th   = Config::lsdDensityTh();
        opts.n_bins       = Config::lsdNBins();
        opts.min_length   = min_line_length;
        lsd->detect( img, lines, Config::lsdScale(), 1, opts);
        // filter lines
        if( lines.size()>Config::lsdNFeatures() && Config::lsdNFeatures()!=0  )
        {
            // sort lines by their response
            sort( lines.begin(), lines.end(), sort_lines_by_response() );
            //sort( lines.begin(), lines.end(), sort_lines_by_length() );
            lines.resize(Config::lsdNFeatures());
            // reassign index
            for( int i = 0; i < Config::lsdNFeatures(); i++  )
                lines[i].class_id = i;
        }
        lbd->compute( img, lines, ldesc);
    }
}

void StereoFrame::matchPointFeatures(BFMatcher* bfm, Mat pdesc_1, Mat pdesc_2, vector<vector<DMatch>> &pmatches_12  )
{
    bfm->knnMatch( pdesc_1, pdesc_2, pmatches_12, 2);
}

void StereoFrame::matchLineFeatures(BFMatcher* bfm, Mat ldesc_1, Mat ldesc_2, vector<vector<DMatch>> &lmatches_12  )
{
    bfm->knnMatch( ldesc_1, ldesc_2, lmatches_12, 2);
}

void StereoFrame::pointDescriptorMAD( const vector<vector<DMatch>> matches, double &nn_mad, double &nn12_mad )
{

    vector<vector<DMatch>> matches_nn, matches_12;
    matches_nn = matches;
    matches_12 = matches;

    // estimate the NN's distance standard deviation
    double nn_dist_median;
    sort( matches_nn.begin(), matches_nn.end(), compare_descriptor_by_NN_dist() );
    nn_mad = matches_nn[int(matches_nn.size()/2)][0].distance;
    for( int j = 0; j < matches_nn.size(); j++)
        matches_nn[j][0].distance = fabsf( matches_nn[j][0].distance - nn_dist_median );
    sort( matches_nn.begin(), matches_nn.end(), compare_descriptor_by_NN_dist() );
    nn_mad = 1.4826 * matches_nn[int(matches_nn.size()/2)][0].distance;

    // estimate the NN's 12 distance standard deviation
    double nn12_dist_median;
    sort( matches_12.begin(), matches_12.end(), compare_descriptor_by_NN12_ratio() );
    nn_dist_median = matches_12[int(matches_12.size()/2)][0].distance / matches_12[int(matches_12.size()/2)][1].distance;
    for( int j = 0; j < matches_12.size(); j++)
        matches_12[j][0].distance = fabsf( matches_12[j][0].distance / matches_12[j][1].distance - nn_dist_median );
    sort( matches_12.begin(), matches_12.end(), compare_descriptor_by_NN_dist() );
    nn12_mad =  1.4826 * matches_12[int(matches_12.size()/2)][0].distance;

}

void StereoFrame::lineDescriptorMAD( const vector<vector<DMatch>> matches, double &nn_mad, double &nn12_mad )
{

    vector<vector<DMatch>> matches_nn, matches_12;
    matches_nn = matches;
    matches_12 = matches;

    // estimate the NN's distance standard deviation
    double nn_dist_median;
    sort( matches_nn.begin(), matches_nn.end(), compare_descriptor_by_NN_dist() );
    nn_mad = matches_nn[int(matches_nn.size()/2)][0].distance;
    for( int j = 0; j < matches_nn.size(); j++)
        matches_nn[j][0].distance = fabsf( matches_nn[j][0].distance - nn_dist_median );
    sort( matches_nn.begin(), matches_nn.end(), compare_descriptor_by_NN_dist() );
    nn_mad = 1.4826 * matches_nn[int(matches_nn.size()/2)][0].distance;

    // estimate the NN's 12 distance standard deviation
    double nn12_dist_median;
    sort( matches_12.begin(), matches_12.end(), compare_descriptor_by_NN12_dist() );
    nn12_mad = matches_12[int(matches_12.size()/2)][1].distance - matches_12[int(matches_12.size()/2)][0].distance;
    for( int j = 0; j < matches_12.size(); j++)
        matches_12[j][0].distance = fabsf( matches_12[j][1].distance - matches_12[j][0].distance - nn_dist_median );
    sort( matches_12.begin(), matches_12.end(), compare_descriptor_by_NN_dist() );
    nn12_mad =  1.4826 * matches_12[int(matches_12.size()/2)][0].distance;

}

double StereoFrame::lineSegmentOverlapStereo( double spl_obs, double epl_obs, double spl_proj, double epl_proj  )
{

    double sln    = min(spl_obs,  epl_obs);
    double eln    = max(spl_obs,  epl_obs);
    double spn    = min(spl_proj, epl_proj);
    double epn    = max(spl_proj, epl_proj);

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

Mat StereoFrame::plotStereoFrame()
{

    // create new image to modify it
    Mat img_l_aux;
    img_l.copyTo( img_l_aux );
    if( img_l_aux.channels() == 1 )
        cvtColor(img_l_aux,img_l_aux,CV_GRAY2BGR);

    // Variables
    unsigned int    r = 0, g, b = 0;
    Point2f         p,q;
    double          thick = 1.5;
    int             k = 0, radius  = 3;

    // plot point features
    for( vector<PointFeature*>::iterator pt_it = stereo_pt.begin(); pt_it != stereo_pt.end(); pt_it++)
    {
        if( (*pt_it)->inlier )
        {
            g = 200;
            p = cv::Point( int((*pt_it)->pl(0)), int((*pt_it)->pl(1)) );
            circle( img_l_aux, p, radius, Scalar(b,g,r), thick);
        }
    }

    // plot line segment features
    for( vector<LineFeature*>::iterator ls_it = stereo_ls.begin(); ls_it != stereo_ls.end(); ls_it++)
    {
        if( (*ls_it)->inlier )
        {
            g = 200;
            p = cv::Point( int((*ls_it)->spl(0)), int((*ls_it)->spl(1)) );
            q = cv::Point( int((*ls_it)->epl(0)), int((*ls_it)->epl(1)) );
            line( img_l_aux, p, q, Scalar(b,g,r), thick);
        }
    }

    return img_l_aux;

}

Mat StereoFrame::plotStereoMatches()
{

    // create new image to modify it
    Mat img_l_aux, img_r_aux, img_aux;
    img_l.copyTo( img_l_aux );
    img_r.copyTo( img_r_aux );
    if( img_l_aux.channels() == 1 )
        cvtColor(img_l_aux,img_l_aux,CV_GRAY2BGR);
    if( img_r_aux.channels() == 1 )
        cvtColor(img_r_aux,img_r_aux,CV_GRAY2BGR);

    // Variables
    unsigned int    r, g, b;
    Point2f         p,q;
    double          thick = 1.5;
    int             k = 0, radius  = 3;

    // Feature detection and description
    vector<KeyPoint> points_l, points_r;
    vector<KeyLine>  lines_l, lines_r;
    Mat pdesc_l_, pdesc_r_, ldesc_l_, ldesc_r_;

    double min_line_length_th = Config::minLineLength() * std::min( cam->getWidth(), cam->getHeight() ) ;
    if( Config::lrInParallel() )
    {
        auto detect_l = async(launch::async, &StereoFrame::detectFeatures, this, img_l, ref(points_l), ref(pdesc_l_), ref(lines_l), ref(ldesc_l_), min_line_length_th, 0 );
        auto detect_r = async(launch::async, &StereoFrame::detectFeatures, this, img_r, ref(points_r), ref(pdesc_r_), ref(lines_r), ref(ldesc_r_), min_line_length_th, 0 );
        detect_l.wait();
        detect_r.wait();
    }
    else
    {
        detectFeatures(img_l,points_l,pdesc_l_,lines_l,ldesc_l_,min_line_length_th,0);
        detectFeatures(img_r,points_r,pdesc_r_,lines_r,ldesc_r_,min_line_length_th,0);
    }

    // Points stereo matching and plot
    if( Config::hasPoints() && !(points_l.size()==0) && !(points_r.size()==0) )
    {
        BFMatcher* bfm = new BFMatcher( NORM_HAMMING, false );
        vector<vector<DMatch>> pmatches_lr, pmatches_rl, pmatches_lr_;
        Mat pdesc_l_;
        stereo_pt.clear();
        // LR and RL matches
        if( Config::bestLRMatches() )
        {
            if( Config::lrInParallel() )
            {
                auto match_l = async( launch::async, &StereoFrame::matchPointFeatures, this, bfm, pdesc_l, pdesc_r, ref(pmatches_lr) );
                auto match_r = async( launch::async, &StereoFrame::matchPointFeatures, this, bfm, pdesc_r, pdesc_l, ref(pmatches_rl) );
                match_l.wait();
                match_r.wait();
            }
            else
            {
                bfm->knnMatch( pdesc_l, pdesc_r, pmatches_lr, 2);
                bfm->knnMatch( pdesc_r, pdesc_l, pmatches_rl, 2);
            }
        }
        else
            bfm->knnMatch( pdesc_l, pdesc_r, pmatches_lr, 2);

        // sort matches by the distance between the best and second best matches
        double nn12_dist_th  = Config::minRatio12P();

        // resort according to the queryIdx
        sort( pmatches_lr.begin(), pmatches_lr.end(), sort_descriptor_by_queryIdx() );
        if(Config::bestLRMatches())
            sort( pmatches_rl.begin(), pmatches_rl.end(), sort_descriptor_by_queryIdx() );

        // bucle around pmatches
        for( int i = 0; i < pmatches_lr.size(); i++ )
        {
            int lr_qdx, lr_tdx, rl_tdx;
            lr_qdx = pmatches_lr[i][0].queryIdx;
            lr_tdx = pmatches_lr[i][0].trainIdx;
            double disp_ = points_l[lr_qdx].pt.x - points_r[lr_tdx].pt.x;
        }

        // plot stereo matches
        cv::drawMatches(img_l_aux, points_l, img_r_aux, points_r, pmatches_lr, img_aux);

    }

    return img_aux;

}

}
