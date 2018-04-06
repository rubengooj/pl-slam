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

#include "keyFrame.h"

#include <iostream>

#include <opencv2/imgproc.hpp>

namespace PLSLAM{

KeyFrame::KeyFrame( const StereoFrame* sf )
{
    kf_idx    = -1;
    x_kf_w    = logmap_se3( T_kf_w );

    T_kf_w    = sf->Tfw ;
    xcov_kf_w = sf->Tfw_cov;

    stereo_frame = new StereoFrame( sf->img_l, sf->img_r, kf_idx, sf->cam );
    stereo_frame->pdesc_l   = sf->pdesc_l;
    stereo_frame->pdesc_r   = sf->pdesc_r;
    stereo_frame->ldesc_l   = sf->ldesc_l;
    stereo_frame->ldesc_r   = sf->ldesc_r;

    std::vector<PointFeature*> &stereo_pt = stereo_frame->stereo_pt;
    stereo_pt.resize(sf->stereo_pt.size());
    for (int idx = 0; idx < stereo_pt.size(); ++idx)
        stereo_pt[idx] = sf->stereo_pt[idx]->safeCopy();

    std::vector<LineFeature*> &stereo_ls = stereo_frame->stereo_ls;
    stereo_ls.resize(sf->stereo_ls.size());
    for (int idx = 0; idx < stereo_ls.size(); ++idx)
        stereo_ls[idx] = sf->stereo_ls[idx]->safeCopy();
}

KeyFrame::KeyFrame( const StereoFrame* sf, int kf_idx_ )
{
    kf_idx    = kf_idx_;
    T_kf_w    = sf->Tfw;
    x_kf_w    = logmap_se3( T_kf_w );
    xcov_kf_w = sf->Tfw_cov;

    stereo_frame = new StereoFrame( sf->img_l, sf->img_r, kf_idx, sf->cam );
    stereo_frame->pdesc_l   = sf->pdesc_l;
    stereo_frame->pdesc_r   = sf->pdesc_r;
    stereo_frame->ldesc_l   = sf->ldesc_l;
    stereo_frame->ldesc_r   = sf->ldesc_r;

    std::vector<PointFeature*> &stereo_pt = stereo_frame->stereo_pt;
    stereo_pt.resize(sf->stereo_pt.size());
    for (int idx = 0; idx < stereo_pt.size(); ++idx)
        stereo_pt[idx] = sf->stereo_pt[idx]->safeCopy();

    std::vector<LineFeature*> &stereo_ls = stereo_frame->stereo_ls;
    stereo_ls.resize(sf->stereo_ls.size());
    for (int idx = 0; idx < stereo_ls.size(); ++idx)
        stereo_ls[idx] = sf->stereo_ls[idx]->safeCopy();
}

KeyFrame::~KeyFrame() {

    delete stereo_frame;
}

Mat KeyFrame::plotKeyFrame()
{
    // create new image to modify it
    Mat img_l_aux;
    stereo_frame->img_l.copyTo( img_l_aux );
    if( img_l_aux.channels() == 1 )
        cvtColor(img_l_aux, img_l_aux, CV_GRAY2BGR, 3);
    else if (img_l_aux.channels() == 4)
        cvtColor(img_l_aux, img_l_aux, CV_BGRA2BGR, 3);
    else if (img_l_aux.channels() != 3)
        throw std::runtime_error(std::string("[KeyFrame->plotKeyFrame] unsupported image format: ") +
                                 std::to_string(img_l_aux.channels()));
    img_l_aux.convertTo(img_l_aux, CV_8UC3);

    // Variables
    unsigned int    r = 0, g, b = 0;
    Point2f         p,q;
    double          thick = 1.5;
    int             k = 0, radius  = 3;
    // plot point features
    for( vector<PointFeature*>::iterator pt_it = stereo_frame->stereo_pt.begin(); pt_it != stereo_frame->stereo_pt.end(); pt_it++)
    {
        if( (*pt_it)->idx != -1 )
        {
            g = 200;
            p = cv::Point( int((*pt_it)->pl(0)), int((*pt_it)->pl(1)) );
            circle( img_l_aux, p, radius, Scalar(b,g,r), thick);
        }
    }
    // plot line segment features
    for( vector<LineFeature*>::iterator ls_it = stereo_frame->stereo_ls.begin(); ls_it != stereo_frame->stereo_ls.end(); ls_it++)
    {
        if( (*ls_it)->idx != -1 )
        {
            g = 200;
            p = cv::Point( int((*ls_it)->spl(0)), int((*ls_it)->spl(1)) );
            q = cv::Point( int((*ls_it)->epl(0)), int((*ls_it)->epl(1)) );
            line( img_l_aux, p, q, Scalar(b,g,r), thick);
        }
    }

    return img_l_aux;
}

}
