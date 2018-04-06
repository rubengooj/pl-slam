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

#include "mapFeatures.h"

namespace PLSLAM{

// Point features

MapPoint::MapPoint(int idx_, Vector3d point3D_, Mat desc_, int kf_obs_, Vector2d obs_, Vector3d dir_, double sigma2_ ) :
    idx(idx_), point3D(point3D_), inlier(true)
{
    desc_list.push_back( desc_ );
    obs_list.push_back( obs_ );
    kf_obs_list.push_back( kf_obs_ );
    dir_list.push_back( dir_ );
    sigma_list.push_back(sigma2_);
    med_obs_dir = dir_;
    med_desc    = desc_;
}

void MapPoint::addMapPointObservation(Mat desc_, int kf_obs_, Vector2d obs_, Vector3d dir_,  double sigma2_ )
{
    desc_list.push_back( desc_ );
    obs_list.push_back( obs_ );
    kf_obs_list.push_back( kf_obs_ );
    dir_list.push_back( dir_ );
    sigma_list.push_back(sigma2_);
    updateAverageDescDir();
}

void MapPoint::updateAverageDescDir()
{

    // descriptor
    // - check distances between all the observed descriptors
    int n = desc_list.size();
    MatrixXf conf_desc(n,n);
    for(int i = 0; i < n; i++ )
    {
        conf_desc(i,i) = 0;
        for(int j = i+1 ; j < n; j++ )
        {
            int d = norm(desc_list[i],desc_list[j],NORM_HAMMING);
            conf_desc(i,j) = d;
            conf_desc(j,i) = d;
        }
    }

    // - select the one with least mean distance to the rest
    int max_dist = 99999;
    int max_idx  = 0;
    for(int i = 0; i < n; i++)
    {
        vector<int> dist_idx;
        for(int j = 0; j < n; j++)
            dist_idx.push_back( conf_desc(i,j) );
        sort( dist_idx.begin(), dist_idx.end() );
        int idx_median = dist_idx[ int(1+0.5*(n-1)) ];
        if( idx_median < max_dist )
        {
            max_dist = idx_median;
            max_idx  = i;
        }
    }
    med_desc = desc_list[max_idx];

    // direction
    Vector3d med_dir;
    for(int i = 0; i < n; i++)
        med_dir += dir_list[i];
    med_obs_dir = med_dir / n;

}

// Line segment features

MapLine::MapLine(int idx_, Vector6d line3D_, Mat desc_, int kf_obs_, Vector3d obs_, Vector3d dir_, Vector4d pts_, double sigma2_) :
    idx(idx_), line3D(line3D_), inlier(true)
{
    desc_list.push_back( desc_ );
    obs_list.push_back( obs_ );
    kf_obs_list.push_back( kf_obs_ );
    dir_list.push_back( dir_ );
    pts_list.push_back(pts_);
    sigma_list.push_back(sigma2_);
    med_obs_dir = dir_;
    med_desc    = desc_;
}

void MapLine::addMapLineObservation(Mat desc_, int kf_obs_, Vector3d obs_, Vector3d dir_, Vector4d pts_, double sigma2_)
{
    desc_list.push_back( desc_ );
    obs_list.push_back( obs_ );
    kf_obs_list.push_back( kf_obs_ );
    dir_list.push_back( dir_ );
    pts_list.push_back(pts_);
    sigma_list.push_back(sigma2_);
    updateAverageDescDir();
}

void MapLine::updateAverageDescDir()
{

    // descriptor
    // - check distances between all the observed descriptors
    int n = desc_list.size();
    MatrixXf conf_desc(n,n);
    for(int i = 0; i < n; i++ )
    {
        conf_desc(i,i) = 0;
        for(int j = i+1 ; j < n; j++ )
        {
            int d = norm(desc_list[i],desc_list[j],NORM_HAMMING);
            conf_desc(i,j) = d;
            conf_desc(j,i) = d;
        }
    }

    // - select the one with least mean distance to the rest
    int max_dist = 99999;
    int max_idx  = 0;
    for(int i = 0; i < n; i++)
    {
        vector<int> dist_idx;
        for(int j = 0; j < n; j++)
            dist_idx.push_back( conf_desc(i,j) );
        sort( dist_idx.begin(), dist_idx.end() );
        int idx_median = dist_idx[ int(1+0.5*(n-1)) ];
        if( idx_median < max_dist )
        {
            max_dist = idx_median;
            max_idx  = i;
        }
    }
    med_desc = desc_list[max_idx];

    // direction
    Vector3d med_dir;
    for(int i = 0; i < n; i++)
        med_dir += dir_list[i];
    med_obs_dir = med_dir / n;

}

}
