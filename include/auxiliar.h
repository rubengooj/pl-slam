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

#include <iostream>

#include <cv.h>
#include <opencv2/features2d/features2d.hpp>
#include <line_descriptor_custom.hpp>
#include <line_descriptor/descriptor_custom.hpp>
using namespace cv;
using namespace line_descriptor;

#include <vector>
using namespace std;

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;
typedef Matrix<double,6,6> Matrix6d;

// Kinematics functions
Matrix3d skew(Vector3d v);
Matrix3d fast_skewexp(Vector3d v);
Vector3d skewcoords(Matrix3d M);
Matrix3d skewlog(Matrix3d M);
MatrixXd kroen_product(MatrixXd A, MatrixXd B);
Matrix3d v_logmap(VectorXd x);
MatrixXd diagonalMatrix(MatrixXd M, unsigned int N);


Matrix4d inverse_se3(Matrix4d T);
Matrix4d expmap_se3(Vector6d x);
Vector6d logmap_se3(Matrix4d T);
Matrix6d adjoint_se3(Matrix4d T);
Matrix6d uncTinv_se3(Matrix4d T, Matrix6d covT );
Matrix6d unccomp_se3(Matrix4d T1, Matrix6d covT1, Matrix6d covTinc );
Vector6d reverse_se3(Vector6d x);


Vector3d logarithm_map_so3(Matrix3d R);
MatrixXd der_logarithm_map(Matrix4d T);
MatrixXd der_logarithm_map_appr(Matrix4d T, double delta);
double diffManifoldError(Matrix4d T1, Matrix4d T2);
bool is_finite(const MatrixXd x);
bool is_nan(const MatrixXd x);
double angDiff(double alpha, double beta);
double angDiff_d(double alpha, double beta);

// Auxiliar functions and structs for vectors
double vector_stdv_mad( VectorXf residues);
double vector_stdv_mad( vector<double> residues);
double vector_stdv_mad_nozero( vector<double> residues);
double vector_mean(vector<double> v);
double vector_stdv(vector<double> v);
double vector_stdv(vector<double> v, double v_mean);

struct compare_descriptor_by_NN_dist
{
    inline bool operator()(const vector<DMatch>& a, const vector<DMatch>& b){
        return ( a[0].distance < b[0].distance );
    }
};

struct compare_descriptor_by_NN12_dist
{
    inline bool operator()(const vector<DMatch>& a, const vector<DMatch>& b){
        return ( a[1].distance - a[0].distance > b[1].distance-b[0].distance );
    }
};

struct compare_descriptor_by_NN12_ratio
{
    inline bool operator()(const vector<DMatch>& a, const vector<DMatch>& b){
        return ( a[0].distance / a[1].distance > b[0].distance / b[1].distance );
    }
};

struct sort_descriptor_by_queryIdx
{
    inline bool operator()(const vector<DMatch>& a, const vector<DMatch>& b){
        return ( a[0].queryIdx < b[0].queryIdx );
    }
};

struct sort_descriptor_by_2nd_queryIdx
{
    inline bool operator()(const vector<DMatch>& a, const vector<DMatch>& b){
        return ( a[1].queryIdx < b[1].queryIdx );
    }
};

struct sort_descriptor_by_trainIdx
{
    inline bool operator()(const vector<DMatch>& a, const vector<DMatch>& b){
        return ( a[0].trainIdx < b[0].trainIdx );
    }
};

struct sort_confmat_by_score
{
    inline bool operator()(const Vector2d& a, const Vector2d& b){
        return ( a(1) > b(1) );
    }
};

struct sort_lines_by_response
{
    inline bool operator()(const KeyLine& a, const KeyLine& b){
        return ( a.response > b.response );
    }
};

struct sort_lines_by_length
{
    inline bool operator()(const KeyLine& a, const KeyLine& b){
        return ( a.lineLength > b.lineLength );
    }
};

