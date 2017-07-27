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

#include <auxiliar.h>

#define PI std::acos(-1.0)

/* Kinematics functions */

Matrix4d inverse_se3(Matrix4d T){
    Matrix4d Tinv = Matrix4d::Identity();
    Matrix3d R;
    Vector3d t;
    t = T.block(0,3,3,1);
    R = T.block(0,0,3,3);
    Tinv.block(0,0,3,3) =  R.transpose();
    Tinv.block(0,3,3,1) = -R.transpose() * t;
    return Tinv;
}

Matrix4d expmap_se3(Vector6d x){
    Matrix3d R, V, s, I = Matrix3d::Identity();
    Vector3d t, w;
    Matrix4d T = Matrix4d::Identity();
    w = x.tail(3);
    t = x.head(3);
    double theta = w.norm();
    if( theta < 0.000001 )
        R = I;
    else{
        s = skew(w)/theta;
        R = I + s * sin(theta) + s * s * (1.0f-cos(theta));
        V = I + s * (1.0f - cos(theta)) / theta + s * s * (theta - sin(theta)) / theta;
        t = V * t;
    }
    T.block(0,0,3,4) << R, t;
    return T;
}

Vector6d logmap_se3(Matrix4d T){
    Matrix3d R, Id3 = Matrix3d::Identity();
    Vector3d Vt, t, w;
    Matrix3d V = Matrix3d::Identity(), w_hat = Matrix3d::Zero();
    Vector6d x;
    Vt << T(0,3), T(1,3), T(2,3);
    w  << 0.f, 0.f, 0.f;
    R = T.block(0,0,3,3);
    double cosine = (R.trace() - 1.f)/2.f;
    if(cosine > 1.f)
        cosine = 1.f;
    else if (cosine < -1.f)
        cosine = -1.f;
    double sine = sqrt(1.0-cosine*cosine);
    if(sine > 1.f)
        sine = 1.f;
    else if (sine < -1.f)
        sine = -1.f;
    double theta  = acos(cosine);
    if( theta > 0.000001 ){
        w_hat = theta*(R-R.transpose())/(2.f*sine);
        w = skewcoords(w_hat);
        Matrix3d s;
        s = skew(w) / theta;
        V = Id3 + s * (1.f-cosine) / theta + s * s * (theta - sine) / theta;
    }
    t = V.inverse() * Vt;
    x.head(3) = t;
    x.tail(3) = w;
    return x;
}

Matrix6d adjoint_se3(Matrix4d T){
    Matrix6d AdjT = Matrix6d::Zero();
    Matrix3d R = T.block(0,0,3,3);
    AdjT.block(0,0,3,3) = R;
    AdjT.block(0,3,3,3) = skew( T.block(0,3,3,1) ) * R ;
    AdjT.block(3,3,3,3) = R;
    return AdjT;
}

Matrix6d uncTinv_se3(Matrix4d T, Matrix6d covT ){
    Matrix6d covTinv = Matrix6d::Zero();
    Matrix6d adjTinv;
    adjTinv = adjoint_se3( inverse_se3(T) );
    covTinv = adjTinv * covT * adjTinv.transpose();
    return covTinv;
}

Matrix6d unccomp_se3(Matrix4d T1, Matrix6d covT1, Matrix6d covTinc ){
    Matrix6d covT2 ; // covariance of T2 = T1 * inverse(Tinc)
    Matrix6d adjT1 = adjoint_se3(T1);
    covT2 = covT1 + adjT1 * covTinc * adjT1.transpose();
    return covT2;
}

Vector6d reverse_se3(Vector6d x){
    Vector6d x_out;
    x_out.head(3) = x.tail(3);
    x_out.tail(3) = x.head(3);
    return x_out;
}


Matrix3d skew(Vector3d v){

    Matrix3d skew;

    skew(0,0) = 0; skew(1,1) = 0; skew(2,2) = 0;

    skew(0,1) = -v(2);
    skew(0,2) =  v(1);
    skew(1,2) = -v(0);

    skew(1,0) =  v(2);
    skew(2,0) = -v(1);
    skew(2,1) =  v(0);

    return skew;
}

Matrix3d fast_skewexp(Vector3d v){
    Matrix3d M, s, I = Matrix3d::Identity();
    double theta = v.norm();
    if(theta==0.f)
        M = I;
    else{
        s = skew(v)/theta;
        M << I + s * sin(theta) + s * s * (1.f-cos(theta));
    }
    return M;
}

Vector3d skewcoords(Matrix3d M){
    Vector3d skew;
    skew << M(2,1), M(0,2), M(1,0);
    return skew;
}

Matrix3d skewlog(Matrix3d M){
    Matrix3d skew;
    double val = (M.trace() - 1.f)/2.f;
    if(val > 1.f)
        val = 1.f;
    else if (val < -1.f)
        val = -1.f;
    double theta = acos(val);
    if(theta == 0.f)
        skew << 0,0,0,0,0,0,0,0,0;
    else
        skew << (M-M.transpose())/(2.f*sin(theta))*theta;
    return skew;
}

MatrixXd kroen_product(MatrixXd A, MatrixXd B){
    unsigned int Ar = A.rows(), Ac = A.cols(), Br = B.rows(), Bc = B.cols();
    MatrixXd AB(Ar*Br,Ac*Bc);
    for (unsigned int i=0; i<Ar; ++i)
        for (unsigned int j=0; j<Ac; ++j)
            AB.block(i*Br,j*Bc,Br,Bc) = A(i,j)*B;
    return AB;
}

Matrix3d v_logmap(VectorXd x){
    Vector3d w;
    double theta, theta2, theta3;
    Matrix3d W, I, V;
    w << x(0), x(1), x(2);
    theta = w.norm();   theta2 = theta*theta; theta3 = theta2*theta;
    W = skew(w);
    I << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    if(theta>0.00001)
        V << I + ((1-cos(theta))/theta2)*W + ((theta-sin(theta))/theta3)*W*W;
    else
        V << I;
    return V;
}

MatrixXd diagonalMatrix(MatrixXd M, unsigned int N){
    MatrixXd A = MatrixXd::Zero(N,N);
    for(unsigned int i = 0; i < N; i++ ){
        A(i,i) = M(i,i);
    }
    return A;
}

Vector3d logarithm_map_so3(Matrix3d R){
    Matrix3d Id3 = Matrix3d::Identity();
    Vector3d w;
    Matrix3d V = Matrix3d::Identity(), w_hat = Matrix3d::Zero();
    w << 0.f, 0.f, 0.f;
    double cosine = (R.trace() - 1.f)/2.f;
    if(cosine > 1.f)
        cosine = 1.f;
    else if (cosine < -1.f)
        cosine = -1.f;
    double sine = sqrt(1.0-cosine*cosine);
    if(sine > 1.f)
        sine = 1.f;
    else if (sine < -1.f)
        sine = -1.f;
    double theta  = acos(cosine);
    if( theta > 0.000001 ){
        w_hat << theta*(R-R.transpose()) / (2.f*sine);
        w = skewcoords(w_hat);
    }
    return w;
}

MatrixXd der_logarithm_map(Matrix4d T)
{

    MatrixXd dlogT_dT = MatrixXd::Zero(6,12);

    // Approximate derivative of the logarithm_map wrt the transformation matrix
    Matrix3d L1 = Matrix3d::Zero();
    Matrix3d L2 = Matrix3d::Zero();
    Matrix3d L3 = Matrix3d::Zero();
    Matrix3d Vinv = Matrix3d::Identity();
    Vector6d x = logmap_se3(T);

    // estimates the cosine, sine, and theta
    double b;
    double cos_ = 0.5 * (T.block(0,0,3,3).trace() - 1.0 );
    if(cos_ > 1.f)
        cos_ = 1.f;
    else if (cos_ < -1.f)
        cos_ = -1.f;
    double theta  = acos(cos_);
    double theta2 = theta*theta;
    double sin_   = sin(theta);
    double cot_   = 1.0 / tan( 0.5*theta );
    double csc2_  = pow( 1.0/sin(0.5*theta) ,2);

    // if the angle is small...
    if( cos_ > 0.9999 )
    {
        b = 0.5;
        L1(1,2) = -b;
        L1(2,1) =  b;
        L2(0,2) =  b;
        L2(2,0) = -b;
        L3(0,1) = -b;
        L3(1,0) =  b;
        // form the full derivative
        dlogT_dT.block(3,0,3,3) = L1;
        dlogT_dT.block(3,3,3,3) = L2;
        dlogT_dT.block(3,6,3,3) = L3;
        dlogT_dT.block(0,9,3,3) = Vinv;
    }
    // if not...
    else
    {
        // rotation part
        double k;
        Vector3d a;
        a(0) = T(2,1) - T(1,2);
        a(1) = T(0,2) - T(2,0);
        a(2) = T(1,0) - T(0,1);
        k = ( theta * cos_ - sin_ ) / ( 4 * pow(sin_,3) );
        a = k * a;
        L1.block(0,0,3,1) = a;
        L2.block(0,1,3,1) = a;
        L3.block(0,2,3,1) = a;
        // translation part
        Matrix3d w_skew = skew( x.tail(3) );
        Vinv += w_skew * (1.f-cos_) / theta2 + w_skew * w_skew * (theta - sin_) / pow(theta,3);
        Vinv  = Vinv.inverse().eval();
        // dVinv_dR
        Vector3d t;
        Matrix3d B, skew_t;
        MatrixXd dVinv_dR(3,9);
        t = T.block(0,3,3,1);
        skew_t = skew( t );
        // - form a
        a =  (theta*cos_-sin_)/(8.0*pow(sin_,3)) * w_skew * t
            + ( (theta*sin_-theta2*cos_)*(0.5*theta*cot_-1.0) - theta*sin_*(0.25*theta*cot_+0.125*theta2*csc2_-1.0))/(4.0*theta2*pow(sin_,4)) * w_skew * w_skew * t;
        // - form B
        Vector3d w;
        Matrix3d dw_dR;
        w = x.tail(3);
        dw_dR.row(0) << -w(1)*t(1)-w(2)*t(2), 2.0*w(1)*t(0)-w(0)*t(1), 2.0*w(2)*t(0)-w(0)*t(2);
        dw_dR.row(1) << -w(1)*t(0)+2.0*w(0)*t(1), -w(0)*t(0)-w(2)*t(2), 2.0*w(2)*t(1)-w(1)*t(2);
        dw_dR.row(2) << -w(2)*t(0)+2.0*w(0)*t(2), -w(2)*t(1)+2.0*w(1)*t(2), -w(0)*t(0)-w(1)*t(1);
        B = -0.5*theta*skew_t/sin_ - (theta*cot_-2.0)*dw_dR/(8.0*pow(sin_,2));
        // - form dVinv_dR
        dVinv_dR.col(0) = a;
        dVinv_dR.col(1) = -B.col(2);
        dVinv_dR.col(2) = B.col(1);
        dVinv_dR.col(3) = B.col(2);
        dVinv_dR.col(4) = a;
        dVinv_dR.col(5) = -B.col(0);
        dVinv_dR.col(6) = -B.col(1);
        dVinv_dR.col(7) = B.col(0);
        dVinv_dR.col(8) = a;
        // form the full derivative
        dlogT_dT.block(3,0,3,3) = L1;
        dlogT_dT.block(3,3,3,3) = L2;
        dlogT_dT.block(3,6,3,3) = L3;
        dlogT_dT.block(0,9,3,3) = Vinv;
        dlogT_dT.block(0,0,3,9) = dVinv_dR;
    }

    return dlogT_dT;

}

MatrixXd der_logarithm_map_appr(Matrix4d T, double delta)
{

    MatrixXd dlogT_dT = MatrixXd::Zero(6,12);
    // Approximate derivative of the logarithm_map wrt the transformation matrix
    int k = 0;
    for( int j = 0; j < 4; j++)
    {
        for(int i = 0; i < 3; i++)
        {
            Matrix4d Taux = T;
            Taux(i,j) += delta;
            dlogT_dT.col(k) = ( logmap_se3(Taux)-logmap_se3(T) ) / delta;
            k++;
        }
    }
    return dlogT_dT;

}

double diffManifoldError(Matrix4d T1, Matrix4d T2){
    return ( logmap_se3(T1)-logmap_se3(T2) ).norm();
}

bool is_finite(const MatrixXd x){
    return ((x - x).array() == (x - x).array()).all();
}

bool is_nan(const MatrixXd x){
    for(unsigned int i = 0; i < x.rows(); i++){
        for(unsigned int j = 0; j < x.cols(); j++){
            if(std::isnan(x(i,j)))
                return true;
        }
    }
    return false;
}

double angDiff(double alpha, double beta){
    double theta = alpha - beta;
    if(theta>PI)
        theta -= 2.f * PI;
    if(theta<-PI)
        theta += 2.f * PI;
    return theta;
}

double angDiff_d(double alpha, double beta){
    double theta = alpha - beta;
    if(theta > 180.0)
        theta -= 360.0;
    if(theta<-PI)
        theta += 360.0;
    return theta;
}

/* Auxiliar functions and structs for vectors */

double vector_stdv_mad( VectorXf residues)
{
    // Return the standard deviation of vector with MAD estimation
    int n_samples = residues.size();
    sort( residues.derived().data(),residues.derived().data()+residues.size());
    double median = residues( n_samples/2 );
    residues << ( residues - VectorXf::Constant(n_samples,median) ).cwiseAbs();
    sort(residues.derived().data(),residues.derived().data()+residues.size());
    double MAD = residues( n_samples/2 );
    return 1.4826 * MAD;
}

double vector_stdv_mad( vector<double> residues)
{
    if( residues.size() != 0 )
    {
        // Return the standard deviation of vector with MAD estimation
        int n_samples = residues.size();
        sort( residues.begin(),residues.end() );
        double median = residues[ n_samples/2 ];
        for( int i = 0; i < n_samples; i++)
            residues[i] = fabsf( residues[i] - median );
        sort( residues.begin(),residues.end() );
        double MAD = residues[ n_samples/2 ];
        return 1.4826 * MAD;
    }
    else
        return 0.0;
}

double vector_mean(vector<double> v)
{
    double sum = 0.0;
    for( int i = 0; i < v.size(); i++ )
        sum += v[i];
    return sum / v.size();
}

double vector_stdv(vector<double> v)
{
    double mean = 0.0, e = 0.0;
    for( int i = 0; i < v.size(); i++ )
        mean += v[i];
    mean /= v.size();
    for( int i = 0; i < v.size(); i++ )
        e += (v[i] - mean)*(v[i] - mean);
    return sqrt(1.0/v.size()*e);
}

double vector_stdv(vector<double> v, double v_mean)
{
    double e = 0.0;
    for( int i = 0; i < v.size(); i++ )
        e += (v[i] - v_mean)*(v[i] - v_mean);
    return sqrt(1.0/v.size()*e);
}

double vector_stdv_mad_nozero( vector<double> residues)
{
    if( residues.size() != 0 )
    {
        // Return the standard deviation of vector with MAD estimation
        int n_samples = residues.size();
        sort( residues.begin(),residues.end() );
        // filter zeros
        vector<double> residues_f;
        for( int i = 0; i < n_samples; i++ )
            if( residues[i] > 0.0f )
                residues_f.push_back( residues[i] );

        // estimate robust stdv
        n_samples = residues_f.size();
        if( n_samples != 0 )
        {
            double median = residues_f[ n_samples/2 ];
            for( int i = 0; i < n_samples; i++)
                residues_f[i] = fabsf( residues_f[i] - median );
            sort( residues_f.begin(),residues_f.end() );
            double MAD = residues_f[ n_samples/2 ];
            return 1.4826 * MAD;
        }
        else
            return 0.0;
    }
    else
        return 0.0;
}
