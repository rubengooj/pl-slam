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

#ifdef HAS_MRPT
#include <voScene.h>
#include <mrpt/utils/CTicTac.h>
#endif

#include <stereoFrame.h>
#include <stereoFrameHandler.h>
#include <ctime>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/filesystem.hpp>
#include <yaml-cpp/yaml.h>

using namespace StVO;

int main(int argc, char **argv)
{

    // read inputs
    int frame_offset = 0, frame_number = 0, frame_step = 1;
    string dataset_name;
    if( argc < 2 || argc == 3 || argc == 5 ||argc == 7 || argc > 8 )
    {
        cout << endl << "Usage: ./dbgStVO <dataset_name>" << endl
             << "Options:\n"
             << "\t-o,Offset (number of frames to skip in the dataset directory\n"
             << "\t-n,Number of frames to process the sequence"
             << "\t-s,Parameter to skip s-1 frames (default 1)"
             << endl;
        return -1;
    }
    dataset_name = argv[1];
    if( argc == 4 || argc == 6 || argc == 8 )
    {
        int nargs = (argc-2)/2;
        for( int i = 0; i < nargs; i++ )
        {
            int j = 2*i + 2;
            if( string(argv[j]) == "-o" )
                frame_offset = atoi(argv[j+1]);
            else if( string(argv[j]) == "-n" )
                frame_number = atoi(argv[j+1]);
            else if( string(argv[j]) == "-s" )
                frame_step = atoi(argv[j+1]);
            else
            {
                cout << endl << "Usage: ./dbgStVO <dataset_name>" << endl
                     << "Options:\n"
                     << "\t-o,Offset (number of frames to skip in the dataset directory\n"
                     << "\t-n,Number of frames to process the sequence"
                     << "\t-s,Parameter to skip s-1 frames (default 1)"
                     << endl;
                return -1;
            }
        }
    }

    // read dataset root dir fron environment variable
    string dataset_dir( string( getenv("DATASETS_DIR") ) + dataset_name );

    // read content of the .yaml dataset configuration file
    YAML::Node dset_config = YAML::LoadFile(dataset_dir+"/dataset_params.yaml");

    // setup camera
    YAML::Node cam_config = dset_config["cam0"];
    string camera_model = cam_config["cam_model"].as<string>();
    PinholeStereoCamera*  cam_pin;
    bool rectify = false;
    if( camera_model == "Pinhole" )
    {
        // if EuRoC or Falcon yaml file
        if( cam_config["Kl"].IsDefined() )
        {
            rectify = true;
            Mat Kl, Kr, Dl, Dr, R, t;
            vector<double> Kl_ = cam_config["Kl"].as<vector<double>>();
            vector<double> Kr_ = cam_config["Kr"].as<vector<double>>();
            vector<double> Dl_ = cam_config["Dl"].as<vector<double>>();
            vector<double> Dr_ = cam_config["Dr"].as<vector<double>>();
            Kl = ( Mat_<float>(3,3) << Kl_[0], 0.0, Kl_[2], 0.0, Kl_[1], Kl_[3], 0.0, 0.0, 1.0 );
            Kr = ( Mat_<float>(3,3) << Kr_[0], 0.0, Kr_[2], 0.0, Kr_[1], Kr_[3], 0.0, 0.0, 1.0 );
            // load rotation and translation
            vector<double> R_ = cam_config["R"].as<vector<double>>();
            vector<double> t_ = cam_config["t"].as<vector<double>>();
            R = Mat::eye(3,3,CV_64F);
            t = Mat::eye(3,1,CV_64F);
            int k = 0;
            for( int i = 0; i < 3; i++ )
            {
                t.at<double>(i,0) = t_[i];
                for( int j = 0; j < 3; j++, k++ )
                    R.at<double>(i,j) = R_[k];
            }
            // load distortion parameters
            int Nd = Dl_.size();
            Dl = Mat::eye(1,Nd,CV_64F);
            Dr = Mat::eye(1,Nd,CV_64F);
            for( int i = 0; i < Nd; i++ )
            {
                Dl.at<double>(0,i) = Dl_[i];
                Dr.at<double>(0,i) = Dr_[i];
            }
            // if dtype is equidistant (now it is default)
            if( cam_config["dtype"].IsDefined() )
            {
                cam_pin = new PinholeStereoCamera(
                    cam_config["cam_width"].as<double>(),
                    cam_config["cam_height"].as<double>(),
                    cam_config["cam_bl"].as<double>(),
                    Kl, Kr, R, t, Dl, Dr, true);

            }
            else
            // create camera object for EuRoC
                cam_pin = new PinholeStereoCamera(
                    cam_config["cam_width"].as<double>(),
                    cam_config["cam_height"].as<double>(),
                    cam_config["cam_bl"].as<double>(),
                    Kl, Kr, R, t, Dl, Dr,false);
        }
        // else
        else
            cam_pin = new PinholeStereoCamera(
                cam_config["cam_width"].as<double>(),
                cam_config["cam_height"].as<double>(),
                fabs(cam_config["cam_fx"].as<double>()),
                fabs(cam_config["cam_fy"].as<double>()),
                cam_config["cam_cx"].as<double>(),
                cam_config["cam_cy"].as<double>(),
                cam_config["cam_bl"].as<double>(),
                cam_config["cam_d0"].as<double>(),
                cam_config["cam_d1"].as<double>(),
                cam_config["cam_d2"].as<double>(),
                cam_config["cam_d3"].as<double>()  );
    }
    else
    {
        cout << endl << "Not implemented yet." << endl;
        return -1;
    }

    // setup image directories
    string img_dir_l = dataset_dir + "/" + dset_config["images_subfolder_l"].as<string>();
    string img_dir_r = dataset_dir + "/" + dset_config["images_subfolder_r"].as<string>();

    // get a sorted list of files in the img directories
    boost::filesystem::path img_dir_path_l(img_dir_l.c_str());
    if (!boost::filesystem::exists(img_dir_path_l))
    {
        cout << endl << "Left image directory does not exist: \t" << img_dir_l << endl;
        return -1;
    }
    boost::filesystem::path img_dir_path_r(img_dir_r.c_str());
    if (!boost::filesystem::exists(img_dir_path_r))
    {
        cout << endl << "Right image directory does not exist: \t" << img_dir_r << endl;
        return -1;
    }

    // get all files in the img directories
    size_t max_len_l = 0;
    std::list<std::string> imgs_l;
    boost::filesystem::directory_iterator end_itr;
    for (boost::filesystem::directory_iterator file(img_dir_path_l); file != end_itr; ++file)
    {
        boost::filesystem::path filename_path = file->path().filename();
        if (boost::filesystem::is_regular_file(file->status()) &&
                (filename_path.extension() == ".png"  ||
                 filename_path.extension() == ".jpg"  ||
                 filename_path.extension() == ".jpeg" ||
                 filename_path.extension() == ".pnm"  ||
                 filename_path.extension() == ".tiff") )
        {
            std::string filename(filename_path.string());
            imgs_l.push_back(filename);
            max_len_l = max(max_len_l, filename.length());
        }
    }
    size_t max_len_r = 0;
    std::list<std::string> imgs_r;
    for (boost::filesystem::directory_iterator file(img_dir_path_r); file != end_itr; ++file)
    {
        boost::filesystem::path filename_path = file->path().filename();
        if (boost::filesystem::is_regular_file(file->status()) &&
                (filename_path.extension() == ".png"  ||
                 filename_path.extension() == ".jpg"  ||
                 filename_path.extension() == ".jpeg" ||
                 filename_path.extension() == ".pnm"  ||
                 filename_path.extension() == ".tiff") )
        {
            std::string filename(filename_path.string());
            imgs_r.push_back(filename);
            max_len_r = max(max_len_r, filename.length());
        }
    }

    // sort them by filename; add leading zeros to make filename-lengths equal if needed
    std::map<std::string, std::string> sorted_imgs_l, sorted_imgs_r;
    int n_imgs_l = 0, n_imgs_r = 0;
    for (std::list<std::string>::iterator img = imgs_l.begin(); img != imgs_l.end(); img++ )
    {
        sorted_imgs_l[std::string(max_len_l - img->length(), '0') + (*img)] = *img;
        n_imgs_l++;
    }
    for (std::list<std::string>::iterator img = imgs_r.begin(); img != imgs_r.end(); img++ )
    {
        sorted_imgs_r[std::string(max_len_r - img->length(), '0') + (*img)] = *img;
        n_imgs_r++;
    }

    // skip #frame_offset frames in the beginning, and grab #frame_number frames
    int k = 0;
    if( frame_number == 0 )
        frame_number = min(n_imgs_l,n_imgs_r);
    for (std::map<std::string, std::string>::iterator it_l = sorted_imgs_l.begin(), it_r = sorted_imgs_r.begin();
         ( it_l != sorted_imgs_l.end() && it_r != sorted_imgs_r.end() ); k++)
    {
        if( (k<frame_offset) || (k>frame_offset+frame_number-1) || ( k % frame_step != 0 ) )
        {
            sorted_imgs_l.erase( it_l++ );
            sorted_imgs_r.erase( it_r++ );
        }
        else
        {
            ++it_l;
            ++it_r;
        }
    }

    // create scene
    #ifdef HAS_MRPT
    string scene_cfg_name;
    if( (dataset_name.find("kitti")!=std::string::npos) ||
        (dataset_name.find("malaga")!=std::string::npos)  )
        scene_cfg_name = "../config/scene_config.ini";
    else
        scene_cfg_name = "../config/scene_config_indoor.ini";
    voScene scene(scene_cfg_name);
    Matrix4d Tcw, Tfw = Matrix4d::Identity(), Tfw_prev = Matrix4d::Identity(), T_inc;
    Vector6d cov_eig;
    Matrix6d cov;
    Tcw = Matrix4d::Identity();
    Tcw << 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1;
    scene.initializeScene(Tfw,false);
    mrpt::utils::CTicTac clock;
    #endif

    // initialize and run PL-StVO
    int frame_counter = 0;
    double t1;
    StereoFrameHandler* StVO = new StereoFrameHandler(cam_pin);
    for( std::map<std::string, std::string>::iterator it_l = sorted_imgs_l.begin(), it_r = sorted_imgs_r.begin();
         it_l != sorted_imgs_l.end(), it_r != sorted_imgs_r.end(); ++it_l, ++it_r, frame_counter++ )
    {
        // load images
        boost::filesystem::path img_path_l = img_dir_path_l / boost::filesystem::path(it_l->second.c_str());
        boost::filesystem::path img_path_r = img_dir_path_r / boost::filesystem::path(it_r->second.c_str());
        Mat img_l( imread(img_path_l.string(), CV_LOAD_IMAGE_UNCHANGED) );  assert(!img_l.empty()); // it depends on the OpenCV version!!!
        Mat img_r( imread(img_path_r.string(), CV_LOAD_IMAGE_UNCHANGED) );  assert(!img_r.empty());        
        // if images are distorted
        if( rectify )
        {
            Mat img_l_rec, img_r_rec;
            cam_pin->rectifyImagesLR(img_l,img_l_rec,img_r,img_r_rec);
            img_l = img_l_rec;
            img_r = img_r_rec;
        }
        // initialize
        if( frame_counter == 0 )
            StVO->initialize(img_l,img_r,0);
        // run
        else
        {
            #ifdef HAS_MRPT
            clock.Tic();
            #endif
            // PL-StVO
            StVO->insertStereoPair( img_l, img_r, frame_counter );
            // solve with robust kernel and IRLS
            StVO->optimizePose(  );
            #ifdef HAS_MRPT
            t1 = 1000 * clock.Tac(); //ms
            #endif
            T_inc   = StVO->curr_frame->DT;
            cov     = StVO->curr_frame->DT_cov;
            cov_eig = StVO->curr_frame->DT_cov_eig;
            // update scene
            #ifdef HAS_MRPT
            scene.setText(frame_counter,t1,StVO->n_inliers_pt,StVO->matched_pt.size(),StVO->n_inliers_ls,StVO->matched_ls.size());
            scene.setCov( cov );
            scene.setPose( T_inc );
            imwrite("../config/aux/img_aux.png",StVO->curr_frame->plotStereoFrame());
            scene.setImage( "../config/aux/img_aux.png" );
            //scene.setImage( img_l );
            scene.updateScene();
            #endif
            // console output
            cout.setf(ios::fixed,ios::floatfield); cout.precision(8);
            cout << "Frame: " << frame_counter << "\tRes.: " << StVO->curr_frame->err_norm;
            cout.setf(ios::fixed,ios::floatfield); cout.precision(3);
            cout << " \t Proc. time: " << t1 << " ms\t ";
            if( Config::adaptativeFAST() )  cout << "\t FAST: "   << StVO->orb_fast_th;
            if( Config::hasPoints())        cout << "\t Points: " << StVO->matched_pt.size() << " (" << StVO->n_inliers_pt << ") " ;
            if( Config::hasLines() )        cout << "\t Lines:  " << StVO->matched_ls.size() << " (" << StVO->n_inliers_ls << ") " ;
            cout << endl;
            // update StVO
            StVO->updateFrame();
        }
    }

    // wait until the scene is closed
    #ifdef HAS_MRPT
    while( scene.isOpen() );
    #endif

    return 0;
}

