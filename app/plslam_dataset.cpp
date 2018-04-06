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

#ifdef HAS_MRPT
#include <slamScene.h>
#endif

#include <stereoFrame.h>
#include <stereoFrameHandler.h>
#include <boost/filesystem.hpp>

#include <mapFeatures.h>
#include <mapHandler.h>

#include <dataset.h>
#include <timer.h>

using namespace StVO;
using namespace PLSLAM;

void showHelp();
bool getInputArgs(int argc, char **argv, std::string &dataset_name, int &frame_offset, int &frame_number, int &frame_step, std::string &config_file);

int main(int argc, char **argv)
{

    // read inputs
    string dataset_name, config_file;
    int frame_offset = 0, frame_number = 0, frame_step = 1;
    if (!getInputArgs(argc, argv, dataset_name, frame_offset, frame_number, frame_step, config_file)) {
        showHelp();
        return -1;
    }

    if (!config_file.empty()) SlamConfig::loadFromFile(config_file);

    if (SlamConfig::hasPoints() &&
            (!boost::filesystem::exists(SlamConfig::dbowVocP()) || !boost::filesystem::is_regular_file(SlamConfig::dbowVocP()))) {
        cout << "Invalid vocabulary for points" << endl;
        return -1;
    }

    if (SlamConfig::hasLines() &&
            (!boost::filesystem::exists(SlamConfig::dbowVocL()) || !boost::filesystem::is_regular_file(SlamConfig::dbowVocL()))) {
        cout << "Invalid vocabulary for lines" << endl;
        return -1;
    }

    // read dataset root dir fron environment variable
    boost::filesystem::path dataset_path(string( getenv("DATASETS_DIR")));
    if (!boost::filesystem::exists(dataset_path) || !boost::filesystem::is_directory(dataset_path)) {
        cout << "Check your DATASETS_DIR environment variable" << endl;
        return -1;
    }

    dataset_path /= dataset_name;
    if (!boost::filesystem::exists(dataset_path) || !boost::filesystem::is_directory(dataset_path)) {
        cout << "Invalid dataset path" << endl;
        return -1;
    }

    cout << endl << "Initializing PL-SLAM...." << flush;

    string dataset_dir = dataset_path.string();
    PinholeStereoCamera*  cam_pin = new PinholeStereoCamera((dataset_path / "dataset_params.yaml").string());
    Dataset dataset(dataset_dir, *cam_pin, frame_offset, frame_number, frame_step);

    // create scene
    string scene_cfg_name;
    if( (dataset_name.find("kitti")!=std::string::npos) ||
        (dataset_name.find("malaga")!=std::string::npos)  )
        scene_cfg_name = "../config/scene_config.ini";
    else
        scene_cfg_name = "../config/scene_config_indoor.ini";
    slamScene scene(scene_cfg_name);
    Matrix4d Tcw, Tfw = Matrix4d::Identity();
    Tcw = Matrix4d::Identity();
    scene.setStereoCalibration( cam_pin->getK(), cam_pin->getB() );
    scene.initializeScene(Tfw);

    // create PLSLAM object
    PLSLAM::MapHandler* map = new PLSLAM::MapHandler(cam_pin);

    cout << " ... done. " << endl;

    Timer timer;

    // initialize and run PL-StVO
    int frame_counter = 0;
    StereoFrameHandler* StVO = new StereoFrameHandler(cam_pin);
    Mat img_l, img_r;
    while (dataset.nextFrame(img_l, img_r))
    {
        if( frame_counter == 0 ) // initialize
        {
            StVO->initialize(img_l,img_r,0);
            PLSLAM::KeyFrame* kf = new PLSLAM::KeyFrame( StVO->prev_frame, 0 );
            map->initialize( kf );
            // update scene
            scene.initViewports( img_l.cols, img_r.rows );
            scene.setImage(StVO->prev_frame->plotStereoFrame());
            scene.updateSceneSafe( map );
        }
        else // run
        {
            // PL-StVO
            timer.start();
            StVO->insertStereoPair( img_l, img_r, frame_counter );
            StVO->optimizePose();
            double t1 = timer.stop(); //ms
            cout << "------------------------------------------   Frame #" << frame_counter
                 << "   ----------------------------------------" << endl;
            cout << endl << "VO Runtime: " << t1 << endl;

            // check if a new keyframe is needed
            if( StVO->needNewKF() )
            {
                cout <<         "#KeyFrame:     " << map->max_kf_idx + 1;
                cout << endl << "#Points:       " << map->map_points.size();
                cout << endl << "#Segments:     " << map->map_lines.size();
                cout << endl << endl;

                // grab StF and update KF in StVO (the StVO thread can continue after this point)
                PLSLAM::KeyFrame* curr_kf = new PLSLAM::KeyFrame( StVO->curr_frame );
                // update KF in StVO
                StVO->currFrameIsKF();
                map->addKeyFrame( curr_kf );
                // update scene
                scene.setImage(StVO->curr_frame->plotStereoFrame());
                scene.updateSceneSafe( map );
            }
            else
            {
                scene.setImage(StVO->curr_frame->plotStereoFrame());
                scene.setPose( StVO->curr_frame->DT );
                scene.updateScene();
            }

            // update StVO
            StVO->updateFrame();
        }

        frame_counter++;
    }

    // finish SLAM
    map->finishSLAM();
    scene.updateScene( map );

    // perform GBA
    cout << endl << "Performing Global Bundle Adjustment..." ;
    map->globalBundleAdjustment();
    cout << " ... done." << endl;
    scene.updateSceneGraphs( map );

    // wait until the scene is closed
    while( scene.isOpen() );

    return 0;
}

void showHelp() {
    cout << endl << "Usage: ./imgPLSLAM <dataset_name> [options]" << endl
         << "Options:" << endl
         << "\t-o Offset (number of frames to skip in the dataset directory" << endl
         << "\t-n Number of frames to process the sequence" << endl
         << "\t-s Parameter to skip s-1 frames (default 1)" << endl
         << "\t-c Config file" << endl
         << endl;
}

bool getInputArgs(int argc, char **argv, std::string &dataset_name, int &frame_offset, int &frame_number, int &frame_step, std::string &config_file) {

    if( argc < 2 || argc > 10 || (argc % 2) == 1 )
        return false;

    dataset_name = argv[1];
    int nargs = argc/2 - 1;
    for( int i = 0; i < nargs; i++ )
    {
        int j = 2*i + 2;
        if( string(argv[j]) == "-o" )
            frame_offset = stoi(argv[j+1]);
        else if( string(argv[j]) == "-n" )
            frame_number = stoi(argv[j+1]);
        else if( string(argv[j]) == "-s" )
            frame_step = stoi(argv[j+1]);
        else if (string(argv[j]) == "-c")
            config_file = string(argv[j+1]);
        else
            return false;
    }

    return true;
}
