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

#include "slamScene.h"

#include <opencv2/imgproc.hpp>

// Auxiliar functions
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

namespace PLSLAM{

// Constructors and destructor

slamScene::slamScene(){

    sbb     = 1.0f;
    sref    = 0.2f;
    srad    = 0.005f;
    sline   = 2.0f;
    saxis   = 0.5f;
    sfreq   = 3.0f;
    szoom   = 200.0f;
    selli   = 2.0f;
    selev   = 0.f;
    sazim   = -90.f;
    sfrust  = 0.2f;
    slinef  = 1.f;
    win     = new CDisplayWindow3D("3D Scene",1920,1080);

    hasText         = true;
    hasAxes         = true;
    hasLegend       = false;
    hasHelp         = false;
    hasGT           = false;
    hasComparison   = false;
    hasImg          = true;
    hasLines        = false;
    hasPoints       = false;
    isKitti         = true;

}

slamScene::slamScene(string configFile){

    CConfigFile config(configFile);
    sbb             = config.read_double("Scene","sbb",1.f);
    sref            = config.read_double("Scene","sref",0.2f);
    srad            = config.read_double("Scene","srad",0.005f);
    sline           = config.read_double("Scene","sline",2.f);
    saxis           = config.read_double("Scene","saxis",0.5f);
    sfreq           = config.read_double("Scene","sfreq",3.f);
    szoom           = config.read_double("Scene","szoom",3.f);
    selli           = config.read_double("Scene","selli",1.f);
    selev           = config.read_double("Scene","selev",30.f);
    sazim           = config.read_double("Scene","sazim",-135.f);
    sfrust          = config.read_double("Scene","sfrust",0.2f);
    slinef          = config.read_double("Scene","slinef",1.f);
    win             = new CDisplayWindow3D("3D Scene",1920,1080);

    hasText         = config.read_bool("Scene","hasText",true);
    hasAxes         = config.read_bool("Scene","hasAxes",true);
    hasLegend       = config.read_bool("Scene","hasLegend",false);
    hasHelp         = config.read_bool("Scene","hasHelp",false);
    hasGT           = config.read_bool("Scene","hasGT",false);
    hasTraj         = config.read_bool("Scene","hasTraj",true);
    hasComparison   = config.read_bool("Scene","hasComparison",false);
    hasImg          = config.read_bool("Scene","hasImg",false);
    hasLines        = config.read_bool("Scene","hasLines",false);
    hasPoints       = config.read_bool("Scene","hasPoints",false);
    hasFrustum      = config.read_bool("Scene","hasFrustum",false);
    isKitti         = config.read_bool("Scene","isKitti",true);

    Matrix4d x_cw;
    x_cw << 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1;
    CPose3D x_aux(getPoseFormat(x_cw));
    pose_ini = x_aux;

}

slamScene::~slamScene(){

}

// Initialize the scene

void slamScene::initializeScene(Matrix4d x_0){

    // Initialize the scene and window
    win->setCameraElevationDeg(selev);
    win->setCameraAzimuthDeg(sazim);
    win->setCameraZoom(szoom);
    theScene = win->get3DSceneAndLock();

    // Initialize poses of different objects
    CPose3D x_aux(getPoseFormat(x_0));
    pose    = x_aux;
    pose1   = x_aux;
    pose_0  = x_aux;
    pose_gt = pose_ini;
    x_ini   = x_0;
    pose.getAsVector(v_aux);
    pose1.getAsVector(v_aux1);
    pose_gt.getAsVector(v_auxgt);

    // Initialize the camera object
    frustumL_.setFromValues(0,0,0,  0, -90.f*CV_PI/180.f, -90.f*CV_PI/180.f);
    frustumR_.setFromValues(sfrust,0,0,  0, -90.f*CV_PI/180.f, -90.f*CV_PI/180.f);
    frustObj  = opengl::CFrustum::Create();
    {
        frustObj->setPose(pose_0+frustumL_);
        frustObj->setLineWidth (slinef);
        frustObj->setScale(sfrust);
        frustObj->setColor_u8( TColor(200,0,0) );
        theScene->insert(frustObj);
    }
    frustObj1 = opengl::CFrustum::Create();
    {
        frustObj1->setPose(pose_0+frustumR_ );
        frustObj1->setLineWidth (slinef);
        frustObj1->setScale(sfrust);
        frustObj1->setColor_u8( TColor(200,0,0) );
        theScene->insert(frustObj1);
    }

    srefObj = opengl::stock_objects::CornerXYZSimple();
    srefObj->setPose(pose_0);
    srefObj->setScale(sref);

    // Initialize the axes
    if(hasAxes){
        axesObj = opengl::CAxis::Create();
        axesObj->setFrequency(sfreq);
        axesObj->enableTickMarks(false);
        axesObj->setAxisLimits(-saxis,-saxis,-saxis, saxis,saxis,saxis);
        axesObj->setColor(0.0,0.0,0.0);
        //theScene->insert( axesObj );
    }

    // Initialize the text [TODO]
    if(hasText){
        string text = "KeyFrame: \t0 \nFrequency: \t0 Hz \nPoints:   \t0 \nLines:    \t0";
        win->addTextMessage(0.85,0.95, text, TColorf(.0,.0,.0), 0, mrpt::opengl::MRPT_GLUT_BITMAP_TIMES_ROMAN_10 );
    }

    // Initialize the lines
    lineObj = opengl::CSetOfLines::Create();
    lineObj->setLineWidth(1.0);
    lineObj->setColor(0,0,0);
    theScene->insert( lineObj );

    lineObj_local = opengl::CSetOfLines::Create();
    lineObj_local->setLineWidth(1.0);
    lineObj_local->setColor(200,0,200);
    theScene->insert( lineObj_local );

    // Initialize the point cloud
    pointObj = opengl::CPointCloud::Create();
    pointObj->setPointSize(2.0);
    pointObj->setColor(0,0,0);
    theScene->insert( pointObj );

    pointObj_local = opengl::CPointCloud::Create();
    pointObj_local->setPointSize(2.0);
    pointObj_local->setColor(200,0,0);
    theScene->insert( pointObj_local );

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

}

void slamScene::initViewports(int W, int H){

    theScene = win->get3DSceneAndLock();

    // Initialize the viewports
    setHelp();
    setLegend();
    image = theScene->createViewport("image");
    image->setViewportPosition(20, 20, W/2, H/2);

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

}

// Update the scene

bool slamScene::updateScene(){

    theScene = win->get3DSceneAndLock();
    bool restart = false;

    theScene->removeObject( frustObj );
    theScene->removeObject( frustObj1 );

    // Update the camera pose
    {
        CPose3D x_aux1(getPoseFormat(x));
        pose = pose + x_aux1;
        v_aux_ = v_aux;
        pose.getAsVector(v_aux);
        frustObj  = opengl::CFrustum::Create();
        {
            frustObj->setPose(pose+frustumL_);
            frustObj->setLineWidth (slinef);
            frustObj->setScale(sfrust);
            frustObj->setColor_u8( TColor(200,0,0) );
            theScene->insert(frustObj);
        }
        frustObj1 = opengl::CFrustum::Create();
        {
            frustObj1->setPose(pose+frustumR_ );
            frustObj1->setLineWidth (slinef);
            frustObj1->setScale(sfrust);
            frustObj1->setColor_u8( TColor(200,0,0) );
            theScene->insert(frustObj1);
        }
    }

    // Update the image
    if(hasImg)
        image->setImageView( img_mrpt_image );

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

    return restart;

}

bool slamScene::updateSceneVO( Matrix4d T_last_kf ){

    theScene = win->get3DSceneAndLock();
    bool restart = false;

    theScene->removeObject( frustObj );
    theScene->removeObject( frustObj1 );

    // Update the camera pose
    pose = CPose3D( getPoseFormat( T_last_kf ) );
    {
        CPose3D x_aux1(getPoseFormat(x));
        pose = pose + x_aux1;
        v_aux_ = v_aux;
        pose.getAsVector(v_aux);
        frustObj  = opengl::CFrustum::Create();
        {
            frustObj->setPose(pose+frustumL_);
            frustObj->setLineWidth (slinef);
            frustObj->setScale(sfrust);
            frustObj->setColor_u8( TColor(200,0,0) );
            theScene->insert(frustObj);
        }
        frustObj1 = opengl::CFrustum::Create();
        {
            frustObj1->setPose(pose+frustumR_ );
            frustObj1->setLineWidth (slinef);
            frustObj1->setScale(sfrust);
            frustObj1->setColor_u8( TColor(200,0,0) );
            theScene->insert(frustObj1);
        }
    }

    // Update the image
    if(hasImg)
        image->setImageView( img_mrpt_image );    

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

    return restart;

}

bool slamScene::updateScene(const MapHandler* map){

    theScene = win->get3DSceneAndLock();
    bool restart = false;

    theScene->removeObject( kfsObj );
    theScene->removeObject( kfsLinesObj );
    theScene->removeObject( frustObj );
    theScene->removeObject( frustObj1 );

    // Represent KFs
    CPose3D kf_pose;
    Vector3d Pi;
    try{
    Pi = map->map_keyframes[0]->T_kf_w.col(3).head(3);
    }catch(...){}
    Vector3d Pj;
    kfsObj      = opengl::CSetOfObjects::Create();
    kfsLinesObj = opengl::CSetOfLines::Create();
    for( vector<KeyFrame*>::const_iterator it = map->map_keyframes.begin(); it != map->map_keyframes.end(); it++ )
    {
        if( (*it)!=NULL )
        {
            // if last keyframe
            if( (*it)->kf_idx == map->map_keyframes.back()->kf_idx )
            {
                kf_pose = CPose3D( getPoseFormat( (*it)->T_kf_w ) );
                opengl::CFrustumPtr frust_ = opengl::CFrustum::Create();
                {
                    frust_->setPose( kf_pose + frustumL_ );
                    frust_->setLineWidth (slinef);
                    frust_->setScale(sfrust / 2.0 );
                    if( (*it)->local )
                        frust_->setColor_u8( TColor(200,0,0) );
                    else
                        frust_->setColor_u8( TColor(0,0,200) );
                    kfsObj->insert( frust_ );
                }
                // represent spanning tree
                Pj = (*it)->T_kf_w.col(3).head(3);
                kfsLinesObj->appendLine( Pi(0),Pi(1),Pi(2), Pj(0),Pj(1),Pj(2) );
                Pi = Pj;
                // represent VO frustum
                frustObj  = opengl::CFrustum::Create();
                {
                    frustObj->setPose( kf_pose + frustumL_ );
                    frustObj->setLineWidth (slinef);
                    frustObj->setScale(sfrust);
                    frustObj->setColor_u8( TColor(200,0,0) );
                    theScene->insert(frustObj);
                }
                frustObj1 = opengl::CFrustum::Create();
                {
                    frustObj1->setPose( kf_pose + frustumR_ );
                    frustObj1->setLineWidth (slinef);
                    frustObj1->setScale(sfrust);
                    frustObj1->setColor_u8( TColor(200,0,0) );
                    theScene->insert(frustObj1);
                }
                pose = kf_pose;
            }
            // if not
            else
            {
                kf_pose = CPose3D( getPoseFormat( (*it)->T_kf_w ) );
                opengl::CFrustumPtr frust_ = opengl::CFrustum::Create();
                {
                    frust_->setPose( kf_pose + frustumL_ );
                    frust_->setLineWidth (slinef);
                    frust_->setScale(sfrust / 2.0);
                    if( (*it)->local )
                        frust_->setColor_u8( TColor(200,0,0) );
                    else
                        frust_->setColor_u8( TColor(0,0,200) );
                    kfsObj->insert( frust_ );
                }
                // represent spanning tree
                Pj = (*it)->T_kf_w.col(3).head(3);
                kfsLinesObj->appendLine( Pi(0),Pi(1),Pi(2), Pj(0),Pj(1),Pj(2) );
                Pi = Pj;
            }
        }
    }
    kfsLinesObj->setLineWidth(0.5f);
    kfsLinesObj->setColor(0,0.5,0);
    theScene->insert( kfsLinesObj );
    theScene->insert( kfsObj );

    // Represent point LMs
    if( hasPoints )
    {
        pointObj->clear();
        pointObj_local->clear();
        for( vector<MapPoint*>::const_iterator it = map->map_points.begin(); it!=map->map_points.end(); it++)
        {
            try{
            if( (*it)!=NULL )
            {
                if( (*it)->local )
                    pointObj_local->insertPoint( (*it)->point3D(0),(*it)->point3D(1),(*it)->point3D(2) );
                else
                    pointObj->insertPoint( (*it)->point3D(0),(*it)->point3D(1),(*it)->point3D(2) );
            }
            }catch(...){}
        }
    }

    // Represent line LMs
    if( hasLines )
    {
        lineObj->clear();
        lineObj_local->clear();
        for( vector<MapLine*>::const_iterator it = map->map_lines.begin(); it!=map->map_lines.end(); it++)
        {
            try{
            if( (*it)!=NULL )
            {
                Vector6d L;
                L = (*it)->line3D;
                if( (*it)->local )
                    lineObj_local->appendLine( L(0),L(1),L(2),L(3),L(4),L(5) );
                else
                    lineObj->appendLine( L(0),L(1),L(2),L(3),L(4),L(5) );
            }
            }catch(...){}
        }
    }

    // Update the text
    if(hasText){
        string text = "KeyFrame: \t" + to_string(frame) + " \nFrequency: \t" + to_string_with_precision(1000.f/time,4) + " Hz \nPoints:   \t" + to_string(nPoints) + "\nLines:    \t" + to_string(nLines);
        //string text = "Frame: \t \t" + to_string(frame) + " \n" + "Frequency: \t" + to_string_with_precision(1000.f/time,4) + " Hz \n" + "Lines:  \t" + to_string(nLines) + " (" + to_string(nLinesH) + ") \nPoints: \t" + to_string(nPoints) + " (" + to_string(nPointsH) + ")";
        win->addTextMessage(0.85,0.95, text, TColorf(.0,.0,.0), 0, mrpt::opengl::MRPT_GLUT_BITMAP_TIMES_ROMAN_10 );
    }

    // Update the image
    if(hasImg)
        image->setImageView( img_mrpt_image );

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

    // Key events
    if(win->keyHit()){
        key = win->getPushedKey(&kmods);
        if(key == MRPTK_SPACE){                     // Space    Reset VO
            theScene->clear();
            initializeScene(x_ini);
        }
        else if (key == MRPTK_ESCAPE){              // Esc      Restart VO
            theScene->clear();
            initializeScene(x_ini);
            restart = true;
        }
        else if ( (key == 104) || (key == 72) ){    // H        help
            hasHelp   = !hasHelp;
            if(!hasHelp)
                help->setViewportPosition(2000, 2000, 300, 376);
            else
                help->setViewportPosition(1600, 20, 300, 376);
        }
        else if ( (key == 103) || (key == 71) ){    // G        legend
            hasLegend   = !hasLegend;
            if(!hasLegend)
                legend->setViewportPosition(2000, 2000, 250, 97);
            else
                legend->setViewportPosition(20, 900, 250, 97);
        }
        else if ( (key ==  97) || (key == 65) ){    // A        axes
            hasAxes   = !hasAxes;
            if(!hasAxes){
                axesObj.clear();
            }
            else{
                axesObj = opengl::CAxis::Create();
                axesObj->setFrequency(sfreq);
                axesObj->enableTickMarks(false);
                axesObj->setAxisLimits(-saxis,-saxis,-saxis, saxis,saxis,saxis);
                theScene->insert( axesObj );
            }
        }
        else if ( (key == 112) || (key == 80) ){    // P        points
            hasPoints = !hasPoints;
            if(!hasPoints){
                elliObjP.clear();
            }
            else{
                elliObjP = opengl::CSetOfObjects::Create();
                elliObjP->setScale(selli);
                elliObjP->setPose(pose);
                theScene->insert(elliObjP);
            }
        }
        else if ( (key == 108) || (key == 76) ){    // L        lines
            hasLines  = !hasLines;
            if(!hasLines){
                elliObjL.clear();
            }
            else{
                elliObjL = opengl::CSetOfObjects::Create();
                elliObjL->setScale(selli);
                elliObjL->setPose(pose);
                theScene->insert(elliObjL);
            }
        }
        else if ( (key == 116) || (key == 84) ){    // T        text
            hasText  = !hasText;
            if(!hasText){
                string text = "";
                win->addTextMessage(0.85,0.95, text, TColorf(1.0,1.0,1.0), 0, mrpt::opengl::MRPT_GLUT_BITMAP_HELVETICA_10 );
            }
        }
        else if ( (key == 105) || (key == 73) ){    // I        image
            hasImg   = !hasImg;
            if(isKitti){
                if(hasImg)
                    image->setViewportPosition(20, 20, 621, 188);
                else
                    image->setViewportPosition(2000, 2000, 621, 188);
            }
            else{
                if(hasImg)
                    image->setViewportPosition(20, 20, 400, 300);
                else
                    image->setViewportPosition(2000, 2000, 400, 300);
            }

        }
    }

    return restart;

}

bool slamScene::updateSceneSafe(const MapHandler* map){

    theScene = win->get3DSceneAndLock();
    bool restart = false;

    theScene->removeObject( kfsObj );
    theScene->removeObject( kfsLinesObj );
    theScene->removeObject( frustObj );
    theScene->removeObject( frustObj1 );

    // grab indices (so we don't have memory access problems)
    int pt_idx = map->map_points.size();
    int ls_idx = map->map_lines.size();

    // Represent KFs
    CPose3D kf_pose;
    Vector3d Pi;
    try{
    Pi = map->map_keyframes[0]->T_kf_w.col(3).head(3);
    }catch(...){}
    Vector3d Pj;
    kfsObj      = opengl::CSetOfObjects::Create();
    kfsLinesObj = opengl::CSetOfLines::Create();
    for( vector<KeyFrame*>::const_iterator it = map->map_keyframes.begin(); it != map->map_keyframes.end(); it++ )
    {
        if( (*it)!=NULL )
        {
            // if last keyframe
            if( (*it)->kf_idx == map->map_keyframes.back()->kf_idx )
            {
                kf_pose = CPose3D( getPoseFormat( (*it)->T_kf_w ) );
                opengl::CFrustumPtr frust_ = opengl::CFrustum::Create();
                {
                    frust_->setPose( kf_pose + frustumL_ );
                    frust_->setLineWidth (slinef);
                    frust_->setScale(sfrust / 2.0 );
                    if( (*it)->local )
                        frust_->setColor_u8( TColor(200,0,0) );
                    else
                        frust_->setColor_u8( TColor(0,0,200) );
                    kfsObj->insert( frust_ );
                }
                // represent spanning tree
                Pj = (*it)->T_kf_w.col(3).head(3);
                kfsLinesObj->appendLine( Pi(0),Pi(1),Pi(2), Pj(0),Pj(1),Pj(2) );
                Pi = Pj;
                // represent VO frustum
                frustObj  = opengl::CFrustum::Create();
                {
                    frustObj->setPose( kf_pose + frustumL_ );
                    frustObj->setLineWidth (slinef);
                    frustObj->setScale(sfrust);
                    frustObj->setColor_u8( TColor(200,0,0) );
                    theScene->insert(frustObj);
                }
                frustObj1 = opengl::CFrustum::Create();
                {
                    frustObj1->setPose( kf_pose + frustumR_ );
                    frustObj1->setLineWidth (slinef);
                    frustObj1->setScale(sfrust);
                    frustObj1->setColor_u8( TColor(200,0,0) );
                    theScene->insert(frustObj1);
                }
                pose = kf_pose;
            }
            // if not
            else
            {
                kf_pose = CPose3D( getPoseFormat( (*it)->T_kf_w ) );
                opengl::CFrustumPtr frust_ = opengl::CFrustum::Create();
                {
                    frust_->setPose( kf_pose + frustumL_ );
                    frust_->setLineWidth (slinef);
                    frust_->setScale(sfrust / 2.0);
                    if( (*it)->local )
                        frust_->setColor_u8( TColor(200,0,0) );
                    else
                        frust_->setColor_u8( TColor(0,0,200) );
                    kfsObj->insert( frust_ );
                }
                // represent spanning tree
                Pj = (*it)->T_kf_w.col(3).head(3);
                kfsLinesObj->appendLine( Pi(0),Pi(1),Pi(2), Pj(0),Pj(1),Pj(2) );
                Pi = Pj;
            }
        }
    }
    kfsLinesObj->setLineWidth(0.5f);
    kfsLinesObj->setColor(0,0.5,0);
    theScene->insert( kfsLinesObj );
    theScene->insert( kfsObj );

    // Represent point LMs
    if( hasPoints )
    {
        pointObj->clear();
        pointObj_local->clear();
        for( int i = 0; i < pt_idx; i++ )
        {
            MapPoint* it = map->map_points[i];
            if( it!=NULL )
            {
                if( it->local && it->inlier )
                    pointObj_local->insertPoint( it->point3D(0),it->point3D(1),it->point3D(2) );
                else
                    pointObj->insertPoint( it->point3D(0),it->point3D(1),it->point3D(2) );
            }
        }
    }

    // Represent line LMs
    if( hasLines )
    {
        lineObj->clear();
        lineObj_local->clear();
        for( int i = 0; i < ls_idx; i++ )
        {
            MapLine* it = map->map_lines[i];
            if( it!=NULL && it->inlier )
            {
                Vector6d L;
                L = it->line3D;
                if( it->local )
                    lineObj_local->appendLine( L(0),L(1),L(2),L(3),L(4),L(5) );
                else
                    lineObj->appendLine( L(0),L(1),L(2),L(3),L(4),L(5) );
            }
        }
    }

    // Update the text
    if(hasText){
        string text = "KeyFrame: \t" + to_string(frame) + " \nFrequency: \t" + to_string_with_precision(1000.f/time,4) + " Hz \nPoints:   \t" + to_string(nPoints) + "\nLines:    \t" + to_string(nLines);
        //string text = "Frame: \t \t" + to_string(frame) + " \n" + "Frequency: \t" + to_string_with_precision(1000.f/time,4) + " Hz \n" + "Lines:  \t" + to_string(nLines) + " (" + to_string(nLinesH) + ") \nPoints: \t" + to_string(nPoints) + " (" + to_string(nPointsH) + ")";
        win->addTextMessage(0.85,0.95, text, TColorf(.0,.0,.0), 0, mrpt::opengl::MRPT_GLUT_BITMAP_TIMES_ROMAN_10 );
    }

    // Update the image
    if(hasImg)
        image->setImageView( img_mrpt_image );

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

    // Key events
    if(win->keyHit()){
        key = win->getPushedKey(&kmods);
        if(key == MRPTK_SPACE){                     // Space    Reset VO
            theScene->clear();
            initializeScene(x_ini);
        }
        else if (key == MRPTK_ESCAPE){              // Esc      Restart VO
            theScene->clear();
            initializeScene(x_ini);
            restart = true;
        }
        else if ( (key == 104) || (key == 72) ){    // H        help
            hasHelp   = !hasHelp;
            if(!hasHelp)
                help->setViewportPosition(2000, 2000, 300, 376);
            else
                help->setViewportPosition(1600, 20, 300, 376);
        }
        else if ( (key == 103) || (key == 71) ){    // G        legend
            hasLegend   = !hasLegend;
            if(!hasLegend)
                legend->setViewportPosition(2000, 2000, 250, 97);
            else
                legend->setViewportPosition(20, 900, 250, 97);
        }
        else if ( (key ==  97) || (key == 65) ){    // A        axes
            hasAxes   = !hasAxes;
            if(!hasAxes){
                axesObj.clear();
            }
            else{
                axesObj = opengl::CAxis::Create();
                axesObj->setFrequency(sfreq);
                axesObj->enableTickMarks(false);
                axesObj->setAxisLimits(-saxis,-saxis,-saxis, saxis,saxis,saxis);
                theScene->insert( axesObj );
            }
        }
        else if ( (key == 112) || (key == 80) ){    // P        points
            hasPoints = !hasPoints;
            if(!hasPoints){
                elliObjP.clear();
            }
            else{
                elliObjP = opengl::CSetOfObjects::Create();
                elliObjP->setScale(selli);
                elliObjP->setPose(pose);
                theScene->insert(elliObjP);
            }
        }
        else if ( (key == 108) || (key == 76) ){    // L        lines
            hasLines  = !hasLines;
            if(!hasLines){
                elliObjL.clear();
            }
            else{
                elliObjL = opengl::CSetOfObjects::Create();
                elliObjL->setScale(selli);
                elliObjL->setPose(pose);
                theScene->insert(elliObjL);
            }
        }
        else if ( (key == 116) || (key == 84) ){    // T        text
            hasText  = !hasText;
            if(!hasText){
                string text = "";
                win->addTextMessage(0.85,0.95, text, TColorf(1.0,1.0,1.0), 0, mrpt::opengl::MRPT_GLUT_BITMAP_HELVETICA_10 );
            }
        }
        else if ( (key == 105) || (key == 73) ){    // I        image
            hasImg   = !hasImg;
            if(isKitti){
                if(hasImg)
                    image->setViewportPosition(20, 20, 621, 188);
                else
                    image->setViewportPosition(2000, 2000, 621, 188);
            }
            else{
                if(hasImg)
                    image->setViewportPosition(20, 20, 400, 300);
                else
                    image->setViewportPosition(2000, 2000, 400, 300);
            }

        }
    }

    return restart;

}

void slamScene::updateSceneGraphs( const MapHandler* map )
{

    // Reset scene
    theScene = win->get3DSceneAndLock();
    theScene->removeObject( kfsObj );
    pointObj->clear();
    pointObj_local->clear();
    lineObj->clear();
    lineObj_local->clear();

    // Represent KFs
    CPose3D kf_pose;
    kfsObj = opengl::CSetOfObjects::Create();
    for( vector<KeyFrame*>::const_iterator it = map->map_keyframes.begin(); it!=map->map_keyframes.end(); it++)
    {
        if( (*it)!=NULL )
        {
            kf_pose = CPose3D( getPoseFormat( (*it)->T_kf_w ) );
            opengl::CFrustumPtr frust_ = opengl::CFrustum::Create();
            {
                frust_->setPose( kf_pose + frustumL_ );
                frust_->setLineWidth (slinef);
                frust_->setScale(sfrust/2.f);
                if( (*it)->local )
                    frust_->setColor_u8( TColor(0,0,200) );
                else
                    frust_->setColor_u8( TColor(0,0,200) );
                kfsObj->insert( frust_ );
            }
        }
    }
    theScene->insert( kfsObj );

    // Represent graph
    int Nkf = map->full_graph.size();
    for( int i = 0; i < Nkf; i++ )
    {
        for( int j = 0; j < Nkf; j++ )
        {
            if( map->map_keyframes[i] != NULL && map->map_keyframes[j] != NULL )
            {
                if( map->full_graph[i][j] >= SlamConfig::minLMCovGraph() )
                {
                    Vector3d Pi = map->map_keyframes[i]->T_kf_w.col(3).head(3);
                    Vector3d Pj = map->map_keyframes[j]->T_kf_w.col(3).head(3);
                    opengl::CSimpleLinePtr obj = opengl::CSimpleLine::Create();
                    obj->setLineCoords(Pi(0),Pi(1),Pi(2), Pj(0),Pj(1),Pj(2));
                    obj->setLineWidth(0.5f);
                    obj->setColor(0,0.5,0);
                    theScene->insert( obj );
                }
            }
        }
    }

    // Represent point LMs
    if( hasPoints )
    {
        pointObj->clear();
        pointObj_local->clear();
        for( vector<MapPoint*>::const_iterator it = map->map_points.begin(); it!=map->map_points.end(); it++)
        {
            if( (*it)!=NULL )
            {
                if( (*it)->local )
                    pointObj_local->insertPoint( (*it)->point3D(0),(*it)->point3D(1),(*it)->point3D(2) );
                else
                    pointObj->insertPoint( (*it)->point3D(0),(*it)->point3D(1),(*it)->point3D(2) );
            }
        }
    }

    // Represent line LMs
    if( hasLines )
    {
        lineObj->clear();
        lineObj_local->clear();
        for( vector<MapLine*>::const_iterator it = map->map_lines.begin(); it!=map->map_lines.end(); it++)
        {
            if( (*it)!=NULL )
            {
                Vector6d L;
                L = (*it)->line3D;
                if( (*it)->local )
                    lineObj_local->appendLine( L(0),L(1),L(2),L(3),L(4),L(5) );
                else
                    lineObj->appendLine( L(0),L(1),L(2),L(3),L(4),L(5) );
            }
        }
    }

    // Update the image
    if(hasImg)
        image->setImageView_fast( img_mrpt_image );

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

}

// Setters

void slamScene::setText(int frame_, float time_, int nPoints_, int nPointsH_, int nLines_, int nLinesH_){
    frame    = frame_;
    time     = time_;
    nPoints  = nPoints_;
    nPointsH = nPointsH_;
    nLines   = nLines_;
    nLinesH  = nLinesH_;

}

void slamScene::setText(int frame_, float time_, int nPoints_, int nLines_ ){
    frame    = frame_;
    time     = time_;
    nPoints  = nPoints_;
    nLines   = nLines_;
}

void slamScene::setGT(Matrix4d xgt_){
    xgt = xgt_;
}

void slamScene::setComparison(Matrix4d xcomp_){
    xcomp = xcomp_;
}

void slamScene::setImage(const Mat &image_){

    Mat aux;
    Size img_sz(0.5*image_.cols, 0.5*image_.rows);
    cv::resize( image_, aux, img_sz );

    bool color;
    if (aux.channels() == 3) {
        aux.convertTo(aux, CV_8UC3);
        color = true;
    }
    else if (aux.channels() == 1) {
        aux.convertTo(aux, CV_8UC1);
        color = false;
    }
    else
        throw std::runtime_error(std::string("[SlamScene->setImage] unsupported image format: ") +
                                 std::to_string(aux.channels()));

    img_mrpt_image.loadFromMemoryBuffer(img_sz.width, img_sz.height, color, aux.data, false);
}

void slamScene::setImage(const string &image_){
    img_mrpt_image.loadFromFile(image_,1);
}

void slamScene::setLegend(){
    // Initialize the legend
    legend = theScene->createViewport("legend");
    if(hasLegend)
        legend->setViewportPosition(20, 900, 250, 97);
    else
        legend->setViewportPosition(2000, 2000, 250, 97);

    if(!hasGT) {
        if(hasComparison){
            img_legend = "../config/aux/legend_comp.png";
            img_mrpt_legend.loadFromFile(img_legend,1);
            legend->setImageView_fast( img_mrpt_legend );
        }
        else
            img_mrpt_legend.loadFromFile("",1);
    }
    else if(hasComparison){
        img_legend = "../config/aux/legend_full.png";
        img_mrpt_legend.loadFromFile(img_legend,1);
        legend->setImageView_fast( img_mrpt_legend );
    }
    else{
        img_legend = "../config/aux/legend.png";
        img_mrpt_legend.loadFromFile(img_legend,1);
        legend->setImageView_fast( img_mrpt_legend );
    }
}

void slamScene::setHelp(){
    // Initialize the legend
    help = theScene->createViewport("help");
    img_help = "../config/aux/help.png";
    img_mrpt_help.loadFromFile(img_help,1);
    help->setImageView_fast( img_mrpt_help );
    if(hasHelp)
        help->setViewportPosition(1600, 20, 300, 376);
    else
        help->setViewportPosition(2000, 2000, 300, 376);
}

void slamScene::setStereoCalibration(Matrix3d K_, float b_){
    f  = K_(0,0);
    cx = K_(0,2);
    cy = K_(1,2);
    b  = b_;
}

void slamScene::setPose(Matrix4d x_){
    x = x_;
}

void slamScene::setPoints(CMatrixFloat pData_){
    pData = pData_;
}

void slamScene::setLines(CMatrixFloat lData_){
    lData = lData_;
}

void slamScene::setKF(Matrix4d Tfw){
    // Initialize the camera object
    opengl::CSetOfObjectsPtr kfbb = opengl::stock_objects::BumblebeeCamera();
    {
        CPose3D pose( getPoseFormat(Tfw) );
        kfbb->setPose( pose );
        kfbb->setScale(sbb*10);
        theScene->insert(kfbb);
    }
}

// Public methods

bool slamScene::waitUntilClose(){
    while(win->isOpen());
    return true;
}

bool slamScene::isOpen(){
    return win->isOpen();
}

// Auxiliar methods

CPose3D slamScene::getPoseXYZ(VectorXd x){
    CPose3D pose(x(0),x(1),x(2),x(3),x(4),x(5));
    return pose;
}

CMatrixDouble slamScene::getPoseFormat(Matrix4d T){
    CMatrixDouble T_(4,4);
    for(unsigned int i = 0; i < 4; i++){
        for(unsigned int j = 0; j < 4; j++){
            T_(i,j) = T(i,j);
        }
    }
    return T_;
}

CMatrixDouble33 slamScene::getCovFormat(MatrixXd cov_){
    CMatrixDouble33 cov3;
    Matrix3d        cov3_eigen = cov_.block(0,0,3,3);

    for(unsigned int i = 0; i < 3; i++){
        for(unsigned int j = 0; j < 3; j++){
            cov3(i,j) = cov3_eigen(i,j);
        }
    }
    return cov3;
}

bool slamScene::getYPR(float &yaw, float &pitch, float &roll){
    double y, p, r;
    pose.getYawPitchRoll(y,p,r);
    yaw   = y;
    pitch = p;
    roll  = r;
}

bool slamScene::getPose(Matrix4d &T){
    CMatrixDouble44 T_;
    pose.getHomogeneousMatrix(T_);
    for(unsigned int i = 0; i < 4; i++){
        for(unsigned int j = 0; j < 4; j++){
            T(i,j) = T_(i,j);
        }
    }
}

}
