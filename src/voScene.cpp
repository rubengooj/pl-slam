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

#include <voScene.h>

// Auxiliar functions
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

// Constructors and destructor

voScene::voScene(){

    sbb     = 1.0f;
    sref    = 0.2f;
    srad    = 0.005f;
    sline   = 2.0f;
    saxis   = 0.5f;
    sfreq   = 3.0f;
    szoom   = 3.0f;
    selli   = 2.0f;
    selev   = 30.f;
    sazim   = -135.f;
    sfrust  = 0.2f;
    slinef  = 0.1f;
    win     = new CDisplayWindow3D("3D Scene",1920,1080);

    hasCamFix       = true;
    hasText         = true;
    hasAxes         = false;
    hasLegend       = false;
    hasHelp         = false;
    hasCov          = true;
    hasTraj         = true;
    hasGT           = false;
    hasChange       = false;
    hasComparison   = false;
    hasImg          = false;
    hasLines        = false;
    hasPoints       = false;
    hasFrustum      = false;

}

voScene::voScene(string configFile){

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
    slinef          = config.read_double("Scene","slinef",0.1f);
    win             = new CDisplayWindow3D("3D Scene",1920,1080);

    hasCamFix       = config.read_bool("Scene","hasCamFix",true);
    hasText         = config.read_bool("Scene","hasText",true);
    hasAxes         = config.read_bool("Scene","hasAxes",true);
    hasLegend       = config.read_bool("Scene","hasLegend",false);
    hasHelp         = config.read_bool("Scene","hasHelp",false);
    hasCov          = config.read_bool("Scene","hasCov",true);
    hasGT           = config.read_bool("Scene","hasGT",false);
    hasTraj         = config.read_bool("Scene","hasTraj",true);
    hasChange       = config.read_bool("Scene","hasChange",false);
    hasComparison   = config.read_bool("Scene","hasComparison",false);
    hasImg          = config.read_bool("Scene","hasImg",false);
    hasLines        = config.read_bool("Scene","hasLines",false);
    hasPoints       = config.read_bool("Scene","hasPoints",false);
    hasFrustum      = config.read_bool("Scene","hasFrustum",false);
    isKitti         = config.read_bool("Scene","isKitti",false);

    Matrix4d x_cw;
    x_cw << 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1;
    CPose3D x_aux(getPoseFormat(x_cw));
    pose_ini = x_aux;

}

voScene::~voScene(){

}

// Initialize the scene

void voScene::initializeScene(Matrix4d x_0, bool has_gt){

    // Initialize the scene and window
    win->setCameraElevationDeg(selev);
    win->setCameraAzimuthDeg(sazim);
    win->setCameraZoom(szoom);
    theScene = win->get3DSceneAndLock();

    hasGT = has_gt;

    // Initialize poses of different objects
    if(hasChange)
        change.setFromValues(0,0,0,  0, 0, -90*CV_PI/180.f);
    else
        change.setFromValues(0,0,0,  0, 0, 0);
    CPose3D x_aux(getPoseFormat(x_0));
    pose    = x_aux;
    pose1   = x_aux;
    pose_0  = x_aux;
    pose_gt = pose_ini - change;
    x_ini   = x_0;
    pose.getAsVector(v_aux);
    pose1.getAsVector(v_aux1);
    pose_gt.getAsVector(v_auxgt);

    // Initialize the camera object
    bbObj = opengl::stock_objects::BumblebeeCamera();
    {
        bbObj->setPose(pose_0);
        bbObj->setScale(sbb);
        theScene->insert(bbObj);
    }
    srefObj = opengl::stock_objects::CornerXYZSimple();
    srefObj->setPose(pose_0);
    srefObj->setScale(sref);
    theScene->insert(srefObj);
    frustumL_.setFromValues(0,0,0,  0, -90.f*CV_PI/180.f, -90.f*CV_PI/180.f);
    frustumR_.setFromValues(0.12f,0,0,  0, -90.f*CV_PI/180.f, -90.f*CV_PI/180.f);
    if(hasFrustum){
        frustObj = opengl::CFrustum::Create();
        {
            frustObj->setPose(pose_0+frustumL_);
            frustObj->setLineWidth (slinef);
            frustObj->setScale(sfrust);
            theScene->insert(frustObj);
        }
        frustObj1 = opengl::CFrustum::Create();
        {
            frustObj1->setPose(pose_0+frustumR_);
            frustObj1->setLineWidth (slinef);
            frustObj1->setScale(sfrust);
            theScene->insert(frustObj1);
        }
    }

    // Initialize the axes
    if(hasAxes){
        axesObj = opengl::CAxis::Create();
        axesObj->setFrequency(sfreq);
        axesObj->enableTickMarks(false);
        axesObj->setAxisLimits(-saxis,-saxis,-saxis, saxis,saxis,saxis);
        theScene->insert( axesObj );
    }

    // Initialize the ground truth camera object
    if( hasGT ){
        gtObj = opengl::stock_objects::BumblebeeCamera();
        {
            gtObj->setPose(pose_ini);
            gtObj->setScale(sbb);
            theScene->insert(gtObj);
        }
        srefObjGT = opengl::stock_objects::CornerXYZSimple();
        {
            srefObjGT->setPose(pose_ini);
            srefObjGT->setScale(sref);
            theScene->insert(srefObjGT);
        }
    }

    // Initialize a second camera to compare when changing parameters
    if(hasComparison){
        bbObj1 = opengl::stock_objects::BumblebeeCamera();
        {
            bbObj1->setPose(pose_0);
            bbObj1->setScale(sbb);
            theScene->insert(bbObj1);
        }
        srefObj1 = opengl::stock_objects::CornerXYZSimple();
        {
            srefObj1->setPose(pose_0);
            srefObj1->setScale(sref);
            theScene->insert(srefObj1);
        }
    }

    // Initialize the text
    if(hasText){
        string text = "Frame: \t \t0 \nFrequency: \t0 Hz \nLines:  \t0 (0)\nPoints: \t0 (0)";
        win->addTextMessage(0.85,0.95, text, TColorf(.0,.0,.0), 0, mrpt::opengl::MRPT_GLUT_BITMAP_HELVETICA_10 );
    }

    // Initialize the covariance ellipse
    ellPose = CPose3D(0.05f,0.f,0.f,0.f,0.f,0.f);
    elliObj = opengl::CEllipsoid::Create();
    elliObj->setPose(pose_0);
    elliObj->setScale(selli);
    elliObj->setQuantiles(2.0);
    elliObj->setColor(0,0.3,0);
    elliObj->enableDrawSolid3D(true);
    theScene->insert(elliObj);

    // Initialize the lines
    lineObj = opengl::CSetOfLines::Create();
    lineObj->setLineWidth(1.0);
    lineObj->setColor(0,0,0);
    theScene->insert( lineObj );

    // Initialize the point cloud
    pointObj = opengl::CPointCloud::Create();
    pointObj->setPointSize(3.0);
    pointObj->setColor(0,0,0);
    theScene->insert( pointObj );

    // Initialize the viewports
    setHelp();
    setLegend();

    image = theScene->createViewport("image");
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

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

}

// Update the scene

bool voScene::updateScene(){

    theScene = win->get3DSceneAndLock();
    bool restart = false;

    // Update camera pose
    CPose3D x_aux(getPoseFormat(x));
    pose = pose + x_aux;
    v_aux_ = v_aux;
    pose.getAsVector(v_aux);
    if(hasTraj){
        opengl::CSimpleLinePtr obj = opengl::CSimpleLine::Create();
        obj->setLineCoords(v_aux_(0),v_aux_(1),v_aux_(2), v_aux(0),v_aux(1),v_aux(2));
        obj->setLineWidth(sline);
        obj->setColor(0,0,0.7);
        theScene->insert( obj );
    }
    bbObj->setPose(pose);
    srefObj->setPose(pose);
    if(hasFrustum){
        frustObj->setPose(pose+frustumL_);
        frustObj1->setPose(pose+frustumR_);
    }

    // Set camera pointing to the current trajectory
    if(hasCamFix)
        win->setCameraPointingToPoint(v_aux(0),v_aux(1),v_aux(2));

    // Update the GT camera pose
    if(hasGT){
        CPose3D x_auxgt(getPoseFormat(xgt));
        //pose_gt = pose_gt + x_auxgt;
        pose_gt = x_auxgt;
        v_auxgt_ = v_auxgt;

        pose_gt.getAsVector(v_auxgt);
        float y_ = v_auxgt(1);
        float z_ = v_auxgt(2);
        float b_ = v_auxgt(4);
        float c_ = v_auxgt(5);
        v_auxgt(1) =  z_;
        v_auxgt(2) = -y_;        
        v_auxgt(4) = -c_;
        v_auxgt(5) =  b_;
        pose_gt = CPose3D(TPose3D(v_auxgt(0),v_auxgt(1),v_auxgt(2),v_auxgt(3),v_auxgt(4),v_auxgt(5)));


        if(hasTraj){
            opengl::CSimpleLinePtr obj = opengl::CSimpleLine::Create();
            obj->setLineCoords(v_auxgt_(0),v_auxgt_(1),v_auxgt_(2), v_auxgt(0),v_auxgt(1),v_auxgt(2));
            obj->setLineWidth(sline);
            obj->setColor(0,0,0);
            theScene->insert( obj );
        }
        gtObj->setPose(pose_gt);
        srefObjGT->setPose(pose_gt);
    }

    // Update the comparison camera pose
    if(hasComparison){
        CPose3D x_aux1(getPoseFormat(xcomp));
        pose1 = pose1 + x_aux1;
        v_aux1_ = v_aux1;
        pose1.getAsVector(v_aux1);
        if(hasTraj){
            opengl::CSimpleLinePtr obj = opengl::CSimpleLine::Create();
            obj->setLineCoords(v_aux1_(0),v_aux1_(1),v_aux1_(2), v_aux1(0),v_aux1(1),v_aux1(2));
            obj->setLineWidth(sline);
            obj->setColor(0,0.7,0);
            theScene->insert( obj );
        }
        bbObj1->setPose(pose1);
        srefObj1->setPose(pose1);
    }

    // Update the text
    if(hasText){
        string text = "Frame: \t \t" + to_string(frame) + " \n" + "Frequency: \t" + to_string_with_precision(1000.f/time,4) + " Hz \n" + "Lines:  \t" + to_string(nLines) + " (" + to_string(nLinesH) + ") \nPoints: \t" + to_string(nPoints) + " (" + to_string(nPointsH) + ")";
        win->addTextMessage(0.85,0.95, text, TColorf(.0,.0,.0), 0, mrpt::opengl::MRPT_GLUT_BITMAP_HELVETICA_10 );
    }

    // Update the covariance
    if(hasCov){
        elliObj->setPose(pose+ellPose);
        elliObj->setCovMatrix(getCovFormat(cov));
    }

    // Update the image
    if(hasImg){
        image->setImageView_fast( img_mrpt_image );
        //image->setImageView( img_mrpt_image );
    }

    // Update the lines
    lineObj->clear();
    if(hasLines){
        plotLinesCovariances();
        for(int i = 0; i < size(lData,2); i++)
            lineObj->appendLine(lData(0,i),lData(1,i),lData(2,i), lData(3,i),lData(4,i),lData(5,i));
        lineObj->setPose(pose);
    }

    // Update the points
    pointObj->clear();
    if(hasPoints){
        plotPointsCovariances();
        for(int i = 0; i < size(pData,2); i++)
            pointObj->insertPoint(pData(0,i),pData(1,i),pData(2,i));
        pointObj->setPose(pose);
    }

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

    // Key events   -       TODO: change the trick to employ viewports
    if(win->keyHit()){       
        key = win->getPushedKey(&kmods);
        if(key == MRPTK_SPACE){                     // Space    Reset VO
            theScene->clear();
            initializeScene(x_ini,false);
        }
        else if (key == MRPTK_ESCAPE){              // Esc      Restart VO
            theScene->clear();
            initializeScene(x_ini,false);
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
        else if ( (key == 102) || (key == 70) ){    // F        frustum
            hasFrustum   = !hasFrustum;
            if(!hasFrustum){
                frustObj.clear();
                frustObj1.clear();
            }
            else{
                frustObj = opengl::CFrustum::Create();
                frustObj->setLineWidth (slinef);
                frustObj->setScale(sfrust);
                frustObj->setPose(pose+frustumL_);
                theScene->insert(frustObj);
                frustObj1 = opengl::CFrustum::Create();
                frustObj1->setLineWidth (slinef);
                frustObj1->setScale(sfrust);
                frustObj1->setPose(pose+frustumR_);
                theScene->insert(frustObj1);
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
        else if ( (key ==  99) || (key == 67) ){    // C        covariance
            hasCov   = !hasCov;
            if(!hasCov){
                elliObj.clear();
            }
            else{
                elliObj = opengl::CEllipsoid::Create();
                elliObj->setScale(selli);
                elliObj->setQuantiles(2.0);
                elliObj->setColor(0,0.3,0);
                elliObj->enableDrawSolid3D(true);
                elliObj->setPose(pose);
                theScene->insert(elliObj);
            }
        }
        else if (  key == 43){                      // +        Increases the scale of the covariance ellipse
            selli += 0.5f;
            if(hasCov)
                elliObj->setScale(selli);
        }
        else if (  key == 45){                      // -        Decreases the scale of the covariance ellipse
            selli -= 0.5f;
            if(hasCov)
                elliObj->setScale(selli);
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

bool voScene::updateScene( list<PointFeature*> matched_pt ){

    theScene = win->get3DSceneAndLock();
    bool restart = false;

    // Update camera pose
    CPose3D x_aux(getPoseFormat(x));
    pose = pose + x_aux;
    v_aux_ = v_aux;
    pose.getAsVector(v_aux);
    if(hasTraj){
        opengl::CSimpleLinePtr obj = opengl::CSimpleLine::Create();
        obj->setLineCoords(v_aux_(0),v_aux_(1),v_aux_(2), v_aux(0),v_aux(1),v_aux(2));
        obj->setLineWidth(sline);
        obj->setColor(0,0,0.7);
        theScene->insert( obj );
    }
    bbObj->setPose(pose);
    srefObj->setPose(pose);
    if(hasFrustum){
        frustObj->setPose(pose+frustumL_);
        frustObj1->setPose(pose+frustumR_);
    }

    // Set camera pointing to the current trajectory
    if(hasCamFix)
        win->setCameraPointingToPoint(v_aux(0),v_aux(1),v_aux(2));

    // Update the GT camera pose
    if(hasGT){
        CPose3D x_auxgt(getPoseFormat(xgt));
        //pose_gt = pose_gt + x_auxgt;
        pose_gt = x_auxgt;
        v_auxgt_ = v_auxgt;

        pose_gt.getAsVector(v_auxgt);
        float y_ = v_auxgt(1);
        float z_ = v_auxgt(2);
        float b_ = v_auxgt(4);
        float c_ = v_auxgt(5);
        v_auxgt(1) =  z_;
        v_auxgt(2) = -y_;
        v_auxgt(4) = -c_;
        v_auxgt(5) =  b_;
        pose_gt = TPose3D(v_auxgt(0),v_auxgt(1),v_auxgt(2),v_auxgt(3),v_auxgt(4),v_auxgt(5));

        if(hasTraj){
            opengl::CSimpleLinePtr obj = opengl::CSimpleLine::Create();
            obj->setLineCoords(v_auxgt_(0),v_auxgt_(1),v_auxgt_(2), v_auxgt(0),v_auxgt(1),v_auxgt(2));
            obj->setLineWidth(sline);
            obj->setColor(0,0,0);
            theScene->insert( obj );
        }
        gtObj->setPose(pose_gt);
        srefObjGT->setPose(pose_gt);
    }

    // Update the comparison camera pose
    if(hasComparison){
        CPose3D x_aux1(getPoseFormat(xcomp));
        pose1 = pose1 + x_aux1;
        v_aux1_ = v_aux1;
        pose1.getAsVector(v_aux1);
        if(hasTraj){
            opengl::CSimpleLinePtr obj = opengl::CSimpleLine::Create();
            obj->setLineCoords(v_aux1_(0),v_aux1_(1),v_aux1_(2), v_aux1(0),v_aux1(1),v_aux1(2));
            obj->setLineWidth(sline);
            obj->setColor(0,0.7,0);
            theScene->insert( obj );
        }
        bbObj1->setPose(pose1);
        srefObj1->setPose(pose1);
    }

    // Update the text
    if(hasText){
        string text = "Frame: \t \t" + to_string(frame) + " \n" + "Frequency: \t" + to_string_with_precision(1000.f/time,4) + " Hz \n" + "Lines:  \t" + to_string(nLines) + " (" + to_string(nLinesH) + ") \nPoints: \t" + to_string(nPoints) + " (" + to_string(nPointsH) + ")";
        win->addTextMessage(0.85,0.95, text, TColorf(.0,.0,.0), 0, mrpt::opengl::MRPT_GLUT_BITMAP_HELVETICA_10 );
    }

    // Update the covariance
    if(hasCov){
        elliObj->setPose(pose+ellPose);
        elliObj->setCovMatrix(getCovFormat(cov));
    }

    // Update the image
    if(hasImg){
        image->setImageView_fast( img_mrpt_image );
        //image->setImageView( img_mrpt_image );
    }

    // Update the lines
    lineObj->clear();
    if(hasLines){
        plotLinesCovariances();
        for(int i = 0; i < size(lData,2); i++)
            lineObj->appendLine(lData(0,i),lData(1,i),lData(2,i), lData(3,i),lData(4,i),lData(5,i));
        lineObj->setPose(pose);
    }

    // Update the points
    pointObj->clear();
    if(hasPoints){
        /*CPointCloudColouredPtr map_points = CPointCloudColoured::Create();
        map_points->setPointSize(2);
        map_points->enablePointSmooth(1);
        map_points->setPose(pose);

        for( auto it = matched_pt.begin(); it != matched_pt.end(); it++ ){
            map_points->setPo
        }
        theScene->insert( map_points );*/
        plotPointsCovariances();
        for(int i = 0; i < size(pData,2); i++)
            pointObj->insertPoint(pData(0,i),pData(1,i),pData(2,i));
        pointObj->setPose(pose);
    }

    // Re-paint the scene
    win->unlockAccess3DScene();
    win->repaint();

    // Key events   -       TODO: change the trick to employ viewports
    if(win->keyHit()){
        key = win->getPushedKey(&kmods);
        if(key == MRPTK_SPACE){                     // Space    Reset VO
            theScene->clear();
            initializeScene(x_ini,false);
        }
        else if (key == MRPTK_ESCAPE){              // Esc      Restart VO
            theScene->clear();
            initializeScene(x_ini,false);
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
        else if ( (key == 102) || (key == 70) ){    // F        frustum
            hasFrustum   = !hasFrustum;
            if(!hasFrustum){
                frustObj.clear();
                frustObj1.clear();
            }
            else{
                frustObj = opengl::CFrustum::Create();
                frustObj->setLineWidth (slinef);
                frustObj->setScale(sfrust);
                frustObj->setPose(pose+frustumL_);
                theScene->insert(frustObj);
                frustObj1 = opengl::CFrustum::Create();
                frustObj1->setLineWidth (slinef);
                frustObj1->setScale(sfrust);
                frustObj1->setPose(pose+frustumR_);
                theScene->insert(frustObj1);
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
        else if ( (key ==  99) || (key == 67) ){    // C        covariance
            hasCov   = !hasCov;
            if(!hasCov){
                elliObj.clear();
            }
            else{
                elliObj = opengl::CEllipsoid::Create();
                elliObj->setScale(selli);
                elliObj->setQuantiles(2.0);
                elliObj->setColor(0,0.3,0);
                elliObj->enableDrawSolid3D(true);
                elliObj->setPose(pose);
                theScene->insert(elliObj);
            }
        }
        else if (  key == 43){                      // +        Increases the scale of the covariance ellipse
            selli += 0.5f;
            if(hasCov)
                elliObj->setScale(selli);
        }
        else if (  key == 45){                      // -        Decreases the scale of the covariance ellipse
            selli -= 0.5f;
            if(hasCov)
                elliObj->setScale(selli);
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

void voScene::plotPointsCovariances(){

    unsigned int nP = size(pData,2);
    Matrix3d covAux;

    // Keypoints covariance ellipsoids
    elliObjP->clear();
    for(int i = 0; i < nP; i++){
        // Estimation of the uncertainty    (TODO: save the landmarks uncertainties into a structure in the Optimization class)
        float px_hat = pData(3,i) - cx;
        float py_hat = pData(4,i) - cy;
        float disp   = pData(3,i) - pData(5,i);
        float disp2  = disp * disp;
        Matrix3d covP_an;
        covP_an(0,0) = disp2+2.f*px_hat*px_hat;
        covP_an(0,1) = 2.f*px_hat*py_hat;
        covP_an(0,2) = 2.f*f*px_hat;
        covP_an(1,1) = disp2+2.f*py_hat*py_hat;
        covP_an(1,2) = 2.f*f*py_hat;
        covP_an(2,2) = 2.f*f*f;
        covP_an(1,0) = covP_an(0,1);
        covP_an(2,0) = covP_an(0,2);
        covP_an(2,1) = covP_an(1,2);
        covP_an << covP_an * bsigmaP / (disp2*disp2);
        // Insertion of the ellipsoids
        CEllipsoidPtr elliAux_ = opengl::CEllipsoid::Create();        
        elliAux_->setQuantiles(1.0);
        elliAux_->enableDrawSolid3D(true);
        elliAux_->setCovMatrix(getCovFormat(covP_an));
        elliAux_->setPose(CPose3D(pData(0,i),pData(1,i),pData(2,i),0,0,0));
        elliAux_->setColor(0.3,0,0);
        elliAux_->setScale(selli/5.f);
        elliObjP->insert(elliAux_);
    }
    elliObjP->setPose(pose);

}

void voScene::plotLinesCovariances(){

    unsigned int nL = size(lData,2);
    Matrix3d covAux;

    // Line-segments covariance ellipsoids
    elliObjL->clear();
    for(int i = 0; i < nL; i++){
        // First endpoint
        // ====================================================================================================================
        // Estimation of the uncertainty    (TODO: save the landmarks uncertainties into a structure in the Optimization class)
        float px_hat = lData(13,i) - cx;
        float py_hat = lData(14,i) - cy;
        float disp   = lData(13,i) - lData(15,i);
        float disp2  = disp * disp;
        Matrix3d covP_an;
        covP_an(0,0) = disp2+2.f*px_hat*px_hat;
        covP_an(0,1) = 2.f*px_hat*py_hat;
        covP_an(0,2) = 2.f*f*px_hat;
        covP_an(1,1) = disp2+2.f*py_hat*py_hat;
        covP_an(1,2) = 2.f*f*py_hat;
        covP_an(2,2) = 2.f*f*f;
        covP_an(1,0) = covP_an(0,1);
        covP_an(2,0) = covP_an(0,2);
        covP_an(2,1) = covP_an(1,2);
        covP_an << covP_an * bsigmaL / (disp2*disp2);
        // Insertion of the ellipsoids
        CEllipsoidPtr elliAux_ = opengl::CEllipsoid::Create();
        elliAux_->setQuantiles(1.0);
        elliAux_->enableDrawSolid3D(true);
        elliAux_->setColor(0,0.3,0);
        elliAux_->setScale(selli/3.f);
        elliAux_->setCovMatrix(getCovFormat(covP_an));
        elliAux_->setPose(CPose3D(lData(0,i),lData(1,i),lData(2,i),0,0,0));
        elliObjL->insert(elliAux_);
        // Second endpoint
        // ====================================================================================================================
        // Estimation of the uncertainty    (TODO: save the landmarks uncertainties into a structure in the Optimization class)
        px_hat = lData(16,i) - cx;
        py_hat = lData(17,i) - cy;
        disp   = lData(16,i) - lData(18,i);
        disp2  = disp * disp;
        Matrix3d covQ_an;
        covQ_an(0,0) = disp2+2.f*px_hat*px_hat;
        covQ_an(0,1) = 2.f*px_hat*py_hat;
        covQ_an(0,2) = 2.f*f*px_hat;
        covQ_an(1,1) = disp2+2.f*py_hat*py_hat;
        covQ_an(1,2) = 2.f*f*py_hat;
        covQ_an(2,2) = 2.f*f*f;
        covQ_an(1,0) = covQ_an(0,1);
        covQ_an(2,0) = covQ_an(0,2);
        covQ_an(2,1) = covQ_an(1,2);
        covQ_an << covQ_an * bsigmaL / (disp2*disp2);
        // Insertion of the ellipsoids
        elliAux_ = opengl::CEllipsoid::Create();
        elliAux_->setQuantiles(1.0);
        elliAux_->enableDrawSolid3D(true);
        elliAux_->setColor(0,0.3,0);
        elliAux_->setScale(selli/3.f);
        elliAux_->setCovMatrix(getCovFormat(covQ_an));
        elliAux_->setPose(CPose3D(lData(3,i),lData(4,i),lData(5,i),0,0,0));
        elliObjL->insert(elliAux_);
    }
    elliObjL->setPose(pose);

}

// Setters

void voScene::setText(int frame_, float time_, int nPoints_, int nPointsH_, int nLines_, int nLinesH_){
    frame    = frame_;
    time     = time_;
    nPoints  = nPoints_;
    nPointsH = nPointsH_;
    nLines   = nLines_;
    nLinesH  = nLinesH_;

}

void voScene::setCov(MatrixXd cov_){
    cov = cov_;
}

void voScene::setPose(Matrix4d x_){
    x = x_;
}

void voScene::setGT(Matrix4d xgt_){
    xgt = xgt_;
}

void voScene::setComparison(Matrix4d xcomp_){
    xcomp = xcomp_;
}

void voScene::setImage(Mat image_){
    IplImage iplimage = image_;
    //IplImage* iplimage = new IplImage(image_);
    img_mrpt_image.loadFromIplImage( &iplimage );
}

void voScene::setImage(string image_){
    img_mrpt_image.loadFromFile(image_,1);
}

void voScene::setLegend(){
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

void voScene::setHelp(){
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

void voScene::setPoints(CMatrixFloat pData_){
    pData = pData_;
}

void voScene::setLines(CMatrixFloat lData_){
    lData = lData_;
}

void voScene::setStereoCalibration(Matrix3d K_, float b_){
    f  = K_(0,0);
    cx = K_(0,2);
    cy = K_(1,2);
    b  = b_;
    sigmaP = 1.f;   // TODO: READ FROM CONFIG FILE
    sigmaL = 1.f;

    bsigmaL = b*b*sigmaL*sigmaL;
    bsigmaP = b*b*sigmaP*sigmaP;
}

void voScene::setKF(){
    // Initialize the camera object (set KF with the current pose)
    opengl::CSetOfObjectsPtr kfbb = opengl::stock_objects::BumblebeeCamera();
    {
        kfbb->setPose( pose );
        kfbb->setScale(sbb*3);
        theScene->insert(kfbb);
    }
}

void voScene::setKF(Matrix4d Tfw){
    // Initialize the camera object
    opengl::CSetOfObjectsPtr kfbb = opengl::stock_objects::BumblebeeCamera();
    {
        //CPose3D pose( getPoseFormat(Tfw) );
        kfbb->setPose( pose );
        kfbb->setScale(sbb*10);
        theScene->insert(kfbb);
    }
}

// Public methods

bool voScene::waitUntilClose(){
    while(win->isOpen());
    return true;
}

bool voScene::isOpen(){
    return win->isOpen();
}

// Auxiliar methods

CPose3D voScene::getPoseXYZ(VectorXd x){
    CPose3D pose(x(0),x(1),x(2),x(3),x(4),x(5));
    return pose;
}

CMatrixDouble voScene::getPoseFormat(Matrix4d T){
    CMatrixDouble T_(4,4);
    for(unsigned int i = 0; i < 4; i++){
        for(unsigned int j = 0; j < 4; j++){
            T_(i,j) = T(i,j);
        }
    }
    return T_;
}

CMatrixDouble33 voScene::getCovFormat(MatrixXd cov_){
    CMatrixDouble33 cov3;
    Matrix3d        cov3_eigen = cov_.block(0,0,3,3);

    for(unsigned int i = 0; i < 3; i++){
        for(unsigned int j = 0; j < 3; j++){
            cov3(i,j) = cov3_eigen(i,j);
        }
    }
    return cov3;
}

bool voScene::getYPR(float &yaw, float &pitch, float &roll){
    double y, p, r;
    pose.getYawPitchRoll(y,p,r);
    yaw   = y;
    pitch = p;
    roll  = r;
}

bool voScene::getPose(Matrix4d &T){
    CMatrixDouble44 T_;
    pose.getHomogeneousMatrix(T_);
    for(unsigned int i = 0; i < 4; i++){
        for(unsigned int j = 0; j < 4; j++){
            T(i,j) = T_(i,j);
        }
    }
}

