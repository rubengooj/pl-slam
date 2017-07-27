# PL-SLAM #

This code contains an algorithm to compute stereo visual SLAM by using both point and line segment features.

**Authors:** [Ruben Gomez-Ojeda](http://mapir.isa.uma.es/mapirwebsite/index.php/people/164-ruben-gomez), [Francisco Angel Moreno](http://mapir.isa.uma.es/mapirwebsite/index.php/people/199-francisco-moreno-due%C3%B1as), [Davide Scaramuzza](http://rpg.ifi.uzh.ch/people_scaramuzza.html), and [Javier Gonzalez-Jimenez](http://mapir.isa.uma.es/mapirwebsite/index.php/people/95-javier-gonzalez-jimenez)

**Related publication:** [*PL-SLAM: a Stereo SLAM System through the Combination of Points and Line Segments*](http://mapir.isa.uma.es/mapirwebsite/index.php/people/164-ruben-gomez)

If you use PL-SLAM in your research work, please cite:

    @article{gomez2017pl,
      title   = {{PL-SLAM: a Stereo SLAM System through the Combination of Points and Line Segments}},
      author  = {Gomez-Ojeda, Ruben and Moreno, Francisco-Angel and Scaramuzza, Davide and Gonzalez-Jimenez, Javier},
      journal = {arXiv preprint arXiv:1705.09479},
      year    = {2017}
}

The pdf file can be found at [https://arxiv.org/abs/1705.09479](https://arxiv.org/abs/1705.09479).

[![PL-SLAM](https://img.youtube.com/vi/-lCTf_tAxhQ/0.jpg)](https://www.youtube.com/watch?v=-lCTf_tAxhQ)

**Previous publications:**

[Gomez-Ojeda, R., Briales, J., & Gonzalez-Jimenez, J. (2016, October). PL-SVO: Semi-direct monocular visual odometry by combining points and line segments. In Intelligent Robots and Systems (IROS), 2016 IEEE/RSJ International Conference on (pp. 4211-4216). IEEE.](http://mapir.isa.uma.es/rgomez/publications/iros16plsvo.pdf)

[Gomez-Ojeda, R., & Gonzalez-Jimenez, J. (2016, May). Robust stereo visual odometry through a probabilistic combination of points and line segments. In Robotics and Automation (ICRA), 2016 IEEE International Conference on (pp. 2521-2526). IEEE.](http://mapir.isa.uma.es/rgomez/publications/icra16plsvo.pdf).

**License:**

The provided code is published under the General Public License Version 3 (GPL v3). More information can be found in the "LICENSE" also included in the repository.

Please do not hesitate to contact the authors if you have any further questions.


## 1. Prerequisites and dependencies

### OpenCV 3.0.0
It can be easily found at http://opencv.org. 
In the case of line segments, we have modified the *line_descriptor* from the *opencv_contrib* 
[repository](https://github.com/Itseez/opencv_contrib), included in the *3rdparty* folder.

### Eigen3
http://eigen.tuxfamily.org

### Boost
Installation on Ubuntu:
```
sudo apt-get install libboost-dev
```

### YAML
Installation on Ubuntu:
```
sudo apt-get install libyaml-dev
```

### MRPT (Optional)
In case of using the provided representation. 
```
sudo apt-get install libmrpt-dev
```

Download and install instructions can be also found at: http://www.mrpt.org/ .

### Line descriptor (in 3rdparty folder)
We have modified the [*line_descriptor*](https://github.com/opencv/opencv_contrib/tree/master/modules/line_descriptor) module from the [OpenCV/contrib](https://github.com/opencv/opencv_contrib) library (both BSD) which is included in the *3rdparty* folder.


## 2. Configuration and generation

Executing the file *build.sh* will configure and generate the *line_descriptor* and *DBoW2* modules, uncompress the vocabulary files, and then will configure and generate the *PL-SLAM* library for which we generate: **libplslam.so** in the lib folder, and the applications **plstvo_dataset** and **plslam_dataset** that works with our dataset format (explained in the next section).


## 3. Dataset format and usage

The **plslam_dataset** (and **plstvo_dataset**)  basic usage is: 
```
./plslam_dataset  <dataset_path>  
```

where *<dataset_path>* refers to the sequence folder relative to the environment variable *${DATASETS_DIR}* that must be previously set. That sequence folder must contain the dataset configuration file named **dataset_params.yaml** following the examples in **pl-slam/config**, where **images_subfolder_{lr}** refers to the left and right image subfolders.

