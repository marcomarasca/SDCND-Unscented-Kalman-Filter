# Unscented Kalman Filter
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

[//]: # (Image References)

[dataset_1_gif]: ./images/dataset1.gif "Tracking on the first dataset"
[position_img]: ./images/position.png "Estimated vs Measured vs Ground Truth Position"
[position_zoomed_img]: ./images/position_zoomed.png "Position Zoomed"
[velocity_img]: ./images/velocity.png "Esitmated vs Ground Truth Velocity"
[yaw_img]: ./images/yaw.png "Esitmated vs Ground Truth Yaw"
[yawrate_img]: ./images/yaw_rate.png "Esitmated vs Ground Truth Yaw Rate"
[nis_radar_img]: ./images/nis_radar.png "Normalized Innovation Squared (NIS) Radar"
[nis_laser_img]: ./images/nis_laser.png "Normalized Innovation Squared (NIS) Laser"
[nis_radar_only_img]: ./images/nis_radar_only.png "Normalized Innovation Squared (NIS) Radar only"
[nis_laser_only_img]: ./images/nis_laser_only.png "Normalized Innovation Squared (NIS) Laser only"

![alt text][dataset_1_gif]

Overview
---

This repository contains a C++ implementation of the (unscented) kalman filter used to estimate the state of a moving object of interest through measurements coming from lidar and radar sensors that are fused together to perform an estimate. The measurements comes from a simulated environment thanks to the [Udacity Simulator](https://github.com/udacity/self-driving-car-sim) and are fed to the program through [WebSockets](https://en.wikipedia.org/wiki/WebSocket) messages. The [main](./src/main.cpp) file process the incoming messages and parse the type of measurement that is then processed by the [UKF](./src/ukf.cpp) class.

Each message contains the type of sensor measurement (LASER or RADAR), the measurement data and a timestamp, plus the ground truth values (e.g. real values of the position of the tracked object). The messages are encoded as JSON objects and the "sensor_measurement" attribute will contain the raw data from the sensor, the data is encoded as a space separated values string, whose first character encodes the type of measurement, **L** for Laser and **R** for Radar. The values following are the data plus ground truth.

#### Radar

```sensor_type, rho_measured, phi_measured, rhodot_measured, timestamp, x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth```

sensor_type is **R**

#### Lidar

```sensor_type, x_measured, y_measured, timestamp, x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth```

sensor_type is **L**

The program process this data and sends back a JSON object containing the estimated position and the root mean squared error values for the position and velocity (computed using the estimation and ground truth):

```
["estimate_x"] <= Estimated position x
["estimate_y"] <= Estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]
```

As we can see the coordinates coming from the 2 sensors are different in the type of representation, while the lidar outputs cartesian coordinates, the radar provides polar coordinates. This is where the unscented (and extended) kalman filter comes into place where the non-linear data is used (for the radar) and needs to deal with this non-linear process and measurement models. While the [extended version](https://github.com/Az4z3l/CarND-Extended-Kalman-Filter) uses a linear approximation for the computation though a multivariate Taylor Series expansions to linearize a model, the unscented version of the kalman filter takes a different approach without trying to linearize the model, but rather generating a set of (*sigma*) points that approximate the distribution and that are then converted to the measurement space for prediction. 

The model used in this implementation is the *Constant Turn Rate and Velocity* magnitude model (CTRV) that, as the name implies assumes a constant rate of turn and velocity. The unscented kalman filter is not only more precise than the extended version, but it's also computationally less expensive than the extended counterpart.

Results
---

The [results](./extra/results.txt) of the estimations can be compared with the ground truth values supplied by the simulator computing the *Root Mean Squared Error* (RMSE), in the following some results for the first dataset:

|                   | Px     | Py     | Vx     | Vy     |
|-------------------|--------|--------|--------|--------|
| **Radar + Lidar** | 0.0656 | 0.0811 | 0.1500 | 0.1578 |
| **Radar Only**    | 0.1452 | 0.2288 | 0.1877 | 0.2117 |
| **Lidar Only**    | 0.0902 | 0.0954 | 0.2085 | 0.2002 |

As expected the results using both sensors together are sensibly better than using a single sensor at the time. Moreover we can compare the results against the [extended version](https://github.com/Az4z3l/CarND-Extended-Kalman-Filter) of the filter:

|         | Px     | Py     | Vx     | Vy     |
|---------|--------|--------|--------|--------|
| **UKF** | 0.0656 | 0.0811 | 0.1500 | 0.1578 |
| **EKF** | 0.0964 | 0.0853 | 0.4154 | 0.4316 |

Note that in this implementation the *laser updates* are performed using the standard kalman filter equations rather than using the UKF ones (configurable, see line 24 of [ukf.cpp](./src/ukf.cpp)). This because the measurements from the laser are linear and the estimation can be efficiently computed avoiding the UKF transformations (note however that the result of the updates should be the same for both the UKF and KF).

During the run, the program also [outputs](./extra/output.txt) the results at the various steps saving to a (tab separated) file the following data:

| Name        | Description                               |
|-------------|-------------------------------------------|
| px_est      | Estimated Px                              |
| py_est      | Estimated Py                              |
| v_est       | Estimated velocity                        |
| vx_est      | Estimated velocity x                      |
| vy_est      | Estimated velocity y                      |
| yaw_est     | Estimated yaw                             |
| yawrate_est | Estimated yaw rate                        |
| px_meas     | Measured Px                               |
| py_meas     | Measured Py                               |
| px_gt       | Ground Truth Px                           |
| py_gt       | Ground Truth Py                           |
| v_gt        | Ground Truth Velocity                     |
| vx_gt       | Ground Truth Velocity x                   |
| vy_gt       | Ground Truth Velocity y                   |
| yaw_gt      | Ground Truth Yaw                          |
| yawrate_gt  | Ground Truth Yaw Rate                     |
| NIS_laser   | Normalized Innovation Squared (NIS) Laser |
| NIS_radar   | Normalized Innovation Squared (NIS) radar |

The data can be visualized using the [provided jupyter notebook](./extra/ukf_visualization.ipynb), in particular the NIS value can be used to perform a consistency check on the process noise standard deviation of the longitudinal acceleration and yaw acceleration. Comparing against the [Chi-squared distribution](https://en.wikipedia.org/wiki/Chi-squared_distribution) using as reference the 95% for the appropriate degrees of freedom (the dimension of the measurement space) we use 7.815 for radar (3 dimensions) and 5.991 for laser (2 dimensions):

![alt text][nis_radar_img]![alt text][nis_laser_img]

As we can see the NIS value gives us an indication that our parameters choice was realistic, in other words in about 5% of the cases our NIS value is over the reference points (7.815 and 5.991).

We can run the kalman filter excluding the radar or lidar sensor measurements (See lines 17-20 of [ukf.cpp](./src/ukf.cpp)) and check again the results:

![alt text][nis_radar_only_img]![alt text][nis_laser_only_img]

To further check our implementation given that we have the ground truth values, we can compare our estimates directly agains them:

![alt text][position_img]![alt text][position_zoomed_img]

We can visually see that the position estimated by the filter is close to the ground truth, we also have information about the velocity, the yaw and yaw rate:

![alt text][velocity_img]
![alt text][yaw_img]
![alt text][yawrate_img]

Getting Started
---

In order to run the program you need the simulator provided by [Udacity](https://www.udacity.com/) which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even better [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. The version compatible with the simulator is the uWebSocketIO branch **e94b6e1**.

Once the install for uWebSocketIO is complete, the main program can be built and run by doing the following from the project top directory.

1. ```mkdir build```
2. ```cd build```
3. ```cmake .. && make```
4. ```./UnscentedKF```

Note that to compile the program with debug symbols you can supply the appropriate flag to cmake: ```cmake -DCMAKE_BUILD_TYPE=Debug .. && make```.

Now the Udacity simulator can be run selecting the EKF/UKF project, after the dataset selection press start and see the application in action.

![alt text][dataset_1_gif]

#### Other Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

Environment Setup
---

This project was developed under windows using the windows subsystem for linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)) with Ubuntu Bash 16.04 together with [Visual Studio Code](https://code.visualstudio.com/).

The steps to setup the environment under mac, linux or windows (WSL) are more or less the same:

- Review the above dependencies
- Clone the repo and run the appropriate script (./install-ubuntu.sh under WSL and linux and ./install-mac.sh under mac), this should install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) from the branch **e94b6e1**

Under windows (WSL) and linux you can make a clean installation as follows:

1. ```sudo apt-get update```
2. ```sudo apt-get install git```
3. ```sudo apt-get install cmake```
4. ```sudo apt-get install openssl```
5. ```sudo apt-get install libssl-dev```
6. ```git clone https://github.com/Az4z3l/CarND-Unscented-Kalman-Filter```
7. ```sudo rm /usr/lib/libuWS.so```
8. ```./install-ubuntu.sh```

#### Debugging with VS Code

Since I developed this project using WSL and Visual Studio Code it was very useful for me to setup a debugging pipeline. VS Code comes with a official Microsoft cpp extension that can be downloaded directly from the marketplace: https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools. After the installation there are a few things to setup in order to make it work with the subsystem for linux, personally I went with the default Ubuntu distribution.

For the following setup I assume that the repository was cloned in **D:/Dev/CarND-Unscented-Kalman-Filter/**.

##### Setup the language server (for IntelliSense)

From the official documentation [https://github.com/Microsoft/vscode-cpptools/blob/master/Documentation/LanguageServer/Windows%20Subsystem%20for%20Linux.md](https://github.com/Microsoft/vscode-cpptools/blob/master/Documentation/LanguageServer/Windows%20Subsystem%20for%20Linux.md): 

Simply Crtl+P and select "C/Cpp: Edit Configurations", this will create a c_cpp_properties.json file that can be configured as follows:

```json
{
    "name": "WSL",
    "intelliSenseMode": "clang-x64",
    "compilerPath": "/usr/bin/gcc",
    "includePath": [
        "${workspaceFolder}"
    ],
    "defines": [],
    "browse": {
        "path": [
            "${workspaceFolder}"
        ],
        "limitSymbolsToIncludedHeaders": true,
        "databaseFilename": ""
    },
    "cStandard": "c11",
    "cppStandard": "c++17"
}
```

##### Setup the Debugger

From the official documenation [https://github.com/Microsoft/vscode-cpptools/blob/master/Documentation/Debugger/gdb/Windows%20Subsystem%20for%20Linux.md](https://github.com/Microsoft/vscode-cpptools/blob/master/Documentation/Debugger/gdb/Windows%20Subsystem%20for%20Linux.md):

First install gdb in the WSL:

```
sudo apt install gdb
```

Then simply create a lunch configuration from VS Code: "Debug" -> "Add Configuration.." and setup the launch.json as follows:

```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "C++ Launch",
            "type": "cppdbg",
            "request": "launch",
            "program": "/mnt/d/Dev/CarND-Unscented-Kalman-Filter/build/UnscentedKF",
            "args": ["-fThreading"],
            "stopAtEntry": false,
            "cwd": "/mnt/d/Dev/CarND-Unscented-Kalman-Filter/build/",
            "environment": [],
            "externalConsole": true,
            "windows": {
                "MIMode": "gdb",
                "setupCommands": [
                    {
                        "description": "Enable pretty-printing for gdb",
                        "text": "-enable-pretty-printing",
                        "ignoreFailures": true
                    }
                ]
            },
            "pipeTransport": {
                "pipeCwd": "",
                "pipeProgram": "c:\\windows\\sysnative\\bash.exe",
                "pipeArgs": ["-c"],
                "debuggerPath": "/usr/bin/gdb"
            },
            "sourceFileMap": {
                "/mnt/d": "d:\\"
            }
        }
    ]
}
```

Note how the program is mapped directly into the file system of the WSL and piped through bash.exe (the paths are relative to the WSL environment).

Now you are ready to debug the application directly from VS Code, simply compile the application from within the WSL with the debug symbols:

```cmake -DCMAKE_BUILD_TYPE=Debug .. && make```

And run the debugger from VS Code (e.g. F5) :)
