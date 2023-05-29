# Otter_USV_Path_Following

This project consist my master thesis research about path following of Otter surface vehicle that is developed by maritime robotic AS. In this study, the trajectory generated using kinodynamic RRT will be tracked using the Model Predictive Control (MPC) with the assistance of the Line of Sight (LOS) algorithm. In autonomous systems, the vehicle's route planning and trajectory tracking systems are controlled by the Guidance, Navigation, and Control (GNC) architecture. The block diagram of this structure is depicted in the figure 1.

|![GNCDiagram](https://github.com/ferhannb/Otter_USV_Path_Following/assets/29739404/a81bd463-34a0-4b07-8ba3-4219bf73ee8d)|
|:--:| 
| **Figure 1** |


In the Guidance block, the vehicle's path planning is performed. In this work, route planning was conducted using the Kinodynamic Rapidly-Exploring Random Tree (RRT) method. Additionally, path tracking was achieved using the Line-of-Sight (LOS) method. The controller structure for path tracking was implemented using Model Predictive Control (MPC). 

MPC structure was implemented using the CasADi optimization tools. CasADi provides a powerful framework for formulating and solving optimization problems, making it suitable for developing the MPC framework in autonomous systems. By leveraging CasADi's capabilities, the MPC controller was designed and optimized to achieve precise and efficient route tracking in the autonomous system.

The flow chart of the kinodynamic RRT algorithm used for path planning is shown in the figure 2.

|<img src=![RRTFlowChart](https://github.com/ferhannb/Otter_USV_Path_Following/assets/29739404/60fde4e9-d44b-49a6-851c-90154b414164) />|
|:--:| 
| **Figure 2** |
