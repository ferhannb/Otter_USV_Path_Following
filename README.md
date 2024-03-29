# Path Following of Otter USV using model predictive control (MPC)

This project consist my master thesis research about path following of Otter surface vehicle that is developed by maritime robotic AS. In this study, the trajectory generated using kinodynamic RRT will be tracked using the Model Predictive Control (MPC) In autonomous systems, the vehicle's route planning and trajectory tracking systems are controlled by the Guidance, Navigation, and Control (GNC) architecture. The block diagram of this structure is depicted in the figure 1.

|![GNCDiagram](https://github.com/ferhannb/Otter_USV_Path_Following/assets/29739404/a81bd463-34a0-4b07-8ba3-4219bf73ee8d )|
|:--:| 
| **Figure 1** |


In the Guidance block, the vehicle's path planning is performed. In this work, path planning was conducted using the Kinodynamic Rapidly-Exploring Random Tree (RRT) method. Additionally, path following was achieved using MPC. 

MPC structure was implemented using the CasADi optimization tools. CasADi provides a powerful framework for formulating and solving optimization problems, making it suitable for developing the MPC framework in autonomous systems. By leveraging CasADi's capabilities, the MPC controller was designed and optimized to achieve precise and efficient route tracking in the autonomous system.

The flow chart of the kinodynamic RRT algorithm used for path planning is shown in the figure 2.

|![RRTFlowChart](https://github.com/ferhannb/Otter_USV_Path_Following/assets/29739404/60fde4e9-d44b-49a6-851c-90154b414164) |
|:--:| 
| **Figure 2** |

In Path following part, Referances values for path following comes from the path is created by path generation algorithm ($x_{ref},y_{ref} \quad and  \quad \psi_{ref}$) and assigned surge velocity and yaw angular rate (u,r)  give to cost function of MPC as target values so generated path is followed. The reference values to be given to the MPC algorithm are the position values of the points determined on the path and the reference speeds at those points.

Casadi sutructure of MPC algorithm is shown figure 3. By establishing this block structure in the casadi, which is shown in the picture, the optimizaton algorithm of the MPC is established. The solution of this structure is solved according to the instant states of the system at the determined time intervals and the optimum input signals are determined.


|![GNCDiagram](https://github.com/ferhannb/Otter_USV_Path_Following/blob/master/CasadiDiagram.png)|
|:--:| 
| **Figure 3** |

## Simulation scanerio
In order to test the performance of this algorithm, a path was generated with kinodynamic RRT in variousenvironments and this path was followed by MPC.


