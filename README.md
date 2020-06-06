# First Order Reliability methods

Here is possible to find the optimization codes I wrote for computing the reliability index and the failure probability for any limite state function. 

There 4 files: 
*  **[FORM_SQP.m](https://github.com/iagolemos1/First-Order-Reliability-methods/blob/master/FORM_SQP.m)**: The most general code, it solves, prior, any limit state function, implicit or explicit, for normal variables. It uses the Sequential Quadratic Programming solver from MATLAB optmization toolbox. In addition, it's quite usual using the [Rosenblatt Transformation](https://github.com/iagolemos1/Rosenblatt-Transformation/blob/master/Rosenblatt_transform.m) in order to transform a non-normal variable to normal and use the Sequential Quadratic Programming solver to find the reliability index and the failure probability. 

* **[FORM_RF_standard](https://github.com/iagolemos1/First-Order-Reliability-methods/blob/master/FORM_RF_standard.m)**: Solver that uses the Rackwitz-Fiessler algorithm to find the reliability index and the failure probability for explicit limit state function, i.e, function in the standard space. Also known as **FORM Method 2**. 

* **[FORM_RF.m](https://github.com/iagolemos1/First-Order-Reliability-methods/blob/master/FORM_RF.m)**: Solver that uses the Rackwitz-Fiessler algorithm with a Newton-Raphson approach for implicit limit state functions. In other words, it's not necessary to transform the variables into the normal standard space. Also known as **FORM Method 2**. For runing this solver, it's necessary to have the [Rosenblatt Transformation](https://github.com/iagolemos1/Rosenblatt-Transformation/blob/master/Rosenblatt_transform.m) file in your computer. 

* **[FORM_R.m](https://github.com/iagolemos1/First-Order-Reliability-methods-MATLAB/blob/master/FORM_R.m)**: Solver that uses the first Rakcwitz algorithm with improvements by Ayyub-Haldar to find the solution for implicit limit state functions. Also known as **FORM Method 1**. It's necessary to have the [Rosenblatt Transformation](https://github.com/iagolemos1/Rosenblatt-Transformation/blob/master/Rosenblatt_transform.m) file in your computer. 

**Observations:** 
* As the Newton-Raphson scheme may fail to converge, also the Rackwitz-Fiessler algorithm may fail too;

