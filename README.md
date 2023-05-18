# State-feedback-control-for-AM-temperature-model
Algorithms solving the Linear Quadratic Tracking (LQT) problem. We approach the problem both through model-based and data-driven techniques.

We study the problem of controlling the behaviour of a linear state-space system. The optimisation problem is formulated as the minimization of a quadratic performance index. 

We first adress the problem through a classic Control Theory perspective, assuming that an underlying model of the process is known and the time horizon for the optimisation is pre-determined and finite. The theory associated with this approach can be found in the book "Optimal Control; Linear Quadratic Methods" by Anderson and Moore.

Next, we explore the case of the infinite horizon LQT problem, and adress it through Reinforcement Learning techniques, specifically Q-learning. We provide the model-based version, where we assume the dynamics of the system are known.

Finally, we study the problem of learning optimal policies strictly from data, with no prior knowledge of the system dynamics. That is done with the help of the least squares algorithm.

All scripts are written in R. The paper associated with the theory behind the code will soon be on ArXiv.
