# Koopman_operator_for_switch_system
This repository implements Koopman representation for several switched systems.

## Detail information
`num_ex`: numerical example of a switched system
1. `num_ex2.m`: generate data of the numerical system
2. `main.m`: learn Koopman operator and estimate the prediction
---
`SLIP_ex`: example of Spring Loaded Inverted Pendulum (SLIP) model, which switch between stance phase and flight phase
`SLIP_ex_trot`: example of dimensionless SLIP model with resonable parameters of trot gait, which switch between stance phase and flight phase
1. `gen_SLIP_model.m`: generate equations of motion (EOM) of the SLIP model
2. `gen_SLIP_fp.m`: generate fixed points of the SLIP model, i.e. gait library
3. `gen_SLIP_fp_data.m`: generate data from the selected fixed points and save the data in a mat file
4. `gen_nominal.m`: generate norminal basis
5. `main.m`: learn the Koopman operator and estimate the prediction. `main2.m` presents a comparison without switching Koopman operator.
---
`LIP`: example of Linear Inverted Pendulum model, which switch between single stance double stance phase
1. `gen_LIP_data.m`: generate data from the selected gait parameters and save the data in a mat file
2. `gen_nominal.m`: generate norminal basis
3. `main.m`: learn the Koopman operator and estimate the prediction
---
`Notes`: notes for derivation 

## Reference
[1] Bruder, Daniel, C. David Remy, and Ram Vasudevan. "Nonlinear system identification of soft robot dynamics using koopman operator theory." In 2019 International Conference on Robotics and Automation (ICRA), pp. 6244-6250. IEEE, 2019.
[2] Bakker, Craig, Arnab Bhattacharya, Samrat Chatterjee, Casey J. Perkins, and Matthew R. Oster. "Learning koopman representations for hybrid systems." arXiv preprint arXiv:2006.12427 (2020).
[3] Kajita, Shuuji, Fumio Kanehiro, Kenji Kaneko, Kazuhito Yokoi, and Hirohisa Hirukawa. "The 3D linear inverted pendulum mode: A simple modeling for a biped walking pattern generation." In Proceedings 2001 IEEE/RSJ International Conference on Intelligent Robots and Systems. Expanding the Societal Role of Robotics in the the Next Millennium (Cat. No. 01CH37180), vol. 1, pp. 239-246. IEEE, 2001.