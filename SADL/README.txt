SADL package Version 1.0 (Nov, 2023) (Matlab code)

------------------------------------

Copyright (c) 2023 Liang Zhao
Emory University
liang.zhao@gmu.edu

Please cite the following paper in any work that uses this material:

Zhao, Liang, Olga Gkountouna, and Dieter Pfoser. "Spatial auto-regressive dependency 
interpretable learning based on spatial topological constraints." ACM Transactions on 
Spatial Algorithms and Systems (TSAS) 5, no. 3 (2019): 1-28.

@article{zhao2019spatial,
  title={Spatial auto-regressive dependency interpretable learning based on spatial topological constraints},
  author={Zhao, Liang and Gkountouna, Olga and Pfoser, Dieter},
  journal={ACM Transactions on Spatial Algorithms and Systems (TSAS)},
  volume={5},
  number={3},
  pages={1--28},
  year={2019},
  publisher={ACM New York, NY, USA}
}

Together with this README file are the codes and data.

---------------
ENVIRONMENT: 
---------------
Matlab 2020 and later (previous version might also work)


---------------
INSTALLATION: 
---------------
1. unzip the downloaded package
2. open Matlab, then set the package folder as the current path of Matlab.

---------------
EXAMPLE RUN:
---------------
1. Reset the current path to the folder of our code.
2. Run the following command lines:
>> load('SADL.mat');
>> [wat_W_tmp,wat_D_tmp] = SADL1(wat_X_tr,wat_Y_tr,wat_d_mat,1,0,true);
>> [wat_Y_pred_tmp,wat_rmse_tmp] = predict_Y(wat_W_tmp,wat_D_tmp,wat_X_te,wat_Y_te,1);


---------------
DATA FORMAT: 
---------------
Please refer to the comments in files 'SADL1.m' and ''SADL2.m'' for the format of the data in SADL.mat



If any problem, please contact via liang.zhao@emory.edu.