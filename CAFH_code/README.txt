CAFH package Version 1.0 (June, 2018) (Matlab code)

------------------------------------

Copyright (c) 2018 Liang Zhao
George Mason University
lzhao9@gmu.edu

Please cite the following paper in any work that uses this material:

Liang Zhao, Amir Alipour-Fanid, Martin Slawski and Kai Zeng. Prediction-time Efficient Classification Using Feature Computational Dependencies. in Proceedings of the 24st ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD 2018), London, United Kingdom, Aug 2018, to appear.

@inproceedings{zhao2018prediction,
  title={Prediction-time Efficient Classification Using Feature Computational Dependencies},
  author={Zhao, Liang and Alipour-Fanid, Amir and Slawski, Martin and Zeng, Kai},
  booktitle={Proceedings of the 24nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining},
  year={2018},
  organization={ACM}
}

Together with this README file are the codes and data.

---------------
ENVIRONMENT: 
---------------
Matlab 2011-2016 (previous version might also work)


---------------
INSTALLATION: 
---------------
1. Before running the codes, please first download Matlab package "L1General" from the link: http://www.cs.ubc.ca/~schmidtm/Software/L1General.zip
2. unzip the downloaded package
3. open Matlab, then set the package folder as the current path of Matlab.
4. add the package path using the following command
>> addpath(genpath(pwd));

---------------
EXAMPLE RUN:
---------------
1. Reset the current path to the folder of our code.
2. Run the following command lines:
>> load('sampledata.mat')
>> example_run


---------------
DATA FORMAT: 
---------------
Please refer to the comments in file 'CAFH.m' for the format of the data in sampledata.mat



If any problem, please contact Dr. Liang Zhao via lzhao9@gmu.edu.