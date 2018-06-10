function [Ws,cMats] = example_run()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main function of Cost-Aware classification using the FCD Heterogeneous hypergraph (CAFH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (c) 2018 Liang Zhao
% George Mason University
% lzhao9@gmu.edu

%{
Please cite the following paper in any work that uses this material:

Liang Zhao, Amir Alipour-Fanid, Martin Slawski and Kai Zeng. Prediction-time 
Efficient Classification Using Feature Computational Dependencies. in 
Proceedings of the 24st ACM SIGKDD Conference on Knowledge Discovery 
and Data Mining (KDD 2018), London, United Kingdom, Aug 2018, to appear.
%}
% load data and FCD Heterogeneous hypergraph, data format will be introduced
% in the function CAFH.
load('sampledata.mat');

% try different regularization parameters.
lambdas = 2.^(-20:-10);

% feature weights. 
Ws = []; 

% choose the value of l_p quasi-norm (0<p<1), when the power = 1/2, it means
% l_{1/2} quasi-norm. When the power = 2/3, it means l_{1/2} quasi-norm.
% power can be only 1/2 or 2/3, which ensures analytical solutions.
power = 1/2; 
% power = 2/3;

for lambda=lambdas
    % call the ADMM-based algorithm.
    [B,~] = CAFH(data_tr,D,H,lambda,power);
    numFeatures = size(B,1)/2;
    
    % the transformation is guaranteed by Theorem 4.1 in the paper.
    W = (B'*[eye(numFeatures);-eye(numFeatures)])'; 
    Ws = [Ws,W];
end
% calculate the confusion matrix.
cMats = get_predict_result(data_te,Ws,lambdas);

% calculate the accuracy and runtime.
[accs,runtimes] = postprocess(Ws,cMats,D,H);
end

function [results,runtimes] = postprocess(ws,cmats,D,H)
ws = ws(2:end,:);
runtimes = [];
for i=1:size(ws,2)
    ws_basic = H*ws(:,i);
    % estimate the feature generation runtime.
    runtimes = [runtimes,sum(D(ws_basic(:,1)~=0,:))];
end
% get the precision, recall, and F-measure.
results = zeros(0,3);
for c=cmats
    c = c{1,1};
    p = c.tp/(c.tp+c.fp);
    r = c.tp/(c.tp+c.fn);
    f = 2*p*r/(p+r);
    results(end+1,:)=[p,r,f];
end
end