function [W,D] = SADL2(X,Y,D0,i_groups,sigma2,lambda3,verbal)
%{ 
Copyright (c) 2023 Liang Zhao
Emory University
liang.zhao@emory.edu

The main function of the parameter optimization algorithm for the model
SADL-2
the detailed algorithmatic introduction is detailed in the paper:

Zhao, Liang, Olga Gkountouna, and Dieter Pfoser. "Spatial auto-regressive dependency 
interpretable learning based on spatial topological constraints." ACM Transactions on 
Spatial Algorithms and Systems (TSAS) 5, no. 3 (2019): 1-28.

Notations:
    - p: original number of features; p+1: the number of features including dummy feature
    - q: number of feature components
    - n: number of observations

Input:
    - X:     (n*t)*k         training data input matrix: (number of locations * number of times) * number of features.
    - Y:        (n*t)*1     training data output vector: (number of locations * number of times) *1
    - D0:        n*n        incidence matrix among spatial locations: (number of locations * number of locaitons
    - i_groups: 1*g         location grouping: g is number of location groups
    - sigma2:   scalar      any real number (default: 1)
    - lambda3:    scalar    any real number (default: 0)
Output: 
    - W:        1*k        learned feature weights: 1* number of features
    - D:        n*n        learned spatial location dependencies: (number of locations * number of locaitons
%}
numTimes = size(Y,2);
numLocs = size(Y,1);
numGroups = size(i_groups,2);
X = [ones(numTimes*numLocs,1),X];
numFeatures = size(X,2);
if numTimes*numLocs ~= size(X,1)
    error('numTimes and numLocs are not equal to the size of X');
end
% initialization
% sigma2 = 1;
ERR = 0.01;
eps = 0.00001;
lambda1 = 1;
lambda2 = 0;
lambda3 = 0.00001;
I1 = speye(numLocs,numLocs);
I2 = speye(numFeatures,numFeatures);
W = sparse(1,numFeatures);
D = sparse(D0);
E = sparse(diag(D*ones(numLocs,1))-D+ones(numLocs,1)*ones(numLocs,1)'-eps*eye(numLocs));
V = sparse(1,numFeatures);
U = E;
Lam1 = sparse(numLocs,numLocs);
Lam2 = sparse(numLocs,numLocs);
Lam3 = sparse(1,numFeatures);
IterMax = 1000;
XX = X'*X;
YY = Y*Y';
rho = 1;
tmp_val = zeros(size(E));
for i=1:IterMax
    if i == 20
        aa = 1;
    end 
    W_old = W;E_old=E;D_old=D;sigma2_old=sigma2;U_old=U;V_old=V;
    YEX = sparse(1,numFeatures);
    for j=1:numTimes
        YEX = YEX + Y(:,j)'*(I1-D)'*X((j-1)*numLocs+1:j*numLocs,:);
    end
    W = (YEX/sigma2/numTimes+rho*(V-Lam3))/(XX/sigma2/numTimes+rho*I2); %%
    XWY = reshape(X*W',numLocs,[])*Y';
       D = update_D(D,YY,rho,U,Lam1,lambda2,lambda3,sigma2*numTimes,XWY,D0,I1,E,Lam2,eps,i_groups);
    [E,tmp_vals] = update_E(D,Lam2,i_groups);
   WLam = W+Lam3;
    V = WLam;
    V(:,2:end)=sign(WLam(:,2:end)).*max(abs(WLam(:,2:end))-lambda1/rho,0); %%
    [P,Q] = eig(full((I1-D-Lam1)*rho));
    U = rho/2*P*(Q+sqrt(Q*Q+4*rho*I1))*P';
    Lam1 = Lam1 + (U-I1+D);
    for j=1:numGroups
        cur_g = i_groups{1,j};
        tmp_val(cur_g,cur_g)=tmp_vals{1,j};
    end
    Lam2 = Lam2 + tmp_val-E;
    Lam3 = Lam3 + W-V;
    p = norm(U-I1+D,'fro')+norm(W-V,'fro')+norm(tmp_val-E,'fro');
    d = rho*(norm(U-U_old-D+D_old,'fro')+norm(V-V_old,'fro')+norm(E-E_old,'fro'));
    if p < ERR && d < ERR
        break;
    end
    if verbal
        fprintf('p:%e\t p1:%e:\t p2:%e\t p3:%e\t d:%e\t d1:%e\t d2:%e\t d3:%e\t d4:%e\t rho:%f\n',...
            p,norm(U-I1+D,'fro'),norm(tmp_val-E,'fro'),norm(W-V,'fro'),d,rho*(norm(U-U_old,'fro')),rho*(norm(D-D_old,'fro')),...
            rho*(norm(V-V_old,'fro')),norm(E-E_old,'fro'),rho);
    end
end
end
function [E,tmp_vals]=update_E(D,Lam2,i_groups)
numGroups = size(i_groups,2);
E = zeros(size(D));
tmp_vals = {};
for i=1:numGroups
    cur_g = i_groups{1,i};
    cur_D = D(cur_g,cur_g);
    cur_Lam2 = Lam2(cur_g,cur_g);
    [cur_E,cur_tmp_val] = update_sub_E(cur_D,cur_Lam2);
    E(cur_g,cur_g)=cur_E;
    tmp_vals = {tmp_vals{:},cur_tmp_val};
end
E = sparse(E);
end
function [E,tmp_val]=update_sub_E(D,Lam2)
numLocs = size(D,1);
tmp_val = diag(D*ones(numLocs,1))-D+ones(numLocs,1)*ones(numLocs,1)'-eps*eye(numLocs);
E = tmp_val+Lam2;
E = sym_proj(E);
end

function new_D=update_D(D,YY,rho,U,Lamb1,lambda2,lambda3,sigma2T,XWY,D0,I1,E,Lamb2,eps,i_groups)
numGroups = size(i_groups,2);
new_D = zeros(size(D));
for i=1:numGroups
    cur_g = i_groups{1,i};
    cur_D = D(cur_g,cur_g);
    cur_YY = YY(cur_g,cur_g);
    cur_U = U(cur_g,cur_g);
    cur_Lamb1 = Lamb1(cur_g,cur_g);
    cur_XWY = XWY(cur_g,cur_g);
    cur_D0 = D0(cur_g,cur_g);
    cur_I1 = I1(cur_g,cur_g);
    cur_E = E(cur_g,cur_g);
    cur_Lamb2 = Lamb2(cur_g,cur_g);
    cur_D = update_sub_D(cur_D,cur_YY,rho,cur_U,cur_Lamb1,lambda2,lambda3,...
                            sigma2T,cur_XWY,cur_D0,cur_I1,cur_E,cur_Lamb2,eps);
    new_D(cur_g,cur_g)=cur_D;
end
new_D = sparse(new_D);
end

function D = update_sub_D(D,YY,rho,U,Lamb1,lambda2,lambda3,sigma2T,XWY,D0,I1,E,Lamb2,eps)
IterMax = 200;
eta = 0.001;
D_old = sparse(zeros(size(D)));
numLocs = size(D,1);
Mat1 = ones(numLocs,1);
OneMat = Mat1*Mat1';
M = OneMat-eps*eye(numLocs)-E+Lamb2;
MMM = rho*diag(diag(M))*OneMat;
for i=1:IterMax
    dD1 = rho*D*OneMat-rho*(diag(diag(D))*OneMat+diag(D*Mat1))+MMM+rho*D-rho*M;
    dD2 = D*YY/sigma2T-YY/sigma2T+XWY/sigma2T+rho*D+rho*(U-I1+Lamb1)+2*lambda2*D;
    dD = dD1+dD2;
    D = D-eta*dD;
    D = (D+D')/2;
    D = max(D,0);
    D = D.*D0;
    D = sign(D).*max(abs(D)-lambda3/rho,0);
    e = norm(D-D_old,'fro');
    if e<0.001
        break
    end
    D_old = D;
end
end

function E = sym_proj(E)
    numLocs = size(E,2);
    [P,Q] = eig(full(E));
    tmp_E = zeros(size(E));
    for j=1:numLocs
        if Q(j,j)<0
            continue
        end
        tmp_E = tmp_E + P(:,j)*Q(j,j)*P(:,j)';
    end
    E = tmp_E;
end

