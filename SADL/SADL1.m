function [W,D,Y_scale] = SADL1(X,Y,D0,sigma2,lambda3,verbal)
%{ 
Copyright (c) 2023 Liang Zhao
Emory University
liang.zhao@emory.edu

The main function of the parameter optimization algorithm for the model
SADL-1
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
    - sigma2:   scalar      any real number (default: 1)
    - lambda3:    scalar    any real number (default: 0)
Output: 
    - W:        1*k        learned feature weights: 1* number of features
    - D:        n*n        learned spatial location dependencies: (number of locations * number of locaitons
%}
Y_scale = max(max(Y));
Y = Y/Y_scale;
numTimes = size(Y,2);
numLocs = size(Y,1);
X = [ones(numTimes*numLocs,1),X];
numFeatures = size(X,2);
if numTimes*numLocs ~= size(X,1)
    error('numTimes and numLocs are not equal to the size of X');
end
% initialization
% sigma2 = 1;
ERR = 0.001;
lambda1 = 1;
lambda2 = 0;
I1 = speye(numLocs,numLocs);
I2 = speye(numFeatures,numFeatures);
W = sparse(1,numFeatures);
E = sparse(I1-D0);
D = sparse(D0);
V = sparse(1,numFeatures);
U = E;
Lam1 = sparse(numLocs,numLocs);
Lam2 = sparse(numLocs,numLocs);
Lam3 = sparse(1,numFeatures);
IterMax = 1000;
XX = X'*X;
YY = Y*Y';
rho = 1;
for i=1:IterMax
    W_old = W;E_old=E;D_old=D;sigma2_old=sigma2;U_old=U;V_old=V;
    YEX = sparse(1,numFeatures);
    for j=1:numTimes
        YEX = YEX + Y(:,j)'*(I1-D)'*X((j-1)*numLocs+1:j*numLocs,:);
    end
    W = (YEX/sigma2+rho*(V-Lam3))/(XX/sigma2+rho*I2); %%
    XWY = reshape(X*W',numLocs,[])*Y';
    D = update_D(D,YY,Y,X,W,rho,U,Lam1,lambda2,lambda3,sigma2,XWY,D0,I1);
    EYXW = reshape(E*Y,[],1)-X*W';

    WLam = W+Lam3;
    V = WLam;
    V(:,2:end)=sign(WLam(:,2:end)).*max(abs(WLam(:,2:end))-lambda1/rho,0); %%

    [P,Q] = eig(full(numTimes*(I1-D-Lam1)/rho));
    U = rho/2/numTimes*P*(Q+sqrt(Q*Q+4*rho*I1/numTimes))*P';

    Lam1 = Lam1 + (U-I1+D);

    Lam3 = Lam3 + W-V;

    p = norm(U-I1+D,'fro')+norm(W-V,'fro');
    d = rho*(norm(U-U_old-D+D_old,'fro')+norm(V-V_old,'fro'));
    if(p>10*d)
        rho = 2*rho;
    else
        if(10*p<d)
            rho = rho/2;
        end
    end
    if p < ERR && d < ERR
        break;
    end
    if verbal
        fprintf('p:%e\t p1:%e:\t p2:%e\t d:%e\t d1:%e\t d2:%e\t d3:%e\t d4:%e\t rho:%f\n',...
            p,norm(U-I1+D,'fro'),norm(W-V,'fro'),d,rho*(norm(U-U_old,'fro')),rho*(norm(D-D_old,'fro')),...
            rho*(norm(V-V_old,'fro')),rho*(norm(W-W_old,'fro')),rho);
    end
end

end
function D=update_D(D,YY,Y,X,W,rho,U,Lam1,lambda2,lambda3,sigma2,XWY,D0,I1)
IterMax = 100;
eta = 0.01;
D_old = sparse(zeros(size(D)));
numLocs = size(D,1);
for i=1:IterMax
    dD = D*YY/sigma2-YY/sigma2+XWY/sigma2+rho*D+rho*(U-I1+Lam1)+2*lambda2*D;
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
