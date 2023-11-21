function [Y_pred,rmse] = predict_Y(W,D,X_te,Y_te,rho)
numObs = size(X_te,1);
X_te = [ones(numObs,1),X_te];
numLocs = size(D,1);
I1 = eye(numLocs);
Y_pred = reshape((I1-rho*D)\reshape(X_te*W',numLocs,[]),[],1);
Y_te = reshape(Y_te,[],1);
rmse = norm(Y_pred-Y_te,'fro')/sqrt(size(Y_pred,1));
end