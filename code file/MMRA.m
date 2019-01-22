function [W, obj] = MMRA(X, Y, m)
% X: training data each row is a data;
% Y: label vector

d = size(X,2);
W0 = orth(rand(d,m));
maxiter = 10;
[Sw, Sb] = calcSwSb_onevsone1(X, Y);
k = length(Sw);
AA=0;BB=0;
for i = 1:k
    AA = AA+Sw{i};
    BB = BB+Sb{i};
end;

for iter = 1:maxiter+1
    
for i = 1:k
    lam(i) = trace(W0'*Sb{i}*W0)/trace(W0'*Sw{i}*W0);
end;
lambda = min(lam);
% lambda = trace(W0'*Sb*W0)/trace(W0'*Sw*W0);

W = mmda_my(Sw,Sb,lambda,m);
if norm(W-W0,2)<=10e-4 || iter>maxiter 
disp(['loss is ', num2str(norm(W-W0,2))]);
break;
else
W0 = W;
end
obj(iter) = lambda;
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
end
