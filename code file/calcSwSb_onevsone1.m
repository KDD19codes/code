function [SW, SB] = calcSwSb_onevsone1(X, Y)
% X: training data each row is a data;
% Y: label vector

classid = unique(Y);
c = length(classid);
for i=1:c
    X1{i} = X(Y==classid(i),:);
end;

k = 1;
for i=1:c-1
    Xi = X1{i};
    ni = size(Xi,1);
    Yi = [ones(ni,1),zeros(ni,1)];
    for j=i+1:c
        Xj = X1{j};
        Xij = [Xi;Xj];
        nj = size(Xj,1);
        Yj = [zeros(nj,1),ones(nj,1)];
        Yij = [Yi;Yj];
        [Sb, Sw] = calculate_SbSw(Xij,Yij);
        %Sw{k} = Sw{k} + 1*Sb{k} + 0.*eye(size(X,2));
        SW{k} = Sw;
        SB{k} = (ni+nj)*Sb;
        k = k+1;
    end;
end;

% AA=0;BB=0;
% for i = 1:length(Sw)
%     AA = AA+Sw{i};
%     BB = BB+Sb{i};
% end;
% [Sb11, Sw11, L_b, L_w] = calculate_L(X,Y);
% Y1 = TransformL(Y,c);
% [Sb1, Sw1] = calculate_SbSw(X,Y1);
% sw=(c-1)*Sw1./AA;
% sb=size(X,1)*Sb1./BB;
% sw1=Sw1./Sw11;
% sb1=Sb1./Sb11;
% 1;

function [Sb, Sw] = calculate_SbSw(X,Y)
% X: training data each row is a data;
% Y: label matrix n*c


n = size(X,1);
H = eye(n) - 1/n*ones(n);

Sb = X'*H*Y*inv(Y'*Y)*Y'*H*X;
Sw = X'*H*X - Sb;
