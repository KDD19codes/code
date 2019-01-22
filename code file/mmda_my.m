function W = mmda_my(Sw,Sb,lambda,dim)
% function W = mmda_ssdq(mu,dim,opt)
% Max-Min Distance Analysis, by sequential sdp relaxation
% input: mu, c by D matrix, mean vectors in the whittened LDA subspace, D<=c-1,          
%        dim, subspace dimensionality
% output: W, projection matrix, s.t. W^TW = I
%         mu, dimension reduced mean vectors

tic
m = length(Sw);
D = size(Sw{1},1);
% D = size(A0,1);
% [c,D] = size(mu); % c is classnumber
% if D>c-1
%     error('dim should be less than class number minus 1 after LDA');
% end
% 
% % stacking 
AA=0;BB=0;
for i = 1:m
    AA = AA+Sw{i};
    BB = BB+Sb{i};
end;
A0 = BB-lambda*AA;
% iter = 1;
% for ii = 1:c-1
%     for jj = ii+1: c
%         aa(iter,:) = mu(ii,:) - mu(jj,:);
%         A0 = A0 + aa(iter,:)'*aa(iter,:);
%         iter = iter + 1;
%     end
% end
% m = size(aa,1);
A0 = A0/m;

%%% global sdp relaxation
% note that a penalty with 0.001 weight of average class scatters is
% added to the objective of mmda; this almost does not change the minimum
% pairwise distance but helps to improve the scatter of larger pairwise
% distances.

cvx_begin sdp
    variable X(D,D) symmetric
    variable t(1)
    minimize -t-0.001*trace(A0*X)
    subject to 
    trace(X) == dim;
    X >= 0;
    X <= eye(D);
%     for ii = 1:m
%         trace((aa(ii,:)'*aa(ii,:)) * X) >= t;
%     end
    for ii = 1:m
        trace((Sb{ii}-lambda*Sw{ii}) * X)>= t;
    end
cvx_end
[W,~] = eigs(X,dim);

% val = sum((aa*W).^2,2); % calcuting pairwise distance
% [val, oo] = sort(val);
% mm = val(1);

X0 = W*W';
% if nargin>2
%     if strcmp(opt,'g');mu = mu*W;return;end
% end
% refinement by local sdp relaxation
eta = logspace(-2,-4,5);
for iter = 1:length(eta)
    Y0 = inv(X0 + eye(D));
    y0 = trace(Y0*X0);
    cvx_begin sdp
        variable X(D,D) symmetric
        variable t(1)
        minimize -t-0.001*trace(A0*X)
        subject to
        trace(X) == dim;
        X >= 0;
        X <= eye(D);
    for ii = 1:m
        trace((Sb{ii}-lambda*Sw{ii}) * X) >= t;
    end
        trace(Y0*X)<=(1+eta(iter))*y0;
        trace(Y0*X)>=(1-eta(iter))*y0;
    cvx_end
    [W_tmp,~] = eigs(X,dim);

%     val = sum((aa*W_tmp).^2,2); % calcuting pairwise distance
%     [val, oo] = sort(val);
%     mm_tmp = val(1);

  
%     if mm_tmp <mm(end)%decreasing
%         break;
%     end
    
    % updating
%     mm = [mm mm_tmp];
    W = W_tmp;
    X0 = W*W';
end
% mu = mu*W;
toc