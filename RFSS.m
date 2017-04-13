function [mae,rmse] = RFSS(R,test_set,d,lambda,mu,iteration)
%R: users x items matrix
%test_set: test data [userid,itemid,rating]
%k: the dimension used in singular value decomposition
%lambda: regularization
%rho: penalty parameter
%mu: learning rate
%Author: Zhi-Lin Zhao
%date: 2016-8-30
%version:1

%the range of rating
rmin = 0;
rmax = 20;
rho = 1;


%[m,n]: #user and #item
[m,n] = size(R);

%initialize parameters
alpha = 0.1*rand(m,d)./sqrt(d);
x = alpha;
u = zeros(m,d);
l = ones(m,1);
beta = 0.1*rand(n,d)./sqrt(d);
%beta = 1/sqrt(d)*ones(n,d);
y = beta;
v = zeros(n,d);
rho_alpha = ones(m,1)*rho;
rho_beta = ones(n,1)*rho;
pr_alpha = zeros(m,1);
dr_alpha = zeros(m,1);
pr_beta = zeros(n,1);
dr_beta = zeros(n,1);

%indicator function
I = R;
I(I > 0) = 1;

mae = zeros(iteration,1);
rmse = zeros(iteration,1);

for t = 1:iteration
    %update alpha
    pr_alpha = sqrt(sum((alpha - x).^2,2));
    rho_alpha = GetRho(rho_alpha,pr_alpha,dr_alpha );
    p = (alpha*beta').*(repmat(l.^2,[1 n])) - (repmat(l,[1 n]) + repmat(rho,[m n])).*R;
    q = alpha - x + u;
    for i = 1:m
        da = sum(repmat((p(i,:).*I(i,:))',[1 d]).*beta)...
            + lambda*sum(repmat(R(i,:)',[1 d]).*repmat(alpha(i,:),[n 1]))...
            + rho_alpha(i)*q(i);
        alpha(i,:) = alpha(i,:) - mu*da;
    end
    %alpha = alpha./repmat(sqrt(sum(alpha.^2,2)),[1 d]);
    
    %update x
    tmp = alpha + u;
    ox = x;
    tmp2 = repmat(sum(tmp.^2,2),[1 d]);
    II = (tmp2 < 1);
    x = tmp.*II + (tmp./sqrt(tmp2))*(~II);
    dr_alpha = sqrt(sum((repmat(rho_alpha,[1 d]).*(ox - x)).^2,2));
    
    %update beta
    pr_beta = sqrt(sum((beta - y).^2,2));
    rho_beta = GetRho(rho_beta,pr_beta,dr_beta );
    p = (alpha*beta').*(repmat(l.^2,[1 n])) - (repmat(l,[1 n]) + repmat(rho,[m n])).*R;
    q = beta - y + v;
    for j = 1:n
        db = sum(repmat(p(:,j).*I(:,j),[1 d]).*alpha)...
            + lambda*sum(repmat(R(:,j),[1 d]).*repmat(beta(j,:),[m 1]))...
            + rho_beta(j)*q(j);
        beta(j,:) = beta(j,:) - mu*db;
    end
    %beta = beta./repmat(sqrt(sum(beta.^2,2)),[1 d]);
    
    %update y
    tmp = beta + v;
    oy = y;
    tmp2 = repmat(sum(tmp.^2,2),[1 d]);
    II = (tmp2 < 1);
    y = tmp.*II + (tmp./sqrt(tmp2))*(~II);
    dr_beta = sqrt(sum((repmat(rho_beta,[1 d]).*(oy - y)).^2,2));
    
    %update l
    tmp = alpha*beta';
    l = sum(tmp.*R,2)./sum((tmp.^2).*I,2);
    %update u
    u = u + alpha - x;
    %update v
    v = v + beta - y;
    
    %predict ratings
    p = (alpha*beta').*repmat(l,[1 n]);
    p(p > rmax) = rmax;
    p(p < rmin) = rmin;

    pr = zeros(size(test_set,1),1);
    for i = 1:size(test_set,1)
        uu = test_set(i,1);
        vv = test_set(i,2);
        pr(i) = p(uu,vv);
    end
    mae(t) = MAE(pr,test_set(:,3));
    rmse(t) = RMSE(pr,test_set(:,3));
end
% ttt = sum(alpha.*alpha,2);
end

