clear
clc

randn('seed',100) %300

%{
load quadoverlindatlccp.mat
p = 51;
d = 50;
n = 200;

x = c(:,1:d)';
y = c(:,d+1)';
z = c(:,d+2:end)';
%}

p = 20;
d = 20;
n = 50;

x = randn(d,n);
y = randn(1,n);
z = abs(randn(p,n));

%CVX
tic
cvx_begin 
        variable alp(p) nonnegative
        variables bet(d) gam(n,1);
        
        sum(alp) == 1;
        
        for i = 1:n,
            [gam(i,1) y(1,i)-x(:,i)'*bet;
             y(1,i)-x(:,i)'*bet z(:,i)'*alp] == semidefinite(2);
        end
        minimize ( sum(gam) )
            
cvx_end
toc 

cvx_optval

opt_val_store = [];
flag_count = [];

%Initialize 
alpha = zeros(p,1);
beta = zeros(d,1);
gamma = zeros(n,1);

rho = 1;
tol = 1e-4;
mu = 1;

%Multiplier for 2x2 constraints
X = ones(2,2,n);

%Multiplier for nonnegative alpha
Y = ones(p+1,1);

%Slack for 2x2 constraints
S = zeros(2,2,n);
for i = 1:n,
    S(:,:,i) = eye(2);
end

%Slack for nonnegative alpha
T = ones(p+1,1);

tic
for iter = 1:1000,
    
    %iter
    
    %update SDP dual variable vector [alpha; beta; gamma]
    b = -Y(1:p,1) + Y(p+1,1)*ones(p,1) + (rho)*T(1:p,1) - rho*(T(p+1,1)-1)*ones(p,1);
    A = rho*eye(p) + rho*ones(p);
    for i = 1:n,
        b = b + (-X(2,2,i)+rho*S(2,2,i))*z(:,i);
        A = A + rho*z(:,i)*z(:,i)';
    end

    alpha = A\b;
    
    b = zeros(d,1);
    A = zeros(d,d);

    for i = 1:n,
        b = b + (2* X(1,2,i) + 2*rho*(y(1,i)-S(1,2,i)) )*x(:,i);
        A = A + 2*rho*x(:,i)*x(:,i)';
    end

    beta = A\b;

    for i = 1:n,
        gamma(i,1) =  -(1/rho)*(1  + X(1,1,i) - rho*S(1,1,i));
    end
       
    %update slack variable on 2x2 sdp constraints
    for i = 1:n,
        temp = -X(:,:,i) + rho*[-gamma(i,1), x(:,i)'*beta - y(1,i); 
                               x(:,i)'*beta - y(1,i), -z(:,i)'*alpha];
        [V,D] = eigs((1/rho^2)*temp^2 + 4*(mu/rho)*eye(2));
        for j = 1:2,
           if D(j,j) < 0,
                D(j,j) = 0;
           end
        end
        S(:,:,i) = (-temp/rho+V*sqrt(D)*V')/2;  
    end
    
    %update slack variable on nonnegative cone for alpha
    temp = -diag(Y) + rho*diag([-alpha; sum(alpha) - 1]);
    [V,D] = eig((1/rho^2)*temp^2 + 4*(mu/rho)*eye(p+1));
    for j = 1:p+1,
        if D(j,j) < 0,
             D(j,j) = 0;
        end
    end
    T = diag((-temp/rho+V*sqrt(D)*V')/2);
   
    %update multiplier for 2x2 constraints
    for i = 1:n,
        temp = zeros(2,2);
        for j = 1:p,
            temp = temp + alpha(j,1)*[0 0; 0 -z(j,i)];
        end

        for j = 1:d,
            temp = temp + beta(j,1)*[0 x(j,i); x(j,i) 0];
        end

        temp = temp + gamma(i,1)*[-1 0; 0 0];

        temp = temp + S(:,:,i);
        temp = temp - [0 y(1,i); y(1,i) 0];
          
        X(:,:,i) = X(:,:,i) - rho*temp;
    end
    
    %update multiplier for nonnegative cone constraint
    temp = zeros(p+1,1);
    temp(1:p,1) = -alpha + T(1:p,1);
    temp(p+1,1) = sum(alpha) + T(p+1,1) - 1;
    Y = Y - rho*temp;
    
    sm = 0;
    for i = 1:n,
        sm = sm + (y(1,i)-x(:,i)'*beta)^2/(z(:,i)'*alpha);
    end
    
    obj = sum(gamma);
    obj2 = sm;
    
    opt_val = obj;
    opt_val_store = [opt_val_store norm(opt_val-cvx_optval)];
    
    mu = 0.8*mu;
    
    %Stopping 
    flag = 0;
    %Check for 2x2 constraints
    for i = 1:n,
        temp = zeros(2,2);
        for j = 1:p,
            temp = temp + alpha(j,1)*[0 0; 0 -z(j,i)];
        end

        for j = 1:d,
            temp = temp + beta(j,1)*[0 x(j,i); x(j,i) 0];
        end

        temp = temp + gamma(i,1)*[-1 0; 0 0];

        temp = temp + S(:,:,i);
        temp = temp - [0 y(1,i); y(1,i) 0];
        if norm(temp,'fro') < 1e-5
            flag = flag + 1;
        end
    end
    
    %check nonnegative cone constraint
    temp = zeros(p+1,1);
    temp(1:p,1) = -alpha + T(1:p,1);
    temp(p+1,1) = sum(alpha) + T(p+1,1) - 1;
    if norm(temp) < tol
        flag = flag + 1;
    end
    
    flag_count = [flag_count n+1-flag];
    
    if flag == n+1,
        break;
    end
    
end
toc

iter
obj = sum(gamma)

figure();
hold on;
plot(1:iter, opt_val_store(1:iter));
plot(1:iter, flag_count, 'r');
xlabel('iter');
legend('||f(x^*) - f(x_k)||', '# unsatisfied');
title('Interior point ADMM');



