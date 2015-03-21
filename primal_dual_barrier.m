clear
clc

p = 51;
d = 50;
n = 200;

c = zeros(n,p+d+1);

c(:,1:d) = randn(n,d);
c(:,d+1) = randn(n,1);
c(:,d+2:end) = abs(randn(n,p));

A = ones(1,p);
b = 1;

x0 = ones(p,1)./p;
z0 = zeros(d,1);
y0 = -4000;

tol=1.e-6; 
gamma=.6;    
alpha=0.3;    

x=x0;
z=z0;
y=y0;
g=gradienf(x,z,c);

s=g(1:p)-A'*y;
mu = x'*s/p;
iter =0; 

tic
while mu >= tol,
    iter = iter + 1;
    
    mu = gamma*mu;
    
    xr = ones(p,1)./x;
    qqq = [zeros(p,1);-g(p+1:p+d);-(A*x-b);mu*ones(p,1)-diag(x)*s]; %Ax-b can be change to 0

    MM = hessianf(x,z,c);
   
    XX=sparse(1:p,1:p,xr.^2);
    RHS = [MM [-A';sparse(d,1)] [eye(p); sparse(d,p)]; A sparse(1,d+p+1); diag(s) sparse(p,d+1) diag(x) ];
    cond(RHS)
    sol = RHS\qqq;
  
    dx= sol(1:p,1);
    dz= sol(p+1:p+d,1);
    dy= sol(p+d+1,1);
    ds= sol(p+d+2:end,1);
    
    beta = 1;
    ss = 0;
    
    while ss <= 0,
        xx = x + beta*dx;
        yy = y + beta*dy;
        zz = z + beta*dz;
        g=gradienf(xx,zz,c);
        s = g(1:p) - A'*yy;
        beta = alpha*beta;
        ss = min([xx;s]);
    end
    
    x=xx;
    y=yy;
    z=zz;
    mu = x'*s/p;
    if iter == 200,
        break;
    end
end
 
alpha = x;
beta = z;
x = c(:,1:d)';
y = c(:,d+1)';
z = c(:,d+2:end)';
toc
iter
sm = 0;
for i = 1:n,
 sm = sm + (y(1,i)-x(:,i)'*beta)^2/(z(:,i)'*alpha);
end
opt_val = sm





