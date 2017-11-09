function [ xconvergence, fconvergence ] = newton_method( func, x )
%% Newton Method with Armijo step size

% Initializing
n=length(x);
sigma=0.25;
beta=0.5;
epsilon=10^-5;
a0=1;
k=1;
[f,g,H]=func(x);
k_max=2000;
fconvergence=zeros(k_max,1);
xconvergence=zeros(k_max,n);

% Newton Method
while(norm(g)>=epsilon && k<k_max)
    d=solve_for_d(g,H);
    %d= d/ norm(d);
    s=armijo_step(a0,x,func,f,g,d,sigma,beta);
    fconvergence(k,1)=f;
    xconvergence(k,:)=x;
    k=k+1;
    x=x+s*d;
    [f,g,H]=func(x);
end

fconvergence=fconvergence(1:k-1,:);
xconvergence=xconvergence(1:k-1,:);

%% 
function [ d ] = solve_for_d( g,H )
    n=length(g);
    % cholesky factorization
    [L,diagonal]=mcholmz(H); % given to us
    %D=diag(diagonal);
    % step 1: foward substitution
    y=zeros(n,1);
    y(1,1)=-g(1)/L(1,1);
    for i=2:n
        y(i)=((-g(i)-L(i,1:i-1)*y(1:i-1))/L(i,i));
    end
    % step 2
    diagonal=1./diagonal;
    m=diagonal.*y;
    % step 3: back substitution
    L_t=L';
    d=zeros(n,1);
    d(n,1)=m(n)/L_t(n,n);
    for i=1:n-1 
        d(n-i)=(m(n-i)-L_t(n-i,n-i+1:n)*d(n-i+1:n))/L_t(n-i,n-i);
    end
end

end

