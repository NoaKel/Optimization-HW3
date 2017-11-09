function [ x_convergence, fconvergence ] = gradint_descent( func, x )
% Gradient Descent with Armijo step size

% Initializing
n=length(x);
sigma=0.25; 
beta=0.5;
epsilon=10^-5;
a0=1;
k=1;
[f,g]=func(x);
max_k=10^5;
fconvergence=zeros(max_k,1);
x_convergence=zeros(max_k,n);

% Gradient Descent
while(norm(g)>=epsilon)
    fconvergence(k,1)=f;
    x_convergence(k,:)=x;
    d= -g;
    %d= d/ norm(d);
    s=armijo_step(a0,x,func,f,g,d,sigma,beta);
    k=k+1;
    x=x+s*d;
    [f,g]=func(x);
end

fconvergence=fconvergence(1:k-1,:);
x_convergence=x_convergence(1:k-1,:);

end





