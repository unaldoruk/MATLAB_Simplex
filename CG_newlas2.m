function [z,x] = CG_newlas2(H,c,A,b)
tic
%nmax=1000;
%This code uses conjugate gradient method in order to solve quadratic
%programming problems. There are some assumptions behind this code.
% 1-)Although CG takes only equality contraints, by thinking the result
% will also try to satisfy equality, this method is used.
%2-) There is no necessary condtion for A to be SPD, in order to use it in
%the code A'*A is used.
%3-) In order to be sure about non-negative x values, x'*x is used on the
%value.
x=A\b;% Initial point;
tol=1.0e-14;%tolerance limits
%According to resources in reference, the matrices which are not SPD, can
%be used as A'*A.
residual=+A'*b-(c+H*x)-A'*A*x;
resnew=residual;
qnew=residual'*residual;
q0=qnew;
numiter=0;
znew=0.5*x'*H*x+c'*x;
zold=znew*2;
while numiter<length(c) && qnew>tol*tol*q0 && znew<zold+tol
    %Different stopping criteria are used to make sure contraints are satisfied.
    value=(H+A'*A+x'*x)*resnew;% this comes from orthogonolity assumption, and used for calculation of alpha.
    alpha=qnew/(resnew'*value);
    x=x+alpha*resnew;
    x0=x-alpha*resnew;
    qold=residual'*residual;
    if mod(numiter,50)==0% in order to eliminate floating round-off effect
        residual=residual+(+A'*b-(c+H*x)-A'*A*x);
    end
    if mod(numiter,50)~=0
        residual=residual-alpha*value;
    end
    qnew=residual'*residual;
    beta=-qnew/qold;
    resnew=residual+beta*resnew;
    numiter=numiter+1;
    zold=znew;
    znew=0.5*x'*H*x+c'*x;
end
x=(x0+x)/2;%this code generally gives the answer in two steps, and the z values at these steps give one greater than optimum and one less than optimum, therefore the average of these approximately gives the optimum. 
z=0.5*x'*H*x+c'*x;
t=toc;
display(t);% in order to see computation time
a=b-A*x;
display(a);
end




