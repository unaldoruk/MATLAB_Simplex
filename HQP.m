function [z,x] = HQP(H,c,A,b)
%-------------------------------------------------------------------------%
%The following function performs Hildreth's algorithm to solve a         
%quadratic program. It iteratively uses Lagrange multipliers to arrive at
%approximately the optimal solution. 
%   
% min z = 1/2x'Hx + cx
% s.t. Ax <= b
%      x >= 0
%
% 1. Compute the general optimal solution
% 2. Substitute the general optimal solution into constraints and determine
% if they are violated. If not violated: break
% 3. Iterate through the calculation of lambda (L) until 
%                       L_i = L_i+1 +/- tolerance
% 4. Compute the x vector given the lambda values
%
% Variations of the Lambda Calculation
% -------------------------------------
% Summation Version of the Lambda calculation
%         sum1 = sum(D(jj,1:(jj-1))'.*L(1:(jj-1),1));
%         sum2 = sum(D(jj,(jj+1):end)'.*Lold((jj+1):end,1));
%         sum3 = sum1 + sum2 + K(jj,1);
%         result = sum3/D(jj,jj);
%         wii_jj = -result;
%        
% Vectorized Version of the Lambda calculation - used to improve speed
%         step1 = D(jj,:)*L(:,1);
%         step2 = step1 - D(jj,jj)*L(jj,1);
%         step3 = step2 + K(jj,1);
%         step4 = step3/D(jj,jj);
%         step5 = -step4;
%         wii_jj = step5;
%
% Condensed Vectorized Version of the Lambda calculation
%         wii_jj = -((D(jj,:)*L(:,1))+K(jj,1)-(D(jj,jj)*L(jj,1)))/D(jj,jj);
%-------------------------------------------------------------------------%

% Variable Defitions
xgen = zeros(size(H,1),1);
varnum = size(H,2);

%Modification of A and b matrix to add the nonnegativity constraints
nonneg = -eye(varnum);
A = [A; nonneg];

nonzero = zeros(varnum,1);
b = [b; nonzero];


%Begin Algorithm

%General Optimal Solution
xgen = -H\c;

%Check if General Optimal Solution violates constraints
C = (A*xgen > b);
if max(C) == 0
    x = xgen;
    z = (x'*H*x)/2 + c'*x;
    return;
end

%Calculate Lagrange Multiplers
D = A*(H\A');
K = b + A*(H\c);
n = size(K,1);
%m = calculateNumIterations(n);
L = zeros(n,1);
Lold = zeros(n,1);
state = true;
counter = 0;

%for ii = 1:250000
while (state == true)
    Lold = L;
    
   for jj = 1:n
       L(jj,1) = computeLambda(D,jj,L,K);   
   end
   
    difference = max(abs(Lold-L));
    if (difference<= 1e-2) 
        state = false;
        %break;
    end
    counter = counter+1;
end

%x and z values after applying Lagrange
x = xgen-H\A'*L;
z = (x'*H*x)/2 + c'*x;  
counter
end


function [lambda] = computeLambda(D,jj,L,K)
% Computes the Condensed Vectorized Version of the Lagrange Multipliers
% 1. Compute the w value
% 2. As the Lagrange multiplier must be positive, take the max of the w
%    calculated and 0.

%%%%%%%%%%%%%%%%
%           ii %
%wwii_jj = w   %
%           jj %
%%%%%%%%%%%%%%%%

w = -((D(jj,:)*L(:,1))+K(jj,1)-(D(jj,jj)*L(jj,1)))/D(jj,jj);
lambda = max(w,0);

end

function [m] = calculateNumIterations(n)
%Based on the number of constraints, the number of iterations is determined

if (n <= 999)
    m = 1000*n;
else
    m = 100*n;
end

end




