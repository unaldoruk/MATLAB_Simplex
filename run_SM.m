function [z, x] = run_SM(H, c, A, b)
%   minimize  (1/2)*transpose(x)*H*x + transpose(c) * x
%   subject to  A * x <= b ,  x >= 0

% set tolerance for min(r)
tol = 0;

% Get dimensions of input matrices
[size_H, ~] = size(H);
[size_C] = size(c);
[size_A_rows, size_A_cols] = size(A);
[size_B] = size(b);

%This matrix tracks what's in the Basis
left = zeros(size_H + size_A_rows,1);

%We create a matrix (Simplex Tableau) in the form of:
%
%          |  |H A^T (-I) 0 (-c) |
%          |  |A 0     0  I   b  |
%          |    reduced costs    |

Simplex = zeros((size_B + size_C + 1), (3*size_H + 3*size_A_rows + 1));

[size_S_rows, size_S_cols] = size(Simplex);

comp = zeros(2*(size_H + size_A_rows),1);

used = zeros(1, size_S_cols-1);

Simplex(1:size_H, 1:size_H) = H;

Simplex(1:size_H, size_H +1 : size_H + size_A_rows) = transpose(A);

Simplex(1:(end-1), size_A_rows + size_H +1: 2*size_H+ 2*size_A_rows) = eye(size_H+ size_A_rows);

Simplex(1:(end-1), size_A_rows + size_H +1: 2*size_H+ size_A_rows) = -1 * Simplex(1:(end-1), size_A_rows + size_H +1: 2*size_H+ size_A_rows);

Simplex(size_H+1: end-1, 1:size_H) = A;

Simplex(1:size_C, end) = -c;

Simplex(size_C +1: end-1 , end) = b;


% Here we make sure the RHS of the Simplex Tableau is positive, if its not,
% then we multiply the corresponding row with (-1)

for ii = 1:size_S_rows-1
    
    if Simplex(ii, end) < 0
        
        Simplex(ii, 1:end) = -Simplex(ii, 1:end);
        
    end
    
end

%We construct the rest of the Simplex Tableau by introducing artificial
%variables:
%
%          |  |H A^T (-I) 0 (-c) I 0 |
%          |  |A 0     0  I   b  0 I |
%          |    reduced costs        |

Simplex(1:end-1, 2*size_H + 2*size_A_rows + 1 : end-1) = eye(size_H + size_A_rows);

% Calculating reduced costs by summing every element in a column and
% multiply?ng it with(-1)


for ii = 1: (2*size_H + 2*size_A_rows)
    
    Simplex(end, ii) = -1 * sum(Simplex(1:end-1, ii));
    
end

%Find smallest reduced cost and the entering variable

[smallest, smallestCol_No] = min(Simplex(end,1:end-1));

% Start optimization

while smallest < tol
    
    % We use for loops to calculate the Min_Ration
    
    smallestRatio = inf;
    
    %Find smallest ratio
    
    for ii = 1:size_S_rows-1
        
        if Simplex(ii, smallestCol_No) > 0
            
            if  Simplex(ii, end)/Simplex(ii, smallestCol_No) < smallestRatio
                
                smallestRatio = Simplex(ii, end)/Simplex(ii, smallestCol_No);
                
                smallestRatio_row = ii;
                
            end
            
        end
        
    end
    
    % Now we determine the position of the complement of the entering
    % variable in the Simplex Tableau. If the entering variable is x or M
    %(mu) then we "go_up", meaning that the index of the complement is higher than the index (m+n)
    
    if smallestCol_No <= size_H + size_A_rows
        
        go_up = true;
        
    else
        
        go_up = false;
        
    end
    
    
    % Here we determine if he complement of the entering variale is already
    % "in" the basis. Remember that "left" is an array that stores the LHS
    % variables
    in = false;
    
    for ii = 1:size(left)
        
        if go_up == true
            
            if left(ii) == smallestCol_No + size_H + size_A_rows
                
                in = true;
                
            end
            
        else
            
            if left(ii) == smallestCol_No - size_H - size_A_rows
                
                in = true;
                
            end
            
        end
        
    end
    
    %Here is where we determine if we do row operations or not
    
    if in == true
        
        if go_up == true
            
            %The next two if statements ensure that if the complement of the
            %entering variable is in the Basis, it is going to be the
            %departing variable in the same iteration the entering variable gets in the
            %Basis
            
            if comp(smallestCol_No) * comp(smallestCol_No + size_H + size_A_rows) == 0
                
                if left(smallestRatio_row) == comp(smallestCol_No + size_H + size_A_rows)
                    
                    % We call the function doSimplex to start row
                    % operations
                    [Simplex, comp, left, used] = doSimplex(Simplex, smallestRatio_row, smallestCol_No, size_S_rows, comp, left, used);
                    performed = true;
                    
                end
                
            end
            
        else
            
            if comp(smallestCol_No) * comp(smallestCol_No - size_H - size_A_rows) == 0
                
                if left(smallestRatio_row) == comp(smallestCol_No - size_H - size_A_rows)
                    
                    [Simplex, comp, left, used] = doSimplex(Simplex, smallestRatio_row, smallestCol_No, size_S_rows, comp, left, used);
                    performed = true;
                    
                end
                
            end
        end
        
        
        % If the complement of the entering variable is not in the Basis
    else
        
        [Simplex, comp, left, used] = doSimplex(Simplex, smallestRatio_row, smallestCol_No, size_S_rows, comp, left, used);
        performed = true;
        
    end
    
    %Here we check if any row operations were done, because if they are not,
    %we need to search for another entering variable. We store the used
    %entering variable in a matrix called "used" so we skip it when we're
    %searching for the min?mum reduced cost
    if performed ~= true
        
        used(1, smallestCol_No) = 1;
        
    end
    
    performed = false;
    
    %Find smallest
    
    temp = Simplex(end,1:end-1);
    
    found = false;
    
    % Here we make sure that we skip the reduced cost of the variable that
    % we determined cannot enter the basis.
    while found == false
        
        [smallest, smallestCol_No] = min(temp);
        
        if used(1, smallestCol_No) == 1
            
            temp(1, smallestCol_No) = inf;
            
        else
            
            found = true;
            
        end
        
    end
    
end

% The solution to the Quadratic programming problem will be the Right Hand
% Side (so the very right column) after all the iterations have been
% completed. We just take the X values from the right hand side

x = zeros(size_H, 1);

for ii = 1:size_S_rows-1
    
    if left(ii) <= size_H && left(ii) ~= 0
        
        x(left(ii)) = Simplex(ii, end);
        
    end
    
end

% We calculate the optimal objective function value
z = (1/2) * transpose(x) * H * x  + transpose(c) * x ;

end


% This function does the row operations, and modifies the Simplex Tableau.
% It makes the element in the column of the entering variable and in the row
% of the minimum ratio equal to 1 while the other elements in the column
% are made equal to zero
function [Simplex, comp, left, used] = doSimplex(Simplex, smallestRatio_row, smallestCol_No, size_S_rows, comp, left, used)

if Simplex(smallestRatio_row, smallestCol_No) ~= 1
    
    Simplex(smallestRatio_row, 1:end) = (1/ Simplex(smallestRatio_row, smallestCol_No)) * Simplex(smallestRatio_row, 1:end);
    
end

for ii = 1:size_S_rows
    
    if ii ~= smallestRatio_row && Simplex(ii, smallestCol_No) ~= 0
        
        Simplex(ii, 1:end) = Simplex(ii, 1:end) + (( -Simplex(ii, smallestCol_No)) * Simplex(smallestRatio_row, 1:end));
        
    end
    
end

comp(smallestCol_No) = 1;

% Reset these matrices for every iteration the row operations are performed
left(smallestRatio_row,1)= smallestCol_No;
used(1, 1:end) = 0;

end