function PlotTimes(number)

n = zeros(number-1,1);

%SM

SM_time = zeros(number-1, 1);

SM_z = zeros(number-1, 1);

SM_nn = zeros(number-1, 1);

SM_lin = zeros(number-1, 1);

SM_minx = zeros(number-1, 1);

SM_max = zeros(number-1, 1);

SM_dev = zeros(number-1, 1);


%HQP


H_time = zeros(number-1 , 1);

H_z = zeros(number-1 , 1);

H_nn = zeros(number-1, 1);

H_lin = zeros(number-1, 1);

H_minx = zeros(number-1, 1);

H_max = zeros(number-1, 1);

H_dev = zeros(number-1, 1);




%CG

CG_time = zeros(number-1 , 1);

CG_z = zeros(number-1 , 1);

CG_nn = zeros(number-1, 1);

CG_lin = zeros(number-1, 1);

CG_minx = zeros(number-1, 1);

CG_max = zeros(number-1, 1);

CG_dev = zeros(number-1, 1);


for ii = 2:number
    
    n(ii-1)=ii;
    
    data = dataGenerator(ii);
    
    [~,fval] = quadprog(data.H, data.c, data.A, data.b, [], [], data.lb, [], []);
    
    
    tic;
    
    [z, x] = run_SM(data.H, data.c, data.A, data.b);
    
    SM_time(ii-1) = toc;
    
    SM_z(ii-1) = z;
    
    SM_nn(ii-1) = findNoNonNegViolation(x);
    
    SM_lin(ii-1) = findNoLinViolation(x, data.A, data.b);
    
    SM_minx(ii-1) = min(x);
    
    SM_max(ii-1) = max((data.A) * x - data.b);
        
    SM_dev(ii-1) = abs((z - fval) / fval * 100);
    
    
    
    tic;
    
    [z, x] = HQP(data.H, data.c, data.A, data.b);
    
    H_time(ii-1) = toc;
    
    H_z(ii-1) = z;
    
    H_nn(ii-1) = findNoNonNegViolation(x);
    
    H_lin(ii-1) = findNoLinViolation(x, data.A, data.b);
    
    H_minx(ii-1) = min(x);
    
    H_max(ii-1) = max((data.A) * x - data.b);
    
    H_dev(ii-1) = abs((z - fval) / fval * 100);
    
    
    tic;
    
    [z, x] = CG_newlas2(data.H, data.c, data.A, data.b);
    
    CG_time(ii-1) = toc;
    
    CG_z(ii-1) = z;
    
    CG_nn(ii-1) = findNoNonNegViolation(x);
    
    CG_lin(ii-1) = findNoLinViolation(x, data.A, data.b);
    
    CG_minx(ii-1) = min(x);
    
    CG_max(ii-1) = max((data.A) * x - data.b);
    
    CG_dev(ii-1) = abs((z - fval) / fval * 100);
    
    
    
    
    ii = ii+4;
    
    
    
end

%time

figure

plot(n(:) , SM_time(:),'-bo')

xlabel('Dimension of H');

ylabel('Time (sec)');

title('Time comparison');

hold all;

plot(n(:) , H_time(:),'-rx')

plot(n(:) , CG_time(:),'-gs')

legend( [  {'SM'} {'HQP'} {'CG'}  ] );

grid;



%z

figure

plot(n(:) , SM_z(:),'-bo')

xlabel('Dimension of H');

ylabel('Z');

title('Optimal Objective Function Value Comparison');

hold all;

plot(n(:) , H_z(:),'-rx')

plot(n(:) , CG_z(:),'-gs')

legend( [  {'SM'} {'HQP'} {'CG'}  ] );

grid;


%Non Neg

figure

plot(n(:) , SM_nn(:),'-bo')

xlabel('Dimension of H');

ylabel('Number of Non-Negativity Constraints Violated');

title('Number of Non-Negativity Constraints Violated Comparison');

hold all;

plot(n(:) , H_nn(:),'-rx')

plot(n(:) , CG_nn(:),'-gs')

legend( [  {'SM'} {'HQP'} {'CG'}  ] );

grid;



%Lin violated

figure

plot(n(:) , SM_lin(:),'-bo');

xlabel('Dimension of H');

ylabel('Number of Linear Constraints Violated');

title('Number of Linear Constraints Violated Comparison');

hold all;

plot(n(:) , H_lin(:),'-rx');

plot(n(:) , CG_lin(:),'-gs');

legend( [  {'SM'} {'HQP'} {'CG'}  ] );

grid;


% min x

figure

plot(n(:) , SM_minx(:),'-bo');

xlabel('Dimension of H');

ylabel('Min(x)');

title('Minimum X Value Comparison');

hold all;

plot(n(:) , H_minx(:),'-rx');

plot(n(:) , CG_minx(:),'-gs');

legend( [  {'SM'} {'HQP'} {'CG'}  ] );

grid;


% Max Ax-b

figure

plot(n(:) , SM_max(:),'-bo');

xlabel('Dimension of H');

ylabel('Max(Ax - b)');

title('Maximum Ax - b Comparison');

hold all;

plot(n(:) , H_max(:),'-rx');

plot(n(:) , CG_max(:),'-gs');

legend( [  {'SM'} {'HQP'} {'CG'}  ] );

grid;


% dev

figure

plot(n(:) , SM_dev(:),'-bo');

xlabel('Dimension of H');

ylabel('Deviation (%)');

title('Deviation from Quadprog optimal obj Comparison');

hold all;

plot(n(:) , H_dev(:),'-rx');

plot(n(:) , CG_dev(:),'-gs');

legend( [  {'SM'} {'HQP'} {'CG'}  ] );

grid;

end

function counter = findNoNonNegViolation(x)

counter = 0;

for it = 1:size(x)
    
    if x(it) < 0
        
        counter = counter + 1;
        
    end
    
end

end

function counter = findNoLinViolation(x, A, b)

counter = 0;

temper = A * x - b;

for it = 1:size(temper)
    
    if temper(it) > 0
        
        counter = counter + 1;
        
    end
    
end

end
