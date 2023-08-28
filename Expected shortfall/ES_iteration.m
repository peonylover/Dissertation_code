% Scenario reduction iteration algorithm for expected shortfall
clear all;
load ('retm_ascii.mat',"-ASCII");
%Choose different size of original scenarios
retm_ascii=retm_ascii(1:10000,:);
%Record beginning time
t1 = clock;
%% Step 2: Setting parameter
%Parameter
B = 10000;
L = 9640;
rng(3);
%The size of initial subset, here we choose 10% of total
initial= 0.1;
%The parameter of Heuristic procedure
b = 1.1;
delta_b= 0.1;
% Iteration counte
v = 1; 
% Reference number
RN = -inf; 
% Creat a matrix to recourd result
all_result = ["x1","x2","x3","x4","x5","x6","x7","x8","opt value","RN","b","s_p","s_p_a","v","number subset","oos"]; 

%% Randomly choose initial scenario subset (Step 1: preparation)
%scenario  with equal probability 1/n
%Calculate the number of scenarios in subset 
[n_total,~] = size(retm_ascii); 
n_subset = ceil(round(initial*n_total,5));
%Select subset scenario randomly
N_AS = randperm(n_total,n_subset);
%Obtain returen matrix for initial subset scenario
ret = retm_ascii(N_AS,:);
%Obtain and recourd the result
[new_x, opt,   s_p, s_p_a, sort_all_sfy, n_subset,oos] = sequential_produce(ret,retm_ascii, B, L);
result = [new_x', opt,  RN, b, s_p, s_p_a, v, n_subset,oos];
all_result = [all_result; result];
%% Other iteration
%convergence checking
while s_p ~= s_p_a
    
    %Heuristic procedure
    if  RN > n_subset
        %Increase b value
        b = b + delta_b;
    end
    
    RN = n_subset;
    
    %% step 3.1
    %Calculate the size of scenario subset in next iteration 
    n_subset = ceil( b * s_p_a );
    %Calculate the upper and lower index to choose new scenario
    index_upper = n_total-n_subset+1;
    index_lower = n_total;
    %Update return matrix for new scenario subset 
    ret = sort_all_sfy(index_upper:index_lower, 2:9);
    
    %%partial step 5, Update iteration counter
    v = v + 1;
    
    %Repeat sequential procedure
    [new_x, opt,  s_p, s_p_a, sort_all_sfy, n_subset,oos]= sequential_produce(ret,retm_ascii, B, L);
    
    %Update result
    result = [new_x', opt, RN,b, s_p, s_p_a, v, n_subset,oos];
    all_result = [all_result; result];

end    

%Calculate total time
t2 = clock;
t = t2-t1;