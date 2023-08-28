% Scenario reduction iteration algorithm for CVaR
clear all;
load ('retm_ascii.mat',"-ASCII");
%Choose different size of original scenarios
retm_ascii = retm_ascii(1:10000,:);
% record start time
t1= clock;
%% Initial parameter 
a = 0.95; % Risk level
b = 2; % Given constant
v = 1; % Iteration counter
cl = Inf; % Reference significance level for CL
delta_b = 0.5; % Delta b value
all_result = ["x1","x2","x3","x4","x5","x6","x7","x8","vAR","opt value","cl","sub_cl_lo","sub_cl_up","true_cl_low","true_cl_up","i_b","i_b_v","b","v"]; % Creat a matrix to recourd result

%% Randomly choose  scenario subset
%%Select N_AS, with equal probability 1/n, 
%%and statisfy condition sum(p)>= b(1-a)
[n_total,~] = size(retm_ascii); 
%Calculate the number of scenarios in subset with equal probability
n_subset = ceil(round((1-a)*b*n_total,5));
%Select subset scenario randomly
N_AS = randperm(n_total,n_subset);
%Obtain returen matrix for subset scenario
ret = retm_ascii(N_AS,:);

%% First iteration
[new_x, t, opt, sub_cl_low, sub_cl_up, true_cl_low, new_cl, i_b, i_b_v, sort_all_sfy]= sequential_produce(b,a,ret,retm_ascii);

%Recourd the result
result = [new_x', t, opt, cl, sub_cl_low, sub_cl_up, true_cl_low, new_cl, i_b, i_b_v, b, v];
all_result = [all_result; result];

%% Other iteration
while i_b ~= i_b_v
   
    if new_cl < cl
        %Update reference significance bound if true bound decrease
        cl = new_cl;
    else
        %Otherwise, increase b value 
        b = b + delta_b;
    end 
    
    %%Choose new scenarios
    %Calculate number of scenario in new subset
    n_subset = ceil(round((1-a)*b*n_total,5));
    %Calculate the upper and lower index to choose new scenario
    index_upper = n_total-n_subset+1;
    index_lower = n_total;
    %Update return matrix for new scenario subset
    ret = sort_all_sfy(index_upper:index_lower, 2:9);
    
    %Update iteration counter
    v = v + 1;
    
    %Repeat sequential produce
    [new_x, t, opt, sub_cl_low, sub_cl_up, true_cl_low, new_cl, i_b, i_b_v, sort_all_sfy]= sequential_produce(b,a,ret,retm_ascii);
    
    %Update result
    result = [new_x', t, opt, cl, sub_cl_low, sub_cl_up, true_cl_low, new_cl, i_b, i_b_v, b, v];
    all_result = [all_result; result];
end    
%% record total time
t2= clock;
t=t2-t1;

