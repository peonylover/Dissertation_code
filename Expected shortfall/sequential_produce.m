% Scenario reduction iteration algorithm for expected shortfall
function [x, cvx_optval, s_p, s_p_a, sort_all_sfy, n_subset, oos] = sequential_produce (ret,retm_ascii, B, L)

%Input
%B: Budget
%L: return liability
%ret: return matrix for selected scenario subset
%retm_ascii: return matrix for all scenario

%Output
%x: optimal solution
%cvx_optval: optimal value
%s_p: the number of scenario with positive loss for scenario subset
%s_p_a: the number of scenario with positive loss for total scenario
%sort_all_sfy: the ordered data matrix which store loss and corresponding return for all scenario
%n_subset: the size of scenario subset
%oos: out-of-sample shortfall corresponding to the optimal solution
%     of subset scenario

%% Calculate the number of scenarios in subset and all scenario
[n_subset,~] = size(ret); 
[n_total,~] = size(retm_ascii);
 
%%  step 3.2 :Solve ES model subproblem
cvx_begin
  variables x(8);
  variables sfy(n_subset);

  minimize (norm(sfy,1)/(n_total));
  subject to
  B*ret*x +sfy >= L;
  x>= 0;
  sfy>=0;
  sum(x) == 1;
cvx_end

%%  step 3.3 Scenarios in subset
%To avoid calculation error in sfy by cvx, recalculate loss for subset
re_sfy = L - B*ret*x;
%Order loss in subset
ordered_sfy = sort(re_sfy);

%find positive shortfall set 
nonzero_index = find(ordered_sfy > 0);

%calculate the number of scenario with positive loss in susbet 
s_c = nonzero_index(1);
s_p = n_subset - s_c +1;

%% step 3.4 and partial step 5 ( for All scenario )
%Calculate loss for all scenario
all_sfy = L - B*retm_ascii*x;
%Creat matrix to store loss and corresponding return for all scenario
data_all_sfy = [all_sfy retm_ascii];
%Sort matrix accoridng to loss
sort_all_sfy = sortrows(data_all_sfy,1);
%Calculate the number of scenario with positive loss in total scenario
[s_c_a, ~] = find(sort_all_sfy > 0);
s_p_a = n_total - s_c_a(1) + 1;

%% Calculate out-of-sample shortfall
oos =mean(max(L-B*retm_ascii*x,0));

end