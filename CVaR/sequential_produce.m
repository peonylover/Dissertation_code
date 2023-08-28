%for CVaR
function [x, t, cvx_optval, sub_cl_low, sub_cl_up, true_cl_low, true_cl_up, i_b, i_b_v, sort_all_sfy] = sequential_produce (b,a,ret,retm_ascii)

%Input
%b: givern constant
%a: risk level,such as 0.95, 0.99
%ret: return matrix for scenario subset
%retm_ascii: return matrix for all scenario

%Output
%x: optimal solution
%t: Var value
%cvx_optval: optimal value
%sub_cl_low: lower bound of 1 minus the true confidence level for subproblem
%sub_cl_up: upper bound of 1 minus the true confidence level for subproblem
%true_cl_low: lower bound of 1 minus the true significance level 
%true_cl_up: upper bound of 1 minus the true significance level
%i_b: index position corresponding to the value-at risk from the ordered loss subset scenario
%i_b_v: the index position from the ordered loss list considering all scenarios
%sort_all_sfy: the ordered data matrix which store loss and corresponding return for all scenario

%% Calculate the number of scenarios in subset and all scenario
[n_subset,~] = size(ret); 
[n_total,~] = size(retm_ascii);

%% Initial parameter
B = 10000;
L = 10100;
 
%% Solve cVar model subproblem
cvx_begin
  variables x(8);
  variables sfy(n_subset);
  variables t(1);

  minimize (t + norm(sfy,1)/(n_total*(1-a)));
  subject to
  B*ret*x +sfy >= L-t;
  x>= 0;
  sfy>=0;
  sum(x) == 1;
cvx_end

%% Scenarios in subset
%calculate index position for subset 
sum_p = 0;
i_a = 0;
while sum_p < a
    sum_p =sum_p +1/n_subset;
    i_a = i_a +1;
end
i_b = n_subset -i_a +1;

%Order loss in subset
%To avoid calculation error in sfy by cvx, recalculate loss for subset
re_sfy = max(L - B*ret*x-t,0);
ordered_sfy = sort(re_sfy);

%Find loss value corresponding i_a
z_ia = ordered_sfy(i_a);

%Calculate bounds for subproblem significance level
sub_cl_up = i_b/n_subset;
sub_cl_low = sub_cl_up- 1/n_subset;

%% All scenario
%Calculate loss for all scenario
all_sfy = max(L - B*retm_ascii*x-t,0);

%Creat matrix to store loss and corresponding return for all scenario
data_all_sfy = [all_sfy retm_ascii];
%Sort all loss
sort_all_sfy = sortrows(data_all_sfy,1);
%Calculate index position
[i_a_v, ~] = find(sort_all_sfy == z_ia);
i_b_v = n_total - i_a_v + 1;

%Calculate bounds for true significance level
true_cl_up = i_b_v/n_total;
true_cl_low = true_cl_up- 1/n_total;
end