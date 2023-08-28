% Importance Sampling algorithm  for CVaR
clear all;
load ('retm_ascii.mat',"-ASCII");
%Choose different size of original scenarios
retm_ascii=retm_ascii(1:10000,:);
% record start time
t1= clock;

%Initial parameter
B = 10000;
L = 10100;
a = 0.95;

%potential scenario
extra = 0.05;
all_result= ["oos"];
%Size of large enough initial scenario subset(10% of total)
initial_size = 0.1;

%% step 1 : Generate adequately large scenarios
[n_total,~] = size(retm_ascii); 
n_large = initial_size*n_total;

%Assume random choose scenario here
%rng(1);
N_AS = randperm(n_total, n_large);
ret_large = retm_ascii(N_AS,:);

%% step 2 solve sub-problem
cvx_begin
  variables x(8);
  variables sfy(n_large);
  variable t(1);

  minimize (t + norm(sfy,1)/(n_large*(1-a)));
  subject to
  B*ret_large*x +sfy >= L -t;
  x>= 0;
  sfy>=0;
  sum(x) == 1;
cvx_end

%% step 3 include potential scenarios 
% calculate and sort (loss-t)+  for large enough set  
re_sfy = max(L - B*ret_large*x -t ,0); 
sort_re_sfy = sort(re_sfy);

%calculate the number of scenario with positive (loss-t)+ 
%on large enough scenario set
sub_nonzero_index = find(sort_re_sfy ~= 0);
critical_index_nonzero = sub_nonzero_index(1);
number_nonzero = n_large -critical_index_nonzero +1;

%calculate the number of needed extra scenarios
add_numb = round( number_nonzero * extra);

%use loop to find the number of scenarios that meet the required number
%and its new (loss-t)+ value
c = 0;
delta_c = 1;
initial_number_nonzero = number_nonzero;
while number_nonzero < ( initial_number_nonzero + add_numb)
    
    c = c + delta_c;
   
    %recalculate the number of scenario with positive (loss-t)+ value
    %at current iteration
    re_sfy = max(L - B*ret_large*x + c - t  , 0);
    sort_re_sfy = sort(re_sfy);
    sub_nonzero_index = find(sort_re_sfy ~= 0);
    critical_index_nonzero = sub_nonzero_index(1);
    number_nonzero = n_large -critical_index_nonzero +1;
end


%% Step 4 and step 5
% Assign pick up probability
cost= max(L - B*ret_large*x + c,0);
average= mean(cost);

p_u = cost/(average*n_large);

% we can use loop to test average performance
% of final reduced scenario subset
test_number =1 ;
for i=1:test_number
    %choose 100 scenario with specified picking up probability
    n_is = 100;
    index = linspace(1, n_large, n_large);
    out = randsrc(n_is, 1, [index;p_u']);
    %return and cost for small scenario set
    ret_is = ret_large(out,:);
    cost_is = cost(out,:);

    % use final reduced scenario susbet to solve opt 
    cvx_begin
        variables x_is(8);
        variables sfy_is(n_is);
        variable t_is(1);

        minimize (t_is + sum(sfy_is./(cost_is*(1-a))));
        subject to
        B*ret_is*x_is + sfy_is >= L-t_is;
        x_is >= 0;
        sfy_is >= 0;
        sum(x_is) ==1;
    cvx_end
    
    %calculate out-of-sample
    oos = t_is + mean(max(L - B*retm_ascii*x_is - t_is, 0))/(1-a);
    
%record result  
result = [oos];
all_result = [all_result;result];
end

%record computing time
t2= clock;
total=t2-t1;























