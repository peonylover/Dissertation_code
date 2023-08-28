% Importance Sampling algorithm  for expected shortfall
clear all;
load ('retm_ascii.mat',"-ASCII");
% record start time
t1 = clock;
%Initial parameter
B = 10000;
L = 9640;
%Choose different size of original scenarios
retm_ascii=retm_ascii(1:10000,:);
%Percentage of extra potential scenarios(10% here)
extra = 0.1;
%Creat a matrix to record result
all_result= ["oos"];
%Size of large enough initial scenario subset(10% of total)
initial_size = 0.1;

%% step1: Generate adequately large scenarios
%Calculate number of large enough scenario subset
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

  minimize (norm(sfy,1)/(n_large));
  subject to
  B*ret_large*x +sfy >= L;
  x>= 0;
  sfy>=0;
  sum(x) == 1;
cvx_end

%% step 3 include potential scenarios
% calculate and sort shortfall for large enough set 
re_sfy = max(L - B*ret_large*x,0); 
sort_re_sfy = sort(re_sfy);

%calculate the number of scenario with positive shortfall 
%on large enough scenario set
sub_nonzero_index = find(sort_re_sfy ~= 0);
critical_index_nonzero = sub_nonzero_index(1);
number_nonzero = n_large -critical_index_nonzero +1;

%calculate the number of needed extra scenarios
add_numb = round( number_nonzero * extra);

%use loop to find the number of scenarios that meet the required number
%and its new shortfall value
c = 0;
delta_c = 1;
initial_number_nonzero = number_nonzero;
while number_nonzero < ( initial_number_nonzero + add_numb)
    
    c = c + delta_c;
    
    %recalculate the number of scenarion with positive shortfall 
    %at current iteration
    re_sfy = max(L - B*ret_large*x + c , 0);
    sort_re_sfy = sort(re_sfy);
    sub_nonzero_index = find(sort_re_sfy ~= 0);
    critical_index_nonzero = sub_nonzero_index(1);
    number_nonzero = n_large -critical_index_nonzero +1;
end


%% Step 4 and step 5
% Assign pick up probability
cost= max(L - B*ret_large*x +c,0);
average= mean(cost);
p_u = cost/(average*n_large);

% we can use loop to test average performance
% of final reduced scenario subset
test_number = 1;
for i=1:test_number
    %choose 100 scenario for final scenario subset with specified picking up probability
    n_is = 100;
    index = linspace(1, n_large, n_large);
    out = randsrc(n_is, 1, [index;p_u']);
    %return and weight for small scenario set
    ret_is = ret_large(out,:);
    cost_is = cost(out,:);

    %use final reduced scenario susbet to solve opt 
    cvx_begin
        variables x_is(8);
        variables sfy_is(n_is);

        minimize (sum(sfy_is./cost_is));
        subject to
        B*ret_is*x_is +sfy_is >= L;
        x_is>= 0;
        sfy_is>=0;
        sum(x_is)==1;
    cvx_end
    
    %calculate out of sample 
    oos = mean(max(L-B*retm_ascii*x_is,0));

%record result
result = [oos];
all_result = [all_result;result];
end

%record computing time 
t2=clock;
t=t2-t1;





















