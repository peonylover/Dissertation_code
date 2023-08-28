%Choose scenario randomly
clear all;
load ('retm_ascii.mat',"-ASCII");
%Choose different size of original scenarios
retm_ascii = retm_ascii(1:100000,:);
%record beginning time
t1= clock;

% initial 
a = 0.95;

%choose scenario subset randomly
[n_total,~] = size(retm_ascii); 
n_subset = 100;
N_AS = randperm(n_total,n_subset);
%Obtain returen matrix for subset scenario
ret = retm_ascii(N_AS,:);
%n_total is number of selected scenario subset
[n,~] = size(ret); 


%parameter
B = 10000;
L = 10100;
 
%Solve opt
cvx_begin
  variables x(8);
  variables sfy(n);
  variables t(1);

  minimize (t + norm(sfy,1)/(n*(1-a)));
  subject to
  B*ret*x +sfy >= L-t;
  x>= 0;
  sfy>=0;
  sum(x) == 1;
cvx_end

%calculate out-of-sample
oos = t + mean(max(L - B*retm_ascii*x-t, 0))/(1-a);
%calculate total time
t2= clock;
t= t2-t1;

%re_sfy = max(round(L - B*retm_ascii*x-t,5),0);
%ordered_sfy = sort(re_sfy);
