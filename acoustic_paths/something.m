a1 = load("g1.mat");
a2 = load("g2.mat");
aa_2 = a2.g2;
aa_1 = a1.g1;
% [aa_sort,index_sorted] = sort(aa_1,"descend","ComparisonMethod","abs");
% L = length(aa_1);
% x = aa_sort(1);
% y=0;
% for i = 1:L
%     y = y+(aa_sort(L-i+1));
%     X = ["i",i,"y",y];
%     disp(X);
%     if(y>=x)
%         break;
%     end
% end
% M = mean(aa_1)
% V = var(aa_1)

obs_length = 100;
X = zeros(aa_1.length,obs_length);

for i=1:obs_length
