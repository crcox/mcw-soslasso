soslasso_results = zeros(10,10); 
lasso_results = zeros(10,10);
univariate_results = zeros(10,1);
for s=1:10
    for l=1:10
        soslasso_results(s,l) = m(s,l).soslasso; 
        lasso_results(s,l) = m(s,l).lasso;
    end
    univariate_results(s,1) = m(s,1).univariate_overall;
end

[X,Y] = meshgrid(1:10,1:10);

%% soslasso plot
figure(1)
mesh(X,Y,soslasso_results)
xlabel(gca,'sigma')
ylabel(gca,'lambda')
zlabel(gca,'dprime')
set(gca,'xtickLabel',sigma)
set(gca,'ytickLabel',lambda)
title('soslasso')

%% lasso plot
figure(2)
mesh(X,Y,lasso_results)
xlabel(gca,'sigma')
ylabel(gca,'lambda')
zlabel(gca,'dprime')
set(gca,'xtickLabel',sigma)
set(gca,'ytickLabel',lambda)
title('lasso')

%% univariate plot
figure(3)
plot(1:10,univariate_results)
xlabel(gca,'sigma')
zlabel(gca,'dprime')
set(gca,'xtickLabel',sigma)
title('univariate')