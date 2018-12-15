iterations = 1000;

% Neighbourhood
neighbours1 = [0 1 0; 1 0 1; 0 1 0];
neighbours2 = [1 1 1; 1 0 1; 1 1 1];
neighbours3 = [0 0 0 1 1; 0 0 1 1 1; 0 1 0 1 0; 1 1 1 0 0; 1 1 0 0 0];
%neighbours_set = {neighbours1, neighbours2, neighbours3};
neighbours_set = {neighbours1};

Zmatrix = zeros(length(y_all), 1, nc, length(neighbours_set));
alpha = zeros(nc, iterations + 1, length(neighbours_set));
beta = zeros(iterations + 1, length(neighbours_set));
Plog = zeros(iterations, length(neighbours_set));
acc = zeros(iterations, length(neighbours_set));

%[theta, ~] = normmix_gibbs(y_all, nc);
theta = cell(nc,1);
for k = 1:nc
    theta{k}.mu = zeros(length(component),1);
    theta{k}.Sigma = eye(length(component));
end


for i_neighbours = 1:length(neighbours_set)
    zmat = Zmatrix(:, :, :, i_neighbours);
    neighbours = neighbours_set{i_neighbours};
    for iter = 1:iterations
        alpha_post = mrf_gaussian_post(alpha(:, iter, i_neighbours)', theta, y_all);
        zmat = mrf_sim(zmat,neighbours, alpha_post, beta(iter, i_neighbours), 1);
        [~,Mz,Mf] = mrf_sim(zmat, neighbours, alpha_post, beta(iter, i_neighbours), 0);
        Plog(iter, i_neighbour) = sum(log(Mz(logical(zmat))));
        for k = 1:nc
            [mu, Sigma] = gibbs_mu_sigma(y_all(logical(zmat(:, :, k))));
            theta{k}.mu = mu;
            theta{k}.Sigma = Sigma;
        end
        [alpha(2:end, iter + 1, i_neighbours), beta(iter + 1, i_neighbours), acc(iter, i_neighbours)] = ...
                      gibbs_alpha_beta(alpha(2:end, iter, i_neighbours)', beta(iter, i_neighbours), zmat, neighbours, 10, 3.5e-2);
    end
    Zmatrix(:,:,:,i_neighbour) = zmat;
    figure();
    imagesc(reshape(Zmatrix(:,1,nc,i_neighbour), sz(1:2)))
    print(['proj3/output/mrf', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta),'_', num2str(i_neighbour), '.png'], '-dpng');
    close;
end
acc = mean(acc)
burn_in = 200;

plot(Plog)
print(['proj3/output/mrf', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta), '.png'], '-dpng');
close;

figure()
plot(Plog)
figure()
hold on
for i_neighbours = 1:length(neighbours_set)
  plot(alpha(:, :, i_neighbours))
end
figure()
plot(beta)
