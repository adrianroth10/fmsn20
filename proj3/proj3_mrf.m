N1 = [0 1 0; 1 0 1; 0 1 0];
N2 = [1 1 1; 1 0 1; 1 1 1];
N3 = [0 0 0 1 1; 0 0 1 1 1; 0 1 0 1 0; 1 1 1 0 0; 1 1 0 0 0];
N = {N1,N2,N3};
alpha = 0;
beta = 0;
iterations = 1000;
%[theta, prior] = normmix_gibbs(y_all, nc);
Zmatrix = zeros(length(y_all),1,nc,length(N));
Plog = zeros(iterations,length(N));
theta = cell(nc,1);
for k = 1:nc
    theta{k}.mu = zeros(length(component),1);
    theta{k}.Sigma = eye(length(component));
end

for i_neighbour = 1:length(N)
    zmat = Zmatrix(:,:,:,i_neighbour);
    neighbours = N{i_neighbour};
    for iter = 1:iterations
        alpha_post=mrf_gaussian_post(alpha,theta,y_all);
        zmat = mrf_sim(zmat,neighbours,alpha_post,beta,1);
        [~,Mz,Mf] = mrf_sim(zmat,neighbours,alpha_post,beta,0);
        Plog(iter, i_neighbour) = sum(log(Mz(logical(zmat))));
        for k = 1:nc
            [mu, Sigma] = gibbs_mu_sigma(y_all(logical(zmat(:,:,k))));
            theta{k}.mu = mu;
            theta{k}.Sigma = Sigma;
        end
        [alpha, beta, acc] = gibbs_alpha_beta(alpha, beta, zmat, Mf, 10);
    end
    Zmatrix(:,:,:,i_neighbour) = zmat;
    figure();
    imagesc(reshape(Zmatrix(:,1,nc,i_neighbour), sz(1:2)))
    print(['proj3/output/mrf', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta),'_', num2str(i_neighbour), '.png'], '-dpng');
    close;
end


plot(Plog)
print(['proj3/output/mrf', num2str(i_component), '_', num2str(nc), '_', num2str(is_beta), '.png'], '-dpng');
close;