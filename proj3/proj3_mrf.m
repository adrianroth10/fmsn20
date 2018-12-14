N1 = [0 1 0; 1 0 1; 0 1 0];
N2 = [1 1 1; 1 0 1; 1 1 1];
N3 = [0 0 0 1 1; 0 0 1 1 1; 0 1 0 1 0; 1 1 1 0 0; 1 1 0 0 0];
N = {N1,N2,N3};
alpha = 0;
beta = 0;
%[theta, prior] = normmix_gibbs(y_all, nc);


for i_neighbour = 1:length(N)
neighbours = N{i_neighbour};
alpha_post=mrf_gaussian_post(alpha,theta,y_all);
z=mrf_sim(z0,N,alpha_post,beta,iter);
[mu, Sigma] = gibbs_mu_sigma(y_all);


x = zeros(128,128,3);
iter = 100;
Plog = zeros(iter,1);
for i = 1:iter
x=mrf_sim(x,N3,0,1,1);
[~,Mz,Mf]=mrf_sim(x,N3,0,0.1,0);
Plog(i) = sum(log(Mz(logical(x))));
imagesc(x)
drawnow
end
end
%mrf_gaussian_post
%mrf_sim
%gibbs_alpha_beta
%gibbs_mu_Sigma
