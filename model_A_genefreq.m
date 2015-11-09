function genefreq=model_A_genefreq(N,M,theta)
% function genefreq=model_A_genefreq(N,M,theta)
% N = sample size (number of genomes in sample)
% M = genome size (number of genes in genome)
% theta = gene transfer parameter

aux = 1:N;
genefreq = M*theta./aux.* ...
    exp(gammaln(N+1)-gammaln(N+1-aux) ...
    +gammaln(theta+N-aux)-gammaln(theta+N));

end