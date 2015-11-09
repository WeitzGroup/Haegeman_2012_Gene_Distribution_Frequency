function genefreq=model_C_genefreq(N,M,frac,theta)
% function genefreq=model_C_genefreq(N,M,frac,theta)
% N = sample size (number of genomes in sample)
% M = genome size (number of genes in genome)
% frac = fraction of rigid core genome
% theta = gene transfer parameter

aux = 1:N;
genefreq = M*theta./aux.* ...
    exp(gammaln(N+1)-gammaln(N+1-aux) ...
    +gammaln(theta+N-aux)-gammaln(theta+N));
genefreq = (1-frac)*genefreq;
genefreq(end) = genefreq(end)+M*frac;

end