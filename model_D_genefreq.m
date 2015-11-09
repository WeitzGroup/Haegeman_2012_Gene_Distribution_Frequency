function genefreq=model_D_genefreq(N,M,frac,theta1,theta2)
% function genefreq=model_D_genefreq(N,M,frac,theta1,theta2)
% N = sample size (number of genomes in sample)
% M = genome size (number of genes in genome)
% frac = fraction of flexible core genome
% theta1 = gene transfer parameter for part outside flexible core
% theta2 = gene transfer parameter for flexible core genome

aux = 1:N;
aux1 = M*theta1./aux.* ...
    exp(gammaln(N+1)-gammaln(N+1-aux) ...
    +gammaln(theta1+N-aux)-gammaln(theta1+N));
aux2 = M*theta2./aux.* ...
    exp(gammaln(N+1)-gammaln(N+1-aux) ...
    +gammaln(theta2+N-aux)-gammaln(theta2+N));
genefreq = (1-frac)*aux1+frac*aux2;

end