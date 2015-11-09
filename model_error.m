function [err,aux]=model_error(data,N,M,model,pars)
% function [err,aux]=model_error(data,N,M,model,pars)
% data = (empirical) gene frequency distribution
% N = sample size (number of genomes in sample)
% M = genome size (number of genes in genome)
% if model==1, model A -- constant population size
%     pars = theta
% if model==2, model B -- exponentially growing population
%     pars = [theta beta]
% if model==3, model C -- rigid core genome
%     pars = [frac theta] (frac = fraction of core genome)
% if model==4, model D -- flexible core genome
%     pars = [frac theta1 theta2] (frac and theta2 for core genome)

if model == 1,
    theta = pars;
    if theta < 0,
        err = Inf;
        return
    end
    aux = model_A_genefreq(N,M,theta);
elseif model == 2,
    th = pars(1);
    be = pars(2);
    if th < 0 || be < 0,
        err = Inf;
        return
    end
    aux = model_B_genefreq(N,M,th,be);
elseif model == 3,
    frac = pars(1);
    theta = pars(2);
    if frac < 0 || frac > 1 || theta < 0,
        err = Inf;
        return
    end
    aux = model_C_genefreq(N,M,frac,theta);
elseif model == 4,
    frac = pars(1);
    theta1 = pars(2);
    theta2 = pars(3);
    if frac < 0 || frac > 1 || ...
       theta1 < 0 || theta2 < 0 || theta1 < theta2,
        err = Inf;
        return
    end
    aux = model_D_genefreq(N,M,frac,theta1,theta2);
else
    disp('not an appropriate model code')
    return
end

err = sum((sqrt(data)-sqrt(aux)).^2)/N;

end
