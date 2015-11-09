function [pars,err,aux]=model_fitparams(data,model)
% function [pars,err,aux]=model_fitparams(data,model)
% data = (empirical) gene frequency distribution
% if model==1, model A -- constant population size
%     pars = theta
% if model==2, model B -- exponentially growing population
%     pars = [theta beta]
% if model==3, model C -- rigid core genome
%     pars = [frac theta] (frac = fraction of core genome)
% if model==4, model D -- flexible core genome
%     pars = [frac theta1 theta2] (frac and theta2 for core genome)

N = length(data);
aux = 1:N;
M = abs(sum(data.*aux))/N;

normdata = data/M;
if model == 1,
    pars0 = .5;
elseif model == 2,
    pars0 = [.5 .5];
elseif model == 3,
    pars0 = [.5 .5];
elseif model == 4,
    pars0 = [.5 2 .5];
else
    disp('not an appropriate model code')
    return
end

opts = optimset('TolFun',1e-8,'TolX',1e-8, ...
    'MaxFunEvals',1e6,'MaxIter',1e6);
pars = fminsearch(@(x) model_error( ...
    normdata,N,1,model,x),pars0,opts);
[err,aux] = model_error(data,N,M,model,pars);

end