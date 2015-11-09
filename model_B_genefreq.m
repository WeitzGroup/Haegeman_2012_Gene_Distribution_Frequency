function popstr=model_B_genefreq(N,M,th,be)
% function popstr=model_B_genefreq(N,M,th,be)
% N = sample size (number of genomes in sample)
% M = genome size (number of genes in genome)
% th = gene transfer parameter
% be = population growth parameter
% if the number of genomes N is large (e.g., N=20),
% computation can take a considerable amount of time.

if N==1,
    popstr=M;
    return
end

n=N; popstr=zeros(1,n);
for k=1:n,
    % construct index vectors
    qq=[]; mm=[]; ll=[];
    for q=1:n,
        for m=max([1 k+q-n]):q,
            aux=max([1 k+q-n]):min([k m]);
            ll=[ll aux];
            mm=[mm m*ones(1,length(aux))];
            qq=[qq q*ones(1,length(aux))];
        end
    end
    % construct transitions
    m1a=zeros(1,sum(ll>1)); i1a=m1a; j1a=m1a; c1a=0;
    m1b=zeros(1,sum((ll==1)&(qq>1))); i1b=m1b; j1b=m1b; c1b=0;
    m2=zeros(1,sum(mm>ll)); i2=m2; j2=m2; c2=0;
    m3=zeros(1,sum(ll>1)); i3=m3; j3=m3; c3=0;
    m4=zeros(1,sum(mm-ll>1)); i4=m4; j4=m4; c4=0;
    m5=zeros(1,sum(mm>ll)); i5=m5; j5=m5; c5=0;
    m6=zeros(1,sum(qq>mm)); i6=m6; j6=m6; c6=0;
    for idx=1:length(ll),
        l=ll(idx); m=mm(idx); q=qq(idx);
        % TRANSITION 1a: prefactor theta
        if l>1,
            c1a=c1a+1;
            m1a(c1a)=l/2;
            i1a(c1a)=idx;
            j1a(c1a)=-1;
        end
        % TRANSITION 1b: prefactor theta
        if (l==1)&&(q>1),
            c1b=c1b+1;
            m1b(c1b)=1/2;
            i1b(c1b)=idx;
            j1b(c1b)=1;
        end
        % TRANSITION 2: prefactor theta
        if m>l,
            c2=c2+1;
            m2(c2)=(m-l)/2;
            i2(c2)=idx;
            j2(c2)=find((ll==l)&(mm==m-1)&(qq==q));
        end
        % TRANSITION 3: prefactor N(0)/N(t)
        if l>1,
            c3=c3+1;
            m3(c3)=l*(l-1)/2;
            i3(c3)=idx;
            j3(c3)=find((ll==l-1)&(mm==m-1)&(qq==q-1));
        end
        % TRANSITION 4: prefactor N(0)/N(t)
        if m-l>1,
            c4=c4+1;
            m4(c4)=(m-l)*(m-l-1)/2;
            i4(c4)=idx;
            j4(c4)=find((ll==l)&(mm==m-1)&(qq==q-1));
        end
        % TRANSITION 5: prefactor N(0)/N(t)
        if m>l,
            c5=c5+1;
            m5(c5)=l*(m-l);
            i5(c5)=idx;
            j5(c5)=-1;
        end
        % TRANSITION 6: prefactor N(0)/N(t)
        if q>m,
            c6=c6+1;
            m6(c6)=(q-m)*(q+m-1)/2;
            i6(c6)=idx;
            j6(c6)=find((ll==l)&(mm==m)&(qq==q-1));
        end
    end
    % solve differential equation
    y0=zeros(1,length(ll));
    y0(end)=1;
    tf=40; if be>0, tf=log(1+be*tf)/be; end
%   Updated with ode45 instead of ode15s to improve numerical
%   robustness, particularly with large datasets.
    [~,y]=ode45(@(t,x) model_B_genefreq_rhs(t,x,th,be, ...
        m1a,i1a,m1b,i1b,m2,i2,j2,m3,i3,j3,m4,i4,j4,m5,i5,m6,i6,j6), ...
        [0 .5 1]*tf,y0);
    % extract population structure
    popstr(k)=exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1))*y(3,1);
end

popstr=M*popstr;

end
