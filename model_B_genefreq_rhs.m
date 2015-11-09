function dx = calc_popstr_varsize_analy_rhs(t,x,th,be,m1a,i1a,m1b,i1b,m2,i2,j2,m3,i3,j3,m4,i4,j4,m5,i5,m6,i6,j6)
% function dx = calc_popstr_varsize_analy_rhs(t,x,th,be,m1a,i1a,m1b,i1b,m2,i2,j2,m3,i3,j3,m4,i4,j4,m5,i5,m6,i6,j6)

dx=zeros(size(x));

dx(1)=dx(1)+th*sum(m1b'.*x(i1b));
dx(j2)=dx(j2)+th*m2'.*x(i2);
dx(j3)=dx(j3)+exp(be*t)*m3'.*x(i3);
dx(j4)=dx(j4)+exp(be*t)*m4'.*x(i4);
dx(j6)=dx(j6)+exp(be*t)*m6'.*x(i6);

dx(i1a)=dx(i1a)-th*m1a'.*x(i1a);
dx(i1b)=dx(i1b)-th*m1b'.*x(i1b);
dx(i2)=dx(i2)-th*m2'.*x(i2);
dx(i3)=dx(i3)-exp(be*t)*m3'.*x(i3);
dx(i4)=dx(i4)-exp(be*t)*m4'.*x(i4);
dx(i5)=dx(i5)-exp(be*t)*m5'.*x(i5);
dx(i6)=dx(i6)-exp(be*t)*m6'.*x(i6);

end

