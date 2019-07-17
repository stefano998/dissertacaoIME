m=20;n=10;   
randn("state",[99]);obsout1=1; obsout2=3;
L_perf=[163854.9;6446.2;57037.0;126209.5;101128.6;296885.8;398014.4;60449.1;173710.4;167264.2;
110227.2;155928.2;52875.0;62904.2;3889.5;42705.7;98891.2;115779.2;113222.5;46428.8];

A=[1	0	0	0	0	0	0	0	0	0;
-1	0	0	0	0	0	0	0	0	1;
0	0	0	0	0	1	0	0	0	-1;
0	0	0	0	0	1	-1	0	0	0;
0	0	0	0	0	0	1	0	0	0;
0	0	0	0	0	0	-1	1	0	0;
0	0	0	0	0	0	0	1	0	0;
0	0	0	0	0	0	0	1	-1	0;
-1	0	0	0	0	0	0	0	1	0;
0	0	0	0	0	0	0	0	1	-1;
0	0	0	0	0	-1	0	0	1	0;
0	0	0	0	1	0	0	0	0	-1;
-1	1	0	0	0	0	0	0	0	0;
0	-1	1	0	0	0	0	0	0	0;
0	0	-1	1	0	0	0	0	0	0;
0	0	0	-1	1	0	0	0	0	0;
0	0	0	0	1	-1	0	0	0	0;
-1	0	1	0	0	0	0	0	0	0;
0	0	0	1	0	0	0	0	0	-1;
0 1 0 0 0 0 0 0 0 -1];

I = eye(m);
A1 = [A -A -I I];
d=[49 41 38 34 22 13 23 48 15 24 62 49 35 43 20 28 19 39 27 21]';
dp=1.*sqrt(d);

  for q=1:m
    do a=randn(1);
    until (a<=3)
    #dp(q)*a
    L(q,1)=L_perf(q)+dp(q)*a;
end
L=round(L.* 10)./10;
L(obsout1)=L(obsout1)+50;    ###############outliers
L(obsout2)=L(obsout2)+200;   ###############outliers

c1=zeros(1,2*n);
p=ones(1,m);
c = cat(2,c1,p,p);

ctype=[];
for i=1:m
  ctype = [ctype, "S"];
 end
vartype=[];
for i=1:(2*(n+m))
  vartype = [vartype, "C"];
 end
param.dual=1; 
param.lpsolver=1; 
[xopt, fopt, erro, extra] = glpk (c, A1, L, lb=[], ub=[], ctype, vartype, s=1, param);
for i=1:m
  v(i)=xopt(2*n+i)-xopt(2*n+i+m);
end
abs(v)
#"erros totais dos outliers"
#L(obsout1)-L_perf(obsout1)
#L(obsout2)-L_perf(obsout2)