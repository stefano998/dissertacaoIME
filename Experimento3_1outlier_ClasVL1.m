fator=3;  # [|v|:1; F:2; F':3] #obs: F=|v|/median(|v|);F'=|v|/mad(|v|);
VC=12.9     #obs: VCL1 calculados: |v|:29.2; F:6.9; F':12.9

for s=1:4
rand("state",[0]);randn("state",[0]);
  if s==1
    ini_int_erro=3; fim_int_erro=6;
  elseif s==2
    ini_int_erro=6; fim_int_erro=12;
  elseif s==3
    ini_int_erro=12; fim_int_erro=25;
  elseif s==4
    ini_int_erro=25; fim_int_erro=100;end
m=20; n=10;
qtd_itr=2000;                   
%qtd_itr se refere aos cenarios sem outliers. para cada um deles, serão 100 cenarios com outliers.
%assim, a qtd de iteracoes serah qtd_itr*100
cont_acerto=0;erro2=0; erro3=0;

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
d=[49 41 38 34 22 13 23 48 15 24 62 49 35 43 20 28 19 39 27 21];
dp=1.*sqrt(d);

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


for q=1:qtd_itr
    #Entrada L (rede perfeita) e simulaçoes em L
  L=[163854.9;6446.2;57037.0;126209.5;101128.6;296885.8;398014.4;60449.1;173710.4;167264.2;
110227.2;155928.2;52875.0;62904.2;3889.5;42705.7;98891.2;115779.2;113222.5;46428.8];
       
    for z=1:m
         do a=randn(1);
         until (a<=3)
       Laleat(z,1)=L(z)+dp(z)*a;
    end
    
    for k=1:100   
    Lgross=Laleat;outidt=[];
    j=randi([1 m]);
    choice=randi([-1 0]);
    if choice==0 
      choice=1;end
    Lgross(j)=Lgross(j)+choice*dp(j)*(ini_int_erro+(fim_int_erro-ini_int_erro)*rand(1));
    
    [xopt, fopt, erro, extra] = glpk (c, A1, Lgross, lb=[], ub=[], ctype, vartype, s=1, param);
    for i=1:m
      v(i)=abs(xopt(2*n+i)-xopt(2*n+i+m));
    end
    vSemZeros=v;
    vSemZeros(find(vSemZeros==0))=[];
    med=median(vSemZeros);
    mad=median(abs(vSemZeros-med));

    for i=1:m
      if fator==1
        Y(i)=v(i);
      elseif fator==2
        Y(i)=v(i)/med;
      elseif fator==3
        Y(i)=v(i)/mad;
        end
      
      if Y(i)>VC
        outidt=[outidt;i];end
    end

    if outidt==[j]
      cont_acerto=cont_acerto+1;
      end        
end
end
ini_int_erro
cont_acerto
end
