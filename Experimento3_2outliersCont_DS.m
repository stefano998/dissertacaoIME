warning ("off");
W=3.29;
if 1
mem=[1	1	2
2	1	3
3	1	6
4	1	9
5	1	10
6	1	11
7	1	12
8	1	13
9	1	14
10	1	15
11	1	16
12	1	17
13	1	18
14	1	19
15	1	20
16	2	3
17	2	4
18	2	5
19	2	6
20	2	7
21	2	8
22	2	9
23	2	10
24	2	11
25	2	12
26	2	13
27	2	14
28	2	15
29	2	16
30	2	17
31	2	18
32	2	19
33	2	20
34	3	4
35	3	5
36	3	6
37	3	7
38	3	8
39	3	9
40	3	10
41	3	11
42	3	12
43	3	13
44	3	14
45	3	15
46	3	16
47	3	17
48	3	18
49	3	19
50	3	20
51	4	7
52	4	9
53	4	10
54	4	11
55	4	12
56	4	13
57	4	14
58	4	15
59	4	16
60	4	17
61	4	18
62	4	19
63	4	20
64	5	8
65	5	9
66	5	10
67	5	11
68	5	12
69	5	13
70	5	14
71	5	15
72	5	16
73	5	17
74	5	18
75	5	19
76	5	20
77	6	9
78	6	10
79	6	11
80	6	12
81	6	13
82	6	14
83	6	15
84	6	16
85	6	17
86	6	18
87	6	19
88	6	20
89	7	8
90	7	9
91	7	10
92	7	11
93	7	12
94	7	13
95	7	14
96	7	15
97	7	16
98	7	17
99	7	18
100	7	19
101	7	20
102	8	9
103	8	10
104	8	11
105	8	12
106	8	13
107	8	14
108	8	15
109	8	16
110	8	17
111	8	18
112	8	19
113	8	20
114	9	10
115	9	11
116	9	12
117	9	13
118	9	14
119	9	15
120	9	16
121	9	17
122	9	18
123	9	19
124	9	20
125	10	11
126	10	12
127	10	13
128	10	14
129	10	15
130	10	16
131	10	17
132	10	18
133	10	19
134	10	20
135	11	12
136	11	13
137	11	14
138	11	15
139	11	16
140	11	17
141	11	18
142	11	19
143	11	20
144	12	13
145	12	14
146	12	15
147	12	18
148	12	19
149	12	20
150	13	16
151	13	17
152	13	19
153	14	16
154	14	17
155	14	19
156	15	17
157	15	20
158	16	18
159	16	20
160	17	18
161	17	19
162	17	20
163	18	19
164	18	20
165	19	20];
  end
for s=1:4
    if s==1
    ini_int_erro=3; fim_int_erro=6;
  elseif s==2
    ini_int_erro=6; fim_int_erro=12;
  elseif s==3
    ini_int_erro=12; fim_int_erro=25;
  elseif s==4
    ini_int_erro=25; fim_int_erro=100;end
rand("state",[2]);randn("state",[2]);
qtd_itr=2000; cont_acerto=0;
m1=20; n=10;m=m1;
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

d=[49 41 38 34 22 13 23 48 15 24 62 49 35 43 20 28 19 39 27 21];
dp=1.*sqrt(d);p=dp.^(-2);
P=diag(p);

for q=1:qtd_itr
    #Entrada L (rede perfeita) e simula�oes em L
L=[163854.9;6446.2;57037.0;126209.5;101128.6;296885.8;398014.4;60449.1;173710.4;167264.2;
110227.2;155928.2;52875.0;62904.2;3889.5;42705.7;98891.2;115779.2;113222.5;46428.8];
    m=m1;   
    for z=1:m
         do a=randn(1);
         until (a<=3)
       Laleat(z,1)=L(z)+dp(z)*a;
    end

    
    for k=1:100  
    Lgross=Laleat;m=m1;   
    
    f=randi([1 165]);
    j=mem(f,2); p=mem(f,3);
    choice=randi([-1 0]);
    if choice==0 
      choice=1;end
    Lgross(j)=Lgross(j)+choice*dp(j)*(ini_int_erro+(fim_int_erro-ini_int_erro)*rand(1));
    
    choice=randi([-1 0]);
    if choice==0 
      choice=1;end
    Lgross(p)=Lgross(p)+choice*dp(p)*(ini_int_erro+(fim_int_erro-ini_int_erro)*rand(1));
    
    x=inv(A'*P*A)*A'*P*Lgross;
    v=A*x-Lgross;
    Ev=inv(P)-A*inv(A'*P*A)*A';
    ww=zeros(m, 1);
    for i=1:m
      ww(i)=abs(v(i)/sqrt(abs(Ev(i,i))));end;
    if max(ww)>W
      if max(ww)==ww(j) | max(ww)==ww(p)
          if max(ww)==ww(j)
            out1=j;out2=p; 
            if p>j
              out2=out2-1;end
          else out1=p; out2=j;
            if j>p
              out2=out2-1; end
          end
                
        m=m-1;
        A2=A;Lgross2=Lgross;P2=P;
        A2(out1,:)=[];Lgross2(out1,:)=[];P2(out1,:)=[];P2(:,out1)=[];
        x=inv(A2'*P2*A2)*A2'*P2*Lgross2;
        v=A2*x-Lgross2;
        Ev=inv(P2)-A2*inv(A2'*P2*A2)*A2';
        ww=zeros(m, 1);
        for i=1:m
          ww(i)=abs(v(i)/sqrt(abs(Ev(i,i))));end;
        if max(ww)>W
         if max(ww)==ww(out2)
       
        m=m-1;
        A3=A2;Lgross3=Lgross2;P3=P2;
        A3(out2,:)=[];Lgross3(out2,:)=[];P3(out2,:)=[];P3(:,out2)=[];
        x=inv(A3'*P3*A3)*A3'*P3*Lgross3;
        v=A3*x-Lgross3;
        Ev=inv(P3)-A3*inv(A3'*P3*A3)*A3';
        ww=zeros(m, 1);
        for i=1:m
          ww(i)=abs(v(i)/sqrt(abs(Ev(i,i))));end;   
          
        if max(ww)<W
          cont_acerto=cont_acerto+1;end;
end
end
end
end
    
end
end
ini_int_erro
cont_acerto
end
