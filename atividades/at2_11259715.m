clear all
format long
% Vanderson da Silva dos Santos
% Email: vanderson.santos@usp.br 
% NUSP: 11259715

%dados enuncido
len = 2500/1000; %[km]
zaa = (0.2 + J*0.6)*len;%[ohms]
zbb = (0.3 + J*0.5)*len;%[ohms]
zcc = (0.1 + J*0.4)*len;%[ohms]
zab = zba = J*0*len;  %[ohms]
zac = zca = J*0*len;  %[ohms]
zbc = zcb = J*0*len;  %[ohms]
vlnom = 9000; %[V] Tensão de linha nominal da carga
scnom = (2+J*1.5)*1000000; %[VA] Potência complexa trifásica nominal da carga

%alfa
alfa = (-0.5)+(J*sqrt(3)/2);

%impedância  e corrente carga
zcarga = conj(vlnom*vlnom/scnom);
ia = conj(scnom)/(vlnom*sqrt(3));

%matrizes
I = ia*[1;(alfa^2);alfa];
Zlin = [[zaa,zab,zac];[zba,zbb,zbc];[zac,zbc,zcc]];
Zc = [[zcarga,0,0];[0,zcarga,0];[0,0,zcarga]];

Vlin = Zlin*I; 
V = (Zlin+Zc)*I;
disp("Q1");
disp("tensão gerador");
disp(abs(V));

disp("Q2");
Slin = Vlin(1)*conj(I(1))+Vlin(2)*conj(I(2))+Vlin(3)*conj(I(3));
disp("perda linha");
disp(Slin);

%%linha é transposta
zp = (zaa+zbb+zcc)/3;
zm = (zab+zac+zbc)/3;

disp("Q3")
disp(real(zp));
disp(imag(zp));
%Zlin = [[zp,zm,zm];[zm,zp,zm];[zm,zm,zp]]
disp("Q4");
disp("as perdas nas linhas não mudam com a transposição da matriz");
disp("perda linha");
disp(Slin);





