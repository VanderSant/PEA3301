clear all
format long
% Vanderson da Silva dos Santos
% Email: vanderson.santos@usp.br 
% NUSP: 11259715

%%condicionais do exercício: Transformador ideal e Carga de impedância constante
%variaveis
Vg = 1000;     %[V] -tensão gerador
Vc = 220;      %[V] -tensão carga nominal
SAc = 20e3;    %[VA] -potência aparente carga
FPc = 0.85;    %-fator de potencia da carga
ZPc = 0.2 + J*0.6; %[ohms/km] - Impedância própria dos condutores
ZPm = J*0.2;   %[ohms/km] - Impedância mútua entre os condutores
VTp = 1000;    %[V] - Tensão nominal do primário do transformador
VTs = 220;     %[V] - Tensão nominal do secundário do transformador
LENp = 340e-3; %[km] - Comprimento do trecho primário
LENs = 45e-3;  %[m] - Comprimento do trecho secundário
f = 60; %[Hz] - frquência

%%bases
sb = 20e3;  %[VA] - Potência de base
%tensão
vb1 = 1000; %[V] - Tensão de base no primário
vb2 = 220;  %[V] - Tensão de base no secundário
%impedânica
zb1 = conj((vb1^2)/sb); %[ohms] - Impedância base 1
zb2 = conj((vb2^2)/sb); %[ohms] - Impedância base 2
%corrente
ib1 = conj(sb/vb1);
ib2 = conj(sb/vb2);

%%impedâncias
%condutores
ZP1 = ZPc*LENp; %[ohms] - Impedância própria dos condutores 1
ZM1 = ZPm*LENp; %[ohms] - Impedância mútua dos condutores 1
ZP2 = ZPc*LENs; %[ohms] - Impedância própria dos condutores 2
ZM2 = ZPm*LENs; %[ohms] - Impedância mútua dos condutores 2
%impedância Carga
modZc = abs(conj((Vc^2)/(SAc)));
Zc = modZc*FPc + J*modZc*sin(acos(FPc)); %[ohms] - Impedância da carga

%%valores em p.u.
Zm1_pu = 2*(ZP1-ZM1)/zb1; %impedânica da linha 1 em pu
Zm2_pu = 2*(ZP2-ZM2)/zb2; %impedânica da linha 2 em pu
Zt_pu = Zm1_pu + Zm2_pu;  %impedânica da linha toda em pu
Zc_pu = Zc/zb2; %impedância da carga em pu
Vg_pu = Vg/vb1; %tensão do gerador em pu 

%%%resultados
I_pu = Vg_pu/(Zt_pu+Zc_pu); %corrente total em pu
Vc_pu=Vg_pu-Zt_pu*I_pu; %tensão da carga em pu
Slin_pu = Vg_pu*conj(I_pu)-Vc_pu*conj(I_pu) ; %potencia nas linhas em pu
pot_ativa = 100*real(Slin_pu); %potência ativa nas linhas em pu
%%Adição de capacitor
fp_new = 0.92; %novo fator de potência
S_c = Vc_pu*vb2*conj(I_pu*ib2); %Potência carga
S_c_r_new = real(S_c)*tan(acos(fp_new)); %Nova potência carga
Q_c = imag(S_c) - S_c_r_new; %Diferença entre a potência reativa antiga e nova 
C = Q_c/(2*pi*f*(abs(Vc_pu)*vb2)^2); %capacitância do capacitor
%%Potência com novo capacitor
Zc_r = 1/(J*2*pi*f*C); %impedância do capacitor
Zc_eq = Zc_r*Zc/(Zc_r+Zc); %impedânica equivalente carga+capacitor
Zc_r_pu = Zc_eq/zb2; %impedânica equivalente carga+capacitor em pu
I_pu_new = Vg_pu/(Zt_pu+Zc_r_pu); %corrente nova em pu
pot_ativa_new=Vg_pu*conj(I_pu_new)-(Zc_r_pu*I_pu_new)*conj(I_pu_new);
%prints
disp(['Q1: V na carga =',num2str(abs(Vc_pu)),' [pu]']);
disp(['Q2: i gerador = ',num2str(abs(I_pu)),' [pu]']);
disp(['Q3: Perda circuito = ',num2str(abs(pot_ativa)),' [%]']);
disp(['Q4: Capacitância = ',num2str((C*1e+6)),' [uF]']);
disp(['Q5: Perda circuito = ',num2str(100*real(pot_ativa_new)),' [%]']);