clear all
format long
% Vanderson da Silva dos Santos
% Email: vanderson.santos@usp.br 
% NUSP: 11259715

%===============================================================
% FUNÇÕES
%===============================================================
function printFasores(V)
    for i = 1:length(V)
        disp([num2str(abs(V(i))), ' \ ', num2str(rad2deg(angle(V(i)))),'°'])
    end 
end

function y = printFasor(x, nome)
    disp([nome, ' = ', num2str(abs(x)), ' \ ', num2str(rad2deg(angle(x))), '°']);
end

function fasor = make_complex(absol,angle)
    e = 2.71828182845;
    angle_pi = angle*pi/180;
    fasor = absol*(e^(1J*angle_pi));
    return
endfunction

function admi = admitancia(re)
    if (re == 0)
        admi = inf;
        return
    else
        admi = 1/re;
        return
    endif
endfunction

function vn_n = tensao_neutro(va,vb,vc,ra,rb,rc,rn)
    #só é valida sem mútuas
    ya = admitancia(ra);
    yb = admitancia(rb);
    yc = admitancia(rc);
    yn = admitancia(rn);
    vn_n = ((va*ya)+(vb*yb)+(vc*yc))/(ya+yb+yc+yn);
    return
endfunction

function imp = triangulo_para_estrela(zab,zbc,zca)
    imp(1) = (zab*zca)/(zab+zbc+zca);
    imp(2) = (zab*zbc)/(zab+zbc+zca);
    imp(3) = (zbc*zca)/(zab+zbc+zca);
    return
endfunction

%===============================================================
% CONSTANTES
%===============================================================
alfa = exp(J*120*pi/180);

%===============================================================
% DADOS ENUNCIADO
%===============================================================
vab = 500; % [V] (Tensao de linha)
za = 34;
zb = 24;
zc = 29;

%===============================================================
% PROBLEMA 1
%===============================================================

disp("Problema 1");
van = vab/make_complex(sqrt(3),30);
vbn = alfa^2*van;
vcn = alfa*van;
vnn = tensao_neutro(van,vbn,vcn,za,zb,zc,inf);
ia = (van - vnn)/za;
ib = (vbn - vnn)/zb;
ic = (vcn - vnn)/zc;
%% Leituras dos wattimetros
w2 = real((van - vbn)*conj(ia))
w1 = real((vcn - vbn)*conj(ic))
w1+w2

%===============================================================
% PROBLEMA 2
%===============================================================

disp("Problema 2");
w2_inv = imag((van - vbn)*conj(ia))
w1_inv = imag((vcn - vbn)*conj(ic))

%===============================================================
% PROBLEMA 3
%===============================================================

disp("Problema 3");
vnn = tensao_neutro(van,0,vcn,za,inf,zc,inf);
ia = (van - vnn)/za;
ic = (vcn - vnn)/zc;
%% Leituras dos wattimetros
w2 = real((van - vbn)*conj(ia))
w1 = real((vcn - vbn)*conj(ic))
s = (van-vnn)*conj(ia) + (vcn-vnn)*conj(ic)
blondel = w1+w2
