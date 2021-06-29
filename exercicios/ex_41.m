clear all
format long
% Vanderson da Silva dos Santos
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
alfa = exp(J*120*pi/180); #operador alfa
grau_mais_30 = exp(J*30*pi/180); #fasor de modulo e fase 30°
grau_menos_30 = exp(J*(-30)*pi/180); #fasor de modulo e fase -30°
T = [1, 1, 1;        #matrix T
     1, alfa^2, alfa;
     1, alfa, alfa^2]; 

%===============================================================
% DADOS ENUNCIADO
%===============================================================
f = 60;
w = 2*pi*f;

vab = 220;
c = 0.15e-6;
za = 1/(J*w*c);
zb = zc = 16e3;

%===============================================================
% PROBLEMA
%===============================================================
vn = vab/make_complex(sqrt(3),30);

vn_array = vn*[1;
              (alfa^2);
               alfa];
               
vnn = tensao_neutro(vn_array(1),vn_array(2),vn_array(3),za,zb,zc,inf);

printFasores(vnn);

z_array = [za,0,0;
           0,zb,0;
           0,0,zc];

i_array = z_array\vn_array;
printFasores(i_array);













