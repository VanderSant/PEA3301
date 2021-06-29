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
ia = make_complex(15,15);
ib = make_complex(16, -120);
i2 = make_complex(2,-150);

%===============================================================
% PROBLEMA
%===============================================================
i1 = (ib - ia + i2*(1 - alfa))/((alfa^2)-1);
i0 = ia - i2 - i1;
ic  = i0 + (alfa*i1) + (alfa*alfa*i2);
in = ia + ib + ic;

printFasor(ia,"ia");
printFasor(ib,"ib");
printFasor(i2,"i2");
disp("\n");
printFasor(i1,"i1");
printFasor(i0,"i0");
printFasor(ic,"ic");
printFasor(in,"in");





               

