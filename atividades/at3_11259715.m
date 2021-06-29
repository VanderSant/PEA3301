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
alfa = exp(J*120*pi/180);

%===============================================================
% DADOS ENUNCIADO
%===============================================================
iab = 6 + 15 + 12;
ibc = 6 + 15*(alfa**2) + 12*(alfa);
ica = 6 + 15*alfa + 12*(alfa**2);

ia = iab - ica;
ib = ibc - iab;
ic = ica - ibc;
abs(ia)
abs(ib)
abs(ic)