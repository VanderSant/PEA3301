clear all
format long
% Vanderson da Silva dos Santos
% Email: vanderson.santos@usp.br 
% NUSP: 11259715

%===============================================================
% CONSTANTES
%===============================================================
alfa = exp(J*120*pi/180); #operador alfa
grau_mais_30 = exp(J*30*pi/180); #fasor de modulo e fase 30°
grau_menos_30 = exp(J*(-30)*pi/180); #fasor de modulo e fase -30°
T = [[1, 1, 1];        #matrix T
     [1, alfa^2, alfa];
     [1, alfa, alfa^2]]; 
     
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

function z_sim = get_z_sim(za,zb,zc,zab,zac,zbc)
    alfa = exp(J*120*pi/180);
    T = [[1, 1, 1];
         [1, alfa^2, alfa];
         [1, alfa, alfa^2]];
    z_sim = [[za,zab,zac];[zab,zb,zbc];[zac,zbc,zc]];
    z_sim = inv(T)*z_sim*T;
    return
endfunction

function matrix = round_matrix(ma,n)
    matrix = ma;
    matrix = round(matrix .* (10**n)) ./ (10**n);
    return
endfunction

%===============================================================
% DADOS ENUNCIADO
%===============================================================


%===============================================================
% PROBLEMA
%===============================================================