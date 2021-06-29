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
grau_30 = exp(J*30*pi/180);
grau_menos_30 = exp(J*(-30)*pi/180);

%===============================================================
% DADOS ENUNCIADO
%===============================================================

% potência constante
% ligação delta estrela
vmt = 15*1000;
vbt = 220;
s_nom = 500*1000;
s_apa_nom_carga = 200*1000;
fp_carga = 0.9; % indutivo
zmt = (0.2 + J*(0.4))*2;
zbt = (0.2 + J*(0.4))*(20/1000);

vmt_base = 15*1000;
vbt_base = 220;

vg_pu = 1;

erro = 1/1000000;
%===============================================================
% PROBLEMA
%===============================================================
%% observações iniciais
% vbt = vbt*(sqrt(3));

phi_carga = acos(fp_carga);
s_carga = s_apa_nom_carga*(cos(phi_carga) + J*sin(phi_carga));

%% valores de base
s_base = s_nom;

zmt_base = conj((vmt_base**2)/(s_base));
zbt_base = conj((vbt_base**2)/(s_base));

imt_base = vmt_base/(zmt_base*sqrt(3));
ibt_base = vbt_base/(zbt_base*sqrt(3));

%% valores em pu

s_carga_pu = s_carga/s_base;

zmt_pu = zmt/zmt_base;
zbt_pu = zbt/zbt_base;

v_carga_pu = vg_pu;
ibt_pu = conj(s_carga_pu/v_carga_pu);
imt_pu = ibt_pu*grau_menos_30;

va = vg_pu - zmt_pu*imt_pu;
vb = va*grau_30;
new_v_carga_pu = vb - zbt_pu*ibt_pu;

while(abs(new_v_carga_pu-v_carga_pu) > erro)
    v_carga_pu = new_v_carga_pu;
    
    ibt_pu = conj(s_carga_pu/v_carga_pu);
    imt_pu = ibt_pu*grau_menos_30;
    
    va = vg_pu - zmt_pu*imt_pu;
    vb = va*grau_30;
    new_v_carga_pu = vb - zbt_pu*ibt_pu;
endwhile;

v_carga_pu = new_v_carga_pu;
printFasor(v_carga_pu,"v_carga_pu");

imt = imt_pu*imt_base;
printFasor(imt,"imt");

ibt = ibt_pu*ibt_base;
printFasor(ibt,"ibt");

s_carga_tri = real(v_carga_pu*conj(ibt_pu)*s_base/1000)

s_gerador_tri = real(vg_pu*conj(imt_pu)*s_base/1000)


%% respostas