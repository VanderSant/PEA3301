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
Van = 127;
Z_ = (8 + J*4);
Zl_ = (1 + j*2);

%===============================================================
% EXERCÍCIO
%===============================================================
disp("a)")
Va = Van * 1;
Vc = Van * alfa;
aux = (1/(Zl_ + Z_))+(1/Zl_);
Ve = tensao_neutro(Va,0,Vc,Z_+Zl_,Inf,Z_+Zl_,Zl_);
disp(Ve)
ia = (Va-Ve)/(Zl_ + Z_);
ic = (Vc-Ve)/(Zl_ + Z_);
in = Ve/Zl_;
printFasor(ia,"ia");
printFasor(ic,"ic");
printFasor(in,"in");

disp("\nb)")
disp("tensões de fase")
Vfac = ia*Z_;
printFasor(Vfac,"Vfac");
Vfcc = ic*Z_;
printFasor(Vfcc,"Vfac");
disp("tensões de linha")
VabC = Vfac - 0;
VbcC = 0 - Vfcc;
VcaC = Vfcc - Vfac;
printFasor(VabC,"Va'b'");
printFasor(VbcC,"Vb'c'");
printFasor(VcaC,"Vc'a'");

disp("\nc)")
s_c = Vfac*conj(ia) + Vfcc*conj(ic);
printFasores(s_c);
s_g = (Va-Ve)*conj(ia) + (Vc-Ve)*conj(ic) + Ve*conj(in);
printFasores(s_g);

disp("\nd)")
Vab = Van*make_complex(sqrt(3),30);
Vbc = Van*(alfa**2)*make_complex(sqrt(3),30);
Vca = Van*alfa*make_complex(sqrt(3),30);
w1 = (Vab)*conj(ia);
w2 = (-Vbc)*conj(ic);
w3 = (VabC)*conj(ia);
w4 = (-VbcC)*conj(ic);
real(w1)
real(w2)
real(w3)
real(w4)


