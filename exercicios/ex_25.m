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
    angle_pi = angle*pi/180
    fasor = absol*(e^(1J*angle_pi));
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
disp("\na)");
Vf = Van*[[1];
      [alfa**2];
      [alfa]];
Vl =  make_complex(sqrt(3),30)*Vf
printFasores(Vl);

disp("\nb)");
If = Vf/(Z_ + Zl_) 
printFasores(If);
disp("corrente de fase e linha são iguais geradores simetricos")

disp("\nc)");
Va_n_f = Z_*If
Va_n_l = Va_n_f*(make_complex(sqrt(3),30)) 
disp("tensão de fase")
printFasores(Va_n_f);
disp("tensão de linha")
printFasores(Va_n_l);

disp("\nd)");
disp("como o sistema é simétrico, n'=0 e como consequência In=0")

disp("\ne)");
pot_g =  (Vf(1)*conj(If(1))) + (Vf(2)*conj(If(2))) + (Vf(3)*conj(If(3)));
pot_c =  (Va_n_f(1)*conj(If(1))) + (Va_n_f(2)*conj(If(2))) + (Va_n_f(3)*conj(If(3)));
printFasor(pot_g,"pot gerador");
printFasor(pot_c,"pot carga");

disp("\nf)")



