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

function matrix = round_matrix(ma,n)
    matrix = ma;
    matrix = round(matrix .* (10**n)) ./ (10**n);
    return
endfunction

%===============================================================
% DADOS ENUNCIADO
%===============================================================
vab = 11*1000;
len_lin = 1.6;
z1 = ((0.2) + (J*0.4))*len_lin;
z0 = ((0.6) + (J*1.2))*len_lin; 

%===============================================================
% PROBLEMA 1 - fase-curto
%===============================================================
z2 = z1;

z_sim = [[z0,0,0];[0,z1,0];[0,0,z2]]; 

van = vab/(make_complex(sqrt(3),30));
v_abc = van*[1;alfa^2;alfa];
v_123 = T\v_abc;
%disp("v_123 = ");
%printFasores(v_123)

i_1 = v_123(2)/(z_sim(1,1) + z_sim(2,2) + z_sim(3,3));
i_123_matrix = i_1*[1;1;1];
i_abc_matrix = T*i_123_matrix;
%disp("i_abc = ");
%printFasores(i_abc_matrix);
disp("1)");
printFasor(i_abc_matrix(1),"icc");

v_1_a = -z_sim (1,1)*i_1;
v_2_a = v_123(2)-z_sim (2,2)*i_1;
v_3_a = -z_sim (3,3)*i_1;
v_123_a = [v_1_a;v_2_a;v_3_a];
v_abc_a = T*v_123_a;
#disp("v_abc_a_matrix = ")
#printFasores(v_abc_a);

vb = vab/sqrt(3);
v_abc_a_pu = v_abc_a/vb;
disp("2)");
printFasor(v_abc_a_pu(2),"v_abc_a_pu");


%===============================================================
% PROBLEMA 2 - dupla-fase
%===============================================================
i1 = 0;
i2 = (v_123(2) - v_123(3))/(z_sim(2,2) + z_sim(3,3));
i3 = -i2;
i_123_matrix = [i1;i2;i3];
i_abc_matrix = T*i_123_matrix;
#printFasores(i_abc_matrix);
disp("3)");
printFasor(i_abc_matrix(2),"i_abc_matrix");
disp("\n");

%===============================================================
% PROBLEMA 3 - sem-mútuas
%===============================================================
z_sim = [[z0,0,0];[0,z1,0];[0,0,z2]];
z_matrix = T*z_sim*(T^(-1));
z_matrix(1,2) = z_matrix(1,3) = 0;
z_matrix(2,1) = z_matrix(2,3) = 0;
z_matrix(3,1) = z_matrix(3,2) = 0;
%disp(z_matrix);
z_sim = (T^(-1))*z_matrix*(T);

van = vab/(make_complex(sqrt(3),30));
v_abc = van*[1;alfa^2;alfa];
v_123 = T\v_abc;
%disp("v_123 = ");
%printFasores(v_123)

i_1 = v_123(2)/(z_sim(1,1) + z_sim(2,2) + z_sim(3,3));
i_123_matrix = i_1*[1;1;1];
i_abc_matrix = T*i_123_matrix;
%disp("i_abc = ");
%printFasores(i_abc_matrix);
disp("4)");
printFasor(i_abc_matrix(1),"icc");

v_1_a = -z_sim (1,1)*i_1;
v_2_a = v_123(2)-z_sim (2,2)*i_1;
v_3_a = -z_sim (3,3)*i_1;
v_123_a = [v_1_a;v_2_a;v_3_a];
v_abc_a = T*v_123_a;
#disp("v_abc_a_matrix = ")
#printFasores(v_abc_a);

vb = vab/sqrt(3);
v_abc_a_pu = v_abc_a/vb;
disp("5)");
printFasor(v_abc_a_pu(2),"v_abc_a_pu");

%===============================================================
% PROBLEMA 4 - dupla-fase-sem-mútua
%===============================================================
i1 = 0;
i2 = (v_123(2) - v_123(3))/(z_sim(2,2) + z_sim(3,3));
i3 = -i2;
i_123_matrix = [i1;i2;i3];
i_abc_matrix = T*i_123_matrix;
#printFasores(i_abc_matrix);
disp("6)");
printFasor(i_abc_matrix(2),"i_abc_matrix");
