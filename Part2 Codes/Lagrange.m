clear;
%Initialize the parameter of the links
%mass, acceleration of gravity and the moment of inertia
syms m1 m2 g i1 i2
%The length of two links
syms l1 l2
%the angle of two links
syms t q1(t)
syms   q2(t)
%The center of mass of link 1 and 2
syms c1 c2

%The (x,y) position of Link 1
x1_c = c1*cos(q1(t));
y1_c = c1*sin(q1(t));
%The (x,y) position of Link 2
x2_c = l1*cos(q1(t))+c2*cos(q1(t)+q2(t));
y2_c = l1*sin(q1(t))+c2*sin(q1(t)+q2(t));

%The velocity and accelerated velocity of Link 1
%In the direction of X
x1_v = diff(x1_c,t);
x1_a = diff(x1_c,t, 2);
%In the direction of Y
y1_v = diff(y1_c,t);
y1_a = diff(y1_c,t, 2);

%The velocity and accelerated velocity of Link 2
%In the direction of X
x2_v = diff(x2_c,t);
x2_a = diff(x2_c,t, 2);
%In the direction of Y
y2_v = diff(y2_c,t);
y2_a = diff(y2_c,t, 2);

%Then calculate the system's kinetic energy
%The translational kinetic energy
KE_trans = 0.5*m1*(x1_v^2+y1_v^2) + 0.5*m2*(x2_v^2+y2_v^2);
%The rotation energy
KE_rot = 0.5*i1*((diff(q1,t))^2) + 0.5*i2*((diff(q1,t)+diff(q2,t))^2);
%The total energy
KE = KE_trans + KE_rot;

%Calculate the system's position energy (Only gravity)
PE = m1*g*y1_c + m2*g*y2_c;

%Operator of Lagrange
L_ORIGINAL = KE - PE;
symvar(L_ORIGINAL)
symvar(KE)

%Map the parameter
syms q_1 q1_d q_2 q2_d
L = subs(L_ORIGINAL, [q1(t) ,diff(q1(t),t), q2(t) ,diff(q2(t),t)], [q_1, q1_d, q_2, q2_d]);
K = subs(KE, [q1(t) ,diff(q1(t),t), q2(t) ,diff(q2(t),t)], [q_1, q1_d, q_2, q2_d]);

%∂K/∂q1_d and ∂L/∂q2_d
dKdq1d = diff(K, q1_d);
dKdq2d = diff(K, q2_d);

%Change the parameter back
dKdq1d = subs(dKdq1d, [q_1, q1_d, q_2, q2_d], [q1(t) ,diff(q1(t),t), q2(t) ,diff(q2(t),t)]);
dKdq2d = subs(dKdq2d, [q_1, q1_d, q_2, q2_d], [q1(t) ,diff(q1(t),t), q2(t) ,diff(q2(t),t)]);

%d/dt (∂K/∂q1_d)
der_of_dKdq1d = diff(dKdq1d, t);
der_of_dKdq2d = diff(dKdq2d, t);

%∂L/∂q1 and ∂L/∂q2
dLdq1 = diff(L, q_1);
dLdq2 = diff(L, q_2);
%Change the parameter back
dLdq1 = subs(dLdq1, [q_1, q1_d, q_2, q2_d], [q1(t) ,diff(q1(t),t), q2(t) ,diff(q2(t),t)]);
dLdq2 = subs(dLdq2, [q_1, q1_d, q_2, q2_d], [q1(t) ,diff(q1(t),t), q2(t) ,diff(q2(t),t)]);

%d/dt (∂L/∂q1_d) - ∂L/∂q1 = Q1
syms Q1
EoM_T1_LeftPart =  der_of_dKdq1d - dLdq1;
EoM_T1_LeftPart = simplify(EoM_T1_LeftPart);
EoM_T1_RightPart = Q1;
EoM_T1 = EoM_T1_LeftPart==EoM_T1_RightPart;
simplify(EoM_T1);
display(EoM_T1);

%d/dt (∂L/∂q2_d) - ∂L/∂q2 = Q2
syms Q2
EoM_T2_LeftPart =  der_of_dKdq2d - dLdq2;
EoM_T2_LeftPart = simplify(EoM_T2_LeftPart);
EoM_T2_RightPart = Q2;
EoM_T2 = EoM_T2_LeftPart==EoM_T2_RightPart;
simplify(EoM_T2);
display(EoM_T2);

%Separate the equation
%Map the parameter
syms q_1 q_2 q1d q1dd q2d q2dd
old_list1 = [q1, diff(q1,1), diff(q1,2)];
new_list1 = [q_1, q1d, q1dd];
old_list2 = [q2, diff(q2,1), diff(q2,2)];
new_list2 = [q_2, q2d, q2dd];
EoM_T1 = subs(EoM_T1, [old_list1, old_list2], [new_list1, new_list2]);
EoM_T2 = subs(EoM_T2, [old_list1, old_list2], [new_list1, new_list2]);
S = solve([EoM_T1, EoM_T2], [q1dd, q2dd]);
Result1 = S.q1dd;
Result2 = S.q2dd;



tf_i_should_create_SL_block = true;
if(true==tf_i_should_create_SL_block)
    MODEL_NAME = 'Lagr_Model';
    if(4==exist(MODEL_NAME))
        close_system(MODEL_NAME, 0);
        delete(MODEL_NAME);
    end
    new_system(MODEL_NAME)
    open_system(MODEL_NAME)
    
    INPUT_VAR_ORDER = { 'i1', 'i2', ...
                        'c1', 'c2', 'l1', 'g', 'm1', 'm2', ...
                        'q1d', 'q_1', ...
                        'q2d', 'q_2', ...
                        'Q1', 'Q2'}; 

    % Put BOTH equations into one block 
    matlabFunctionBlock( [MODEL_NAME,'/The_Lagrange_Model_for_Two_Links_Robot'], Result1, Result2, ...
                         'Optimize', false, ...
                         'Vars',     INPUT_VAR_ORDER, ...
                         'Outputs', {'q1dd_', 'q2dd_'}   );     
end


