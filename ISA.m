clearvars
clc
close all

g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

h_G0_row=[0,11,25,47,53,79,90,105]*1e3;
T_0_row=[288.16,216.66,216.66,282.66,282.66,165.66,165.66];
p_0_row=[101330,22633,2488.7,120.45,58.323,1.0095,0.10444];
a_0_row=[-.0065,.003,-.0045,.004];

N_layer=50;
h_G1_col=linspace(h_G0_row(1),h_G0_row(2),N_layer).';
h_G2_col=linspace(h_G0_row(2),h_G0_row(3),N_layer).';
h_G3_col=linspace(h_G0_row(3),h_G0_row(4),N_layer).';
h_G4_col=linspace(h_G0_row(4),h_G0_row(5),N_layer).';
h_G5_col=linspace(h_G0_row(5),h_G0_row(6),N_layer).';
h_G6_col=linspace(h_G0_row(6),h_G0_row(7),N_layer).';
h_G7_col=linspace(h_G0_row(7),h_G0_row(8),N_layer).';

T1_col=T_0_row(1)+a_0_row(1).*(h_G1_col-h_G0_row(1));
T2_col=repmat(T_0_row(2),N_layer,1);
T3_col=T_0_row(3)+a_0_row(2).*(h_G3_col-h_G0_row(3));
T4_col=repmat(T_0_row(4),N_layer,1);
T5_col=T_0_row(5)+a_0_row(3).*(h_G5_col-h_G0_row(5));
T6_col=repmat(T_0_row(6),N_layer,1);
T7_col=T_0_row(7)+a_0_row(4).*(h_G7_col-h_G0_row(7));

p1_col=p_0_row(1).*(T1_col./T_0_row(1)).^(-g_0./a_0_row(1)/R);
p2_col=p_0_row(2).*exp(-g_0.*(h_G2_col-h_G0_row(2))./R./T_0_row(2));
p3_col=p_0_row(3).*(T3_col./T_0_row(3)).^(-g_0./a_0_row(2)/R);
p4_col=p_0_row(4).*exp(-g_0.*(h_G4_col-h_G0_row(4))./R./T_0_row(4));
p5_col=p_0_row(5).*(T5_col./T_0_row(5)).^(-g_0./a_0_row(3)/R);
p6_col=p_0_row(6).*exp(-g_0.*(h_G6_col-h_G0_row(6))./R./T_0_row(6));
p7_col=p_0_row(7).*(T7_col./T_0_row(7)).^(-g_0./a_0_row(4)/R);

h_G_col=[h_G1_col;h_G2_col;h_G3_col;h_G4_col;h_G5_col;h_G6_col;h_G7_col];
T_col=[T1_col;T2_col;T3_col;T4_col;T5_col;T6_col;T7_col];
p_col=[p1_col;p2_col;p3_col;p4_col;p5_col;p6_col;p7_col];

h_col=r.*h_G_col./(r+h_G_col);
rho_col=p_col./R./T_col;
a_col=sqrt(gamma.*R.*T_col);

%% plot
plot(h_col./1000,h_G_col./1000)
xlabel('h (km)')
ylabel('h_G (km)')

figure
tiledlayout(1,2)
nexttile
plot(T_col-273,h_G_col./1000)
xlabel('T (C)')
ylabel('h_G (km)')

nexttile
plot(a_col./1000.*60.*60,h_G_col./1000)
xlabel('a (km/h)')
ylabel('h_G (km)')


figure
tiledlayout(1,2)
nexttile
plot(p_col./1e5,h_G_col./1000)
xlabel('p (bar)')
ylabel('h_G (km)')


nexttile
plot(rho_col,h_G_col./1000)
xlabel('rho (kg/m^{3})')
ylabel('h_G (km)')