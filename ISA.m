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
h_G_col=nan(7*N_layer,1);
T_col=nan(7*N_layer,1);
p_col=nan(7*N_layer,1);
for n_layer=1:7
    ind_vec=(n_layer-1)*N_layer+(1:N_layer);
    h_G_col(ind_vec)=linspace(h_G0_row(n_layer),h_G0_row(n_layer+1),N_layer).';
    if rem(n_layer,2)~=0 % n_layer is odd
        n_a=(n_layer+1)/2;
        T_col(ind_vec)=T_0_row(n_layer)+a_0_row(n_a).*(h_G_col(ind_vec)-h_G0_row(n_layer));
        p_col(ind_vec)=p_0_row(n_layer).*(T_col(ind_vec)./T_0_row(n_layer)).^(-g_0./a_0_row(n_a)/R);
    else % n_layer is even
        T_col(ind_vec)=repmat(T_0_row(n_layer),N_layer,1);
        p_col(ind_vec)=p_0_row(n_layer).*exp(-g_0.*(h_G_col(ind_vec)-h_G0_row(n_layer))./R./T_0_row(n_layer));
    end
end

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