function [h_vec,T_vec,p_vec,rho_vec,a_vec]=isa_prop(h_G_vec)

g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

h_G0_row=[0,11,25,47,53,79,90,105]*1e3;
T_0_row=[288.16,216.66,216.66,282.66,282.66,165.66,165.66];
p_0_row=[101330,22633,2488.7,120.45,58.323,1.0095,0.10444];
a_0_row=[-.0065,.003,-.0045,.004];

h_vec=r.*h_G_vec./(r+h_G_vec);

sz=size(h_G_vec);
T_vec=nan(sz);
p_vec=nan(sz);
for n=1:length(h_G_vec)
    if h_G_vec(n)>=h_G0_row(1) && h_G_vec(n)<=h_G0_row(2)
        T_vec(n)=T_0_row(1)+a_0_row(1).*(h_G_vec(n)-h_G0_row(1));
        p_vec(n)=p_0_row(1).*(T_vec(n)./T_0_row(1)).^(-g_0./a_0_row(1)/R);
    elseif h_G_vec(n)<=h_G0_row(3)
        T_vec(n)=T_0_row(2);
        p_vec(n)=p_0_row(2).*exp(-g_0.*(h_G_vec(n)-h_G0_row(2))./R./T_0_row(2));
    elseif h_G_vec(n)<=h_G0_row(4)
        T_vec(n)=T_0_row(3)+a_0_row(2).*(h_G_vec(n)-h_G0_row(3));
        p_vec(n)=p_0_row(3).*(T_vec(n)./T_0_row(3)).^(-g_0./a_0_row(2)/R);
    elseif h_G_vec(n)<=h_G0_row(5)
        T_vec(n)=T_0_row(4);
        p_vec(n)=p_0_row(4).*exp(-g_0.*(h_G_vec(n)-h_G0_row(4))./R./T_0_row(4));
    elseif h_G_vec(n)<=h_G0_row(6)
        T_vec(n)=T_0_row(5)+a_0_row(3).*(h_G_vec(n)-h_G0_row(5));
        p_vec(n)=p_0_row(5).*(T_vec(n)./T_0_row(5)).^(-g_0./a_0_row(3)/R);
    elseif h_G_vec(n)<=h_G0_row(7)
        T_vec(n)=T_0_row(6);
        p_vec(n)=p_0_row(6).*exp(-g_0.*(h_G_vec(n)-h_G0_row(6))./R./T_0_row(6));
    elseif h_G_vec(n)<=h_G0_row(8)
        T_vec(n)=T_0_row(7)+a_0_row(4).*(h_G_vec(n)-h_G0_row(7));
        p_vec(n)=p_0_row(7).*(T_vec(n)./T_0_row(7)).^(-g_0./a_0_row(4)/R);
    else
        error('Input heights are out of range!')
    end

end
rho_vec=p_vec./R./T_vec;
a_vec=sqrt(gamma.*R.*T_vec);