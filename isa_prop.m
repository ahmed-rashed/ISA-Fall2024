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
    if h_G_vec(n)<h_G0_row(1) || h_G_vec(n)>h_G0_row(end)
        warning("Invalid value in the input. Values must lie in the region ["+h_G0_row(1)+","+h_G0_row(end)+"]");
    end

    for n_layer=1:7
        if h_G_vec(n)<=h_G0_row(n_layer+1)
            if rem(n_layer,2)~=0 % n_layer is odd
                n_a=(n_layer+1)/2;
                T_vec(n)=T_0_row(n_layer)+a_0_row(n_a).*(h_G_vec(n)-h_G0_row(n_layer));
                p_vec(n)=p_0_row(n_layer).*(T_vec(n)./T_0_row(n_layer)).^(-g_0./a_0_row(n_a)/R);
            else % n_layer is even
                T_vec(n)=T_0_row(n_layer);
                p_vec(n)=p_0_row(n_layer).*exp(-g_0.*(h_G_vec(n)-h_G0_row(n_layer))./R./T_0_row(n_layer));
            end
            break
        end
    end
end
rho_vec=p_vec./R./T_vec;
a_vec=sqrt(gamma.*R.*T_vec);