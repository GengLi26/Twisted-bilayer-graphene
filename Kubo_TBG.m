%% Tight bonding method calculation: band structure of bilayer graphene
% Ref: Rep. Prog. Phys. 76 (2013) 056503 (28pp)
clear all

%% calculate moire lattice constant
theta=1.05; % magic angle (degree)
theta=theta*pi/180; % magic angle (rad)
a_mono=2.46e-10; % monolayer graphene lattice constance (angstrom)
a_moire=a_mono/(2*sin(theta/2)); % moire lattice constant(angstrom)

e = 1.602177e-19; % electron charge (C)
hbar = 1.054572e-34; % reduced Planck constant ( J . s )

Kx=2*pi/(a_moire*sqrt(3)); %% x axis of K points
points=20; % avoid kx=0; ky =0;  20 points works for eta=0.5e-3*e; 

Ky=Kx/sqrt(3); % Y axis of K points

Y_max=Kx/(sqrt(3)/2);  %% ky for max FBZ
Ky_array=linspace(-(Kx/200),(Kx/200),points);
dky=abs(Ky_array(2)-Ky_array(1));
Kx_array=Ky_array;
dkx=dky;

Ef=0*e;% Fermi level of graphene unit /s
E=linspace(0.001,0.040,400).*e; % Energy unit J % 350 
omega=E./hbar; % omega unit s^-1

sigma0=e^2/4/hbar; % e^2/4/hbar, 
spin=2; % Spin degeneracy 2
valey=1; % Graphene at Gamma points has valley degeneracy of 1;
scale=spin*valey/(2*pi)^2;  % here I figured out the unit problem. 

eta=0.5e-3*e;  % phenomenological broadening

alpha=1:1:10;  % how many bands?  Ten band model 

% Generate k points in FBZ
for omega_count=1:1:length(omega)
    tic % timer
    for alpha_count=1:1:length(alpha)
        alpha_loop=alpha(alpha_count);
        beta=1:1:length(alpha);
        beta(alpha_loop)=[];
        
        for beta_count=1:1:length(beta)
            beta_loop=beta(beta_count);
            
            %%Integration over a rectangular enclosing Gamma point
            for count_ky=1:1:(length(Ky_array)) %%
                ky=Ky_array(count_ky);

                for count_kx=1:1:length(Kx_array)
                    
                    kx=Kx_array(count_kx);
                    Haml_10_A=Hamil_10(kx,ky,a_moire).*e*1e-3; % conver unit from meV into J 
                    Haml_10_B=Hamil_10(kx+Kx/10000,ky,a_moire).*e*1e-3; % conver unit  from meV  into J 
                    [Vector,E_eig]=eig(Haml_10_A);
                    Band=real(diag(E_eig));
                    vx=(Haml_10_B-Haml_10_A)./(Kx/10000)./hbar; % Don't use kx/10000, or NAN results appears
                    
                    Mean_alpha_beta=Vector(:,alpha_loop)'*vx*Vector(:,beta_loop);
                    Mean_beta_alpha=Vector(:,beta_loop)'*vx*Vector(:,alpha_loop);
                    % S not imputed for simplicity
                    cond_kx(count_kx)=e^2*hbar/i*(fermi(Band(alpha_loop),Ef)-fermi(Band(beta_loop),Ef))...
                        /(Band(alpha_loop)-Band(beta_loop))...
                        *(Mean_alpha_beta*Mean_beta_alpha)...
                        ./(Band(alpha_loop)-Band(beta_loop)+hbar*omega(omega_count)+i*eta);
                   
                end
                cond_kx_ky(count_ky)=sum(cond_kx).*dkx;
            end
            cond_beta(beta_count)=sum(cond_kx_ky).*dky;
        end
        cond_alpha_beta(alpha_count)=sum(cond_beta);
    end
    toc % timer 
    cond_omega(omega_count)=sum(cond_alpha_beta);
    
    omega_count./length(omega)*100
end


%% Plot band structure
% subplot(2,1,1)
plot(E./e,real(cond_omega)./sigma0*scale,'color','k','linewidth',2)
hold on 
% subplot(2,1,2)
% plot(E./e,imag(cond_omega)./sigma0*scale,'color','k','linewidth',2)
% hold on 


