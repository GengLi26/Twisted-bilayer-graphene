%% 21. Faithful tight-binding models and fragile topology of magic-angle bilayer graphene.pdf
% Six band model 
clear all
clc

%% calculate moire lattice constant
theta=1.05; % magic angle (degree)
theta=theta*pi/180; % magic angle (rad)
a_mono=2.46; % monolayer graphene lattice constance (angstrom)
a_moire=a_mono/(2*sin(theta/2)); % moire lattice constant(angstrom)

% Hopping parameters % meV
t_kappa=27; 
delta_pz=0*t_kappa;
delta_p_pm=-0.23*t_kappa;
delta_kappa=0.25*t_kappa;
t_pz=0.17*t_kappa;
t_p_pm=-0.017*t_kappa;
t_pp_pos=-0.065*t_kappa;
t_pp_neg=-0.055*t_kappa;
t_kappa_prime=0.25*t_kappa;
t_ppz_pos=0.095*t_kappa;
t_ppz_neg=0.055*t_kappa;
t_kappap_pm_pos=0.6*t_kappa;
t_kappap_pm_neg=0.2*t_kappa;

u_pz=-6*t_pz+delta_pz;
u_p_pm=3*t_p_pm+delta_p_pm;
u_kappa=-4*(t_kappa+t_kappa_prime)+delta_kappa;

omega=exp(i*2*pi/3);



Kx=2*pi/(a_moire*sqrt(3)); % kx of K points.
points=1000; % mesh points


% From K point to -K point
k_K_nK_x=[Kx:-Kx/points:-Kx];
k_K_nK_y=k_K_nK_x./sqrt(3);
% From -K point to -M point
k_nK_nM_x=-Kx.*ones(1,points+1);
k_nK_nM_y=[-Kx/sqrt(3):Kx/sqrt(3)/points:0];
% From -M point to M point
k_nM_M_x=[-Kx:Kx/points:Kx];
k_nM_M_y=k_nM_M_x.*0;
% From M point to K point
k_M_K_y=[0:Kx/sqrt(3)/points:Kx/sqrt(3)];
k_M_K_x=Kx.*ones(1,length(k_M_K_y));
% Generate kx and ky array 
Kx_array=[k_K_nK_x k_nK_nM_x  k_nM_M_x k_M_K_x];
Ky_array=[k_K_nK_y k_nK_nM_y k_nM_M_y k_M_K_y];
% plot kx and ky for checking the scan trace of k vector. 
%plot(Kx_array,Ky_array)

% x axis for ploting band structure
Kx_plot=[linspace(0,4,length(k_K_nK_x)) ...
         linspace(4,5,length(k_nK_nM_x))...
         linspace(5,5+2*sqrt(3),length(k_nM_M_x))...
         linspace(5+2*sqrt(3),6+2*sqrt(3),length(k_M_K_y))];

% plot(Kx_array,Ky_array)


count=1;

for kx=Kx_array
    %%
    ky=Ky_array(count);
    
    %% Construct 6*6 Hamiltonian
    p011110hc=(b6(kx,ky,0,1,a_moire)+b6(kx,ky,1,1,a_moire)+b6(kx,ky,1,0,a_moire)...
        +b6(kx,ky,0,1,a_moire)'+b6(kx,ky,1,1,a_moire)'+b6(kx,ky,1,0,a_moire)');
    
    H_pz=t_pz*p011110hc;


    C_pp=t_pp_pos*(b6(kx,ky,0,1,a_moire)+b6(kx,ky,-1,-1,a_moire)*omega+b6(kx,ky,1,0,a_moire)*omega')+...
        t_pp_neg*(b6(kx,ky,0,-1,a_moire)+b6(kx,ky,1,1,a_moire)*omega+b6(kx,ky,-1,0,a_moire)*omega');

    H_p_pm=t_p_pm*p011110hc*eye(2)+[0 C_pp';C_pp 0];

    C_ppz=i*t_ppz_pos*...
        [b6(kx,ky,0,1,a_moire)+b6(kx,ky,-1,-1,a_moire)*omega+b6(kx,ky,1,0,a_moire)*omega';...
        -(b6(kx,ky,0,-1,a_moire)+b6(kx,ky,1,1,a_moire)*omega'+b6(kx,ky,-1,0,a_moire)*omega)]...
        -i*t_ppz_neg*...
        [b6(kx,ky,0,-1,a_moire)+b6(kx,ky,1,1,a_moire)*omega+b6(kx,ky,-1,0,a_moire)*omega';...
        -(b6(kx,ky,0,1,a_moire)+b6(kx,ky,-1,-1,a_moire)*omega'+b6(kx,ky,1,0,a_moire)*omega)];

    C_kappap_pm=t_kappap_pm_pos*...
        [b6(kx,ky,-1,0,a_moire) b6(kx,ky,-1,-1,a_moire);...
        b6(kx,ky,-1,-1,a_moire)*omega' omega;...
        omega b6(kx,ky,-1,0,a_moire)*omega']-...
        t_kappap_pm_neg*...
        [b6(kx,ky,-1,-1,a_moire) b6(kx,ky,-1,0,a_moire);...
        omega' b6(kx,ky,-1,-1,a_moire)*omega;...
        b6(kx,ky,-1,0,a_moire)*omega omega'];

    H_kappa_1=t_kappa*...
        [0 b6(kx,ky,-1,0,a_moire) 1;...
        1 0 b6(kx,ky,0,-1,a_moire);...
        b6(kx,ky,1,1,a_moire) 1 0];
    H_kappa_2=t_kappa_prime*...
        [0 b6(kx,ky,-1,-1,a_moire) b6(kx,ky,-1,0,a_moire);...
        b6(kx,ky,0,-1,a_moire) 0 b6(kx,ky,1,0,a_moire);...
        b6(kx,ky,0,1,a_moire) b6(kx,ky,1,1,a_moire) 0];

    H_kappa=H_kappa_1+H_kappa_2+H_kappa_1'+H_kappa_2';

    H_k6=[H_pz+u_pz C_ppz' zeros(1,3);
        C_ppz H_p_pm+eye(2)*u_p_pm C_kappap_pm';
        zeros(3,1)  C_kappap_pm H_kappa+u_kappa*eye(3)];


     % get eigenvalue
     Band(count,:)=sort(real(eig(H_k6)),'ascend');
     count=count+1;
end

plot(Kx_plot,Band,'color','k','linewidth',2)
hold on 
set(gca,'fontsize',28)
set(gca,'XTick',[])
ylabel(['E (meV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])
xlim([0,max(Kx_plot)])
y_l=250;
ylim([-y_l,50])

% Gamma points
x=2.*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
% -K points
x=4.*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
hold on 
% -M points
x=5.*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
% Gamma points
x=(5+1*sqrt(3)).*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 
% M points
x=(5+2*sqrt(3)).*ones(1,points+1);
y=[-y_l:2*y_l./points:y_l];
plot(x,y,'LineWidth',1.5,'color',[0.2 0.2 0.2],'linestyle','--')
hold on 


% Create textbox
annotation('textbox',...
    [0.12 0.06 0.03 0.04],'String','K',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.28 0.06 0.03 0.04],'String','\Gamma',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');


annotation('textbox',...
    [0.43 0.06 0.03 0.04],'String','-K',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');


annotation('textbox',...
    [0.51 0.06 0.03 0.04],'String','-M',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.67 0.06 0.03 0.04],'String','\Gamma',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.8 0.06 0.03 0.04],'String','M',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.89 0.06 0.03 0.04],'String','K',...
    'LineWidth',4,...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');

set(gcf,'PaperOrientation','landscape')
print(gcf, 'PR 6band.pdf', '-dpdf','-r0','-bestfit')

