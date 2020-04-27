%% 21. Faithful tight-binding models and fragile topology of magic-angle bilayer graphene.pdf
% Ten band model 
clear all
clc

%% calculate moire lattice constant
theta=1.05; % magic angle (degree)
theta=theta*pi/180; % magic angle (rad)
a_mono=2.46; % monolayer graphene lattice constance (angstrom)
a_moire=a_mono/(2*sin(theta/2)); % moire lattice constant(angstrom)


%% K line scan 

Kx=2*pi/(a_moire*sqrt(3));
points=1000; % k mesh points

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

    % get eigenvalue
    Band(count,:)=sort(real(eig(Hamil_10(kx,ky,a_moire))),'ascend');
    count=count+1;
end

%% Plot band structure

plot(Kx_plot,Band,'color','k','linewidth',2)
hold on 
set(gca,'fontsize',28)
set(gca,'XTick',[])
ylabel(['E (meV)'],'FontSize',28)
set(gcf,'Position',[500 300 800 600])
xlim([0,max(Kx_plot)])
grid on 
grid minor
y_l=25;
ylim([-y_l,y_l])

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

% save figure as 'pdf'
set(gcf,'PaperOrientation','landscape')
print(gcf, 'PR 10band delta1_mag.pdf', '-dpdf','-r0','-bestfit')

% set(gcf,'PaperOrientation','landscape')
% print(gcf, 'PR 10band delta0.pdf', '-dpdf','-r0','-bestfit')

