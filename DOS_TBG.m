%% DOS of monolayer graphene 

clc
clear all

%% calculate moire lattice constant
theta=1.05; % magic angle (degree)
theta=theta*pi/180; % magic angle (rad)
a_mono=2.46; % monolayer graphene lattice constance (angstrom)
a_moire=a_mono/(2*sin(theta/2)); % moire lattice constant(angstrom)

%%

sigma=0.002; %% Gaussian broadening 

Kx=2*pi/(a_moire*sqrt(3)); % kx of K point


Ky=Kx/sqrt(3); 
points=300; % 300 is OK for sigma=0.002;

ky=linspace(-Kx/(sqrt(3)/2),Kx/(sqrt(3)/2),points);
dky=abs(ky(2)-ky(1));

E=linspace(-4,4,600); % Energy in meV 600
% changeng E step doesn't change final DOS. 

 
for count=1:1:length(E)
    
    
    % Generate k points in FBZ
    for count_ky=1:1:length(ky)
        
        if ky(count_ky)<-Ky|ky(count_ky)>Ky
            Kx_BZ=abs(abs(ky(count_ky))-Kx/(sqrt(3)/2))*sqrt(3);
        else
            Kx_BZ=Kx;
        end
 
        kx=linspace(-Kx_BZ,Kx_BZ,points);
        dkx=abs(kx(2)-kx(1));

        for count_kx=1:1:length(kx)
            % Get eigenvalue of ten band Hamiltonian
            Eigenvalue=real(eig(Hamil_10(kx(count_kx),ky(count_ky),a_moire)));
            Ek=sort(Eigenvalue,'ascend');
            delta=exp(-(E(count)-Ek).^2./(2*sigma^2))./(sqrt(2*pi*sigma^2));

            % 
            delta_kx(count_kx)=sum(delta);
        end
        
        
        DOS_ky(count_ky)=sum(delta_kx).*dkx;

    end
    % Integration over FBZ
    DOS(count)=sum(DOS_ky).*dky;
    
    count/length(E)*100 % percents of code 
end

% for simplicity, DOS is a.u. 

plot(DOS*10e6,E,'color','k','linewidth',2)
%semilogx(DOS*10e6,E,'color','k','linewidth',2)
hold on 

%scatter(DOS*10e6,E)
set(gca,'fontsize',34)
xlabel(['DOS (a.u.)'],'FontSize',34)
%ylabel(['E (meV)'],'FontSize',28)

set(gcf,'Position',[500 300 300 600])
 xlim([0,2])
% ylim([-10,10])

set(gcf,'PaperOrientation','landscape')
print(gcf, 'DOS_TBG.pdf', '-dpdf','-r0','-bestfit')



