
% Stuff for jack

 %Global Params
Mu_0=4*pi;
D=90*0.72527*1e08; %65.27*10^8;      %      
%phi=4.767*10^10; %Phi_0
M= 1.98e33;
G= 6.674e-8;
R_s = 6.958e10;

%Inputted params
U= 30;
a=0.3;
LB = (a/2) + 0.1;
%B_phot = 4.64937
BPhot = 4;
m = 10000; %2538e04;


%things to calc.
phi=(pi*D*BPhot)/2;
%% Shit we need
h = chebfun('h', [LB U]); %We're letting [0.01 2] be the domain because h=0 will break I_eq and the log function, because you've got a 0/0 and I can't be arsed to talk limits.

BExt = chebfun(@(x) (2*phi)./(D*pi*( 1+ h(x).^2)), [LB U] );

g= chebfun(@(x) (G*M)./((R_s + D*h(x)).^2), [LB U]);

I_eq = chebfun(@(x) ((pi*h(x)*D)/(2*phi)).*(BExt(x) + ((BExt(x).^2 + ((4.*m.*g(x))./(D.*h(x))))).^(0.5)), [LB U]);

%[I_Max, hmax] = max(I_eq); %Gives value of Max(I_eq) and its location, which should be hmax.
 I_Max = 1;
 hmax =1;
 
F = chebfun(@(x) 2*I_Max * log((2/a)*hmax) - 2*atan(hmax - a) , [LB U] );   %a/d = 0.1

I_evol2 = chebfun2(@(x,y) ((F(x)./y + (2.*atan(h(x) - a)))./(2.*log((2/a).*h(x)))), [LB U LB U]); 

%savefig(f1, 'IEqvsIEvol.fig')

%% Rootfind

Distance = chebfun2(@(x,y) I_eq(x)-I_evol2(x,y), [LB U LB U]);
%DDiff = diff(Distance,1,1);    % one of these is the x derivative
%DDiff2 = diff(Distance,1,2);    % one of these is the x derivative
r=roots(Distance);

%Instead of doing this by inspection, we can definitely do better.

k=imag(r);
[PhiMod b] = min(k);

Intersection = real(r(b));

I_evol = chebfun(@(x) ((F(x)./PhiMod + (2.*atan(h(x) - a)))./(2.*log((2/a).*h(x)))), [LB U]);
%I_evol3 = chebfun(@(x) ((F(x)./(0.976) + (2.*atan(h(x) - 0.1)))./(2.*log(20.*h(x)))), [0.1 U]);
f1=figure; 
hold on
   plot(I_eq,'-r','Linewidth',2)
plot(I_evol,'-b','Linewidth',2)
 plot(Intersection,I_evol(real(r(b))),'x-k','Linewidth',2);
%plot(I_evol3,'-g')
xlabel('h/D');
ylabel('I');
legend('I_{eq}', 'I_{evol}')
hold off    
        

f2=figure;
plot(r) 
hold on
 plot(r(b),'x-k','Linewidth',4);
 hold off
xlabel('h values at intersection of I_{eq} and I_{evol}')
ylabel('\phi_{mod}')

%% Other stuff

% f2=figure;
% hold on
%  plot(F,'-r')
% plot(I_evol,'-b')
% xlabel('h/D');
% ylabel('I');
% legend('F', 'I_{evol}')
% hold off
% savefig(f2, 'FvsIEvol.fig')
% 
% f3=figure; 
% hold on
%    plot(I_eq,'-r')
% plot(I_EvolSimplified,'-b')
% xlabel('h/D');
% ylabel('I');
% legend('I_{eq}', 'I_{evol}')
% hold off    
%         

% figure(02)
% hold on
% plot(NI_eq)
% plot(NI_evol)
% xlabel('h/D');
% ylabel('I/I_{max}');

% figure(03)
% hold on
% plot(I_evol,'-b')
% plot(I_evolb,'-r')
% xlabel('h/D');
% ylabel('I/I_{max}');
% hold off
% 
% figure(04)
% hold on
% plot(I_evol,'-b')
% plot(I_evolc,'-r')
% plot(I_eq,'-g')
% xlabel('h/D');
% ylabel('I/I_{max}');
% hold off


%Single intersection exists at h= hmax = 0.8705.


%% Finding the intersection.




% u = chebfun(@(x) I_eq(x)-I_evol(x), [0.1 2]);
% 
% L = chebop(0.1, 2);
% L.op = @(x,u) diff(u,1);
% du = L*u;
% % Few ideas: 
% First we need to put phimod in the above as a variable, which is fine.
%
% 1. Formulate this as an differential equation, construct some sensible
% boundary condition and solve. 
%
% 2. Simply rootfind u and du/dx and compare the two, take the point that
% satisfies u = u' = 0, if it exists.



