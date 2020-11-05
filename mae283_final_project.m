clear all;
load('spk-mic10khz.mat');
p=1;
q=1;
dt=1/(25600);
Nyf=pi/dt;
%Splitting the data for training and validation
utr=u(1:4096,1);
ytr=y(1:4096,1);
uval=u(4097:8192,1);
yval=y(4097:8192,1);

%auto correlation estimate 
N=4096;
Ru=xcorr(u,'biased');
tau=-N+1:N-1;
%41 lags?

%cross correlation for lags
[Ryu lag]=xcorr(ytr,utr);


figure(1);
plot(tau,Ru);
xlabel('$\tau$','Interpreter','latex');
ylabel('$\hat Ru$','Interpreter','latex');
title('$Input$ $covariance$','Interpreter','latex');

%input is white noise with lambda=1.1
% reorder the R_u via the variable 'index'
index=[N:2*N-1 2*N-1 1:N-1];
figure(2)
plot(Ru(index));
xlabel('index');ylabel('R_u(\tau)');
title(['Reordered (symmetric) auto covariance for FFT']);
%Phi estimate

P_u=fft(Ru(index));
% disp('Computation of Phi_u on the basis of re-ordered R_u: notice how Phi_u(w) is now real valued!');
 disp(P_u(1:10))

% computation of spectrum via periodogram
U=fft(u);
P_up=U.*conj(U)/N;

% creation of a frequency vector
Deltaf=2*pi/N;
w=[0:Deltaf:pi];

% plot results

% notice how we plot P_u (computed via correlation function) for every 2nd
% value, as it was based on 2N (symmetric) data points in R_u instead of the N
% data points of the actual u or the U=fft(u), allowing us to compare the
% spectrum estimated via correlation functions and estimated via the
% periodogram

figure(3)plot loglog
plot(log(w),real(P_u(1:2:N+1)))
xlabel('$log(\omega)$','Interpreter','latex');
ylabel('$\hat \phi u$','Interpreter','latex');
title('$Spectrum$ $of$ $the$ $input$','Interpreter','latex');

 legend('spectral estimate via correlation function')
% %this signal covers the whole range of pi/del_t
% 
% 
% 
% 
% % %non-parametric estimation
% % %etfe
Getf=etfe([ytr utr]);
% % %Spectral analysis
%Hanning window lag size 30
Gspa=spa([ytr utr]);
[mag phi]=bode(Gspa);


gk2=ifft(bode(Gspa),'symmetric');
% gk=ifft(mag);
% 
% 
% % %estimate is 6 pole 4 zero
% %remove ytr from FIR
% %run from for loop 
% %impulse plot
% %
% %  FIR model for least squares (6th order model)
% gtheta=[utr(7:N) utr(6:N-1) utr(5:N-2) utr(4:N-3) utr(3:N-4) utr(2:N-5) utr(1:N-6) -ytr(6:N-1) -ytr(5:N-2) -ytr(4:N-3) -ytr(3:N-4) -ytr(2:N-5) -ytr(1:N-6)]\ytr(7:N);
% Gtheta=tf(gtheta(1:7)',[1 gtheta(8:13)'],dt);
gtheta=[utr(7:N) utr(6:N-1) utr(5:N-2) utr(4:N-3) utr(3:N-4) utr(2:N-5) utr(1:N-6)]\ytr(7:N);
Gtheta=tf(gtheta',[0 0 0 0 0 1 0],dt);
kbar=13;
U=toeplitz([utr(1) zeros(1,N-1)],utr);
U1=U(1:kbar,1:end);U2=U(kbar+1:end,1:end);
%U2p=eye(N)-U2'*inv(U2*U2')*U2;
U2p=eye(N)-U2'*(U2'\eye(N));
U1m=[U1*U2p]';
thetafir=U1m\[ytr'*U2p]';

% 
 figure(5)
kpl=[zeros(20,1);thetafir;zeros((kbar+10)-length(thetafir),1)];
sgm=std([zeros(20,1);thetafir;zeros((kbar+10)-length(thetafir),1)]);
plot(-20:(kbar+9),kpl,'b');
hold on
plot(-20:(kbar+9),kpl-3*sgm,'r');
hold on
plot(-20:(kbar+9),kpl+3*sgm,'g');
ylabel('$Impulse$ $Response$ $Coefficients$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
title('$Impulse$ $Response$ $Coefficients$','Interpreter','latex');
legend('$g(k)$','$g(k)-\alpha (k)$','$g(k)+\alpha (k)$','Interpreter','latex');

% FIR error
for i=1:N
    yfir(i)=0;
        if(i==1)
            yfir(i)=y(i);
        else
            for j=1:13
                if (j>=i)
                    break;
                else
                    yfir(i)=yfir(i)+thetafir(j)*utr(i-j);
                end
            end
        end
        efir(i)=ytr(i)-yfir(i);
end
 figure(20)
Reufir=xcorr(efir,utr);
sgm1=std(Reufir/1000);
plot(-length(Reufir)/2:length(Reufir)/2-1,Reufir/1000,'b');
hold on
plot(-length(Reufir)/2:length(Reufir)/2-1,Reufir/1000+3*sgm1,'r');
hold on
plot(-length(Reufir)/2:length(Reufir)/2-1,Reufir/1000-3*sgm1,'g');
ylabel('$\hat R_eu$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
title('$Error$ $Input$ $Covariance$','Interpreter','latex');
legend('$\hat R_eu$','$\hat R_eu+\alpha (k)$','$\hat R_eu-\alpha (k)$','Interpreter','latex');





plot(yest,'og');
hold on
plot(ytr);
plot(gtheta,'g');
hold on
plot(0.99*gtheta,'r');
hold on
plot(1.01*gtheta,'r');
xlabel('$k$','Interpreter','latex');
ylabel('$Impulse$ $response$','Interpreter','latex');
title('$FIR$ $estimate$','Interpreter','latex');

%question 3b
tauu=1:2*N-1;

%Reu
[Reu flag]=xcorr(ytr-yest,utr,'biased');
plot(xcorr(yest,utr));

plot(tauu,Reu);
xlabel('$\tau$','Interpreter','latex');
ylabel('$\hat Reu$','Interpreter','latex');
title('$estimated$ $cross-correlation$ $function$','Interpreter','latex');

%hankel matrix for realization
gk2=real(gk2);
H=hankel(gk2(2:128/2),gk2(128/2:128-1));
Hbar=hankel(gk2(3:128/2+1),gk2(128/2+1:128));
[U,S,V]=svd(H);
figure(6);
plot(diag(S),'*')
xlabel('$N$','Interpreter','latex');
ylabel('$Singular$ $Values$','Interpreter','latex');
title('$Singular$ $Values$ $of$ $the$ $Hankel$ $Matrix$','Interpreter','latex');

%Kung's algortihm for realization
n1=10;
Dhat=gk2(1);
Sig=S(1:n1,1:n1);
Un=U(:,1:n1);
Vn=V(1:n1,:);

H1=Un*Sig^(0.5);
H2=Sig^(0.5)*Vn;

Bhat=H2(:,1);
Chat=H1(1,:);

H1p=Sig^(-0.5)*Un';
H2p=Vn'*Sig^(-0.5);

Ahat=H1p*Hbar*H2p;

Gre=ss(Ahat,Bhat,Chat,Dhat,dt);
yre=lsim(Gre,uval);

for i=1:128
 if(i==1)
     gknew(i)=Dhat;
 else
     gknew(i)=Chat*(Ahat^(i-1))*Bhat;
 end
end




plot(-20:19,[zeros(20,1);gknew(1:20)'],'b');
xlabel('$N$','Interpreter','latex');
ylabel('$Impulse$ $Response$ $Coefficients$','Interpreter','latex');
title('$Impulse$ $Response$ $Coefficients$','Interpreter','latex');
 
hold on;
plot(kpl,'b');
plot(gk2);

essr=yval-yre;
Reur=xcorr(essr,uval,'biased');
sgmr=std(Reur);
plot(-length(Reur)/2:length(Reur)/2-1,Reur,'b');
hold on
plot(-length(Reur)/2:length(Reur)/2-1,Reur+3*sgmr,'r');
hold on
plot(-length(Reur)/2:length(Reur)/2-1,Reur-3*sgmr,'g');
xlabel('$N$','Interpreter','latex');
ylabel('$\hat Reu$','Interpreter','latex');
title('$Covariance$ $of$ $error$ $with$ $input$','Interpreter','latex');
legend('$\hat Reu$','$\hat Reu+\alpha$','$\hat Reu-\alpha $','Interpreter','latex');

%Model structure question5
sysmod=oe([ytr utr],[10 8 1 ]);
Goe=tf(sysmod);
bode(Goe);
sysmod=armax([ytr utr],[10 10 10 1]);
Garx=tf(sysmod);
bode(Garx);

yarx=lsim(Garx,uval);
earx=yval-yarx;
Reuarx=xcorr(earx,uval,'biased');
sgmarx=std(Reuarx);
plot(-length(Reuarx)/2:length(Reuarx)/2-1,Reuarx,'b');
hold on
plot(-length(Reuarx)/2:length(Reuarx)/2-1,Reuarx+3*sgmarx,'r');
hold on
plot(-length(Reuarx)/2:length(Reuarx)/2-1,Reuarx-3*sgmarx,'g');
xlabel('$N$','Interpreter','latex');
ylabel('$\hat Reu$','Interpreter','latex');
title('$Covariance$ $of$ $error$ $with$ $input$','Interpreter','latex');
legend('$\hat Reu$','$\hat Reu+\alpha$','$\hat Reu-\alpha $','Interpreter','latex');

%predictor
opt = n4sidOptions('Focus','simulation');
init_sys = n4sid([ytr(1:400) utr(1:400)],10,opt);
sysp=pem([ytr(1:400) utr(1:400)],init_sys);
compare([ytr(1:400) utr(1:400)],sysp,init_sys);



