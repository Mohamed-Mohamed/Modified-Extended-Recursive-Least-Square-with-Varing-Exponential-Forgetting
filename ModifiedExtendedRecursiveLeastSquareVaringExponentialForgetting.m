function [ Guz, Gez ] = ModifiedExtendedRecursiveLeastSquareVaringExponentialForgetting ( u, y, na, nb, d, nc, Ts, Plot )
% this function is used to estimate the parameter of the system
% by using Modified Extended RLS with Varing Exponential Forgetting  method
% estimated     y(z) = Guz u(z) + Gez epslon(z)
%                      z^(-d) B(z^-1)                                                C(z^-1)
% Guz  =  ------------------------------     ,            Gez  =  ------------------------------
%                          A(z^-1)                                                       A(z^-1)
% B(z^-1) = b_0 + b_1 z^-1 + ... + b_nb z^-nb
% A(z^-1) = 1 + a_1 z^-1 + ... + a_na z^-na
% C(z^-1) = 1 + c_1 z^-1 + ... + c_nc z^-nc
% Y=phi'*theta_hat
% where :
% phi'=[y(t-1) ... y(t-na) | u(t-1-d) ... u(t-nb-d) | epslon_p(t-1) ... epslon_p(t-nc)]
% theta_hat =[a_1 ... a_na | b_0 ... b_nb | c_1 ... c_nc]
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs:
% uc    : input required signal
% y      : output signal to the system
% na   : order of thr Den. of the transfer function of the system
% ab   : order of thr Num. of the transfer function of the system
% d     : order of delay of the system
% nc   : order of thr Num. of the transfer function of the system noise
% Ts     : sampling time
% Plot          : used to get the plot of system parameter estimation
%                   1 - if Plot = 0  --> no plot needed
%                   2- if Plot = 1  -->  plot needed and Plot(2:5) == figuires number
%% outputs:
% Gz              : discreate transfer function of the system
% Gez            : discreate transfer function of the system noise

%% function body
N=max(na,nb+d+1);
P{1}=eye(na+nb+nc+1)*10^(6);
theta_hat(:,1)=zeros(na+nb+nc+1,1);
epslon_n(1)=0;
epslon_p(1)=0;
y1(1)=y(1);
lambda(1)=0.9;
for i=2:length(y)
        for j=1:na
                if i-j <=0
                        phiT(i,j)=0;
                else
                        phiT(i,j)=[-y(i-j)];
                end
        end
        for j=0:nb
                if i-j-d <= 0
                        phiT(i,j+1+na)=0;
                else
                        phiT(i,j+1+na)=[u(i-j-d)];
                end
        end
        for j=1:nc
                if i-j <=0
                        phiT(i,j+na+nb+1)=0;
                else
                        phiT(i,j+na+nb+1)=[epslon_p(i-j)];
                end
        end
        K1{i}=P{i-1}*phiT(i,:)'*inv(lambda(i-1)+phiT(i,:)*P{i-1}*phiT(i,:)');
        epslon_n(i)=y(i)-phiT(i,:)*theta_hat(:,i-1);
        theta_hat(:,i)=theta_hat(:,i-1)+K1{i}*epslon_n(i);
        epslon_p(i)=y(i)-phiT(i,:)*theta_hat(:,i);
        P{i}=(eye(na+nb+nc+1)*lambda(i-1)-K1{i}*phiT(i,:))*P{i-1};
        if i > N
                lambdan=1-(1-phiT(i,:)*K1{i})*epslon_n(i)^2/(std(epslon_n)^2*mean(epslon_n));
                lambdap=1-(1-phiT(i,:)*K1{i})*epslon_p(i)^2/(std(epslon_p)^2*mean(epslon_p));
                Sn=std(epslon_n)^2/(std(epslon_n)^2+std(epslon_p)^2);
                Sp=std(epslon_p)^2/(std(epslon_n)^2+std(epslon_p)^2);
                lambda(i)=Sn*lambdan+Sp*lambdap;
        else
                lambda(i)=lambda(i-1);
        end
end

Theta_hat=theta_hat(:,end);
Guz=tf([Theta_hat(na+1:end-nc)'],[1,Theta_hat(1:na)'],Ts);
Gez=tf([1, Theta_hat(na+nb+2:end)'],[1,Theta_hat(1:na)'],Ts);

for l=1:length(phiT(:,1))
        y1(l)=phiT(l,:)*theta_hat(:,l);
end
Y=size(y);
Y1=size(y1);
if Y(1) ~= Y1(1)
        y1=y1';
end
%% plotting
if Plot(1)==1
        figure(Plot(2));
        set(gcf,'color','w')
        hold all;
        for k=1:na+nb+1
                if k <= na
                        H='-';
                else
                        H='--';
                end
                plot((0:length(y1)-1)*Ts,theta_hat(k,:),H,'linewidth',2);
        end
        grid on;
        for m=1:na
                Ylabel{m}=['a_' num2str(m) ', '];
                Leg{m}=['a_' num2str(m) ];
        end
        for n=1:nb+1
                if n<nb+1
                        Ylabel{n+na}=['b_' num2str(n-1) ', '];
                        Leg{n+na}=['b_' num2str(n-1) ];
                else
                        Ylabel{n+na}=['b_' num2str(n-1)];
                        Leg{n+na}=['b_' num2str(n-1) ];
                end
        end
        xlabel('Time (s)','fontsize',18);
        Ylabel=cell2mat(Ylabel);
        ylabel(Ylabel,'fontsize',18);
        legend(Leg)
        title('System parameter estimation with time','fontsize',18)
        
        
        
        figure(Plot(3));
        set(gcf,'color','w')
        hold all;
        for k=[1:na, na+nb+2:na+nb+nc+1]
                if k <= na
                        H='-';
                else
                        H='--';
                end
                plot((0:length(y1)-1)*Ts,theta_hat(k,:),H,'linewidth',2);
        end
        grid on;
        for m=1:na
                Ylabel2{m}=['a_' num2str(m) ', '];
                Leg2{m}=['a_' num2str(m) ];
        end
        for c=1:nc
                if c<nc
                        Ylabel2{c+na}=['c_' num2str(c) ', '];
                        Leg2{c+na}=['c_' num2str(c) ];
                else
                        Ylabel2{c+na}=['c_' num2str(c)];
                        Leg2{c+na}=['c_' num2str(c) ];
                end
        end
        xlabel('Time (s)','fontsize',18);
        Ylabel2=cell2mat(Ylabel2);
        ylabel(Ylabel2,'fontsize',18);
        legend(Leg2)
        title('System noise parameter estimation with time','fontsize',18)
        
        figure(Plot(4));
        set(gcf,'color','w')
        plot((0:length(lambda)-1)*Ts,lambda,'linewidth',2);
        grid on;
        ylabel('\lambda (t)','fontsize',18);
        xlabel('t(s)','fontsize',18);
        
        figure(Plot(5));
        subplot(3,1,1:2)
        set(gcf,'color','w')
        plot((0:length(y1)-1)*Ts,y,(0:length(y1)-1)*Ts,y1,'-o','linewidth',2);
        grid on;
        ylabel('y, y_e_s_t','fontsize',18);
        legend('y','y_e_s_t')
        subplot(3,1,3)
        plot((0:length(y1)-1)*Ts,abs(y-y1))
        grid on;
        xlabel('t(s)','fontsize',18);
        ylabel('y-y_e_s_t','fontsize',18);
        legend('error')
end
end
