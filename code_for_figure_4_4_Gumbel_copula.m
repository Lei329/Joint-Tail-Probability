%Figure 4.3 Gumbel copula when \alpha=1
clear all
clc 

alpha=1.5;%heavy tails
gamma=1;%paremater for pareto distribution
b=1;%constant; see Theorem 1
lambda=1;%paremater for N(t)
x1=100;x2=200;N9=100;dx=(x2-x1)/N9;XX=[x1:dx:x2]; Tail=zeros(1,N9);%paremater for x
kappa=1;%see Remark 1 for FGM case
delta2=2.2;delta3=1;delta4=1;%paremater for Brownian motion with jumps
mu=1;%paremater for Y
phio1= -delta2*alpha+delta3^2*alpha^2/2-delta4*alpha/(mu+alpha);%\phi(\alpha)
phio2= -delta2*2*alpha+delta3^2*(2*alpha)^2/2-delta4*2*alpha/(mu+2*alpha);%\phi(2\alpha)
phio3= -delta2*alpha*kappa+delta3^2*(kappa*alpha)^2/2-delta4*alpha*kappa/(mu+kappa*alpha);%\phi(\kappa*\alpha)


t=1;%time
vartheta=2;%paremater for Gumbel copula
tau=1+b-(1^(vartheta)+b^(vartheta))^(1/vartheta);%see Remark 1 for Gumbel case
K1=(2*lambda^2/phio1)*((exp(phio2*t)-exp(phio1*t))/(phio2-phio1)-(exp(phio2*t)-1)/phio2);%K_1; see Remark 4(i)
K2=lambda*(exp(phio3*t)-1)/phio3;%K_2; see Remark 4(i)
for i=1:N9
    Tail(i)=b*K1*(gamma/(gamma+XX(i)))^(alpha*2)+tau*K2*(gamma/(gamma+XX(i)))^(alpha*kappa);%See Theorem 1
end
plot(XX(1:end-1),Tail,'r');hold on

vartheta=3;%paremater for Gumbel copula
tau=1+b-(1^(vartheta)+b^(vartheta))^(1/vartheta);%see Remark 1 for Gumbel case
for i=1:N9
    Tail(i)=b*K1*(gamma/(gamma+XX(i)))^(alpha*2)+tau*K2*(gamma/(gamma+XX(i)))^(alpha*kappa);%See Theorem 1
end
plot(XX(1:end-1),Tail,'green');hold on

vartheta=4;%paremater for Gumbel copula
tau=1+b-(1^(vartheta)+b^(vartheta))^(1/vartheta);%see Remark 1 for Gumbel case
for i=1:N9
    Tail(i)=b*K1*(gamma/(gamma+XX(i)))^(alpha*2)+tau*K2*(gamma/(gamma+XX(i)))^(alpha*kappa);%See Theorem 1
end
plot(XX(1:end-1),Tail,'blue');hold on

t=2;%time
vartheta=2;%paremater for Gumbel copula
tau=1+b-(1^(vartheta)+b^(vartheta))^(1/vartheta);%see Remark 1 for Gumbel case
K1=(2*lambda^2/phio1)*((exp(phio2*t)-exp(phio1*t))/(phio2-phio1)-(exp(phio2*t)-1)/phio2);%K_1; see Remark 4(i)
K2=lambda*(exp(phio3*t)-1)/phio3;%K_2; see Remark 4(i)
for i=1:N9
    Tail(i)=b*K1*(gamma/(gamma+XX(i)))^(alpha*2)+tau*K2*(gamma/(gamma+XX(i)))^(alpha*kappa);%See Theorem 1
end
plot(XX(1:end-1),Tail,'r--');hold on

vartheta=3;%paremater for Gumbel copula
tau=1+b-(1^(vartheta)+b^(vartheta))^(1/vartheta);%see Remark 1 for Gumbel case
for i=1:N9
    Tail(i)=b*K1*(gamma/(gamma+XX(i)))^(alpha*2)+tau*K2*(gamma/(gamma+XX(i)))^(alpha*kappa);%See Theorem 1
end
plot(XX(1:end-1),Tail,'green--');hold on

vartheta=4;%paremater for Gumbel copula
tau=1+b-(1^(vartheta)+b^(vartheta))^(1/vartheta);%see Remark 1 for Gumbel case
for i=1:N9
    Tail(i)=b*K1*(gamma/(gamma+XX(i)))^(alpha*2)+tau*K2*(gamma/(gamma+XX(i)))^(alpha*kappa);%See Theorem 1
end
plot(XX(1:end-1),Tail,'blue--');hold on

xlabel('x')
ylabel('Joint Tail Probability')
legend('t=1,\vartheta=2','t=1,\vartheta=3','t=1,\vartheta=4','t=2,\vartheta=2','t=2,\vartheta=3','t=2,\vartheta=4','Location','northeast')

