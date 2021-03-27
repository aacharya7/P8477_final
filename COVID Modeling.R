library(deSolve)

#Documented cases only 
COVID=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # rate of change for DOCUMENTED cases
    dS =  -S/N * (beta.I * I + beta.H * H) + mu*N - mu*S;
    dE = S/N * (beta.I * I + beta.H * H)  -  alpha*E -mu*E;
    dI = alpha*E - I*( (gamma_h*theta1*(1-delta2)) + (gamma_f*(1-theta1)*(1-delta2)) + (gamma_i*(1-theta1)*delta2) ) - mu*I;
    dH = gamma_h*theta1*(1-delta2)*I - H*( (gamma_d*delta1)  +  (gamma_r*(1-delta1)) ) - mu*H;
    dR = gamma_r*(1-delta1)*H + (gamma_f*(1-theta1)*(1-delta2))*I - mu*R;
    dDead = (gamma_i*(1-theta1)*delta2)*I + H*(gamma_d*delta1) - mu*Dead;
    dcumInci = alpha * E; # cumulative incidence
    
    # return the rate of change
    list(c(dS,dE, dI,dH, dR, dDead, dcumInci))
  }) # end with(as.list...)
}

N=8500000;
I0=21;
E0=20*I0;
cumInci0=E0+I0;
H0=R0=Dead0=0;
S0=N-I0-E0;
state=c(S=S0,E=E0,I=I0,H=H0,R=R0,Dead=Dead0,cumInci=cumInci0)

## parameters for the COVID outbreak 2019 for DOCUMENTED cases with no intervention
# no intervention
R_0 = 2.7 #from (Wang et al, 2020)
alpha = 1/6.4 #latent period
gamma_f = 1/20.5 #from infected to recovered without hospitalization
gamma_i = 1/16 #from infection to death
gamma_h = 1/5.8 #from infected to hospitalized
gamma_d = 1/(16-5.8) #from hospitalization to death
gamma_r = 1/7 #rate of recovery
mu = 1/(78.7*365) #birth/death rate = 1/life expectancy in days

theta1 = 0.21 #proportion of infected who are hospitalized
delta1 = 0.26 #case fatality rate for hospitalized infections
delta2 = 0.0139 #case fatality rate for non-hospitalized infections

beta.I = ((gamma_h*theta1*(1-delta2)) + (gamma_f*(1-theta1)*(1-delta2)) + (gamma_i*(1-theta1)*delta2)+mu)*R_0*((alpha+mu)/alpha) #transmission in the community
beta.H = 1/9.5 #transmission rate in hospital

docparms.NoCtrl= c(alpha=alpha,
                   beta.I=beta.I,
                   beta.H=beta.H,
                   gamma_f=gamma_f,
                   gamma_i=gamma_i,
                   mu=mu,
                   gamma_h=gamma_h,
                   gamma_d=gamma_d,
                   gamma_r=gamma_r,
                   theta1=theta1,
                   delta1=delta1,
                   delta2=delta2)

times=seq(0,365,by=1) # run for 10 years
COVID.sim=ode(y=state,times=times,func=COVID,parms=docparms.NoCtrl);
inci=COVID.sim[seq(7,nrow(COVID.sim),by=7),'cumInci']-c(0,COVID.sim[seq(7,nrow(COVID.sim)-7,by=7),'cumInci']) #weekly incidence
  #plot 
  plot(inci,ylab='Weekly Incidence',xlab='Week',type='l',lwd=1, main="Weekly Incidence of COVID over 1 Year",col='blue')

death_week=COVID.sim[seq(7,nrow(COVID.sim),by=7),'Dead']-c(0,COVID.sim[seq(7,nrow(COVID.sim)-7,by=7),'Dead']) #weekly incidence of death


i.sim=COVID.sim[,'I']/N  # %I
s.sim=COVID.sim[,'S']/N #%S
d.sim=COVID.sim[,'Dead']/N # %Dead
all.sim=cbind(i.sim,s.sim,d.sim)

#Plot 
par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.8,.5,0))
matplot(COVID.sim[,'time'], COVID.sim[,c('I','Dead')],type='l', 
        lwd=1,col=c('blue','red'),lty=1, main='Infection Prevalence Over 1 Year', cex.main=1,
        ylab='Number of People',xlab='Time (days)')
legend('topright',c('Infected','Dead'),col=c('blue','red'),
       lty=1, cex=1, lwd=1, bty='n')

matplot(all.sim[,c('i.sim','s.sim','d.sim')],type='l', 
     lwd=1,col=c('blue','red','green'),lty=1, main='% of Population Affected by COVID Over 1 Year', cex.main=1,
     ylab='% of N',xlab='Time (days)')
legend('topright',c('%I','%S','%Dead'),col=c('blue','red','green'),
       lty=1, cex=1, lwd=1, bty='n')

par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.8,.5,0))
plot(s.sim,ylab='%S',xlab='Days',type='l',lwd=1, main="%S of COVID Over 1 Year",col='blue')
plot(i.sim,ylab='%I',xlab='Days',type='l',lwd=1, main="%I of COVID Over 1 Year",col='green')
plot(d.sim,ylab='%Dead',xlab='Days',type='l',lwd=1, main="%Dead of COVID Over 1 Year",col='red')


##################### Risk Structured Model with Documented and Undocumented Cases #####################
COVID_risk=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # rate of change for DOCUMENTED cases
    dSD =  -SD/N * (beta.IDD*ID + beta.IDU*ID + beta.H*HD + beta.H*HD ) + mu*ND - mu*SD;
    dED = SD/N * (beta.IDD* ID + beta.IDU*ID + beta.H*HD + beta.H*HD )  -  alpha*ED -mu*ED;
    dID = alpha*ED - ID*( (gamma_h*theta1*(1-delta2)) + (gamma_f*(1-theta1)*(1-delta2)) + (gamma_i*(1-theta1)*delta2) ) - mu*ID;
    dHD = gamma_h*theta1*(1-delta2)*ID - HD*( (gamma_d*delta1)  +  (gamma_r*(1-delta1)) ) - mu*HD;
    dRD = gamma_r*(1-delta1)*HD + (gamma_f*(1-theta1)*(1-delta2))*ID -mu*RD;
    dDeadD = (gamma_i*(1-theta1)*delta2)*ID+ HD*(gamma_d*delta1) -mu*DeadD;
    dcumInciD = alpha * ED; # cumulative incidence
    
    # rate of change for UNDOCUMENTED cases
    dSU =  -SU/N * (beta.IUU*IU + beta.IUD*IU) + mu*NU - mu*SU;
    dEU = SU/N * (beta.IUU*IU + beta.IUD*IU)  -  alpha*EU -mu*EU;
    dIU = alpha*EU - IU*( (gamma_f*(1-delta2)) + (gamma_i*delta2) ) - mu*IU;
    dRU = (gamma_f*(1-delta2))*IU + mu*RU ;
    dDeadU = (gamma_i*(1-theta2)*delta2)*IU - mu*DeadU;
    dcumInciU = alpha * EU; # cumulative incidence
    
    # return the rate of change
    list(c(dSD, dSU, dED, dEU, dID, dIU, dHD, dRD, dRU, dDeadD, dDeadU, dcumInciD, dcumInciU))
  }) # end with(as.list...)
}

N=8500000;
ND=N*(1-0.86);
ID_0=21;
ED_0=20*ID_0;
SD_0=ND-ID_0-ED_0;
cumInciD_0=ED_0+ID_0;
HD_0=RD_0=DeadD_0=0;

NU=N*0.86;
IU_0=21;
EU_0=20*IU_0;
SU_0=NU-IU_0-EU_0;
cumInciU_0=EU_0+IU_0;
RU_0=DeadU_0=0;

state_risk=c(SD=SD_0, SU=SU_0, ED=ED_0, EU=EU_0, ID=ID_0, IU=IU_0, 
             HD=HD_0,RD=RD_0, RU=RU_0, DeadD=DeadD_0, DeadU=DeadU_0, cumInciD=cumInciD_0, cumInciU=cumInciU_0)

## parameters for the COVID outbreak 2019 for risk structured model with no intervention
# no intervention
R_0 = 2.7 #from (Wang et al, 2020)
RD_0= 1.34 #simulated from (Li et al, 2020)
beta.IDD = ((gamma_h*theta1*(1-delta2)) + (gamma_f*(1-theta1)*(1-delta2)) + (gamma_i*(1-theta1)*delta2)+mu)*RU_0*((alpha+mu)/alpha) #transmission from D to U
beta.IUU = ((gamma_f*(1-theta1)*(1-delta2)) + (gamma_i*(1-theta1)*delta2)+mu)*R_0*((alpha+mu)/alpha) #transmission in community
beta.IUD = beta.IDU = 0.52
beta.H = 1/9.5 #transmission rate in hospital

beta.H = 1/9.5 #transmission rate in hospital
alpha = 1/6.4 #latent period
gamma_f = 1/20.5 #from infected to recovered without hospitalization
gamma_i = 1/16 #from infection to death
gamma_h = 1/5.8 #from infected to hospitalized
gamma_d = 1/(16-5.8) #from hospitalization to death
gamma_r = 1/7 #rate of recovery
mu = 1/(78.7*365) #birth/death rate = 1/life expectancy in days
theta1 = 0.21 #proportion of infected who are hospitalized
delta1 = 0.26 #case fatality rate for hospitalized infections
delta2 = 0.0139 #case fatality rate for non-hospitalized infections

risk.parms= c(alpha=alpha,
                   beta.IDD=beta.IDD,
                   beta.IUU = beta.IUU,
                   beta.IUD = beta.IUD,
                   beta.IDU = beta.IDU,
                   beta.H=beta.H,
                   gamma_f=gamma_f,
                   gamma_i=gamma_i,
                   mu=mu,
                   gamma_h=gamma_h,
                   gamma_d=gamma_d,
                   gamma_r=gamma_r,
                   theta1=theta1,
                   delta1=delta1,
                   delta2=delta2)

times=seq(0,365,by=1) # run for 10 years
COVID_risk.sim=ode(y=state_risk,times=times,func=COVID_risk,parms=risk.parms);
inci_D=COVID_risk.sim[seq(7,nrow(COVID_risk.sim),by=7),'cumInciD']-c(0,COVID_risk.sim[seq(7,nrow(COVID_risk.sim)-7,by=7),'cumInciD']) #weekly incidence
inci_U=COVID_risk.sim[seq(7,nrow(COVID_risk.sim),by=7),'cumInciU']-c(0,COVID_risk.sim[seq(7,nrow(COVID_risk.sim)-7,by=7),'cumInciU']) #weekly incidence
#plot
plot(inci_D,ylab='Weekly Incidence',xlab='Week',type='l',lwd=1, main="Weekly Documented Incidence of COVID over 1 Year",col='blue')
plot(inci_U,ylab='Weekly Incidence',xlab='Week',type='l',lwd=1, main="Weekly Undocumented Incidence of COVID over 1 Year",col='red')

death_weekD=COVID_risk.sim[seq(7,nrow(COVID_risk.sim),by=7),'DeadD']-c(0,COVID_risk.sim[seq(7,nrow(COVID_risk.sim)-7,by=7),'DeadD']) #weekly incidence of death
death_weekU=COVID_risk.sim[seq(7,nrow(COVID_risk.sim),by=7),'DeadU']-c(0,COVID_risk.sim[seq(7,nrow(COVID_risk.sim)-7,by=7),'DeadU']) #weekly incidence of death
death.week.inci = cbind(inci_D,inci_U,death_weekD,death_weekU)
death.week.inci_all = cbind(inci_D+inci_U,death_weekD+death_weekU)

i.sim_risk=(COVID_risk.sim[,'ID']+COVID_risk.sim[,'IU'])/N  # %I
s.sim_risk=(COVID_risk.sim[,'SD']+COVID_risk.sim[,'SU'])/N #%S
d.sim_risk=(COVID_risk.sim[,'DeadD']+COVID_risk.sim[,'DeadU'])/N # %Dead
all.sim_risk=cbind(i.sim_risk,s.sim_risk,d.sim_risk)

combined_sim=cbind(COVID_risk.sim[,'ID']+COVID_risk.sim[,'IU'],COVID_risk.sim[,'SD']+COVID_risk.sim[,'SU'],COVID_risk.sim[,'DeadD']+COVID_risk.sim[,'DeadU'],
                   COVID_risk.sim[,'cumInciD']+COVID_risk.sim[,'cumInciU'], COVID_risk.sim[,'time'])

s.simD=COVID_risk.sim[,'SD']/N #%S
s.simU=COVID_risk.sim[,'SU']/N #%S
s.U =cbind(s.simD,s.simU)

#Plot 
par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.8,.5,0))
matplot(COVID_risk.sim[,'time'], COVID_risk.sim[,c('cumInciD','DeadD')],type='l', 
        lwd=1,col=c('blue','red'),lty=1, main='Infection Prevalence Over 1 Year with Risk Structured Model', cex.main=1,
        ylab='Number of People',xlab='Time (days)')
legend('right',c('Documented Infected','Documented Dead'),col=c('blue','red','green','purple'),
       lty=1, cex=1, lwd=1, bty='n')
matplot(COVID_risk.sim[,'time'], COVID_risk.sim[,c('cumInciU','DeadU')],type='l', 
        lwd=1,col=c('blue','red'),lty=1, main='Infection Prevalence Over 1 Year with Risk Structured Model', cex.main=1,
        ylab='Number of People',xlab='Time (days)')
legend('topright',c('Undocumented Infected', 'Undocumented Dead'),col=c('blue','red'),
       lty=1, cex=1, lwd=1, bty='n')

matplot(all.sim_risk[,c('i.sim_risk','s.sim_risk','d.sim_risk')],type='l', 
        lwd=1,col=c('blue','red','green'),lty=1, main='% of Population Affected by COVID Over 1 Year', cex.main=1,
        ylab='% of N',xlab='Time (days)')
legend('topright',c('%I','%S','%Dead'),col=c('blue','red','green'),
       lty=1, cex=1, lwd=1, bty='n')

matplot(s.U[,c('s.simD','s.simU')],type='l', 
        lwd=1,col=c('blue','red','green'),lty=1, main='%S of Population Affected by COVID Over 1 Year', cex.main=1,
        ylab='% of N',xlab='Time (days)')
legend('topright',c('Documented','Undocumented'),col=c('blue','red','green'),
       lty=1, cex=1, lwd=1, bty='n')

matplot(death.week.inci[,c('inci_D','death_weekD')],type='l', 
        lwd=1,col=c('green','purple'),lty=1, main='Documented Weekly Incidence vs Death', cex.main=1,
        ylab='% of N',xlab='Time (weeks)')
legend('topleft',c('Incidence','Deaths'),col=c('green','purple'),
       lty=1, cex=1, lwd=1, bty='n')
matplot(death.week.inci[,c('inci_U','death_weekU')],type='l', 
        lwd=1,col=c('green','purple'),lty=1, main='Undocumented Weekly Incidence vs Death', cex.main=1,
        ylab='% of N',xlab='Time (weeks)')
legend('topleft',c('Incidence','Deaths'),col=c('green','purple'),
       lty=1, cex=1, lwd=1, bty='n')

matplot(combined_sim[,c(1,2,3)],type='l', 
        lwd=1,col=c('red','blue','green'),lty=1, main='Undocumented Weekly Incidence vs Death', cex.main=1,
        ylab='% of N',xlab='Time (weeks)')
legend('topright',c('Infected','Susceptible','Dead'),col=c('red','blue','green'),
       lty=1, cex=1, lwd=1, bty='n')
