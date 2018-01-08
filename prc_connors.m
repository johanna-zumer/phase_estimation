% This model is the Connors-Stevens model, similar to Hodgkin-Huxley, but
% more like neurons in the cortex, being type-I. 
% See Dayan and Abbott pp 166-172 then pp.196-198 and p.224.
clear
dt = 0.000005;
tmax=1;

iclamp_flag = 1; % if this is 1, run under current clamp conditions
vclamp_flag = 0; % otherwise this should be 1, for voltage clamp conditions

istart = 0; % time applied current starts
ilength=1;   % length of applied current pulse
%Ie=2.405e-7;     % magnitude of applied current pulse
Ie=2.5e-7;     % magnitude of applied current pulse
%Ie = 0.78e-7; % threshold for constant spiking with no A-current

vstart = 0.25;  % time to step voltage
vlength = 0.5;  % length of voltage step
V0 = -0.080     % initial voltage before step
Ve = -0.040;    % value of votage stepped to

V_L = -0.070;   % leak reversal potential
E_Na = 0.055;   % reversal for sodium channels
E_K = -0.072;   % reversal for potassium channels
E_A = -0.075;   % reversal for A-type potassium channels

g_L = 3e-6;     % specific leak conductance
g_Na = 1.2e-3;  % specific sodium conductance
g_K = 2e-4;     % specific potassium conductance
g_A = 4.77e-4;  % specific A-tpe potassium conductance
%g_A = 0.0;     % if g_A is zero it switches off the A-current

cm = 10e-9;     % specific membrane capacitance

t=0:dt:tmax;    % time vector
V=zeros(size(t)); % voltage vector
spikes = zeros(size(t)); % vector of spike times
spikenow = 0;   % variable that is 1 during a spike
Vspike = 0.0; % voltage at which a spike is registered

if ( iclamp_flag ) % i.e. if in current-clamp mode 
    V(1) = V_L;    % set the inititial value of voltage     
end

n=zeros(size(t));   % n: potassium activation gating variable
n(1) = 0.0;         % start off at zero
m=zeros(size(t));   % m: sodium activation gating variable
m(1) = 0.0;         % start off at zero
h=zeros(size(t));   % h: sodim inactivation gating variable
h(1) = 0.0;         % start off at zero

a=zeros(size(t));   % A-current activation gating variable
a(1) = 0.0;         % start off at zero
b=zeros(size(t));   % A-current inactivation gating variable
b(1) = 0.0;         % start off at zero

Iapp=zeros(size(t)); % Applied current, relevant in current-clamp mode
if ( iclamp_flag )   % i.e. if in current-clamp mode 
    for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
        Iapp(i) = Ie;
    end
end

Vapp=zeros(size(t)); % Applied voltage, relevant in voltage-clamp mode
if ( vclamp_flag )
    for i = 1:round(vstart/dt)          % % make V0 before pulse
        Vapp(i) = V0;
    end
    for i=round(vstart/dt)+1:round((vstart+vlength)/dt) % make Ve for duration of voltage pulse
        Vapp(i) = Ve;
    end
    for i=round((vstart+vlength)/dt):length(Vapp) % make V0 following pulse
        Vapp(i) = V0;
    end
end

Itot=zeros(size(t)); % in case we want to plot and look at the total current

for i = 2:length(t); % now see how things change through time
    I_L = g_L*(V_L-V(i-1));
    
    Vm = V(i-1)*1000; % converts voltages to mV as needed in the equations on p.224 of Dayan/Abbott
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta. Written out from Dayan/Abbott, units are 1/ms.
    
    
    if ( Vm == -29.7 ) 
        alpha_m = 0.38/0.1;
    else 
        alpha_m = 0.38*(Vm+29.7)/(1-exp(-0.1*(Vm+29.7)));
    end
    beta_m = 15.2*exp(-0.0556*(Vm+54.7));

    alpha_h = 0.266*exp(-0.05*(Vm+48));
    beta_h = 3.8/(1+exp(-0.1*(Vm+18)));
    
    if ( Vm == -45.7 ) 
       alpha_n = 0.02/0.1;
    else
        alpha_n = 0.02*(Vm+45.7)/(1-exp(-0.1*(Vm+45.7)));
    end
    beta_n = 0.25*exp(-0.0125*(Vm+55.7));
     
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1e-3/(alpha_m+beta_m);      % time constant converted from ms to sec
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1e-3/(alpha_h+beta_h);      % time constant converted from ms to sec
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1e-3/(alpha_n+beta_n);      % time constant converted from ms to sec
    n_inf = alpha_n/(alpha_n+beta_n);   
    
    m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
    
    h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
    
    n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
    
    % For the A-type current gating variables, instead of using alpha and
    % beta, we just use the steady-state values a_inf and b_inf along with 
    % the time constants tau_a and tau_b that are found empirically
    % (Dayan-Abbott, p. 224)
    
    a_inf = (0.0761*exp(0.0314*(Vm+94.22))/(1+exp(0.0346*(Vm+1.17))))^(1/3.0);
    tau_a = 0.3632*1e-3 + 1.158e-3/(1+exp(0.0497*(Vm+55.96)));
    
    b_inf = (1/(1+exp(0.0688*(Vm+53.3))))^4;
    tau_b = 1.24e-3 + 2.678e-3/(1+exp(0.0624*(Vm+50)));
    
    a(i) = a(i-1) + (a_inf-a(i-1))*dt/tau_a;    % Update a
    b(i) = b(i-1) + (b_inf-b(i-1))*dt/tau_b;    % Update b
    
    I_Na = g_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i-1)); % total sodium current
    
    I_K = g_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i-1)); % total potassium current
    
    I_A = g_A*a(i)*a(i)*a(i)*b(i)*(E_A-V(i-1)); % total A-type current
    
    Itot(i-1) = I_L+I_Na+I_K+I_A+Iapp(i-1); % total current is sum of leak + active channels + applied current
    
    V(i) = V(i-1) + Itot(i-1)*dt/cm;        % Update the membrane potential, V.

    if ( spikenow == 0 )                % If not already in a spike
        if ( V(i) > Vspike )            % and if the voltage is above a spikeing level
            spikenow = 1;               % then let code know it is in a spike
            spikes(i) = 1;              % and record the time of the spike
        end
    else
        if ( V (i) < Vspike-0.010 )     % If the voltage drops below spiking level
            spikenow = 0;               % then let the code know it is not in a spike
        end
    end
    
    if ( vclamp_flag )      % if we are using voltage clamp
        V(i) = Vapp(i);     % ignore the voltage integration and set V to be the applied voltage
    end
    
end

subplot(3,1,1)
plot(t,V);

spiketimes = find(spikes)

if ( sum(spikes) < 2 ) 
    print "ERROR: NOT ENOUGH SPIKES! "
else
    ISI = dt*diff(spiketimes);
end


period = ISI(2)

firstspiketime = dt*spiketimes(1);      % time of first spike
numpoints = 60;                         % number of points for phase response curve
phaseshift = zeros(numpoints);          % initially set phase shifts to zero
deltaI = 1e-8;                          % amplitude of current pulse
lenpulse = 0.001;                       % length of current pulse on each trial
phase = zeros(numpoints);               % phase in oscillation at which pulse is applied
dshift = period/numpoints;              % difference in phase between neighboring data points
npulse = lenpulse/dt;                   % number of time points to maintain pulse
totalt = 4*period;                      % total reduced time for simulation
totali = totalt/dt;                     % total number of time points

t2=0:dt:totalt;                         % time vector for a shorter total time

for npt=1:numpoints                     % loop through different times to generate a phase shift
    tshift = firstspiketime + dshift*npt;   % time of the current pulse
    tshift
    ishift = tshift/dt;                 % time step number of pulse
    phase(npt) = 2.*3.14159*npt/numpoints;  % corresponding phase in the cycle
    V2=zeros(size(t2));                 % voltage vector set to zero
    V2(1) = V_L;                        % first value intitialized
    
    Itot2=zeros(size(t2));              % set total currents to zero
    spikenow = 0;                       % assume not initially in a spike
    spikes2 = zeros(size(t2));          % set a spike array to zero

    for ( i = 2:totali )                % now loop through time 
        
       I_L = g_L*(V_L-V2(i-1));         % set leak current
    
        Vm = V2(i-1)*1000; % converts voltages to mV as needed in the equations on p.224 of Dayan/Abbott
        
        if ( Vm == -29.7 ) 
            alpha_m = 0.38/0.1;
        else 
            alpha_m = 0.38*(Vm+29.7)/(1-exp(-0.1*(Vm+29.7)));
        end
       beta_m = 15.2*exp(-0.0556*(Vm+54.7));

       alpha_h = 0.266*exp(-0.05*(Vm+48));
        beta_h = 3.8/(1+exp(-0.1*(Vm+18)));
    
        if ( Vm == -45.7 ) 
            alpha_n = 0.02/0.1;
        else
            alpha_n = 0.02*(Vm+45.7)/(1-exp(-0.1*(Vm+45.7)));
        end
        beta_n = 0.25*exp(-0.0125*(Vm+55.7));

        tau_m = 1e-3/(alpha_m+beta_m);      % time constant converted from ms to sec
        m_inf = alpha_m/(alpha_m+beta_m);
    
        tau_h = 1e-3/(alpha_h+beta_h);      % time constant converted from ms to sec
        h_inf = alpha_h/(alpha_h+beta_h);
    
        tau_n = 1e-3/(alpha_n+beta_n);      % time constant converted from ms to sec
        n_inf = alpha_n/(alpha_n+beta_n);   
    
        m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
    
        h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
    
        n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
     
        a_inf = (0.0761*exp(0.0314*(Vm+94.22))/(1+exp(0.0346*(Vm+1.17))))^(1/3.0);
        tau_a = 0.3632*1e-3 + 1.158e-3/(1+exp(0.0497*(Vm+55.96)));
    
        b_inf = (1/(1+exp(0.0688*(Vm+53.3))))^4;
        tau_b = 1.24e-3 + 2.678e-3/(1+exp(0.0624*(Vm+50)));
    
        a(i) = a(i-1) + (a_inf-a(i-1))*dt/tau_a;    % Update a
        b(i) = b(i-1) + (b_inf-b(i-1))*dt/tau_b;    % Update b
    
        I_Na = g_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V2(i-1)); % total sodium current
    
        I_K = g_K*n(i)*n(i)*n(i)*n(i)*(E_K-V2(i-1)); % total potassium current
    
        I_A = g_A*a(i)*a(i)*a(i)*b(i)*(E_A-V2(i-1)); % total A-type current
    
        Itot2(i-1) = I_L+I_Na+I_K+I_A+Iapp(i-1); % total current is sum of leak + active channels + applied current
    
        if ( i > ishift ) && ( i < ishift+npulse+1 ) % this is the time of the current pulse
            Itot2(i-1) = Itot2(i-1) + deltaI;        % then add an extra current
        end
        
        V2(i) = V2(i-1) + Itot2(i-1)*dt/cm;        % Update the membrane potential, V.

        if ( spikenow == 0 )                % If not already in a spike
           if ( V2(i) > Vspike )            % and if the voltage is above a spikeing level
                spikenow = 1;               % then let code know it is in a spike
                spikes2(i) = 1;             % and record the time of the spike
           end
        else
           if ( V2(i) < Vspike-0.010 )      % If the voltage drops below spiking level
                spikenow = 0;               % then let the code know it is not in a spike
           end
        end

    end
 
%    subplot(numpoints,1,npt)
%    plot(t2,V2)
%    hold on
%    plot(t2,1e4*Itot2,'r')
    ISI = dt*diff(find(spikes2));           % ISI holds the set of all inter-spike intervals
    if ( tshift <  firstspiketime+period )  % if the kick was in the second ISI
        phaseshift(npt) = 2.*3.14159*(period-ISI(2))/period;
    else
        phaseshift(npt) = 2.*3.14159*(period-ISI(3))/period; % otherwise assume in the 3rd ISI
    end
    phaseshift(npt)             % write out the phaseshift
    ISI                         % and write out the ISI
    clear spikes2;
    clear ISI;
end
figure(2)
%subplot(3,1,2)
hold on
plot(phase,phaseshift,'.');