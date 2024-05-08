# load libraries 
#install.packages(pacman)
pacman::p_load(deSolve,ggplot2,plyr)

# initiate state variables: host population abundances
init <- c(
  S=500, # susceptible subpopulation
  I1=100, #infected subpop by virulence class 1
  I2=100, #infected subpop by virulence class 2
  I3=100, #infected subpop by virulence class 3
  I4=100, #infected subpop by virulence class 4
  I5=100, #infected subpop by virulence class 5
  I6=100, #infected subpop by virulence class 6
  I7=100, #infected subpop by virulence class 7
  I8=100, #infected subpop by virulence class 8
  I9=100, #infected subpop by virulence class 9
  I10=100 #infected subpop by virulence class 10
)

# initate time of simulation 
times <- seq(0, 500, by = 1)

# simulation 1: vary exclusively host resistance

parameters_resis <- expand.grid(m_t = 0,
                           m_r = seq(0,1,0.05),
                           beta = 0.4,
                           gamma = 0.01,
                           miu_s = 0.1,
                           fecundity = 1
)
#Simulation ODE

SIM_MODEL <- function(times, init, parameters) {
  with(as.list(c(init, parameters)),{
  
  # extract state variables from init vector
    S = as.vector(init[1]) # susceptible host subpop
    I = matrix(init[2:length(init)], nrow=10, ncol=1) # infected host subpop
    
   # can it run without this? below b/c you have it already defined globally
  # modified parameters i.e. "realized" rates   
  #  beta_realized = matrix(data=NA, nrow=10,ncol=1) # realized transmission rate
  #  alpha = seq(0.1,1,0.1) # virulence class
  #  for (z in 1:10){
  #    beta_realized[z] <- alpha[z]*beta*(1-m_r) # realized transmission rate
  #  }
   # miu_i = (sqrt(miu_s)+ (1-m_t)*(beta_realized/beta))**2 # realized mortality rate
    
  # system of differential equations for SI model
    dS <- -1*S*sum(I*beta_realized) - miu_s * S + sum(gamma * I) + fecundity * (S+sum(I))
    dI <- I * beta_realized * S - (gamma+miu_i) * I
    return(list(c(dS,dI)))
  })
}
rootfun <- function(times, init, parameters) {
  dstate <- unlist(SIM_MODEL(times, init, parameters))
  sum(abs(dstate)) - 1e-1 #this function isn't working...
}

# run simulation!
for (i in 1:nrow(parameters_resis)){
  out <- lsodar(func = SIM_MODEL, y = init, parms = parameters_resis[i,],
                times = times, rootfun = rootfun)
  out_df = data.frame(out)
  for (t in 1:length(times)){
    out_df[t,'avg_class'] <- (sum(out_df[t,3:12]*1:10))/(sum(out_df[t,3:12]))
  }
  parameters_resis[i,'avg_class'] = out_df[nrow(out_df),length(out_df)]
}
# generate figure
ggplot(parameters_resis, aes(x=m_r,y=avg_class)) + geom_point() + stat_smooth(method = "loess")
