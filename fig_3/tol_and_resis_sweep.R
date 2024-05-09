# load libraries 
#install.packages("pacman")
pacman::p_load(deSolve,ggplot2,dplyr)
pacman::p_library() #verify packages are loaded

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

# initate simulation duration

times <- seq(0, 500, by = 1)


#initate simulation parameters 

alpha = seq(0.1,1,0.1) # virulence factor

#modify m_t and m_r increments for shorter run time
m_t = seq(0,1,0.1) # tolerance
m_r = seq(0,1,0.1) # resistance
parameter_sweep <- expand.grid(
  m_t = m_t, # tolerance
  m_r = m_r, # resistance
  beta = 0.4, # max transmission rate
  gamma = 0.01, # recovery rate
  miu_s = 0.1, # max mortality rate
  fecundity = 1 # fecundity rate
)


# host-pathogen interaction; defensive host adaptations reduce effect of pathogen transmission and infected mortality rate 

# realized transmission rate
beta_realized = matrix(data=NA, nrow=length(alpha),ncol=length(m_r)) 
for (i in 1:length(m_r)){
  for (z in 1:length(alpha)){
    beta_realized[z,i] <- alpha[z]*parameter_sweep$beta[i]*(1-m_r[i]) 
  }
}
# realized mortality rate
miu_i = matrix(data=NA, nrow=length(alpha),ncol=length(m_t)) 
for (i in 1:length(m_t)){
  for (z in 1:length(alpha)){
    miu_i[z,i] = (sqrt(parameter_sweep$miu_s[i])+ (1-m_t[i])*(beta_realized[z,i]/parameter_sweep$beta[i]))**2
  }
  
}  

#ordinary differential equation (ODE)

SIM_MODEL <- function(times, init, parameters) {
  with(as.list(c(init, parameters)),{
    # extract state variables from init vector
    S = as.vector(init[1]) # susceptible host subpop
    I = matrix(init[2:length(init)], nrow=10, ncol=1) # infected host subpop
    beta_realized = matrix(data=NA, nrow=10,ncol=1) # realized transmission rate
    alpha = seq(0.1,1,0.1) # virulence class
    for (z in 1:10){ # beta transmission is a function of alpha virulence factor
      beta_realized[z] <- alpha[z]*beta*(1-m_r) # realized transmission rate
      }
    miu_i = (sqrt(miu_s)+ (1-m_t)*(beta_realized/beta))**2 # realized mortality rate
    # system of differential equations for SI model
    dS <- -1*S*sum(I*beta_realized) - miu_s * S + sum(gamma * I) + fecundity * (S+sum(I))
    dI <- I * beta_realized * S - (gamma+miu_i) * I
    return(list(c(dS,dI)))
  })
}
rootfun <- function(times, init, parameters) {
  dstate <- unlist(SIM_MODEL(times, init, parameters))
  sum(abs(dstate)) - 1e-1 
}
# run simulation!
for (i in 1:nrow(parameter_sweep)){
  out <- lsodar(func = SIM_MODEL, y = init, parms = parameter_sweep[i,],
                times = times, rootfun = rootfun)
  out_df = data.frame(out)
  for (t in 1:length(times)){
    out_df[t,'avg_class'] <- (sum(out_df[t,3:12]*1:10))/(sum(out_df[t,3:12]))
  }
  parameter_sweep[i,'avg_class'] = out_df[nrow(out_df),length(out_df)]
}

#generate figure of simulation output

sp2 <- ggplot(parameter_sweep, aes(x=m_t,y=m_r, fill = avg_class)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF", #low investment in virulence=low instantaneous transmission rate BUT high life-time transmission (because decreased infection mortality rate)
                       mid = "#FFFFCC", #equal pathogen investment in virulence and life-time transmission
                       high = "#FF0000") #high investment in virulence=high instantaneous transmission rate BUT low life-time transmission rate (because increased infection mortality rate) 
  + coord_fixed()

sp2
