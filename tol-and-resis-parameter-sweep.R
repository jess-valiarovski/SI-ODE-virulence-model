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
