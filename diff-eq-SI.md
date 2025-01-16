ODE Analysis of Pathogen Fitness and R0
——updating markdown in progress!
================
Jess Valiarovski
2025-01-16

## Background

Many environments within ecological communities are hotbeds for evolving
resistant and dangerous pathogens. It has been observed that host
population traits influence the outcome of an emerging disease outbreak.
For example, in the avian flu, there are species of migratory birds that
are asymptomatic carriers who are tolerant, but spread it to other
species that are highly vulnerable, causing rapid population decline and
ecological instability.

The two primary traits that characterize different strains on an
epidemiological scale are transmission and virulence. Our focus will be
on micro parasites as macro parasites exhibit survival and reproduction
outside the host. All pathogenic strains have a level of virulence which
is generalized here defined by infection mortality rate. Virulence is a
phenomenon that arises as a consequence of a pathogen multiplying,
spreading within a host, and transmitting to new host contacts. In the
case of avian flu, a higher virulent strain is prevalent because
pathogen population is maintained as a reservoir in a tolerant host
species, but is a detriment to other species.

## Research Question

The research question I explored in Allison Shaw’s lab with Martha Torstenson and Allison Shaw at the University
of Minnesota Department of Ecology, Evolution, and Behavior is whether
the host adaptations tolerance and resistance have an effect on the
optimal pathogen virulence and how effective is the virulence strategy
in maximizing infection spread. By using a theoretical model, the
systemic dynamics and interplay of disease transmission become apparent.
These insights can be used to inform ecological host population strategy
effectiveness in minimizing pathogen harm. This model extends off the work of Martha Torestenson's in the below citation:

Torstenson, M., & Shaw, A. K. (2024). Pathogen evolution following spillover from a resident to a migrant host population depends on interactions between host pace of life and tolerance to infection. Journal of Animal Ecology, 93, 475–487. https://doi.org/10.1111/1365-2656.14075

## Model Overview

My theoretical model is based on a common mathematical model used the in
the epidemiology field called a susceptible infected ordinary
differential equation (SI ODE) model which depicts how biological
processes such as transmission, infection recovery, and mortality,
affect a host population that is divided into two subgroups: susceptible
(S) or currently infected (I). There are near infinite amount of
ecological and evolutionary strategy combinations that host species use
to combat off disease. The two notable strategies I will model are host
tolerance to the infection (by reducing infection-related mortality
rate) and host resistance (reduce the rate that an individual becomes
infected upon exposure to a pathogen). In my model, host and pathogen
traits mutually contribute to the transmission and mortality processes
which shift population abundances over time. My model assumes that
virulence increases transmission rate, where more virulent pathogens
have the potential to be evolutionary favored long term, which is not
the case in the classic SI epidemiology mathematical models. I have
encoded a transmission mortality trade off in my model, so it is not
immediately evident what virulence phenotype is most competitive to
exist in a host population with different defensive adaptations.

## Realized Transmission and Mortality Rates

To capture the host and pathogen contributions to disease transmission
and mortality, we denote a new rate called realized transmission rate
$\beta(\alpha_z)$ and realized mortality rate $\mu_I$:

$$
\beta(\alpha_z) = \alpha_z \beta_z \cdot (1 - m_R)
$$

Realized transmission is a function of virulence $\alpha_z$ and scalar
term $m_R$ as the susceptible host’s resistance to new infection.

$$
\mu_I = (\sqrt{\mu_S} + (1 - m_T) \left(\frac{\beta(\alpha_z)}{\beta_{\text{max}}}\right))^2
$$

Realized mortality rate of infected hosts by a pathogen with virulence
strategy $z$. Host tolerance $m_T$ decreases mortality rate by
infection. Mortality is parabolic, where increasing transmission is
traded off with higher infection-mortality rates. 

## ODE Model

The ODE model functions to predict host population abundances of 11
cohorts (susceptibles and 10 infected classes by virulence alpha) over
time until stability is reached.

### Rate of Change of Susceptible Individuals

$$
\frac{dS}{dt} = fN - \mu_S S + \gamma I - \sum_{z=0}^{10} \frac{\beta(\alpha_z) SI}{N^q}
$$

### Rate of Change of Infected Individuals

$$
\frac{dI_z}{dt} = \frac{\beta(\alpha_z) SI}{N^q} - \mu_I I - \gamma I
$$

### Average Pathogen Strategy

$$
V(X) = \frac{\sum_{z=0}^{10} \alpha_z I_z}{\sum_{z=0}^{10} I_z}
$$

The average pathogen strategy is calculated as the sum of individuals
infected by each pathogen class and the alpha value associated with it,
divided by the total infected individuals in the population.

## Symbol Definitions

| Symbol               | Meaning                                        | Numerical Range                          |
|----------------------|------------------------------------------------|------------------------------------------|
| $S$                  | Host variable: Susceptible population count    | Any integer                              |
| $I_z$                | Host variable: Infected population by $\alpha$ | Any integer                              |
| $t$                  | Time in days                                   | $t \geq 0$                               |
| $N$                  | Total number of individuals in the population  | Varies                                   |
| $q$                  | Host population density effect on transmission | $0 \leq q \leq 1$                        |
| $\mu_S$              | Host parameter: Natural mortality rate         | —                                        |
| $\alpha_z$           | Virulence factor                               | $0.1 \leq \alpha \leq 1$ in steps of 0.1 |
| $z$                  | Pathogen strategy location in parameter space  | Integers $0$ to $10$                     |
| $f$                  | Per capita fecundity rate                      | —                                        |
| $m_T$                | Host adaptive tolerance to infection mortality | $m_T \leq 1$                             |
| $m_R$                | Host adaptive resistance to new infections     | $m_R \leq 1$                             |
| $\mu_I$              | Realized mortality rate induced by $\alpha$    | Quadratic                                |
| $\beta$              | Base pathogen transmission rate                | Individual $\cdot$ time$^{-1}$           |
| $\beta_{\text{max}}$ | Maximum instantaneous transmission rate        | —                                        |
| $\gamma$             | Recovery rate                                  | $\gamma > 0$                             |
| $V(X)$               | Mean pathogen strategy                         | Evolves                                  |

### ODE Analysis
Figures 1-3 in the github repository are parameter sweep simulations of the equations above.

The $R_0$ equation divides the favorable pathogen rates (that allow it
to proliferate) by the unfavorable/trade-off rates, giving you a
standarized pathogen fitness value.

If $R_0$ is greater than 1 ($R_0$\>1) then the strain will proliferate
and induce an epidemic state in the host population; if $R_0$ is one
($R_0$=1), then it’s at steady state; if $R_0$ is less than 1 then the
strain will eventually go extinct

$$ R_0 = \frac{\beta(\alpha)}{\mu_I + \gamma} \to R_0 = \frac{\alpha_Z \beta^*(1-m_R) S}{\sqrt{\mu_S + (1-m_T)(\frac{\alpha_Z \beta^*(1-m_R)}{\beta_{\text{max}}})^2} + \gamma} $$

Where: - $\beta(\alpha)$: Transmission rate as a function of $\alpha$

- $\mu_I$: Infected mortality rate

- $\gamma$: Recovery rate

## ODE analysis Goals

1.  Investigate the conditions under which $R_0 > 1$, $R_0 = 0$, or
    $R_0 < 1$.
2.  Determine the optimal virulence $\alpha$ value by finding the
    derivative of $R_0$ and locating the global maximum where
    $\frac{d}{d\alpha} = 0$.

## Methodology

### Partial Derivative of R0 with Respect to Alpha

We independently vary host parameters (e.g., tolerance and resistance)
to examine the derivative of alpha (represented as A below) and host
trait constants (represented as H below) to determine \## Partial
Derivative of $R_0$ with Respect to $\alpha$

We independently vary host parameters tolerance T and resistance R to
examine the derivative of $\alpha$ (represented as $A$) and host trait
constants R and T.

1.  **When Tolerance is Zero ($m_r = R$):**

``` r
library(Deriv)

# Define R0 for tolerance
R0_tolerance <- expression(
  alpha * 0.4 * (1 - R) / 
  (sqrt(0.1 + alpha * 0.4 * (1 - R) / 0.4)^2 + 0.01)
)

# Compute the partial derivative
partial_tolerance <- Deriv(R0_tolerance, "alpha")
print(partial_tolerance)
```

    ## expression({
    ##     .e1 <- 1 - R
    ##     .e2 <- alpha * .e1
    ##     .e3 <- 0.11 + .e2
    ##     (0.4 - 0.4 * (.e2/.e3)) * .e1/.e3
    ## })

$$
\left( 0.4 - 0.4 \cdot \frac{A \cdot (1 - R)}{0.11 + A \cdot (1 - R)} \right) \cdot \frac{1 - R}{0.11 + A \cdot (1 - R)}
$$ 2. **When Resistance is Zero ($m_t = T$):**

``` r
# Define R0 for resistance
R0_resistance <- expression(
  alpha * 0.4 / 
  (sqrt(0.1 + (1 - T) * alpha * 0.4 / 0.4)^2 + 0.01)
)

# Compute the partial derivative
partial_resistance <- Deriv(R0_resistance, "alpha")
print(partial_resistance)
```

    ## expression({
    ##     .e2 <- alpha * (1 - T)
    ##     .e3 <- 0.11 + .e2
    ##     (0.4 - 0.4 * (.e2/.e3))/.e3
    ## })

$$
\left( 0.4 - 0.4 \cdot \frac{\alpha \ \cdot (1 - T)}{0.11 + \alpha \ \cdot (1 - T)} \right) \cdot \frac{1 - T}{0.11 + \alpha \ \cdot (1 - T)}
$$ They’re both the same result! The derivatives of resistance and
tolerance are identical equations because the host traits in the model
are constants, but the basic reproductive $R_0$ equations when virulence
level is varied will be different because those same constants affect
the trajectory curve differently when virulence level alpha is varied.
