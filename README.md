# 🦠 **Virulence and Host Strategies: ODE Analysis of Pathogen Evolution**  

📅 **Project Timeline:** September 2023 - May 2024
📍 **Author:** **Jess Valiarovski**  
🎓 **Mentorship: **Martha Torstenson** and **Allison Shaw** 
---

## 📌 **1. Project Overview**  

How do host adaptations impact pathogen evolution? This project explores the relationship between **host resistance, tolerance, and pathogen virulence** using **Ordinary Differential Equations (ODEs)** to model infectious disease dynamics.  

### **🔍 Objectives**  
✅ Develop a **Susceptible-Infected (SI) model** that incorporates host adaptations (tolerance & resistance).  
✅ Use **differential equations** to simulate how virulence evolves over time.  
✅ Analyze how **pathogen reproductive number (R₀)** changes under different host conditions.  
✅ Compare theoretical predictions to **real-world epidemiological trends**.  

📊 **Primary Hypothesis:**  
If hosts invest in **tolerance**, pathogens evolve toward **higher virulence** to compensate for reduced host mortality. If hosts invest in **resistance**, pathogens face **stronger selective pressure**, leading to a **trade-off in virulence and transmission rates**.  

---

## 🛠 **2. Methods & Data Processing**  

### **📊 Epidemiological Modeling (ODEs)**  
This project implements a **Susceptible-Infected (SI) model**, where:  
- **Host resistance (mR)** reduces transmission probability.  
- **Host tolerance (mT)** reduces disease-induced mortality.  
- **Pathogen virulence (α)** determines trade-offs in transmission and host survival.  
- **Evolutionary outcomes** are analyzed by computing **R₀ as a function of virulence**.  

#### **Key Equations Modeled**  
- **Realized Transmission Rate:** β(α) = function of virulence  
- **Realized Mortality Rate:** μ(α) = quadratic trade-off function  
- **Evolutionary Fitness Criterion:** max R₀(α)  

### **📊 Simulation & Analysis in R**  
✅ **Numerical Simulations:** Simulated infection dynamics over time using **deSolve** in R.  
✅ **Analytical Approach:** Derived optimal virulence by solving **dR₀/dα = 0**.  
✅ **Parameter Sweeps:** Explored effects of **different tolerance & resistance values**.  
✅ **Visualization:** Used **ggplot2, deSolve, and Desmos** to analyze model behavior.  

---

## 📌 **3. Key Findings & Insights**  

📌 **Higher tolerance (mT)** leads to the evolution of **more virulent pathogens**.  
📌 **Higher resistance (mR)** leads to **lower transmission rates but does not lower virulence**.  
📌 The **optimal pathogen strategy** (α*) is the same under resistance and tolerance, but the effect on **R₀ differs**.  
📌 **Trade-offs:** Pathogens evolve to **maximize R₀**, balancing transmission benefits vs. host mortality costs.  
📌 These results challenge classical **trade-off theory** and suggest **different evolutionary pressures in resistant vs. tolerant populations**.  

### **📊 Key Figures**
1️⃣ **ODE Simulations:** Show how infection spreads in different host populations.  
2️⃣ **R₀ vs. Virulence Plots:** Identify the optimal virulence strategy under different host conditions.  
3️⃣ **Parameter Sensitivity Analysis:** Demonstrates how increasing tolerance or resistance shifts evolutionary outcomes.  

---

## 📌 **4. Technical Skills Demonstrated**  

✅ **Mathematical Modeling (ODEs, R₀ Analysis, Dynamical Systems)**  
✅ **Numerical Simulations (deSolve, parameter sweeps, R programming)**  
✅ **Data Visualization (ggplot2, Desmos, simulation outputs)**  
✅ **Evolutionary Epidemiology (Host-Pathogen Interactions, Virulence Evolution)**  
---




