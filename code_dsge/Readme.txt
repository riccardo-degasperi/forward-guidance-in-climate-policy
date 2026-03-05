Replication of Forward Guidance in Climate Policy
European Economic Review
Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara

LIST OF FILES

m.files
console.ss: sets the parameters and computes the initial steady state
console_ss_clean: sets the options of the simulation and computes the final steady state
plot_transition: plots the IRFs
saving_irf1: saves the IRFs of the baseline scenario
saving_irf2: saves the IRFs of the stance shock
saving_irf3: saves the IRFs of the path shock

m.functions
find_steady: sets the initial steady state system
find_steady_clean: sets the final steady state system

mat.files
x0: initial values for the initial steady state system
x1: initial values for the final steady state system

mod.files
transition_shock: the model

FOLDER
figures: where figures are stored


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGURES 1 and 9 
- Run console.ss


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGURES D1

- In console_ss_clean.m set:

SPEC=1; 
FAST=1;

- Run console_ss


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGURES D2

- In console.ss set

price=130;

- In console_ss_clean.m set:

SPEC=2; 
FAST=0;

- Run console:ss


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGURES D3

- In console.ss set

price=65;


- In console_ss_clean.m set:
SPEC=3; 
mu=0.5;

- Run console_ss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGURES D4

- In console_ss_clean.m set:

SPEC=4;
Zend=30; 
mu=1;

- Run console_ss