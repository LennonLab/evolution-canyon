#ODD Protocol
The ODD protocol is standard for describing an individual-based model (IBMs). 
We descibe the IBM of evolution canyon as close as reasonable according to the 
ODD protocol.

Grimm, V. *et al*. (2006) A standard protocol for describing individual-based 
and agent-based models. Ecological Modeling. **198,** 115-126.ODD-based model descriptions consist of the following seven elements. In most 
cases it is necessary to have a simulation experiment or model analysis section 
following the model description.## PurposeThe purpose of the evolution canyon IBM is to simulate life history and dispersal 
of individual microorganisms and the assembly of microbial communities within a 
spatially explicit landscape that is divided into two equal-sized and distinct 
halves. Each half has its own unique and heterogeneous habitat.
## Entities & their state variables
**Suggestions of the ODD protocol worth noting:** 
One way to define entities and state variables is the following: if you want to 
stop the model and save it in its current state, so it can be re-started later in 
exactly the same state, what kinds of information must you save? State variables 
should be low level or elementary in the sense that they cannot be calculated from 
other state variables.  

**Individual organisms** Individuals are represented as lists of elements. Likewise,
each individual is contained within a list corresponding to the side of the landscape
the individual is within. For example:

	sCOM[i] = [12, 3, 7]
	nCOM[i] = [12, 32, 11]
	
	The ith individual of the south-facing community belongs to species 12 and has an 
	environmental optima of 0.1 and is located at position x=3, y=7. Meanwhile, the 
	ith individual of the north-facing community belongs to species 12 and is located 
	at the x-y coordinates 32-11.

**Species** Species are characterized by individuals that share a common numeric ID 
(e.g., 12). The species ID then corresponds to specific environmental optima and 
dispersal kernel. This species-level information is stored separate from the list of
individuals (e.g. in python dictionaries) to decrease computational overhead. 

	Example:
	
	oDict1 = {1: 0.2; 12: 0.8}
	oDict2 = {1: 0.1; 12: 0.5}
	dDict  = {1: 0.4; 12: 0.9}
	
	Checking the species optima for the north-facing habitat (oDict1) we see
	that species 1 has an optima of 0.2, while species 12 has an optima of 0.8.
	
##System level state variables
The EC-IBM initiates with ...
	 

##Spatial and temporal scales
The two general aspects of scale are grain (a.k.a. 
resolution) and extent (e.g. total area).

###Space
The environment of EC-IBM is two dimensional and can vary along each axis from 
5 to 100 discrete units. This makes for a potential total **extent** of 25 to 
10,000 discrete patches, each with a **grain** of 1 square unit. 

Note, that all particles move in decimal units the limit of which is determined 
by Python's decimal precision. This means that individual particles can occupy 
practically infinite locations within patches and likewise, squeeze through barriers 
(which only occupy integer x-y coordinates).

###Time
**Extent** of time in EC-IBM models refers to residence time, i.e., the average 
amount of time that individual particles spend in the system. Residence time for 
inert tracer particles can vary across five orders of magnitude. 

**Grain** is the smallest unit over which change can happen. For example, as per 
the original ODD documentation: “One time step represents one year and simulations 
were run for 100 years. One grid cell represents 1 ha and the model landscape 
comprised 1,000 x 1,000 ha; i.e., 10,000 square kilometers”. In contrast, the 
smallest grain achievable by EC-IBM is determined by slowest rate at which individuals 
can undergo BIDE (birth, immigration, death, emigration) processes. For example, under 
high residence times, an individual can move across 0.00001% of the x or y axis in one 
time step. Under low residence times, an individual can move across 10% or greater of 
the x or y axis in one time step.
## Process overview and scheduling### Assembly
The user runs a program that chooses random values for system-level state variables 
including whether disturbance, immigration, speciation, fluid dynamics, etc. will occur 
and at what rates.

### Core simulation process
The EC-IBM begins simulation immediately after assembly from random combinations of state-variables and processes. Instead of operating by definitive time steps (i.e. days, generations), EC-IBM models advance turnover of the environmental matrix according to the initial rate of flow. If the initial rate of flow is 1.0, then the environmental matrix and inert particles would flow 1.0 units of distance. After each iteration of flow, each individual is given the chance to consume, grow, reproduce, starve, die, and to disperse towards resources 
and environmental optima.

### Duration: A run to mean reversion
Once assembled, a EC-IBM model simulates ecological processes (birth, death, dispersal, growth, consumption, etc.) until the system reaches a fluctuaing equilibrium determined by a point of mean reversion (quantified by Hurst's exponent). Mean reversion captures the tendency of a system to repeatedly reverse a directional change in, say, total abundance. The system examines whether a point of mean reversion has occured by recording the total abundance of the system each time a tracer particle exits the system. This ensures that, at the least, enough time 
has passed for an inert particle to enter and exit the system. 
### Fluid dynamics
EC-IBM uses an efficient and powerful method for simulating fluid flow, i.e., a Lattice-Boltzmann Method (LBM). An LBM discretizes the environment into a lattice and attaches to each position in the lattice a number of particle densities and velocities for each of nine directions of movement possible in a 2D environment (N, S, E, W, NE, NW, SE, SW, current position).###Active dispersalEC-IBM models allow individuals to move towards their environmental optima. Rather than a single environmental optima resulting from a single environmental gradient, EC-IBM allows environmental optima to occur as intersections among environmental gradients. Hence, 
individuals potentially have multiple optima resulting from unique and equally optimal intersection of up to 10 environmental gradients.###Simulated life history
EC-IBM models simulate growth, reproduction, and death via weighted random sampling. This simulate the partly probabilistic and partly deterministic nature of environmental filtering and individual-level interactions. 

**Inflow/Entrance:** Resources and individuals enter from any point in the environment. Species identities of inflowing propagules are chosen at random from a log-series distribution, which often approximates the distribution of abundance among species in ecological communities 
(see Hubbell 2001). Along with a species ID and species-specific maximum rates of resource uptake, active dispersal, resource efficiencies, cell maintenance, and environmental optima, each propagule was given a unique ID and a multi-resource cell quota that represent the state of internal resources. The average of these cell quotas determine the probability of reproduction. 

**Dispersal:** Individuals are allowed to actively move along environmental gradients (sometimes against the direction of environmenta flow) and towards their optima. A better match to one's environmental optima increases the chance of reproduction and the individual's ability 
to perform (consume, grow). Individuals are considered to have left or to have flowed out when they pass beyond edges of the environment.

**Reproduction:** Reproduction occurs clonally. Individuals reproduce with a probability determined by the combination of mean cell quota and the proportional match to the environmental optima. The cell quota of each resource is evenly divided between the individual and its daughter. The daughter is given a unique individual ID and the species ID of its mother, unless in the case of speciation, but is allowed small mutations in individual-level state variables.

**Death:** Individuals sampled at random will die if their smallest resource specific cell quota (i.e., N, C, P) is equal to or less than 0. 
## Design conceptsThe ODD protocol calls for eleven design concepts. The IBM of Evolution Canyon adheres to most of these.###Basic principles. 
**Concepts**  

* Stochastic and environmentally induced dormancy
* Environmental filtering
* Process-based stochasticity

**Theories**

* Metacommunity theory:
	*   
* Ecological neutral theory:  
  	* EC-IBM operates via random sampling and can vary from being completely neutral (all individuals having equal vital rates) to completely idiosyncratic (all individual and species are as different as possible). The one aspect of neutral theory that EC-IBM adopts without question is the importance of stochastic life history processes (i.e. weighted or unweighted random fluctuations in population sizes). 

**Hypotheses**


**Modeling approaches**

* Random sampling
* Stochastic switching
###Emergence
EC-IBM uses random sampling and random assembly
its models to avoid imposing strong constraints on
the properties that emerge and to allow unanticipated
combinations of traits and ecological structure to
emerge.

* Abundance
	* Total community abundance
	* Trait-related population size
	* Effective population size
	* Abundance-biomass relation
* Community assembly
	* Species turnover
	* Extinction rate
	* Succession
* Community structure
	* Species richness
	* Species evenness
	* Trait diversity
* Population structure
	* Demography
	* Subpopulation trait variation
###Adaptation
Individuals can move towards their environmental optima.
Populations can become aggregated in areas that provide 
favorable intersections of species optima. Species can 
evolve by the action of the environmental filter on 
subpopulation variation in state variables.###Objectives
Individual seek conditions that match them to the 
environment (e.g., positions along environmental 
gradients). Individuals also seek to acquire resources 
through active searching. In the future, individuals 
will seek to avoid predation.
###Learning 
There is no aspect of individual-based learning in 
EC-IBM###Prediction 
Individuals in EC-IBM do not have the ability to 
anticipate conditions.###Sensing 
Individuals only sense in the sense that they can become dormant or active in response to an en. Otherwise, all encounters are the result of random walks or fluid flow.###Interaction 
At the moment, individuals only interact indirectly through excluding each other from resources (e.g. preemption). 
In the future, individuals will interact as predator-prey, mutualists, resource-dependents, etc. 
Likewise, there is currently no communication, though quorum sensing would be cool.###Stochasticity 
The occurrence of nearly all processes of birth, death, life, immigration, dispersal, emigration, consumption, etc. are conducted via random sampling. In this way, population and community dynamics result, in part, from demographic stochasticity. Likewise, the emergence of life history 
traits proceeds from initially random combinations of traits.### Collectives
Individuals belong to species. Species belong to communities. In the future, EC-IBM will allow 
communities to belong to trophic levels.###Observation
Many EC-IBM models should be run to examine trends in the variation generated. The following is recorded for each EC-IBM model:

* Values of randomly chosen state variables
* Total abundance, $N$
* Species richness, $S$
* Compositional turnover
	* Bray-Curtis
	* Sorensen's
* Species turnover
	* Whittaker's $\beta$ 
* Species evenness
	* Smith and Wilson's evenness, $E_{var}$
	* Simpson's evenness, $E_{1/D}$ 
* Species diversity
	* Shannon's diversity, $H'$ 
	* Simpson's diversity, $D_{1/D}$
* Dominance
	* Absolute, $N_{max}$
	* Relative, $N_{max}/N$

These data are stored in file as R-formatted data.frames.
These files can be directly imported into an R or Python 
environment.## Initialization
The model initiates with a random set of values for state-
variables, 100 to 10,000 randomly drawn individuals from a 
theoretical log-series metacommunity.
These values are saved, so that a EC-IBM model could 
be programmed to replicate an analysis.
## Input dataThe EC-IBM requires no input data.## Submodels & Equations

### Log-series metacommunity
EC-IBM models draw immigrating individuals from a 
theoretical log-series distribution.
Hubbell (2001) states that the regional community 
(i.e., metacommmunity) often has a log-series species 
abundance distribution. The probability density 
function (pdf) of the log-series is:

$$f(k) = -1/ln(1 - p) * p^{k}/k$$

Hubbell (2001) provides explicit detail of the log-series, 
which is also covered in most ecological diversity texts
and even on Wikipedia: https://en.wikipedia.org/wiki/Logarithmic_distribution

**Reference:** Hubbell, S.P. (2001) The unified neutral 
theory of biodiversity and biogeography. Princeton 
University Press, Princeton, NJ.