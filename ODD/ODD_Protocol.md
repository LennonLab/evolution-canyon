#ODD Protocol
The ODD protocol is standard for describing an individual-based model (IBMs). 
We descibe the IBM of evolution canyon as close as reasonable according to the 
ODD protocol.

Grimm, V. *et al*. (2006) A standard protocol for describing individual-based 
and agent-based models. Ecological Modeling. **198,** 115-126.ODD-based model descriptions consist of the following seven elements. In most 
cases it is necessary to have a simulation experiment or model analysis section 
following the model description.## PurposeThe purpose of the evolution canyon IBM is to simulate life history and dispersal 
of individual microorganisms and the assembly of microbial communities within a 
spatially explicit landscape that is divided into two equal-sized halves. 
Each half has its own largely unique and heterogeneous habitat.
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
	

##Spatial and temporal scales
The two general aspects of scale are grain (a.k.a. 
resolution) and extent (e.g. total area).

###Space
The landscape of EC-IBM is a two dimensional quadrangle with an arbitary
extent. Each habitat is square in shape and occupies one side of the landscape. 
While the extent of the landscape is arbitrary, it is overlaid with a lattice
of 1,000 equally spaced points. Each of the 1,000 points is assigned two
values representing two environmental conditions.

Individuals move in decimal units the limit of which is determined by Python's 
decimal precision. This means that individual particles can occupy 
practically infinite locations within the landscape, and hence, a minimal
spatial grain that is, for all intents and purposes, infinitely small.

###Time
EC-IBM models run for a temporal extent of 10 million time steps. Each time 
step (temporal grain) is equivalent to one complete set of birth, death, 
immigration, dispersal, and activity/dormancy transition events. This yields
a total of 50 million life history events.
## Process overview and scheduling### Initiation
The user runs a program (ECmodel.py or ECmodel_MultiProc.py) that initializes
the landscape, including the values of environmental variables for each of 
the 1,000 points in the landscape lattice. The script also initializes the
community with one million individuals chosen at random from a log-series 
distributed metacommunity.

### Core simulation processes
The EC-IBM then simulates the life history events of immigration, death, 
reproduction, dispersal, and the transition to or from dormancy across 
10 million time steps. For each event, an individual is chosen at random 
from the community.  
**Inflow/Entrance:** Resources and individuals enter from any point in the 
environment. Species identities of inflowing propagules are chosen at random 
from a log-series distribution, which often approximates the distribution of 
abundance among species in ecological communities (see Hubbell 2001). Along 
with a species ID and species-specific maximum rates of resource uptake, active 
dispersal, resource efficiencies, cell maintenance, and environmental optima, 
each propagule was given a unique ID and a multi-resource cell quota that 
represent the state of internal resources. The average of these cell quotas 
determine the probability of reproduction. 

**Dispersal:** The direction of dispersal in the EC-IBM is random. However, the distance 
traveled along the x and y axes are random draws from a Gaussian distribution. 
This distance is chosen at random from a species-neutral dispersal kernel. 
This dispersal kernel is a two-dimensional representation of the Guassian 
(normal) distribution, which has two parameters: location and scale. The location 
parameter is simply the individual's x-y coordinate and the scale parameter is an 
arbitrarily chosen value representing the standard deviation around the location. 
The scale parameter determines how far an individual is likely to disperse, and 
is the same for all individuals.  
**Reproduction:** Reproduction occurs clonally. Individuals reproduce with a 
probability determined by their match to the environmental optima. In this way,
the match can range between 0.0 and 1.0, and the relationship of the match (m)
to the probability of reproducing (p) is: m = p.  

**Death:** Individuals die and reproduce with a probability either equal to 0.5
or with a probability determined by the match of their environmental optima to 
the environmental conditions at their location. In this way, the match can range 
between 0.0 and 1.0, and the relationship of the match (m) to the probability
of dying is: 1-m = p.

**transition to and from dormancy:** Depending the on model, individuals go 
dormant or become active with a probability either equal to 0.5 or with a 
probability determined by the match of their environmental optima to the
environmental conditions at their location. The relationship of the match (m) 
to the probability (p) of going active is: p = m, while the relationship of m to
the probability (p) of going dormant is: p = 1 - m.
## Design conceptsThe ODD protocol calls for eleven design concepts. The IBM of Evolution Canyon 
adheres to most of these.### Basic principles. 
**Concepts**  

* Stochastic and environmentally induced dormancy
* Environmental filtering
* Process-based stochasticity

**Theories**
EC-IBM does not incorporate an ecological theory *per se* but 
draws from processes such as dispersal, birth, and death that 
are central to ecological theories such as metacommunity theory
and ecological neutral theory.

**Hypotheses**
EC-IBM was built to test hypotheses relating to how dormancy
and environmental filtering could have led to community compositional
differences observed in the microbial communities at Evolution Canyon
in Israel. These hypotheses are laid out in the main manuscript.

**Modeling approaches**
* Random sampling
* Probabilistic life history
* Spatially explicit and heterogeneous landscape
* Stochastic switching between activity and dormancy
### Emergence
EC-IBM uses random sampling to avoid imposing deterministic 
outcomes community assembly and structure
###Adaptation
EC-IBM does not include individual-level adaptation###Objectives
Individual are not coded to have explicit objectives, 
e.g., searching.
###Learning 
There is no aspect of individual-based learning in 
EC-IBM###Prediction 
Individuals in EC-IBM do not have the ability to 
anticipate conditions.###Sensing 
Individuals only sense in the sense that they can 
become dormant or active in response to the environment.
###Interaction 
EC-IBM includes no individual-level interactions.###Stochasticity 
The occurrence of nearly all processes of birth, death, life, 
immigration, dispersal, emigration, consumption, etc. are 
conducted via random sampling. In this way, population and 
community dynamics result, in part, from demographic stochasticity. 
### Collectives
Individuals belong to species. Species belong to communities.## Initialization
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