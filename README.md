# Malaria Var-Gene Evolution-Transmission Model

## Requirements

* git
* Python 2.7.x (for build scripts)
* C++11-compliant compiler (Intel C++ compiler, LLVM/Clang, or GCC)
* Doxygen (to build code documentation)

## Getting and updating the soure code/Git primer

The most straightforward way to get the code is at the command line:

```sh
git clone git@bitbucket.org:pascualgroup/malariamodel.git malariamodel
cd malariamodel
git submodule init
git submodule update
```

The `submodule` commands load code from some other git repositories that this code depends on.

To update to a new version of the code, do the following in the project directory:

```sh
git pull origin master
git submodule update
```

Always pull before working on changes to the code, or before pushing new changes.

If you have modified the code and want to commit locally, first use `git add` to tell git which files should be committed:

```sh
git add [file1] [file2] [file3] [...]
git add [file4]
```

etc. Then, commit like this:

```sh
git commit -m "Short description of changes"
```

You can have git automatically commit all changed files that are already in the repository using `-a`, but be careful with this as it might do something you don't intend:

```sh
git commit -a -m "Short description of changes"
```

You can always check on the status of which files have been added, modified, etc. using

```sh
git status
```

and you can view a history of commits using

```sh
git log
```

Bitbucket also has a graphical view of commit history.

## Building the executable

To build the executable, simply run

```{bash}
./make.py [varmodelName]
```

in the root directory of the repository. This will create the executable in

```
bin/[varmodelName]
```

The build script chooses a compiler based on availability, preferring the Intel C++ compiler if available, then LLVM/Clang, and finally GCC.

## Building code documentation

There is some code documentation, mostly for parameters and database tables, written in the Doxygen format. To build the code documentation, you need to install [Doxygen](http://www.stack.nl/~dimitri/doxygen/). On Mac OS X, you can install the Doxygen app (binary) in the Applications folder and the script will find it. Otherwise, e.g. on Linux, make sure the `doxygen` command is in your path--if you install using your package
manager, then this should happen automatically.

You can build code documentation by running

```{bash}
./makedoc.py
```

in the root directory. You can then view `doxygen/index.html` in your web browser.
Descriptions of the parameter and database types are linked from the index page.

## Running the model

The model requires a parameters file, in JSON format (see the next section). The recommended way to run the model is to create a working directory containing just, e.g., `parameters.json`, and running the executable in that directory:

```{bash}
cd [path-to]/output
[path-to]/bin/[varmodelName] parameters.json
```

It's often useful to redirect standard error and standard output into files for later examination:

```{bash}
cd [path-to]/output
[path-to]/bin/malariamodel parameters.json 1> stdout.txt 2> stderr.txt
```

This works well for parameter sweeps if a different directory is created for each model
run using a custom script, and also works in a straightforward way with [runm](https://github.com/edbaskerville/runm).

Each run will create a single SQLite database as output, whose filename is specified in the parameters file (see "Output database" below).

## Parameters file

The parameters file contains all output control and model parameters in [JSON](http://json.org) format. In order to prevent confusion, the code contains no default values whatsoever: every parameter used in the simulation must be present in the parameters file, or the program will throw an exception.

Parameters are defined and documented in the file `SimParameters.h`, using a strange-looking but simple macro format that allows symbol names in the defined type to be used directly as keys in the JSON file. E.g.,

```{c++}
ZPPJSON_DEFINE_TYPE(
	MySubJsonType,
	( (Double)(a) )
	( (Double)(b) )
)

ZPPJSON_DEFINE_TYPE(
	MyJsonType,
	( (Double)(x) )
	( (Double)(y) )
	( (MySubJsonType)(z) )
)
```

defines a JSON object type that can be directly mapped from this JSON:

```{javascript}
{
	"x" : 5.3,
	"y" : 3.8,
	"z" : {
		"a" : 100,
		"b" : 88.1
	}
}
```

Fields can be declared as type `Double`, `Int64`, `String`, `Array<T>`, where `T` is any allowed type; `Map<T>`, which maps strings to any allowed type; or can be declared using a custom type (e.g., `MySubJsonType` above).

This mechanism is implemented in the [zppjson](https://github.com/edbaskerville/zppjson) library.

The file `example/example_parameters.json` should currently contain all parameters defined in `SimParameters.h`, and should be kept in sync with that file to be a valid, working example. `SimParameters.h` also includes documentation on parameters in Doxygen format.

## Output database

The output database is in SQLite 3 format, which can be easily accessed from R using the `RSQLite` library or from Python using the built-in `sqlite3` library. In Matlab, the
[`mksqlite`](http://sourceforge.net/projects/mksqlite/) package does the trick: 

I also recommend using this graphical SQLite browser, especially while testing, to
visually see what's going on: 
(http://sqlitebrowser.org)

Output database tables are defined using a similar macro format to parameters, e.g.,

```{c++}
ZPPDB_DEFINE_ROW_TYPE(
	MyDbRowType,
	( (Integer)(penguinId) )
	( (Text)(name) )
	( (Real)(height) )
	( (Real)(weight) )
)
```

Database fields can be only `Integer` (64-bit integer), `Text` (string), or `Real` (double), with names taken from the SQLite conventions for database types.

Current database tables include:

* `genes`: a table of all var genes, including characteristics that affect dynamics, transmissibility, immunityLossRate, source (initial pool = 0, mutation = 1, recombination = 2), functionality (a gene can be non-functional if it's formed through recombination) 
* `loci`: a table of all allele compositions of a specific geneid
* `hosts`: a table of hosts, including birth and death time
* `microsats`: a table of microsatellites allelic compositions
* `strains`: a table of all strains generated in the simulation
* `InfectionDuration`: a table of sampled duration of infection, sampled every 1000 infections
* `recordEIR`: a list of sampled mosquito bites, infectious = 0 means the donor host doesn't have infections, infectious = 1, means it's an infectious bites
* `sampledHosts`: number of hosts sampled at different sampling times to achieve enough hosts with infections, set by parameter: populations.sampleSize. 
* `sampledHostInfections`: a list of all the infections of sampled hosts, including progression of var expression
* `sampledHostImmunity`: a list of immunity held by all sampled hosts (turned off by current implementation)
* `sampledTransmissions`: a list of sampled transmission events (turned off by current implementation)
* `sampledTransmissionInfections`: infection state of hosts involved in sampled transmission events (turned off by current implementation)
* `sampledTransmissionImmunity`: immunity state of hosts involved in sampled transmission events (turned off by current implementation)

(and corresponding tables for "clinical immunity").

Hosts are sampled every `sampleHostsEvery` units of simulation time.

Transmission events are sampled every `sampleTransmissionEventEvery` *transmission events*: so if `sampleTransmissionEventEvery = 100`, the 100th, 200th, 300th, ...
transmission events will be sampled.

## Code organization

The code is object-oriented and mostly hierarchical. At the top level is an instance of the `Simulation` class, which contains collections of `Gene`, `Strain`, and `Population` objects. Each `Population` object contains a collection of `Host` objects, each of which contains current infections and immune history through `Infection` and `ImmuneHistory` objects. In other words, the simulation hierarchy is organized like this:

```
Simulation
	Gene objects
	Strain objects
	Population objects
		Host objects
			Infection objects
			ImmuneHistory objects
```

The relationships are not strictly hierarchical, since `Infection` and `ImmuneHistory` contain references to strains and genes, and during the simulation, interactions happen at multiple levels.

Interactions are generally mediated using member function calls, but in some cases internal state may be accessed directly. That is, data hiding/encapsulation is not strictly enforced.

Raw pointers are used in a number of places in the code. The code would use a bit more
memory, but would be safer going forward, if they were entirely replaced with
`shared_ptr` or `unique_ptr` in all cases.

## Event-based simulation architecture

The simulation is implemented using a global event queue, which contains lightweight event objects that generally do nothing but call a member function on a core simulation object (e.g., perform some action on an `Infection`, `Host`, or `Population`). These events include contact events (paired biting events); immigration; loss of immunity; transitions between infection states; host death; periodic updates to contact rates based on a seasonal function; and meta-operations such as host sampling.

The underlying event queue implementation is similar to the ["next-reaction method" of Gibson and Bruck](http://pubs.acs.org/doi/abs/10.1021/jp993732q), but dependencies among events are specified explicitly by asking the queue to update rates or add/remove events as part of the function that performs the event.

The event queue consists of a large number of events, each of which is associated with a "putative time". In the case of events scheduled to happen at a specific time (e.g., periodic updates of contact rates) or whose times are precalculated by drawing from a random distribution but do not depend on subsequent simulation events (e.g., death times), these times canot change. For events that are Poisson processes whose rates may change as a result of other simulation events, these times really are "putative" and may be recalculated.

In order to facilitate fast updating of changing-rate events, the underlying data structure is an indexed binary heap, so that all operations are `O(log(N))`, where `N` is the number of events on the heap. That is, the number of events that can be executed per second will go down as the size of the simulation goes up, but it will go up only logarithmically, and thus total simulation time should be `O(N log(N))`.

Abstractly, the simulation proceeds by repeating these two steps:
1. Remove next event from event queue.
2. Perform event.

The first step is handled by the underlying event queue implementation: the priority heap always has the lowest-time event ready to be accessed. The second step may include requests to the event queue to update rates for certain events, add events, or remove events.

## Simulation Details

### Simulation loop

A simulation follows these steps:

* Load all parameters from parameters file.
* Initialize population size with `initialPopulationSize` hosts. Sample `nInitialStrains` strains and assign each of them to `nInitialInfections` hosts.
* Set `t = 0`.
* Repeat until `t >= tEnd`:
	* Choose and perform the next event via the Gillespie algorithm (or probabilistic equivalent, e.g. next-reaction method), either:
		* Biting events, with rate `currentBitingRate * currentPopulationSize`, where `currentBitingRate` is a sinusoidal function of time:
			* Choose a random transmitting host.
			* Calculate which strains should be transmitted from the transmitting host (see description below).
			* Modify strains to be transmitted via sexual recombination, possibly with extra circulating strains.
			* Choose destination host, and transmit strains into host
		* Immigration (introduction) events, with rate `immigrationRate`, 
		* Birth events (if birth-death processes are uncoupled), currently assuming a demography that is a exponential distribution, with mean life expectancy = 30; 
		* Death events (if birth-death processes are uncoupled)
		* Within-host events (see below)

Individual var genes may be assigned different attributes:

* functionality
* transmissibility
* mean duration of infection
* duration of immunity

### Transmission process

* Calculate P(transmission) of each active strain based on currently expressed var. To begin with, the probability of transmission will be inversely proportional to the number of concurrent infections.
* Select strains to be transmitted based on their probability of transmission.
* Reassort var genes between strains:
	* If `n` strains were picked up from source host, transmit `n` strains randomly sampled from the original strains and all generated daughter strains.
    * percentage of recombined strain is not determined by a parameter, but rather determined by the number of co-infected strains:
```
    pRecombinant = 1.0-(1.0/double(originalStrains.size()))
```
    * Select pairs of strains to recombine, so that each pair has a constant probability of recombining. For each chosen pair, generate a daughter strain as a random reassortment of parent strains' genes.
    * Shuffle genomes of microsatellites with the same fashion if microsattellites are also simulated. There is a big difference though: each microsatellite is a gene with a fixed location, when two genomes exchange, they exchange by position.
* each host has a `maxMOI` parameter that sets the upper bound of how many co-infection each host can take.

### Within-host dynamics

The within-host dynamics for a single infection look like this:

```
[LIVER STAGE] -> [VAR A INACTIVE] -> [VAR A ACTIVE] -> [VAR B INACTIVE] -> [VAR B ACTIVE] -> ...
```

### Selection Mode `selectionMode`

Three type of models (i.e., within-host dynamics) can be run in this implementation:
1. immune selection (hosts build immunity towards specific epitope alleles or genes)
2. generalized immunity (i.e., clearance is faster if the host is exposed to more times of infection)
3. complete neutrality (i.e., clearance is determined by a clearance rate, and is not influenced by infection histories)

#### Liver Stage

The liver stage is always a fixed time period. During that period, the infection cannot be cleared. This stage actually also include the time required in the mosquito midgut to develop to the asexual transmission stage, which is about 7 days. Therefore usually this parameter is set at 14 days.

#### Activation

During inactive stages, the activation rate is a function of the number of active infections in the host:

activationRate = activationRateConstant * nActiveInfections^activationRatePower

where activationRatePower <= 0.

If activationRatePower == 0, then the activation rate is a constant.

If activationRatePower < 0, then the activation rate is a monotonically decreasing function of nActiveInfections, with activationRate = infinity (immediate activation) if there are no active infections.

#### Deactivation

During active (expressed) stages, the deactivation rate is similarly a function of the number of active infections in the host, including this one:

```
deactivationRate = deactivationRateConstant * nActiveInfections^deactivationRatePower
```

where deactivationRatePower may be positive, negative, or zero.

Not understanding the biology fully, I suspect `deactivationRatePower = 0` (constant deactivation rate) is probably the most sensible?


If selection mode is 1 (immunity selection):
    deactivationRateConstant governs the rate when hosts are not immune to the gene
    if `useAlleleImmunity = true`, hosts build immunity to specific alleles within genes,
        deactivationRate per gene is determined by function `Infection::immuneClearRate`,increases with the number of alleles in the gene hosts are immuned to
    if `useAlleleImmunity = false`, hosts build immunity to a specific gene,
        deactivationRate per gene is 1 if hosts are immune to the gene

if selection mode is 2 or 3, deactivationRate doesn't change

#### Clearance

Currently, clearance--which ends the infection course completely--can happen only when a var gene is active (expressed). The clearance rate varies with the selection mode.

If selection mode is 1 (immunity selection) or 3 (complete neutrality):
    there is no global clearance rate, the time to clearance is determined by cumulative time of expression of all the genes in the genome
    therefore `clearanceRateConstant` is set to 0

If selection mode is 2 (generalized immunity):
    clearanceRateConstant is determined by the function `immunity.checkGeneralImmunity`, which has parameters: `clearanceRateConstantImmune`,`infectionTimesToImmune`
    and generalize immunity fitting function parameters


```
clearanceRate = clearanceRateConstant * nActiveInfections^clearanceRatePower
```
#### Mutation
During the course of infection in the host, strains can mutate, and is governed by the rate 
```
pMutation * locusNumber * genesPerStrain
```
if a gene is selected to mutate, one of its epitope allele will mutate and have a new ID associated. this will create new ID for a gene, and a new ID for a strain,
infectionItr will be redirected to the new strain pointer. related function: Simulation::mutateStrain, Host::hstMutateStrain

Similar microsatellite mutation events is also included. 

#### ectopic recombination
During the course of infection in the host, genes within a strain can exchange their alleles, which is called "ectopic recombination", and is governed by the rate 
```
pIntraRecomb*numGenes * (numGenes -1)/2
```
The more nubmer of genes per genome, the more frequent such recombination occurs. functions: `Host::RecombineStrain`, `Simulation::ectopicRecStrain`
Ectopic recombination is to create a breakpoint between Var gene A and gene B, for example:
Gene A: a1 b1 c1
Gene B: a2 b2 c2
if the breakpoint is between a and b, then resulting daughter genes are
a1 b2 c2 and a2 b1 c1
if any of the daughter is non-functional with a probability proportional to their distance to the parental genome, then the strain will still keep the parental genes
Otherwise, the parental genes will be replaced by daughter genes, and strain ID will be changed, infectionItr will be redirected to the new strain pointer.

These events do not occur in microsatellites

Both mutation and ectopic recombination occur during the infection cycle, therefore the longer they stay in a host, the more likely they'll mutate into a different strain.

### interventions
Currently, interventions can by implemented as `interventions` (refers to IRS) or `MDA`
IRS is a one time event, set with an initial time and a duration (which is another one time event to remove the IRS implementation).
During IRS, biting rate and/or immigration rate can be changed through an mulplitication with the amplitude

MDA is a period event with an interval of `interval` and a maximum count of events = "totalNumber", i.e., after `totalNumber` of implementation, MDA will be removed.
During MDA, hosts' infections will be cleared, unless the host fails to take the drugs effectively (`hostFailRate`) or the strains are resistant (`strainFailRate`)
The effectiveness of drugs are set by `drufEffDuration`.
During MDA, migration rate can also be rescaled to be `immigrationRate * MDAMRateChange`
Note that currently, all hosts whose age is under 3 months are not given MDA.

## Modifying transmission and within-host dynamics

The events related to within-host infection dynamics are handled by the `Infection` class, implemented in the `Infection.cpp` file.

In order to modify the assumptions of the dynamics, it should be sufficient to modify
the rate calculations in:
* `Infection::transitionRate()`
* `Infection::activationRate()`
* `Infection::deactivationRate()`
* `Infection::clearanceRate()`

There are commented-out lines of code in the `activationRate()` function showing how to
retrieve various state variables (host age, number of active infections, etc.) that might
be useful for different assumptions of the dynamics.

If significant changes are desired, it will probably be necessary to modify
`SimParameters.h` as well.

The calculation for transmission probability can similarly be modified:
* `Infection::transmissionProbability()`

## following hosts

Specific hosts can be followed and all of their infections will be recored for a period of time.
This is governed by following.HostNumber and folloinwg.followDuration

## save database after burn-in
`burnIn` governs when the database will start sampling hosts, and writing genes, loci and strains, etc.
*`genes` and `strains` classes have a bool variable to record whether this gene or strain has been written to the sqlite


