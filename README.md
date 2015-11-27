# Malaria Var-Gene Evolution-Transmission Model

## Requirements

* git
* Python 2.7.x (for build scripts)
* C++11-compliant compiler (Intel C++ compiler, LLVM/Clang, or GCC)
* Doxygen (to build code documentation)

## Getting and updating the soure code/Git primer

The most straightforward way to get the code is at the command line:

```sh
git clone git@bitbucket.org:pascualgroup/malariamodel.git varmodel
cd varmodel
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
./make.py
```

in the root directory of the repository. This will create the executable in

```
bin/malariamodel
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
[path-to]/bin/malariamodel parameters.json
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
visually see what's going on: [http://sourceforge.net/projects/mksqlite/](http://sqlitebrowser.org)

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

* `genes`: a table of all var genes, including characteristics that affect dynamics
* `hosts`: a table of hosts, including birth and death time
* `strains`: a table of all strains generated in the simulation
* `sampledHosts`: a list of hosts sampled at different sampling times
* `sampledHostInfections`: a list of all the infections of sampled hosts, including progression of var expression
* `sampledHostImmunity`: a list of immunity held by all sampled hosts
* `sampledTransmissions`: a list of sampled transmission events
* `sampledTransmissionInfections`: infection state of hosts involved in sampled transmission events
* `sampledTransmissionImmunity`: immunity state of hosts involved in sampled transmission events

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
			* Modify strains to be transmitted via mutation and/or recombination, possibly with extra circulating strains.
			* Choose destination host, and transmit strains into host
		* Immigration (introduction) events, with rate `currentIntroductionRate * currentPopulationSize`, where `currentIntroductionRate` is a sinusoidal function of time.
		* Birth events (if birth-death processes are uncoupled; demography details to be worked out)
		* Death events (if birth-death processes are uncoupled)
		* Within-host events (see below)

Individual var genes may be assigned different attributes:

* transmissibility
* mean duration of infection
* duration of immunity

### Transmission process

* Calculate P(transmission) of each active strain based on currently expressed var. To begin with, the probability of transmission will be inversely proportional to the number of concurrent infections.
* Select strains to be transmitted based on their probability of transmission.
* Randomly mutate transmitted strains: each gene in each strain has a constant probability of being sampled from the pool.
* Reassort var genes between strains:
	* Select pairs of strains to recombine, so that each pair has a constant probability of recombining. For each chosen pair, generate a daughter strain as a random reassortment of parent strains' genes.
	* If `n` strains were picked up from source host, transmit `n` strains randomly sampled from the original strains and all generated daughter strains.


### Within-host dynamics

The within-host dynamics for a single infection look like this:

```
[LIVER STAGE] -> [VAR A INACTIVE] -> [VAR A ACTIVE] -> [VAR B INACTIVE] -> [VAR B ACTIVE] -> ...
```


#### Liver Stage

The liver stage is always a fixed time period. During that period, the infection cannot be cleared.

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

#### Clearance

Currently, clearance--which ends the infection course completely--can happen only when a var gene is active (expressed). The clearance rate varies with immunity to the var gene being expressed.

If immune:

```
clearanceRateImmune = clearanceRateConstantImmune * nActiveInfections^clearanceRatePower
```

If not immune:

```
clearanceRateNotImmune = clearanceRateConstantNotImmune * nActiveInfections^clearanceRatePower
```

Presumably it would make sense for `clearanceRateConstantImmune` to be much larger than `clearanceRateConstantNotImmune`.


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