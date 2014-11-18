# Malaria Var-Gene Evolution-Transmission Model

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

## Building the code

TODO

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

## Event-based simulation architecture

The simulation is implemented using a global event queue, which contains lightweight event objects that generally do nothing but call a member function on a core simulation object (e.g., perform some action on an `Infection`, `Host`, or `Population`). These events include contact events (paired biting events); immigration; loss of immunity; transitions between infection states; host death; periodic updates to contact rates based on a seasonal function; and meta-operations such as host sampling.

The underlying event queue implementation is similar to the ["next-reaction method" of Gibson and Bruck](http://pubs.acs.org/doi/abs/10.1021/jp993732q), but dependencies among events are specified explicitly by asking the queue to update rates or add/remove events as part of the function that performs the event.

The event queue consists of a large number of events, each of which is associated with a "putative time". In the case of events scheduled to happen at a specific time (e.g., periodic updates of contact rates) or whose times are precalculated by drawing from a random distribution but do not depend on subsequent simulation events (e.g., death times), these times canot change. For events that are Poisson processes whose rates may change as a result of other simulation events, these times really are "putative" and may be recalculated.

In order to facilitate fast updating of changing-rate events, the underlying data structure is an indexed binary heap, so that all operations are `O(log(N))`, where `N` is the number of events on the heap. That is, the number of events that can be executed per second will go down as the size of the simulation goes up, but it will go up only logarithmically, and thus total simulation time should be `O(N log(N))`.

Abstractly, the simulation proceeds by repeating these two steps:
1. Remove next event from event queue.
2. Perform event.

The first step is handled by the underlying event queue implementation: the priority heap always has the lowest-time event ready to be accessed. The second step may include requests to the event queue to update rates for certain events, add events, or remove events.

## Detailed description of simulation

## Summary of C++ features used

