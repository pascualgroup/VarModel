# Malaria Model

## Getting the soure code

The most straightforward way to do this is at the command line:

	git clone git@bitbucket.org:pascualgroup/malariamodel.git malariamodel
	cd malariamodel
	git submodule init
	git submodule update

The `submodule` commands load three external git repostiaries into this one: the
[Catch](https://github.com/philsquared/Catch) test framework; 
[zppdata](https://bitbucket.org/edbaskerville/zppdata), which contains some JSON and SQLite data management utilities; and
[zppsim](https://bitbucket.org/edbaskerville/zppsim), which contains random number utilities and a general-purpose event queue
for stochastic simulations.
