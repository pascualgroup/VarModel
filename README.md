# Malaria Model

## Getting the soure code

The most straightforward way to do this is at the command line:

<code><pre>
git clone git@bitbucket.org:pascualgroup/malariamodel.git malariamodel
cd malariamodel
git submodule init
git submodule update
</pre></code>

The <code>submodule</code> commands load three external git repostiaries into this one: the
<a href="https://github.com/philsquared/Catch">Catch</a> test framework; 
<a href="https://bitbucket.org/edbaskerville/zppdata">zppdata</a>, which contains some JSON and SQLite data management utilities; and
<a href="https://bitbucket.org/edbaskerville/zppsim">zppsim</a>, which contains random number utilities and a general-purpose event queue
for stochastic simulations.
