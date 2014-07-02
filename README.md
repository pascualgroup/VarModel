# Malaria Var-Gene Evolution-Transmission Model

## Getting and updating the soure code/Git primer

The most straightforward way to get the code is at the command line:

```sh
git clone git@bitbucket.org:pascualgroup/malariamodel.git malariamodel
cd malariamodel
git submodule init
git submodule update
```

The `submodule` commands load three external git repositories into this one: the
[Catch](https://github.com/philsquared/Catch) test framework; 
[zppdata](https://bitbucket.org/edbaskerville/zppdata), which contains some JSON and SQLite data management utilities; and
[zppsim](https://bitbucket.org/edbaskerville/zppsim), which contains random number utilities and a general-purpose event queue
for stochastic simulations.

To get a new version of the code, use `git pull` and then `git submodule update` to make sure submodules are updated if necessary:

```sh
cd [path-to]/malariamodel
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

(Of course, Bitbucket also has a graphical view of commit history.)
