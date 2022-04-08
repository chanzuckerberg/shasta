# shasta/conf directory

Most Shasta users don't need to use the contents of this directory. Instead, use the following Shasta commands:

* `shasta --command listConfigurations` to get a list of available built-in configuration names.
* `shasta --command listConfiguration --config configurationName` to list the details of a specific configuration. Its output will include comments describing the intended usage of the specified configuration.

See [here](../docs/Configurations.html) for more information about Shasta assembly configurations. See below for more detailed information about the contents of this directory. Most of this information is not necessary or useful to most Shasta users.



## Shasta configuration files

The Shasta `--config` option accepts either a built-in configuration name (one of the names listed by `shasta --command listConfigurations`) or a Shasta configuration file. In most cases, a built-in configuration name is used.


A Shasta configuration file can be created by a user for a specific application or can selected from this directory. This directory includes a collection of Shasta configuration file, including ones that are now obsolete. To get a list of currently relevant assembly configurations, use `shasta --command listConfigurations`.

Each of the built-in configurations also has a corresponding configuration file in this directory. There is no reason for any user ever to use one of them.


## Bayesian models for repeat counts

Files named `SimpleBayesianConsensusCaller-X.csv`, where `X` is a number, describe Bayesian models used by Shasta to make probabilistic calls on the repeat count at each base position, using as input the observed repeat counts. The first two lines of each file contain a succinct description of the Bayesian model contained in that file. 

In normal circumstances, these files are not used. Instead, the user specifies one of the built-in Bayesian models - for example `--Assembly.consensusCaller Bayesian:guppy-5.0.7-b`.




