# Real Data for MMHP model

This repository contains the data of ten cohorts mice study from the paper Williamson(2016). Data was collected by the [Curley lab](https://labs.la.utexas.edu/curley/). All procedures were conducted with approval from the Columbia University Institutional Animal Care and Use Committee. If you use the dataset in published research, please cite the following reference:

Williamson, Cait M., Won Lee, and James P. Curley. "Temporal dynamics of social hierarchy formation and maintenance in male mice." Animal behaviour 115 (2016): 259-272.

`mice.RData` is constructed by a list of R dataframe. Each dataframe records the observations for each cohort. 

In each cohort, twelve male mice were placed in a large vivarium at the age of nine weeks. 
These mice fight each other to establish a social hierarchy. For twenty-one consecutive days, 
observations were taken during the dark phase of the light cycle when mice are most active. 
During each observation interval, trained observers recorded the location, time stamp and
action of pair-wise mice interactions. For each interaction, the actor (i.e., winner) and 
recipient (i.e., loser) were recorded, denoted as `Actor` and `Recipient` columns in the dataframe. 
