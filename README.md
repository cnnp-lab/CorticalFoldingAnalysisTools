# CorticalFoldingAnalysisTools
Matlab scripts to extract features of cortical folding from Freesurfer folders. Current version 1.0

## What can this do?

Given a list of Freesurfer subjects, the provided code will automatically extract among other measures: the average cortical thickness (slightly different to FS's own algorithm), the total pial surface area (different to FS's code), and the exposed surface area. With these three quantities, [an analysis to test for the universal scaling of cortical folding](https://doi.org/10.1073/pnas.1610175113) can be performed. The main difference to FS's own code is that this pipeline accounts for the parts on the surface mesh that are not cortex.

Note that the files ?h.pial-outer-smoothed have to exist. See the [Freesurfer website](https://surfer.nmr.mgh.harvard.edu/fswiki/LGI) for more details.



## How to use?

The lib folder must be included on the path in all folders and subfolders.

### Hemisphere-based analysis 
For extraction of the whole-hemisphere features, Hemisphere/extract_FreeSurferHemi_features.m is the main function. Demo code is provided in Hemisphere/demo_extract_FreeSurferHemi_features.m with more details on the usage & functionalities.
This is the code used for the [2016 PNAS publication](https://doi.org/10.1073/pnas.1610175113).

### Lobe-based analysis
For extraction of the lobe-based features, Lobes/extract_FreeSurferLobes_features.m is the main function. Demo code is provided in Lobes/demo_extract_FreeSurferLobe_features.m with more details on the usage & functionalities.
This is the code used for the [2019 Commun Biol publication](https://www.nature.com/articles/s42003-019-0421-7)

### Dependencies
For the lobe-based analysis, you will need the [iso2mesh](http://iso2mesh.sourceforge.net) matlab library on your path.


## How to cite?

If you use this code for your publications, please cite the zenodo entry for this pipeline: [https://doi.org/10.5281/zenodo.3608675](https://doi.org/10.5281/zenodo.3608675)


## Related publications

[Universality in human cortical folding in health and disease, PNAS 2016](https://doi.org/10.1073/pnas.1610175113)

[Human cortical folding across regions within individual brains follows universal scaling law, Commun Biol 2019](https://www.nature.com/articles/s42003-019-0421-7)


## Contributors & thanks

In chronological order:
* Yujiang Wang (2015 - now): developing main algorithms & code
* Andre Muricy (2016): developing lib/find_smooth_labels_adjp.m and associated scripts
* Joe Necus (2016-2019): code checking
* Kathryn Garside (2017-2018): code checking
* Tobias Ludwig (2019-2020): code checking, adding table and merging functionality Hemisphere/extract_FreeSurferHemi_features.m and Lobes/extract_FreeSurferLobes_features.m
