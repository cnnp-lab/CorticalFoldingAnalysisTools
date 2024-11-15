# CorticalFoldingAnalysisTools
Matlab toolbox to extract features of cortical folding from Freesurfer folders. Current version 2.0

## What can this toolbox do?
The MATLAB toolbox is designed to extract various features of cortical folding derived/based on the universal scaling law of cortical morphology (see "Related publications").
The toolbox assumes Freesurfer folder structures as input, and can provide outputs (features of cortical folding) in a variety of ways. See section below on "Further statistical analysis", if you want a starting point. 

In this version, we have now also added a Graphical User Interface (GUI) if you are not familiar with MATLAB.
Our collaborators also have a Python version of this toolbox, **contact us (cnnplab@ncl.ac.uk)** if you want early access, or have general questions.

## How to use?


### Before you start
The lib folder must be included on the path in all folders and subfolders.
You also will need the following files in the surf folder of your Freesurfer subject:
?h.pial, ?h.white, ?h.thickness, and ?h.pial-outer-smoothed
If ?h.pial-outer-smoothed does not exist, follow the [Freesurfer LGI](https://surfer.nmr.mgh.harvard.edu/fswiki/LGI) pipeline. You can just run the first few commands until you get the ?h.pial-outer-smoothed file (just a few minutes), but you can also run the whole LGI pipeline if you want (more time consuming).


### Hemisphere-based analysis 
For extraction of the whole-hemisphere features, Hemisphere/extract_FreeSurferHemi_features.m is the main function. Demo code is provided in Hemisphere/demo_extract_FreeSurferHemi_features.m with more details on the usage & functionalities.
This is the code used for the [2016 PNAS publication](https://doi.org/10.1073/pnas.1610175113).

### Lobe-based analysis
For extraction of the lobe-based features, Lobes/extract_FreeSurferLobes_features.m is the main function. Demo code is provided in Lobes/demo_extract_FreeSurferLobe_features.m with more details on the usage & functionalities.
This is the code used for the [2019 Commun Biol publication](https://www.nature.com/articles/s42003-019-0421-7)

### Dependencies
For the lobe-based analysis, you will need the [iso2mesh](http://iso2mesh.sourceforge.net) matlab library on your path.


## Further statistical analysis

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
