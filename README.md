# CorticalFoldingAnalysisTools
Matlab toolbox to extract features of cortical folding from Freesurfer folders. Current version 3.0

## What can this toolbox do?
The MATLAB toolbox is designed to extract various features of cortical folding derived/based on the universal scaling law of cortical morphology (see "Related publications").
The toolbox assumes Freesurfer folder structures as input, and can provide outputs (features of cortical folding) in a variety of ways. See section below on "Further statistical analysis", if you want a starting point. 

In this version, we have now also added a Graphical User Interface (GUI) if you are not familiar with MATLAB.
Our collaborators also have a Python version of this toolbox, **contact us (cnnplab@ncl.ac.uk)** if you want early access, or have general questions.

## How to use the toolbox?

### Before you start

You need to consider if you wish to use the GUI or prefer to use and adapt the core scripts. The former is recommended as a first step, as it brings together all the key functionalities in one place. A **manual** can be found in the "Graphic UI" folder, and is the easiest place to start.
We generally recommend a Linux/Unix or MacOS environment, which is what we have tested most extensively.

#### Required external library 
[iso2mesh](https://iso2mesh.sourceforge.net/cgi-bin/index.cgi) is needed for the Lobes and Multiscale functionalities. Please make sure this is on your path.

#### Using and adapting the core scripts
If you decide to only use particular core functionalities without the GUI, make sure the points below are followed. (The GUI takes care of these points without need for further user intervention.)
The "lib" folder must be included on the path in all folders and subfolders.
You also will need the following files in the surf folder of your Freesurfer subject:
?h.pial, ?h.white, ?h.thickness, and ?h.pial-outer-smoothed
If ?h.pial-outer-smoothed does not exist, follow the [Freesurfer LGI](https://surfer.nmr.mgh.harvard.edu/fswiki/LGI) pipeline. You can just run the first few commands until you get the ?h.pial-outer-smoothed file (just a few minutes), but you can also run the whole LGI pipeline if you want (more time consuming).

### Hemisphere-based analysis 
If you wish to perform a similar analysis as our [2016 PNAS publication](https://doi.org/10.1073/pnas.1610175113), you can follow the hemisphere-based analysis stream. Our toolbox extracts, for each cortical hemisphere in each subject, the average cortical thickness, total surface area, and exposed surface area. From these quantities you can calculate the "independent components" of cortical morphology, which are suggested as alternative and more robust cortical morphometrics (see [2021 NeuroImage publication](https://doi.org/10.1016/j.neuroimage.2020.117546) for more details).

#### Core scripts
For extraction of the whole-hemisphere features, Hemisphere/extract_FreeSurferHemi_features.m is the main function. Demo code is provided in Hemisphere/demo_extract_FreeSurferHemi_features.m with more details on the usage & functionalities.


### Lobe-based analysis
If you wish to perform a similar analysis as our [2019 Commun Biol publication](https://www.nature.com/articles/s42003-019-0421-7), you can follow the lobe-based analysis stream. The stream is very similar to the hemisphere-based analysis, but instead of analysing an entire cortical hemisphere, we subdivide it into four lobes (frontal, temporal, parietal, and occipital). As before, you can also obtain the "independent components" for each lobe.

#### Core scripts
For extraction of the lobe-based features, Lobes/extract_FreeSurferLobes_features.m is the main function. Demo code is provided in Lobes/demo_extract_FreeSurferLobe_features.m with more details on the usage & functionalities.


### Multiscale analysis
If you wish to perform a similar analysis as our [2024 eLife publication](https://doi.org/10.7554/eLife.92080.4), you can follow the multiscale analysis stream. This stream is effectively a multiscale version of the hemisphere-based analysis, but instead of analysing an entire cortical hemisphere as is, we disaggregate its shape at various spatial scales. 

#### Core scripts
For extraction of the multiscale features, Scales/fastEstimateScale.m is the main function. Demo code is provided in Scales/example_run.m with more details on the usage & functionalities.


## Further statistical analysis
If you used the GUI to obtain outputs, you can upload these directly onto our [Brain MoNoCle](https://cnnplab.shinyapps.io/BrainMoNoCle/) platform for further statistical analysis. More details can be found on the [website](https://cnnplab.shinyapps.io/BrainMoNoCle/), and in our [paper](https://arxiv.org/abs/2406.01107).
If you used the core scripts instead, please take a look at the [Brain MoNoCle](https://cnnplab.shinyapps.io/BrainMoNoCle/) example data to understand how to re-format your outputs.
Of course you can also use any number of other toolboxes and platforms to analyse these outputs further.

## How to cite?

If you use this code for your publications, please cite this GitHub repo and the zenodo entry for this pipeline: [https://doi.org/10.5281/zenodo.3608675](https://doi.org/10.5281/zenodo.3608675)


## Related publications

[Cortical folding scales universally with surface area and thickness, not number of neurons, Science 2015](https://doi.org/10.1126/science.aaa9101)

[Universality in human cortical folding in health and disease, PNAS 2016](https://doi.org/10.1073/pnas.1610175113)

[Human cortical folding across regions within individual brains follows universal scaling law, Commun Biol 2019](https://www.nature.com/articles/s42003-019-0421-7)

[Independent components of human brain morphology, NeuroImage 2021](https://doi.org/10.1016/j.neuroimage.2020.117546)

[Neuro-evolutionary evidence for a universal fractal primate brain shape, eLife 2024](https://doi.org/10.7554/eLife.92080.4)


## Roadmap for future development

If you are interested in a more regional analysis (e.g. DK ROIs), or a surface-based analysis such as in our [MICCAI 2021 publication](https://link.springer.com/chapter/10.1007/978-3-030-87234-2_65), please drop us a message at cnnplab@ncl.ac.uk. We are actively developing these and would love to hear from you.

Our regional analysis stream of multiscale brain morphology is also under development, with first results shown in (this preprint)[https://arxiv.org/abs/2311.13501]. We envisage this functionality to become available in 2025 as a core script.


## Contributors & thanks

In chronological order:
* Yujiang Wang (2015 - now): developing main algorithms & code
* Andre Muricy (2016): developing lib/find_smooth_labels_adjp.m and associated scripts
* Joe Necus (2016-2019): code checking
* Kathryn Garside (2017-2018): code checking
* Tobias Ludwig (2019-2020): code checking, adding table and merging functionality Hemisphere/extract_FreeSurferHemi_features.m and Lobes/extract_FreeSurferLobes_features.m
* Karoline Leiberg (2020 - now): developing main algorithms & code
* Guillermo Besne (2023 - now): GUI development
