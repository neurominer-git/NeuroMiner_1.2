Getting Started
===============

Compatibility
-------------

NeuroMiner has been tested using Windows, Linux (Cent OS), and Mac OS
operating systems and works with MATLAB2020b and above. We recommend the
\"Statistics and Machine Learning\" and the \"Optimization\" toolboxes
in order to use advanced training options in the matLearn toolbox and
for imputation functions.

For non-supported operating systems, NeuroMiner has a number of programs
installed that may require independent compilation using a C compiler
(e.g., gcc for linux or Xcode for mac).

For some functionalities, NeuroMiner resorts to recognized Python
libraries. Thus, a Python installation (tested with Python 3.8) may be
required. If you have never called a Python script from Matlab, you
initially need to set the correct Python environment in Matlab. You can
find instructions on how to set the Python version in the Matlab Help
Center
([link](https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)).
You only have to configure the path once.

Downloading and Initialising NeuroMiner
---------------------------------------

NeuroMiner works in two main modes: matrix and neuroimaging (i.e., NIFTI
or FreeSurfer files). For neuroimaging using NIFTI data format, it is
necessary to have a working copy of the Statistical Parametric Mapping
(SPM) toolbox installed
([download](http://www.fil.ion.ucl.ac.uk/spm/software/download/)). We
also recommend downloading the WFU Pickatlas toolbox
([link](http://fmri.wfubmc.edu/software/pickatlas)) to take advantage of
this selection tool when inputting data. For surface-based analyses, it
is necessary to have a FreeSurfer distribution installed
([download](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)).

Once the program is downloaded, it needs to be added to the MATLAB path
by using the **addpath** command from the command line or by using the
\"Set Path\" button on the toolbar (e.g., \"addpath
/home/NM/NeuroMiner\"). NeuroMiner can then be launched by typing **nm**
on the command line. To load the program faster, you can type **nm
nosplash** and to enter expert mode type **nm nosplash expert**. The
first time that NeuroMiner is run, a file selector box will appear first
asking to define the SPM directory and then the FreeSurfer directory. If
these directories are not available, then simply press cancel and
NeuroMiner will be launched in non-imaging mode--i.e., a mode that
restricts options only to matrix data (see Fig.
[1](#fig:initialisation); upper). Otherwise, NeuroMiner will be
launched in matrix and neuroimaging mode (Fig.
[1](#fig:initialisation); lower).

![Display screen of NeuroMiner in different modes. NeuroMiner will be
loaded in non-imaging mode (upper) when the SPM or FreeSurfer
directories are not provided during start-up (upper). If the directories
have been provided, it will be started in standard imaging and matrix
data mode
(lower).](Images/initialisation_imaging-nonimaging.png)<a name="fig:initialisation"></a>

**User note**: once NeuroMiner is initialised and added to the path, there
can be incompatibilities between some MATLAB functions on some systems.
For example, there can be a incompatibility with the function
**ttest2**. In these cases, NeuroMiner needs to be removed from the path
in order to use the original function.

Nomenclature
------------

Machine-learning uses different nomenclature compared to classical
statistics. The following terms and their definitions are used
throughout the manual (Rob Tibshirani, Statistics 315a, Stanford):

  | **Machine Learning**     | **Statistics**                 |
|--------------------------|--------------------------------|
| Features                 | Predictor variables            |
| Target variable/labels   | Outcome variables              |
| Weights                  | Parameters                     |
| Hyperparameters          | Parameters                     |
| Learning                 | Fitting                        |
| Generalization           | Test set performance           |
| Supervised learning      | Regression/ Classification     |
| Unsupervised learning    | Clustering, density estimation |
| Large grant = $1,000,000 | Large grant = $50,000          |
