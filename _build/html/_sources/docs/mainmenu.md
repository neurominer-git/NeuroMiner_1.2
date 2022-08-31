Main Interface
==============

NeuroMiner is software that has a text-based interface and works from
the MATLAB command line. This implies that all analyses are implemented
by selecting options from text menus and entering parameters when asked
to do so. Each of the menus share a similar format as described in the
following example.

After the neuroimaging paths have been established, the first time that
NeuroMiner is run the user will see the **MAIN INTERFACE** (see
Fig.[1](#fig:NM_start-screen){reference-type="ref"
reference="fig:NM_start-screen"}) that consists of 4 different options
that are numbered 1-4. Below these numerical options there is also an
option to go **Back/Quit (Q)**.Below the dividing line (i.e.,
==========), the user can enter their choice of option in the space
after **Menu choice (1-4/Q) (Default: Q)?**.

![The main starting display of
NeuroMiner](Images/NM_start-screen){#fig:NM_start-screen
width="0.8\\linewidth"}

In this example, the user can do one of three things:

-   1\) to select one of the 4 numerals corresponding to each option and
    press the enter key on the keyboard;

-   2\) to enter \"Q\" to go back to the previous menu, or if they are
    in the **MAIN INTERFACE** to quit the program;

-   3\) to simply press the enter key and the program will activate the
    default option--in this case, to quit the program.

For example, a user could press \"2\" followed by the enter key to
\"Load NeuroMiner structure\". The NeuroMiner MATLAB structure (\"NM\")
is the primary variable that NeuroMiner uses to store data, settings,
analyses, and results. For example, when a user inputs data into
NeuroMiner it is stored in the NM structure as a data matrix. Then when
the user chooses settings or performs an analysis, all of these inputs
and outputs are stored in the structure as well. This means that the NM
structure forms the core of any analysis within NeuroMiner.

It is important to note that the NeuroMiner structure is not
automatically saved during analyses in order to give the user more
control. This means that the user must save the structure periodically
during the analysis using the menu item **10: Save NeuroMiner
structure**.

NeuroMiner will **automatically change** the items in all menus based on
two main criteria:

-   1\) Whether a NeuroMiner MATLAB structure storing all the
    information from a previous analysis has been loaded to the MATLAB
    workspace;

-   2\) The individual options that are selected during a NeuroMiner
    session.

For example, if a NM structure is loaded prior to the start of a
NeuroMiner session then the program will automatically detect this and
the **MAIN INTERFACE** will specifically relate to the loaded analysis.
Similarly, if specific options are chosen during the establishment of
the analysis then the menus will adapt to the chosen settings and
options that do not apply to the analysis will be removed.

**The Main Interface After Data Has Been Loaded and Processed**

Once the user has entered and processed data (see section
[\[mainmenu_3.1_input_data\]](#mainmenu_3.1_input_data){reference-type="ref"
reference="mainmenu_3.1_input_data"}), the **MAIN INTERFACE** will
change to the following:

1: Inspect data used for model discovery and cross-validation\
2: Set up NM parameter workspace\
3: Initialize & manage analyses\
4: Preprocess features\
5: Train supervised classifiers\
6: Visualise classifiers\
7: Display training results\
8: Generate master files\
9: Load NeuroMiner structure\
10: Save NeuroMiner structure\
11: Clear NeuroMiner structure\
12: Change working directory\
13: Open Manual\
14: Utilities\
Back/Quit\[Q\]

These options are major functional domains related to the methods and
stages of machine learning analysis using NeuroMiner. The contents of
each of these options are described throughout the manual, but it is
important to understand the overall analytic workflow that the program
uses by briefly describing each of the options above.

= \[diamond, draw, fill=blue!20, text width=5.5em, text badly centered,
node distance=3cm, inner sep=0pt\] = \[rectangle, draw, fill=blue!20,
text width=8em, text centered, rounded corners, minimum height=4em\] =
\[rectangle, draw, fill=red!20, text width=10em, text centered, rounded
corners, minimum height=4em\] = \[draw, -latex'\] = \[draw,
ellipse,fill=red!20, text centered, text width=6em, node distance=3cm,
minimum height=3em\]

**1: Input data into NM** The entry of data is the most important aspect
of any analysis. When data is entered into NeuroMiner, it is stored
within a MATLAB structure called \"NM\" (i.e., \"NeuroMiner\") in the
workspace. All other analyses and results are then added to this NM
structure as the user progresses. The contents of the structure can also
be viewed from the MATLAB command line.

**2: Define parameter template** Once the data have been entered, the
first step of any analysis is to define machine learning options. Within
this menu, the user can establish the cross-validation framework (e.g.,
nested cross-validation), perform essential preprocessing operations on
the data (e.g., scaling), choose the machine learning algorithm (e.g.,
SVM), define the machine learning algorithm settings (e.g., C values),
and perform other operations on the data.

**3: Initialize/Delete analyses** After the entry of all the settings
for an analysis, the user can then 'initialize' these settings in
preparation for processing.This option saves the parameters in an
analysis structure, which the user can then select when they go to
process the data. This is a great feature because it means that the user
can define multiple different analyses for the same dataset by changing
parameters in the 'parameter template' and then initializing separate
analyses.

**4: Preprocess Features** In NeuroMiner, the analysis can be split into
two parts: preprocessing features and training supervised classifiers.
Preprocessing in NeuroMiner can also be called 'feature preparation and
selection', and involves preparing the data for machine learning
analysis (e.g., a support vector machine). For example, the researcher
may want to covary for age and sex, perform a Principal Components
Analysis, and scale neuroimaging data prior to analysis. Neurominer
performs all preprocessing steps within the cross-validation structure,
which is recognised as the most valid approach. Allowing the user to
preprocess their data first within the cross-validation structure can
save a lot of computational resources and time because it means that the
researcher can run multiple analyses from the same preprocessed data.

**5: Train supervised classifiers** This option allows the user to
select and run the machine learning analysis. For example, the user can
run a support vector machine across a range of C parameters. The results
of this analysis are then stored in the NM structure.

**6: Visualise classifiers** This step is optional. In a standard
machine learning analysis, it is often unclear how each of the features
(e.g., imaging voxels or variables) contribute to the prediction
accuracy and they are often not in the same space as the original
images. For this reason, the user has the option within NeuroMiner
calculate statistics for each feature with the 'visualise classifiers'
option. When using neuroimaging data, this option also projects (i.e.,
transforms) the data into the original space so that it can be viewed in
a similar manner as other neuroimaging results in nii files.

**7: Display training results** Based on the above essential steps, and
the optional step to visualise the classifiers, the results of the
analysis will be ready to be viewed with this step. For example, the
user can view the accuracy, the frequency of model selection, the
classifier weights that are applied to each variable, and any additional
visualisation routines that have been run.

**8: Generate master files** As mentioned above, NeuroMiner stores
settings and analysis results in the NM structure within MATLAB, which
is visible in the workspace. In addition to this, NeuroMiner can also
store individual files related to steps **4: Preprocess Features** and
**5: Train supervised classifiers** above. When NeuroMiner does this, it
stores a single file for each permutation within the cross-validation
framework--as described later. These individual files can then be saved
as \"master files\" using this option.

**9 - 11: Load, save, and clear NeuroMiner structure** These options
allow the user to load the NM structure (where all the settings and
analyses are stored) to the MATLAB workspace, to save the structure as a
file, and to clear the structure from the workspace.

**12: Change working directory** In similarity to most neuroimaging
software, NeuroMiner operates from a \"working directory\" where
preprocessing or analysis files are saved during any analysis. The
working directory is automatically the directory where NeuroMiner is
started, but can be changed using this item.

**13: Open manual** This opens the manual that you are reading right
now.

**14: Utilities** Contains functions related to FreeSurfer and ability
to add SPM/FreeSurfer paths to the NeuroMiner environment.

Load NeuroMiner structure
-------------------------

As outlined in section [1](#mainmenu){reference-type="ref"
reference="mainmenu"}, the NM structure is the core structure that
NeuroMiner uses to store and process data. This option allows the user
to load it from a .mat file.

Save NeuroMiner structure
-------------------------

This option allows the user to save the NeuroMiner structure (NM) as a
.mat file in a folder of the users choosing. Notably, any changes that
are made to the settings of NeuroMiner during an active session are
stored within the NM structure in the workspace. However, program
crashes can lead to the loss of information. We recommend to save the NM
structure at critical points of the analysis: 1) After data entry; 2)
after structure initialization; and 3) after preprocessing or analysis.
This ensures that no data is lost.

Clear NeuroMiner structure
--------------------------

This option will clear the NeuroMiner structure **from the workspace**.
This will allow the user to select another structure or enter different
data.

Change Working Directory
------------------------

NeuroMiner will save files in the directory that it is loaded from. This
option allows the user to change the working directory.

Open Manual
-----------

This feature is currently not functional. Please open the manual
normally.

Utilities
---------

#### 1: Set root paths of neuroimaging tools (SPM/FreeSurfer)

If you want to process structural MRI (nii) or surface files from
FreeSurfer (mgz) then the paths of these tools need to be established.
This option will establish the paths. If the user does not establish the
paths, then the neuroimaging functions will not be displayed.

#### 2: Create PreprocData path master

The option allows the user to update the existing path for only the
pre-processed data.

#### 3: Create CVdatamat path master

This option facilitates the user to update path for the cross-validated
data.

#### 4: Create CVresults path master

This enables to user to configure the cross-validated results in a new
directory.

#### 5: Create VISdatamat path master

This option sets up the user defined-updated path after running the
visualisation.

#### 6: Create OptPreprocParam path master

This is to create a master file that can be used for batch scripts or
when repetitive actions are carried out on a set of files the user
creates the file list once and does not waste to search specific file
sets

#### 7: Create OptModel path master

this is to create a master file that can be used for batch scripts or
when repetitive actions are carried out on a set of files the user
creates the file list once and does not waste to search specific file
sets.

#### 8: Update analyses' paths to new root directory

This is an important utility option which allows the users to update the
NM structure path, if the user wants to continue the analysis in another
directory. By selecting this option, all the paths in the NM structure
are updated according to the new/selected location.
