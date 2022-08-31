# NeuroMiner

NeuroMiner is continuously extended and updated by Nikolaos Koutsouleris, Clara Vetter & Ariane Wiegand. Nikolaos Koutsouleris first developed NeuroMiner
back in 2009 to provide clinical researchers with cutting-edge machine
learning methods for the analysis of heterogeneous data domains, such as
clinical and neurocognitive read-outs, structural and functional
neuroimaging data, and genetic information. The program is an interface
to a large variety of unsupervised and supervised pattern recognition
algorithms that have been developed in the machine learning field over
the last decades. Furthermore, the application implements different
strategies for preprocessing, filtering and fusing of heterogeneous
data, training ensembles of predictors, and for the visualization and
testing of the significance for the computed predictive patterns. The
current release candidate of NeuroMiner has been tested in the Section
for Precision Psychiatry on a variety of datasets containing healthy
controls and patients with psychiatric disorders, and it was designed
specifically to create robust models with a high probability of
generalization to new datasets.

### Publications using NeuroMiner
For reference, we include [a list of
papers below](intro_nm_papers) or see [publications listed in
PubMed](https://www.ncbi.nlm.nih.gov/pubmed/?term=koutsouleris+n),
which were all based on previous versions of the program.

(intro_nm_features)=
### NeuroMiner Feautures
More specifically, using a light-weight and interactive text-based menu
system, NeuroMiner allows the user to:

1. load data easily (e.g., using spreadsheets, NifTi images, or SPM
structures);
2. build a variety of cross-validation frameworks for classification and
regression problems that have become a gold standard in the field (e.g.,
repeated nested cross-validation, leave-site-out cross-validation);
3. apply a number of preprocessing strategies (e.g., scaling, filtering,
many forms of dimensionality reduction, etc.);
4. choose and combine cutting-edge supervised algorithms (e.g., support
vector machine, elastic net, random forest, etc.);
5. apply feature selection procedures (e.g., wrappers), data fusion
techniques, and stacked generalization;
6. apply learned models to new data (external validation).

To assist in selecting and analyzing data, the user can visualize the
data during input, monitor accuracy during learning, and understand the
results of complex analyses using multiple display options. These allow
the user to accurately report the data and also to understand the
underlying machine learning analyses. Furthermore, the ability to apply
the learned models to completely unseen data is important because it is
quickly becoming a standard requirement for all machine-learning
studies. In total, NeuroMiner gives the user the ability to design,
implement, understand, and report machine learning analyses.

#### DISCLAIMER

Please note that NeuroMiner is supplied as is and no formal maintenance
is provided or implied. In no event shall the author of the software
(heretofore known as the Author) be liable to any party for direct,
indirect, special, incidental, or consequential damages, including lost
profits, arising out of the use of this software and its documentation,
even if the Author has been advised of the possibility of such damage.
The Author specifically disclaims any warranties, including, but not
limited to, the implied warranties of merchantability and fitness for a
particular purpose. The software and accompanying documentation provided
hereunder is provided 'as is'. The Author has no obligation to provide
maintenance, support, updates, enhancements, or modifications (but we
plan to).

This is the beta release version of the software and the software is
undergoing regular updates. Please send any comments, questions, or bug
reports to email.neurominer\@gmail.com.

(intro_nm_papers)=
### Papers that have used used NM:

1. [Koutsouleris N, Meisenzahl EM, Davatzikos C, Bottlender R, Frodl T,
Scheuerecker J, Schmitt G, Zetzsche T, Decker P, Reiser M, Moller HJ,
Gaser C. Use of neuroanatomical pattern classification to identify
subjects in at-risk mental states of psychosis and predict disease
transition. Archives of General Psychiatry. 2009; 66(7):700-12.](https://10.1001/archgenpsychiatry.2009.62)
2. Koutsouleris N, Meisenzahl EM, Borgwardt S, Riecher-Rossler A, Frodl
T, Kambeitz J, Kohler Y, Falkai P, Moller H.-J., Reiser M, Davatzikos C.
Individualized differential diagnosis of schizophrenia and mood
disorders using neuroanatomical biomarkers. Brain. 2015 Jul;138(Pt
7):2059-73.
3. Koutsouleris N, Kahn RS, Chekroud AM, Leucht S, Falkai P, Wobrock T,
Derks EM,Fleischhacker WW, Hasan A. Multisite prediction of 4-week and
52-week treatment outcomes in patients with first-episode psychosis: a
machine learning approach. Lancet Psychiatry. 2016 Oct;3(10):935-946.
doi: 10.1016/S2215-0366(16)30171-7.
4. Cabral C, Kambeitz-Ilankovic L, Kambeitz J, Calhoun VD, Dwyer DB, von
Saldern S, Urquijo MF, Falkai P, Koutsouleris N. Classifying
Schizophrenia Using Multimodal Multivariate Pattern Recognition
Analysis: Evaluating the Impact of Individual Clinical Profiles on the
Neurodiagnostic Performance. Schizophrenia Buletinl. 2016 Jul;42 Suppl
1:S110-7. doi: 10.1093/schbul/sbw053.
5. Koutsouleris N, Borgwardt S, Meisenzahl EM, Bottlender R, Moller HJ,
Riecher-Rossler A. Disease prediction in the at-risk mental state for
psychosis using neuroanatomical biomarkers: results from the
FePsy-study. Schizophrenia Bulletin. 2012; 38(6):1234-46
6. Borgwardt SJ, Koutsouleris N, Aston J, Studerus E, Smieskova R,
Riecher-Rossler A, Meisenzahl EM. Distinguishing prodromal from
first-episode psychosis using neuroanatomical pattern recognition:
Evidence from single-subject structural MRI. Schizophrenia Bulletin.
2013; 39(5):1105-14. doi: 10.1093/schbul/sbs095
7. Koutsouleris N, Davatzikos C, Bottlender R, Patschurek-Kliche K,
Scheuerecker J, Decker P, Gaser C, Moller HJ; Meisenzahl E. Early
recognition and disease prediction in the at-risk mental states for
psychosis using neurocognitive pattern classification. Schizophrenia
Bulletin. 2012; 38(6):1200-15
8. Koutsouleris N, Riecher-Rossler A, Meisenzahl E, Smieskova R,
Studerus E, Kambeitz-Ilankovic L, von Saldern S, Cabral C, Reiser M,
Falkai P, Borgwardt S. Detecting the psychosis prodrome across high-risk
populations using neuroanatomical biomarkers. Schizophrenia Bulletin.
2014, 41(2):471-82.
9. Koutsouleris N, Davatzikos C, Borgwardt S, Gaser C, Bottlender R,
Frodl T, Falkai P, Riecher-Rossler A, Moller HJ, Reiser M, Pantelis C,
Meisenzahl E. Accelerated Brain Aging in Schizophrenia and Beyond: A
Neuroanatomical Marker of Psychiatric Disorders. Schizophrenia Bulletin.
2014 Sep;40(5):1140-53
10. Koutsouleris N, Gaser C, Bottlender R, Davatzikos C, Decker P, Jager
M, Schmitt G, Reiser M, Moller HJ, Meisenzahl EM, Use of Neuroanatomical
Pattern Regression to Predict the Structural Brain Dynamics of
Vulnerability and Transition to Psychosis. Schizophrenia Research.
2010;123(2-3):175-187
11. Kambeitz-Ilankovic L, Meisenzahl EM, Cabral C, von Saldern S,
Kambeitz J, Falkai P, Moller HJ, Reiser M, Koutsouleris N. Prediction of
outcome in the psychosis prodrome using neuroanatomical pattern
classification. Schizophrenia Research. 2015;173(3):159-65.
12. Koutsouleris et al., JAMA Psychiatry. Prediction models of
functional outcomes for individuals in the clinical high-risk state for
psychosis or with recent-onset depression: a multimodal, multisite
machine learning analysis. 2018.
