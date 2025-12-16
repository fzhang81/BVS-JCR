{\rtf1\ansi\ansicpg1252\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Bold;\f1\froman\fcharset0 Times-Roman;\f2\fswiss\fcharset0 Helvetica-Bold;
\f3\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
{\*\listtable{\list\listtemplateid1\listhybrid{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{disc\}}{\leveltext\leveltemplateid1\'01\uc0\u8226 ;}{\levelnumbers;}\fi-360\li720\lin720 }{\listname ;}\listid1}
{\list\listtemplateid2\listhybrid{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{disc\}}{\leveltext\leveltemplateid101\'01\uc0\u8226 ;}{\levelnumbers;}\fi-360\li720\lin720 }{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{circle\}}{\leveltext\leveltemplateid102\'01\uc0\u9702 ;}{\levelnumbers;}\fi-360\li1440\lin1440 }{\listname ;}\listid2}
{\list\listtemplateid3\listhybrid{\listlevel\levelnfc0\levelnfcn0\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{decimal\}}{\leveltext\leveltemplateid201\'01\'00;}{\levelnumbers\'01;}\fi-360\li720\lin720 }{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{circle\}}{\leveltext\leveltemplateid202\'01\uc0\u9702 ;}{\levelnumbers;}\fi-360\li1440\lin1440 }{\listname ;}\listid3}}
{\*\listoverridetable{\listoverride\listid1\listoverridecount0\ls1}{\listoverride\listid2\listoverridecount0\ls2}{\listoverride\listid3\listoverridecount0\ls3}}
\margl1440\margr1440\vieww18440\viewh11420\viewkind0
\deftab720
\pard\pardeftab720\sl360\slmult1\partightenfactor0

\f0\b\fs48 \cf0 \expnd0\expndtw0\kerning0
Bayesian Variable Selection via Joint Credible Regions\
\pard\pardeftab720\sl360\slmult1\partightenfactor0

\f1\b0\fs24 \cf0 This repository contains the R implementation for a Bayesian variable selection method applied to computer experiments (specifically demonstrated on the Borehole function).\
The algorithm identifies active variables by constructing a 
\f0\b Joint Credible Region (JCR)
\f1\b0  based on a generated "inert" (dummy) variable. Variables whose posterior parameters fall outside this reference region are classified as active.\
\pard\pardeftab720\sl360\slmult1\partightenfactor0

\f2\b\fs36 \cf0 \outl0\strokewidth0 \strokec2 Key Features\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sl360\slmult1\partightenfactor0
\ls1\ilvl0
\fs24 \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Inert Variable Construction:
\f3\b0  Generates an orthogonal inert variable to serve as a baseline for noise.\
\ls1\ilvl0
\f2\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Bayesian Inference:
\f3\b0  Uses rstan to sample from the posterior distribution of the model parameters (Gaussian Process with ARD kernel).\
\ls1\ilvl0
\f2\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Selection Mechanism:
\f3\b0  Constructs a 95% Joint Credible Region (ellipse) using the posterior samples of the inert variable's parameters ($\\beta$ and transformed $\\theta$).\
\ls1\ilvl0
\f2\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Visualization:
\f3\b0  visually maps the selection boundary and variable positions.\
\pard\pardeftab720\sl360\slmult1\partightenfactor0

\f2\b\fs36 \cf0 Dependencies\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sl360\slmult1\partightenfactor0
\ls2\ilvl0
\fs24 \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 R
\f3\b0  (>= 4.0.0)\
\ls2\ilvl0
\f2\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 R Packages:
\f3\b0 \
\pard\tx940\tx1440\pardeftab720\li1440\fi-1440\sl360\slmult1\partightenfactor0
\ls2\ilvl1\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u9702 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 rstan (for MCMC sampling)\
\ls2\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u9702 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 MASS (for statistical functions)\
\ls2\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u9702 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 car (for data ellipses)\
\pard\pardeftab720\sl360\slmult1\partightenfactor0

\f2\b\fs36 \cf0 Usage\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sl360\slmult1\partightenfactor0
\ls3\ilvl0
\f3\b0\fs24 \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	1	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Ensure your data files (X.csv, y.csv) and the Stan model file (L_cov_exp_quad_ARD5_samplesigma.stan) are in the working directory.\
\ls3\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	2	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Configure the file paths and simulation parameters in the script:\uc0\u8232 R\
\ls3\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	3	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 X_path <- 'path/to/your/design_matrix.csv'\
\ls3\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	4	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 y_path <- 'path/to/your/response.csv'\
\ls3\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	1	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 stan_file <- 'your_model.stan'\uc0\u8232 \
\ls3\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	2	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Run the main R script. The code will:\
\pard\tx940\tx1440\pardeftab720\li1440\fi-1440\sl360\slmult1\partightenfactor0
\ls3\ilvl1\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u9702 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Perform MCMC sampling.\
\ls3\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u9702 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Compute the median $\\beta$ and length-scale parameters.\
\ls3\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u9702 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Generate a selection plot.\
\ls3\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u9702 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Output the selection status (Selected/Not Selected) for each variable.\
}