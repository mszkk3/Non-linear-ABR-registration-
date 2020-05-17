{\rtf1\ansi\ansicpg1252\cocoartf1671\cocoasubrtf600
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww32920\viewh14160\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 This project contains three matlab functions required for extracting wave latencies and amplitudes from auditory brainstem responses (ABRs) usiing non-linear curve registration through continuous monotone registration.\
(i) preproc.m:\
This function preprocesses the ABRs by normalizing them by their root-mean-square amplitude after optionally pre-aligning them either to a separately entered target response (tgt) or to the subejct-average response for a prespecified condition (TGT). The pre-alignment is performed by linear time shifting. \
(ii) nlcurvereg.m:\
This function non-linearly aligns the pre-processed ABRs either with the prespecified target (tgt) or the subject-average ABR for the prespecified condition (TGT), or, if neither are defined, to the garnd-average of all to-be-aligned ABRs. \
(iii) xtractlatamp.m:\
This function requires the user to pick the structural peaks and troughs of the required waves and then uses these to derive the waves\'92 latencies and amplitudes for each ABR. \
}