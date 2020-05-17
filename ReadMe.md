This project contains three matlab functions required for extracting wave latencies and amplitudes from auditory brainstem responses (ABRs) usiing non-linear curve registration through continuous monotone registration.

(i) preproc.m:
This function preprocesses a set if individual ABRs by normalizing them by their root-mean-square amplitude after optionally pre-aligning them either to a separately defined target response (tgt), or to the subject-average response for a prespecified condition (TGT; if the data set comprises more than one condition). The pre-alignment is performed by linear time shifting. 
(ii) nlcurvereg.m:
This function non-linearly aligns the pre-processed ABRs either with the predefined target (tgt) or with the subject-average ABR for the prespecified target condition (TGT), or - if neither are defined - to the garnd-average of all to-be-aligned ABRs.
(iii) xtractlatamp.m:
This function requires the user to pick the structural peaks and troughs of the required waves and then uses these to derive the waves' latencies and amplitudes for each ABR.

In addition, the project also contains a matlab script, examples.m, which demonstrates how the above functions should be used in different situations. example.m uses the data contained in the .mat file chirdata.mat. 
