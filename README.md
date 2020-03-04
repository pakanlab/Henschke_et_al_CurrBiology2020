# Henschke*, Dylda*, Katsanevaki*, Dupuy, Currie, Amvrosiadis, Pakan* and Rochefort*, Curr. Biology 2020
### *Reward association enhances stimulus-specific representation in primary visual cortex*

MATLAB code from the Pakan lab and Rochefort lab used for analysis in this paper.

## General notes
- The functions used for analysis are defined in **getFuncHandleRGroup.m**. 
These files can be found in the folder **\analysis**.
The analysis functions may call other specific functions (e.g. compute_lmi.m), which can be found in **\analysis\projectfunctions**, 
as well as utilities tools (e.g. for data handling) which are in the folder **\analysis\utilities**.
- Code downloaded from external sources can be found in the folder **\external** (licence in corresponding folder).
- The documentation for each function can be found in the corresponding matlab files.

Below are the functions used for each figure in the paper. 
Letters on the left indicate the corresponding panel within each figure.

## Figure 1
    D-F) OriSelective
	G) ML_GTdecoder
## Figure 2
    B-D,G-I) OriSelective
    E,J) ML_GTdecoder
## Figure 3
    C,E) OriSelective
    D) OriSelective, OriVar
## Figure 4
    C) successRateSMI, CorrResponsive
    D) CorrResponsive, CorrSelective
    E) tempMatch_taskOri
    E) OriSelective
## Figure 5
    A) CorrResponsive, CorrSelective
    B) CorrSelective, OriSelective
    C-D) stimDiscrim, stimDiscrimVR
## Figure 6
    A-B) rwdResponsive, OriSelective, rwdOnsetData
## Figure S1
    B) successRateSMI
    C) successRateSMI, OriSelective
## Figure S2 
    OriSelective
## Figure S3
    A-B,E,F) runTime 
    C-D) runTime, OriSelective
    E) LMI
	
	
