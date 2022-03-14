# Threat-memory-in-human-olfactory-cortex

*Analysis scripts
  Behavior_analysis.m: for analyses of odor detection tasks (pre-conditioning/Precond, post-conditioning Day 1/Postcond, post-conditioning Day 9/Postcond2) and risk ratings (post-conditioning Day 1/Post, post-conditioning Day 9/Post2)
  
  RSA_analysis.m: for RSA analysis
    Phases - PrecEPI1: pre-conditioning;  PostEPI3: post-conditioning Day 1; P2EPI4: post-conditioning Day 9
    ROIs - r/lapc: right/left anterior piriform cortex; ppc: posterior piriform cortex; ofc: orbitofrontal cortex; amg: amygdala; hip: hippocampus
    Odor conditions - cond1: odor1; cond2: odor2; cond3: odor3; cond4: odor4; cond5: odor5
  
  Tuning_shift_analysis.m: for tuning analysis
  
  
*Behavioral data
  /Odor Detection Task contains three subfolders that correspond to data during three different phases: /Pre-conditioning, /Post-Day1-conditioning, /Post-Day9-conditioning
  
  /Odor Risk Rating contains post-conditioning risk ratings
  risk1: post-conditioning Day 1; risk2: post-conditioning Day 9
  
  /Respiration contains extracted respiration data (original files too large)
  
  /SCR contains preprocessed SCR data during conditioning


*Stimuli
  Contains images and sounds used as UCS during conditioning
  Images: Neutral (7021 - 7224); Disgust (7301-9354, vomit)
  Sounds: Neutral (N1 - N4); Disgust (D6 - D12)
  

*fMRI data
Within each individual subject folder, 
  /beta_matrices_LSA contains .mat data files (e.g. sub1_PrecEPI1_rapc_cond1.mat), which has a "betaMatrix" (15 odor trials X No. voxels) variable of each ROI and odor condition
  
  /coords contains .mat data files (e.g. sub1_lapc_maskidx.mat), which has a "coords" variable that records the MNI coordinates of the selected voxels in the ROI after applying the functional constraint


