# rcbFLOW: rCBF fitting, lineup, and observation wizard

Passning av perfusionsdata enligt Kober at al. NMR Biomed 21(8):781-92. Calculates regional cerebral blood flow (rCBF) from flow-sensitive alternating inversion recovery echo planar imaging (FAIR-EPI) acquisitions over a user-defined region of interest (ROI) mask in one brain slice per batch data set and automatically registers the resultant masked rCBF image to the average over a user-defined stack of corresponding anatomical reference image slices. Updated for MATLAB 2023b. Parallel Processing Toolbox required. 

Add to MATLAB path before calling rcbFLOW: 

rcbFLOW.m
rcbFLOW_invrecexp.m
rcbFLOW_DMTI2.m
rcbFLOW_DMTIb.m

### Inputs (by message prompt)

* batch_processing (string; '0' or '1'): '0' processes a single dataset; '1' enables batch analysis of a list of inputs defined as in the example CSV below, each additionally defined in the case of batch_processing = 0

Example batch input CSV: 

FAIR_directory_fullpath,RARE_directory_fullpath,first_RARE_slice_filename, (optional: FAIR_mask_file_fullpath, RARE_mask_file_fullpath)
C:\Users\kswanberg\CBFDir1,C:\Users\kswanberg\RAREDir1,MRIm12.dcm,C:\Users\kswanberg\FAIRmasks\1.dcm,C:\Users\kswanberg\RAREmasks\1.dcm,
C:\Users\kswanberg\CBFDir2,C:\Users\kswanberg\RAREDir2,MRIm14.dcm,C:\Users\kswanberg\FAIRmasks\2.dcm,C:\Users\kswanberg\RAREmasks\2.dcm,
C:\Users\kswanberg\CBFDir3,C:\Users\kswanberg\RAREDir3,MRIm10.dcm,C:\Users\kswanberg\FAIRmasks\3.dcm,C:\Users\kswanberg\RAREmasks\3.dcm,
C:\Users\kswanberg\CBFDir75,C:\Users\kswanberg\RAREDir75,MRIm12.dcm,C:\Users\kswanberg\FAIRmasks\75.dcm,C:\Users\kswanberg\RAREmasks\75.dcm,

* auto_masking for FAIR-EPI (string; '0', '1', or '2'): '0' requires the user to manually select the analysis ROI on a sample FAIR-EPI image, while '1' uses a multi-step segmentation process to do so without manual user input. '2' enables loading of previously traced masks saved as dcms. 

* auto_masking for RARE (string; '0', '1', or '2'): '0' requires the user to manually select the analysis ROI on a sample RARE image, while '1' uses a multi-step segmentation process to do so withoutmanual user input. '2' enables loading of previously traced masks saved as dcms. 

* anti_aliasing (string; '0' or '1'): '0' will not implement anti-aliasing correction on rCBF images, while '1' will. 

### Outputs (per every dataset analyzed, by file export to single new directory): 
     
* xxx_rCBF.dcm (DCM): Masked rCBF slice calculated from FAIR-EPI
* xxx_T1c.dcm (DCM): Masked T1c slice calculated from FAIR-EPI
* xxx_RARE_uncorrected.dcm (DCM): Averaged anatomical reference slice
* xxx_RARE_flipped_resized.dcm (DCM): Flipped and resized anatomical reference slice
* xxx_RARE_flipped_resized_registered.dcm (DCM): Flipped, resized, and registered anatomical reference slice
* FAIR_EPI_mask_origsize_xxx.dcm (DCM): Mask for FAIR-EPI (rCBF and T1c) slice
* RARE_mask_origsize_xxx.dcm (DCM): Mask for RARE reference slice
* RARE_mask_resized_xxx.dcm (DCM): Mask for RARE reference slice resized for overlay on rCBF and T1c slices 
* xxx_Fig1 (PNG): Figure 1. Difference between selective and global inversion conditions for each inversion time (TI) 
* xxx_Fig2 (PNG): Figure 2. FAIR-EPI ROI mask (either manual, automatic, or loaded) for rCBF calculation 
* xxx_Fig3 (PNG): Figure 3. RARE ROI mask (either manual, automatic, or loaded) for registration
* xxx_Fig4 (PNG): Figure 4. Masked calculated R1c image
* xxx_Fig5 (PNG): Figure 5. Masked calculated rCBF image 
* xxx_Fig6 (PNG): Figure 6. Masked calculated T1c image
* xxx_Fig7 (PNG): Figure 7. Averaged anatomical reference image 
* xxx_Fig8 (PNG): Figure 8. Flipped anatomical reference image 
* xxx_Fig9 (PNG): Figure 9. Masked calculated rCBF image overlaid on averaged anatomical reference image 
* xxx_Fig10 (PNG): Figure 10. Masked rCBF image overlaid on flipped, registered, and resized averaged anatomical reference image        

### Citation 

Work that employed code from rcbFLOW can cite it as follows: 

Gottschalk, M. and Swanberg, K.M. (2024). rcbFLOW: rCBF fitting, lineup, and observation wizard v. 1.0. Source code. https://github.com/kswanberg/rcbFLOW.


### Developer

Please send comments and questions to [Kelley Swanberg](mailto:kelley.swanberg@med.lu.se). 