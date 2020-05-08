## Low Complexity Adaptive Sonar Imaging: Window width calculation, oversampling test and MVDR comparison

| Author  | Jo Inge Buskenes  |
|:----|:----|
| Affiliation | University of Oslo / The Norwegian Defense Research Establishment (FFI) |
| License | None, but credits are appreciated |

### Documentation
| Article | Low-Complexity Adaptive Sonar Imaging |
|:----|:----|
| Authors | *J.I. Buskenes, R.E. Hansen and A. Austeng* |
| Journal | IEEE Journal of Oceanic Engineering |
| DOI     | https://doi.org/10.1109/JOE.2016.2565038 |
| Date    | January 2017 |

### Source file

The source file calcWindow3dSimple.py is a Python function for generating array tapers/windows for LCA

### Media 1: Lateral oversampling test

This animation demonstrates the effect of changing lateral oversampling factor on a set of LCA images. Upon display the
images were all bilinearly upinterpolated to the same size.

![alt text](media1_oversampling_test.gif "Oversampling test")

### Media 2: "Histogram" of typical MVDR window responses

This animation ilustrates typical spatial window responses (far-field, amplitude and phase) of the minimum variance distorionless response beamformer (MVDR, also known av Capon) in shadow, highlight and speckle regions in an image, shown as a function of subarray length L and temporal averaging K. The shipwreck is of Holmengraa located outside Horten, Norway.

At L = 1 the window responses are rectangular in all areas in the image. At L = 2 we observe window responses with slight amplitude variations in the highlight region, but with the phase being 0◦ in the illuminated seafloor sector. Already at L = 2 the MVDR is able to greatly suppress noise. This is a common observation; adding a little flexibility to adapt to the scene has a dramatic effect, but allowing full flexibility is much less significant. This shows that MVDR can be run with subarray sizes L ∈ [M/2, 5M/8], but only with temporal averaging K = 1 or above.

![alt text](media2_typical_mvdr_windows.gif "Typical MVDR windows")

