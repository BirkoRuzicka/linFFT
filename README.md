# linFFT
Algorithm for automatized detection and analysis of lineament structures from images using Fast Fourier Transform (FFT)

The MATLAB version of this function (linFFT) is part of the electronic supplementary files for

Ruzicka, B.R., Schroeter, M., Pack, A., and Boehnhardt, H.: Detecting and analysing geomorphological structures in images of comet 67P/Churyumov-Gerasimenko

published in the Monthly Notices of the Royal Astronomical Society, Volume 503, Issue 3, Pages 3449-3459, in May 2021 (https://doi.org/10.1093/mnras/stab618)

The entire method, particularly a step-by-step guide for using the function to analyse a larger image, is explained in detail in the article. The function is made available under a Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) license.

Technical requirements:
- linFFT was designed to run on MATLAB R2019b and requires the Image Processing Toolbox and the Signal Processing Toolbox.
- the graphical outputs (figures) generated by the function are designed to be viewed with an undocked, maximised figure-window. For optimal display, the minimum required screen resolution is 1920 x 1080.

To use linFFT for detecting and analysing linear features in an image, run these commands in MATLAB:

```
>> image = imread('yourfilename.png');
>> linFFT(image);
```

The function has several input parameters ("settings") which are explained in the function header. The function was designed with defaults for all these settings, but the user will likely need to adjust at least the parameters 'framesize' and 'threshold' according to the input image. Suitable parameters need to be found by trial and error.

The function is able to generate several different figures to visualise the workflow and/or the results. The input parameter 'dispfig' is used to decide which (if any) figure is displayed. <br>
&nbsp;&nbsp;&nbsp;&nbsp; dispfig = -1	displays a map of backtransformations along the detected lineaments (cf. Figure 8B) <br>
&nbsp;&nbsp;&nbsp;&nbsp; dispfig >= 1	displays an overview of the steps of the workflow (Figure 1) and results (see example below)<br>
&nbsp;&nbsp;&nbsp;&nbsp; dispfig =  0	no figures are generated

<br>
An example view for dispfig >= 1:<br>

![example of dispfig >= 1](https://raw.githubusercontent.com/BirkoRuzicka/linFFT/main/frame_Gh.png)

Further notes:
- linFFT is compatible with any image file type that can be read into MATLAB using imread.
- To use linFFT on a single image without dividing it into frames, set 'framesize' == image width. In this case, the input image must be square and image width must end in '01' (e.g. 101, 301, 1201).

