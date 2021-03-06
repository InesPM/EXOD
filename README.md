# EXOD
## EPIC-pn XMM-Newton Outburst Detector

The aim of this project is to provide a tool to spot variable sources in EPIC-pn XMM-Newton observations.
This is done by computing the variability of each pixel of the detector instead of generating lightcurves of individual sources. It allows to find faint transients for which the lightcurves are not generated by XMM-Newton's pipeline since they are too faint.

We encourage the potential users to read the users guide (EXOD_users_guide.pdf), and especially follow the tutorial presented in section 6.

If you use EXOD for your research, please acknowledge it by citing I. Pastor-Marazuela, N. A. Webb, D. T. Wojtowicz and J. van Leeuwen, 2020, A&A, 640, A124 https://arxiv.org/abs/2005.08673.

Original project "Variabilitectron" created by Damien Wojtowicz. See previous versions:
https://framagit.org/DWojtowicz/Variabilitectron
https://framagit.org/InesPM/Variabilitectron

## Tutorial

Let's set some useful parameters; the path to where the scripts are located and where we want to store our data, as well as the observation ID we want to analyse:
```
SCRIPTS=/path/EXOD/scripts
FOLDER=/path/data
obs=0652250701
```
Next we will download the raw data we need for the variability analysis, and then filter the observation:
```
bash $SCRIPTS/download_observation.sh $FOLDER $obs
bash $SCRIPTS/filtering.sh -f $FOLDER -o $obs
```
Now we are ready to perform the variability analysis:
```
python3 -W"ignore" $SCRIPTS/detector.py -path $FOLDER/$obs --render --ds9 --ds9 -tw 100 -dl 8 -bs 3
```

The whole process should take a few minutes, and it depends on the duration of the observation.
An example of the output of these commands can be found in the folder `examples`.

As the render_all.py output shows, one variable source has been detected with TW=100 s, DL=8, bs=3.



![variability](../master/example/variability_whole.png)

If at least a variable source has been detected, we can get its lightcurve and probability of constancy as follows:
```
bash $SCRIPTS/lightcurve.sh -f $FOLDER -s $SCRIPTS -o $obs -dl 8 -tw 100 -gtr 1.0 -bs 3 -id 1
```
