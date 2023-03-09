# mucho-galfit

This repository has the code for running galfit in parallel on WISE images.  We are starting with a sample drawn from the Virgo Filament Survey, and we will then expand to run on a larger survey.

The galfit models will be used to estimate the size of the star-forming disk from the W3 images.  We will compare with the size of the stellar disks, measured from legacy r-band images.


# Installation
These instructions are based on these tutorials https://carpentries-incubator.github.io/python-intermediate-development/12-virtual-environments/index.html
## create a virtual environment

```
cd github/mucho-galfit
python3 -m venv venv
```

then activate environment
```
source venv/bin/activate
```

## Install python requirements

```
pip3 install -r requirements.txt
```

