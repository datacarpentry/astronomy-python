---
title: Setup
---

## Overview

This workshop is designed to be run on your local machine. First, you will need to download the data we
use in the workshop. Then, you need to set up your machine to use the required software. Lastly, you will
download and run a Jupyter Notebook that contains test code to check that your installation was 
successful.

## Data





## Software

You will need to install Python, Jupyter, and some additional libraries.
[Python](http://python.org) is a popular language for
scientific computing, and great for general-purpose programming as
well. For this workshop we use Python version 3.x. 
Installing all of its scientific packages individually can be
a bit difficult, so we recommend an all-in-one installer.
We will use Anaconda.

### Anaconda
Download and install [Anaconda](https://www.anaconda.com/products/individual#anaconda-installers).

To create a new Conda environment, which includes the additional packages we will be using
in this workshop, you will need to download an environment file from the lesson 
repository.  On Mac or Linux, you can download it using `wget` on the command line: 

```
wget https://raw.githubusercontent.com/datacarpentry/astronomy-python/gh-pages/code/environment.yml
```

Or you can download the environment file by right-clicking 
[this link](https://raw.githubusercontent.com/datacarpentry/astronomy-python/gh-pages/code/environment.yml) and 
selecting "Save Link As". Make sure the file saves as `environment.yml`, not `environment.yml.txt`.

In a Terminal or Jupyter Prompt, make sure you are in the folder where `environment.yml` is stored, and run:

```
conda env create -f environment.yml
```

Then, to activate the environment you just created, run:

```
conda activate AstronomicalData
```

### Jupyter

We will test our environment setup using a test notebook. Before you launch Jupyter, download this 
notebook. On Mac or Linux, you can download it using `wget` on the command line: 

```
wget https://raw.githubusercontent.com/datacarpentry/astronomy-python/gh-pages/code/test_setup.ipynb
```

Or you can download the test notebook by right-clicking 
[this link](https://raw.githubusercontent.com/datacarpentry/astronomy-python/gh-pages/code/test_setup.ipynb) and 
selecting "Save Link As". Make sure the file saves as `test_setup.ipynb`, not `test_setup.ipynb.txt`.

To start Jupyter, make sure you have activated your new conda environment, then run:

```
jupyter notebook
```
The notebook should open automatically in your browser. If it does not or you wish to use a different 
browser, open this link: [http://localhost:8888](http://localhost:8888).

Now open the notebook you downloaded, `test_setup.ipynb`, and read through the instructions there. 
Make sure to run the cells that contain `import` statements.
If they work and you get no error messages, **you are ready for the workshop**.

Please contact your instructors if you experience any problems with these installation instructions. If 
you are working through these materials independently, let us know about any problems you encounter by 
[filing an issue on the lesson's GitHub repository](https://github.com/datacarpentry/astronomy-python/issues) 
or emailing team@carpentries.org.

{% include links.md %}
