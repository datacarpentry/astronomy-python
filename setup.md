---
title: Setup
---

## Overview

This workshop is designed to be run on your local machine. First, you will need to download the data we
use in the workshop. Then, you need to set up your machine to use the required software. Lastly, you will
download and run a Jupyter Notebook that contains test code to check that your installation was 
successful.

## Data

To start your installation process download this [zip file](https://zenodo.org/record/5637441/files/student_download.zip?download=1).
Make sure to note where in your file system are are downloading the file to.
If your file does not automatically unzip into a directory called `student_download`, you can unzip it with the following steps:
* Mac: navigate in a finder window to where you downloaded the file and double click `student_download.zip`
* Windows: navigate in a file explorer window to where you downloaded the file, right click the `student_download.zip` file, and select `Extract All`
* Linux: open a terminal and navigate to where you downloaded the zip file. Type `unzip student_download.zip`

You should now have a directory called `student_download`.
In this directory you will find files that you will use during the workshop as well as files that you will use in the remainder of the set up process.

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
in this workshop, you will need the environment file (`environment.yml`) you downloaded in the data section.

In a Terminal or Jupyter Prompt, make sure you are in the `student_download` directory, where `environment.yml` is stored, and run:

```
conda env create -f environment.yml
```

Then, to activate the environment you just created, run:

```
conda activate AstronomicalData
```

### Jupyter

We will test our environment setup using a test notebook (`test_setup.ipynb`) that you downloaded in the data section.

In a Terminal or Jupyter Prompt, make sure you are in the `student_download` directory. 
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
