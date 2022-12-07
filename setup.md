---
title: Setup
---

## Overview

This workshop is designed to be run on your local machine. First, you will need to download the data we
use in the workshop. Then, you need to set up your machine to use the required software. Lastly, you will
download and run a Jupyter Notebook that contains test code to check that your installation was 
successful.

**Optional** - You may also be interested in reading the journal article that we will explore during the workshop - 
[Off the Beaten Path: *Gaia* Reveals GD-1 Stars outside of the Main Stream](https://iopscience.iop.org/article/10.3847/2041-8213/aad7b5) by 
Adrian M. Price-Whelan and Ana Bonaca.

## Data

To start your installation process download this [zip file](https://figshare.com/ndownloader/files/35777540). Move the downloaded file to your Desktop.
If your file does not automatically unzip into a directory called `student_download`, you can unzip it with the following steps:
* Mac: navigate in a finder window to your Desktop and double click `student_download.zip`
* Windows: navigate in a file explorer window to your Desktop, right click the `student_download.zip` file, and select `Extract All`
* Linux: open a terminal and navigate to your Desktop. Type `unzip student_download.zip`

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

In a Terminal, Jupyter Prompt or Anaconda Prompt, make sure you are in the `student_download` directory. 
To start Jupyter, make sure you have activated your new conda environment, then run:

```
jupyter notebook
```
The notebook should open automatically in your browser. If it does not or you wish to use a different 
browser, open this link: [http://localhost:8888](http://localhost:8888).

Now open the notebook you downloaded, `test_setup.ipynb`, and read through the instructions there. 
Make sure to run the cells that contain `import` statements.
If they work and you get no error messages, **you are ready for the workshop**.

> ## Why didn't the imports work? 
> Occasionally learners will need to take one additional step to make Jupyter run within the environment we have created. 
> If your imports fail, close Jupyter by closing its terminal, and try running the following from your Anaconda prompt (
> Terminal or otherwise):
> 
> ~~~
> python -m ipykernel install --user --name=AstronomicalData
> ~~~
> Then start Jupyter up again:
> ~~~
> jupyter notebook
> ~~~
> This time, when you open your notebook, navigate to the Kernel menu --> Change Kernel --> select **AstronomicalData** . 
> This will ensure that the relevant packages are all available. 
> You can seek installation help if this looks confusing!
{: .callout}


Please contact your instructors if you experience any problems with these installation instructions. If 
you are working through these materials independently, let us know about any problems you encounter by 
[filing an issue on the lesson's GitHub repository](https://github.com/datacarpentry/astronomy-python/issues) 
or emailing team@carpentries.org.

{% include links.md %}
