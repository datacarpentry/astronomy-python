---
layout: lesson
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---
The Foundations of Astronomical Data Science curriculum covers a range of core concepts aimed at efficient study of the 
growing datasets facing the astronomical community. In particular, this curriculum teaches learners to work with 
astronomy-specific packages (Astropy, Astroquery) and more general packages (Pandas) to perform database operations
(SQL queries, joins, filtering) to create publication-quality data visualisations. These concepts are demonstrated on the
large, all-sky, multi-dimensional dataset from the [Gaia satellite](https://sci.esa.int/web/gaia), which measures the 
positions, motions, and distances of approximately a billion stars in our Milky Way galaxy with unprecedented accuracy and
precision, and the [Pan-STARRS](https://panstarrs.stsci.edu/) photometric survey, which precisely measures light output and 
distribution from many stars. Together, these datasets are used to reproduce part of the analysis from the article 
[“Off the beaten path: Gaia reveals GD-1 stars outside of the main stream”](https://arxiv.org/abs/1805.00425) by Drs. Adrian
M. Price-Whelan and Ana Bonaca. The lessons show how to identify and visualize the GD-1 stellar stream, which is a globular
cluster that has been tidally stretched by the Milky Way.

This lesson can be taught in approximately 10 hours and covers the following topics:
* Incremential creation of complex ADQL and SQL queries, including asynchronous queries.
* Transforming coordinates between common coordinate systems.
* Working with common astronomical file formats, including FITS, HDF5, and CSV.
* Writing functions to make your work less error-prone and more reproducible.
* Customising all elements of a plot and creating complex, multi-panel, publication-quality graphics.
* AZALEE WHAT IS MISSING HERE?

<!-- this is an html comment -->

{% comment %} This is a comment in Liquid {% endcomment %}

> ## Prerequisites
> 
> This lesson assumes you have a working knowledge of Python and some previous exposure to the Bash shell. 
> These requirements can be fulfilled by:  
> a) completing a Software Carpentry Python workshop **or**  
> b) completing a Data Carpentry Ecology workshop (with Python) **and** a Data Carpentry Genomics workshop **or**  
> c) independent exposure to both Python and the Bash shell. 
> 
> If you're unsure whether you have enough experience to participate in this workshop, please read over
> [this detailed list]({{ page.root }}/prereqs/index.html), which gives all of the functions, operators, and other concepts you will need
> to be familiar with.
> 
> In addition, this lesson assumes that learners have some familiarity with astronomical concepts, including 
> SOME THINGS THA AZALEE WILL ADD HERE. Participants should bring their own laptops and plan to participate actively.
{: .prereq}

{% include links.md %}
