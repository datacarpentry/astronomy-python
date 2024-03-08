---
permalink: index.html
site: sandpaper::sandpaper_site
---

The Foundations of Astronomical Data Science curriculum covers a range of core concepts necessary to efficiently study the ever-growing datasets developed in modern astronomy. In particular, this curriculum teaches learners to perform database operations (SQL queries, joins, filtering) and to create publication-quality data visualisations. Learners will use software packages common to the general and astronomy-specific data science communities ([Pandas](https://pandas.pydata.org), [Astropy](https://www.astropy.org), [Astroquery](https://astroquery.readthedocs.io/en/latest/) combined with two astronomical datasets: the large, all-sky, multi-dimensional dataset from the [Gaia satellite](https://sci.esa.int/web/gaia), which measures the positions, motions, and distances of approximately a billion stars in our Milky Way galaxy with unprecedented accuracy and precision; and the [Pan-STARRS photometric survey](https://panstarrs.stsci.edu/), which precisely measures light output and distribution from many stars. Together, the software and datasets are used to reproduce part of the analysis from the article ["Off the beaten path: Gaia reveals GD-1 stars outside of the main stream"](https://arxiv.org/abs/1805.00425) by Drs. Adrian M. Price-Whelan and Ana Bonaca. This lesson shows how to identify and visualize the GD-1 stellar stream, which is a globular cluster that has been tidally stretched by the Milky Way.

GD-1 is a stellar stream around the Milky Way. This means it is a collection of stars that we believe was once part of a bound clump, but the gravitational influence of the Milky Way has torn it apart and spread it over an arc that traces out its orbit on the sky.  This is interesting, because if the original bound clump was a dwarf galaxy, understanding its orbit with sufficient precision allows us to measure the mass of the Milky Way, which is very important for understanding the future and past of the Milky Way as a whole. But that is much easier to do if we have a coordinate system aligned with the stream because that makes fitting the location of the stars much easier mathematically - it becomes more linear instead of some complicated curve.  Additionally, this stream is especially interesting because it has "gaps", which have a natural interpretation as being caused by the influence of small clumps of dark matter passing near the stream. Knowing the typical rate of these gaps tells you about the typical size and density of these clumps, which turns out to be one of the best probes we have of the fine structure of dark matter.

This lesson can be taught in approximately 10 hours and covers the following topics:

- Incremental creation of complex ADQL and SQL queries.
- Using Astroquery to query a remote server in Python.
- Transforming coordinates between common coordinate systems using Astropy units and coordinates.
- Working with common astronomical file formats, including FITS, HDF5, and CSV.
- Managing your data with Pandas DataFrames and Astropy Tables.
- Writing functions to make your work less error-prone and more reproducible.
- Creating a reproducible workflow that brings the computation to the data.
- Customising all elements of a plot and creating complex, multi-panel, publication-quality graphics.

Interested in teaching these materials? We have an [onboarding video](https://www.youtube.com/watch?v=gfaNFaKIOrY) and accompanying
[slides](https://docs.google.com/presentation/d/1YosDXx1gBGpBxf6fCEaazFQwZ2dYTWgtYSdPEeD09yo/edit#slide=id.p) available to prepare Instructors to
teach this lesson. After watching this video, please contact [team@carpentries.org](mailto:team@carpentries.org) so that we can record your status as an onboarded Instructor. Instructors who have completed onboarding will be given priority status for teaching at Centrally-Organised
Data Carpentry Foundations of Astronomical Data Science workshops.

<!-- this is an html comment -->

::::::::::::::::::::::::::::::::::::::::::  prereq

## Prerequisites

This lesson assumes you have a working knowledge of Python and some previous exposure to the Bash shell.
These requirements can be fulfilled by:  
a) completing a Software Carpentry Python workshop **or**  
b) completing a Data Carpentry Ecology workshop (with Python) **and** a Data Carpentry Genomics workshop **or**  
c) independent exposure to both Python and the Bash shell.

If you're unsure whether you have enough experience to participate in this workshop, please read over
[this detailed list](instructors/prereqs.md), which gives all of the functions, operators, and other concepts you will need
to be familiar with.

In addition, this lesson assumes that learners have some familiarity with astronomical concepts, including
reference frames, proper motion, color-magnitude diagrams, globular clusters, and isochrones. Participants should bring their own laptops and plan to participate actively.


::::::::::::::::::::::::::::::::::::::::::::::::::




