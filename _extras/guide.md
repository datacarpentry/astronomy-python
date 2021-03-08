---
title: "Instructor Notes"
---

## Instructor notes

### Overview
This lesson guides students through analyzing data from a large database. Scientifically, we are identifying stars in GD-1, a stellar stream in the Milky Way (creating Figure 1 in TODO: link to paper). The first part of this lesson (1-5) shows students how to prototype a query, starting by querying a subset of the sky we ultimately want and then building up stronger and stronger filters. With our filters in place, lesson 6 performs the full query giving us a dataset to visualize in lesson 7. Lesson 7 demonstrates best practices and tips and tricks to efficiently and effectively visualize data. 

Because this lesson follows a single dataset throughout, its easy for students (and instructors) to lose sight of the bigger picture and focus instead on the scientific goals. At the beginning of each lesson it is recommended that the instructor discuss both the scientific goal of the lesson (with frequent references to Figure 1) and highlight the big picture skills that we hope each student takes away from the lesson, beyond the specific science case. At the end of the lesson the instructors should recap the same information, highlighting the best practices covered. TODO: reference slide show.

### Introduction
* Astronomical datasets are getting larger and larger. As a result, the days of downloading a raw image and performing the end to end analysis are numbered. Instead, already reduced data (and simulations) are being stored in databases. If you know how to ask these databases the right questions, they are an incredibly powerful tool. The ability to interact with these large databases is fundamental to the future of astronomy as well as a skill that is highly valued outside of astronomy.
* Motivate this statement with your science whether it is transients, exoplanets, or weak lensing (or something else)
* Introduce the tools we will be teaching:
    - SQL
    - Python
    - Astropy
    - Astroquery 
    - Pandas
    - Matplotlib
* Note that this is a balance of astronomy specific tools which allow us to efficiently interact with our specific data and broadly used tools (including industry)
* Introduce the skills (and best practices surrounding those skills) we will be teaching:
    - Developing and testing SQL and Python
    - Working with remote databases
    - Storing data locally
    - Validating data and analysis pipelines
    - Visualizing results
* Introduce the dataset:
    - Science
        - GD-1 is a stream of stars in the Milky Way that is believed to have been a globular cluster that was stretched out by the gravitational potential of the Milky Way
        - Identify stars in GD-1 using Gaia and Pan-STARRs data which is stored in the ESA Gaia database
        - Give a one sentence description of the Gaia and Pan-STARRs surveys
        - Reproduce Figure 1 of Price-Whelan & Bonaca, 20??
    - Why this dataset?

{% include links.md %}