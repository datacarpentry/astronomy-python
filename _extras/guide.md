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
* Introduce the tools we will be teaching and why you think these are important for astronomers to know:
    - SQL: general purpose language to interact with databases - used widely in astronomy and industry
    - Python: general purpose programming language - used widely in astronomy and industry
    - Astropy: package with common astronomical functions and objects to facilitate analysis of astronomical data
    - Astroquery: python interface to SQL-flavored astronomical databases (this might be a good place to show all of the catalogs that astroquery supports) 
    - Pandas: python package to deal with arrays and databases - used widely by astronomy and industry
    - Matplotlib: python visualization package - used widely by astronomy and industry
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
        - Uses data from two broadly used surveys which are accessible through astroquery (which can be used to access a large number of astronomical surveys)
        - Two large catalogs make this too big to be done locally
        - Compelling science case: this analysis extends the known structure of GD-1, discovers new substructure, and identifying a candidate progenitor
        - Compelling use case: the time from data release to publication on the arXiv was less than a week. This models how efficiently this workflow can be and also how one could set up the tools before a data release to quickly turn around a paper
* Show Figure 1 from Price-Whelan & Bonaca and describe the workflow that we're going to go through
    - First we are going to select stars purely based on their location in the sky near where we know GD-1 exists. Until we figure out our full query and have a more robust filter, we're going to limit ourselves to a very small spatial subset of these stars so that we are downloading a manageable dataset. We'll also do a rough filter on color (to avoid nearby faint foreground stars which are likely to be red) and on parallax (again to avoid foreground stars)
    - This is a good time to introduce the GD-1 frame (coordinates phi1 and phi2) where phi1 is parallel to the stream direction and phi2 is perpendicular to the stream direction. This is the easiest frame to visualize the stream in and often to set query constraints.
    - The first filter we apply will be a filter in proper motion. Because GD-1 is a stretched out globular cluster, the stars should all be moving along the stream (parallel to the x-axis). If we plot the proper motion of all of the stars in the spatial region we selected, we should see a clump of stars with no motion in phi2 and a non-zero motion in phi1. These stars are GD-1. The first filter we apply will be to select for stars with a proper motion that falls in this region of proper motion space (top right plot). The stars selected are shown on the top left plot. This is all done with the Gaia catalog
    - This allows us to see GD-1, but there is still a lot of contamination. To further refine our filter, we are going to use the fact that we believe GD-1 was a globular cluster and therefore all of the stars formed at the same time. This means we can build a g-i color magnitude diagram (CMD) and select stars that fall on a single isochrone corresponding to the age of GD-1. To build this color magnitude diagram we will need the g and i magnitudes from the Pan-STARRS catalog. The CMD of all of the stars selected based on proper motion is shows in the bottom right panel along with the best-fit isochrone. We will then select stars only around the main sequence of the isochrone. The results of this selection can be seen in the bottom left panel. You can see here that GD-1 stands out much more prominently.
    - Once we have prototyped this query on a small subset of the spatial extent of GD-1, we will re-run the query on the full extent and recreate this figure.
* For each lesson start with an overview - often its good to show Figure 1 and talk about what you just did and what you are going to do. Be sure to highlight the learning objectives and discuss how these tools and skills are applicable beyond this dataset and science case.
* End each lesson with a summary of what you just did and going over the best-practices. Take this as an opportunity to connect the skills used for this specific science case to the broader skill set that can be applied to any science case.
* Share links to the lesson with students. There is a lot of typing and if a student gets stuck debugging a typo, having the lesson can allow them to catch up. Encourage them to type along and only use the lesson as a back up.
* It is recommended that each day be a new notebook (they get really slow)
* It is recommended that you do not repeat variable names and that you do not modify queries in place. If a student misses a step then they don't have an easy way to determine that they don't have the same code that you have.
* Although students should have seen functions in the SWC python lesson, its always good to reiterate the process of prototyping what you want to do. Then if you find yourself wanting to do it over and over, make it into a function (a good rule of thumb is if you're copy and pasting, then think about whether you should be using a function)
* Highlight throughout the lesson what commands are ADQL specific and what are general SQL commands
### Lesson 1: Queries
* remind students: if you get exactly 2000 results, check that you did an async query and not a synchronous query
### Lesson 2: Coordinates and units
* Highlight how the WHERE and CONTAINS statements interact - that CONTAINS returns a 1 when something is contained within the defined shape and you're checking WHERE 1=CONTAINS and this is true when something is located in a region
* Students struggle with what we're transforming between and why. Repeat this as many times as possible. The Gaia catalog is in a universal frame (ICRS) but its easier to visualize and create filters specific to GD-1 when the axes of our reference frame are aligned with the GD-1 stream direction. So we will be switching between them often. To the GD-1 frame to build our filter then to the ICRS frame to actually query the database, then back to the GD-1 frame to visualize our results
* Emphasize the basic structure of the SkyCoord object: location along axis 1, location along axis 2, frame. For example in the ICRS frame, axis 1 is ra and axis 2 is dec. In the GD-1 frame, axis 1 is phi1 and axis 2 is phi2
* Mention that this is a benefit of using a unified framework like astropy. You can build a custom frame object and then have access to all of the other astropy tools that deal with coordinate transformations
* When you define phi1_min, phi1_max, etc go back to Figure 1 and show learners the region you are defining, connecting the min and max values to the coordinates in the GD-1 frame. This is another place you can mention the benefit of using the GD-1 frame is that we can define a rectangle around the stream.
### Lesson 3: Proper motion 
### Lesson 4: Coordinate transformation and selection
### Lesson 5: Joining tables
### Lesson 6: Photometry
### Lesson 7: Visualization
{% include links.md %}