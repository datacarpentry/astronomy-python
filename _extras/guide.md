---
title: "Instructor Notes"
---

## Instructor notes

### Overview
This lesson guides learners through analyzing data from a large database. Scientifically, we are identifying stars in GD-1, a stellar stream in the Milky Way (creating Figure 1 in "[Off the beaten path: Gaia reveals GD-1 stars outside of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian Price-Whelan and Ana Bonaca.). The first part of this lesson (1-5) shows learners how to prototype a query, starting by querying a subset of the sky we ultimately want and then building up stronger and stronger filters locally. With our filters in place, lesson 6 performs the full query remotely, giving us a dataset to visualize in lesson 7. Lesson 7 demonstrates best practices and tips and tricks to efficiently and effectively visualize data.

Because this lesson follows a single dataset throughout, its easy for learners (and instructors) to lose sight of the bigger picture and focus instead on the scientific goals. At the beginning of each lesson it is recommended that the instructor discuss both the scientific goal of the lesson (with frequent references to Figure 1) and highlight the big picture skills that we hope each student takes away from the lesson, beyond the specific science case. At the end of the lesson the instructors should recap the same information, highlighting the best practices covered. TODO: reference slide show.
### Decisions Made
* We explicitly recommend using Jupyter notebooks rather than Jupyter lab as we've found them more reliable across different platforms
### Introduction

* Astronomical datasets are getting larger and larger. As a result, the days of downloading a raw image and performing the end to end analysis are numbered. Instead, already reduced data (and simulations) are being stored in databases. If you know how to ask these databases the right questions, they are an incredibly powerful tool. The ability to interact with these large databases is fundamental to the future of astronomy as well as a skill that is highly valued outside of astronomy.

* Motivate this statement with your science whether it is transients, exoplanets, or weak lensing (or something else).

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
* Share links to the lesson with learners. There is a lot of typing and if a student gets stuck debugging a typo, having the lesson can allow them to catch up. Encourage them to type along and only use the lesson as a back up.
* It is recommended that each day be a new notebook (Note from ABD: I think there's a tradeoff here.  Breaking a complex analysis into multiple notebooks is generally good practice, but for the workshop it might cause more problems than it is worth.)
* It is recommended that you do not repeat variable names and that you do not modify queries in place. If a student misses a step then they don't have an easy way to determine that they don't have the same code that you have. (Note from ABD: I might soften this recommendation and suggest that modifying a query in place is a good way to develop a query incrementally, but that they might want to name a few milestones to refer back to)
* Although learners should have seen functions in the SWC python lesson, its always good to reiterate the process of prototyping what you want to do. Then if you find yourself wanting to do it over and over, make it into a function (a good rule of thumb is if you're copy and pasting, then think about whether you should be using a function)
* Highlight throughout the lesson what commands are ADQL specific and what are general SQL commands

* Errors will happen throughout the workshop. In an online setting it can be especially useful for learners to share their screen with the class why you troubleshoot together. This gives other learners the opportunity to see how you troubleshoot and to see common errors. However, be aware of the time and don't derail the entire lesson if someone is having an edge case issue.

### Lesson 1: Queries

* At the beginning of this lesson, possible points of confusion are (1) there are two functions with similar named, `load_tables` and `load_table`, and neither of them actually downloads a table; rather, they download metadata about the tables.

* Note that even though we don't use the parallax column for our analysis, we include it in the results for use in the data exploration section

* For `query2` (and most subsequent queries) if you get exactly 2000 results, it's probably because you used `launch_job` rather than `launch_job_async`. This is a good opportunity to remind learners: if you get exactly 2000 results, check that you did an async query and not a synchronous query

* One of the challenges of debugging queries is that `astroquery` provides basically no debugging information other than a generic error message.  So it is important to emphasize careful checking of queries and incremental development: starting with a known-good query and making small changes.

* throughout the lessons we check the type of our output. This is a good best practice to make sure that you're getting what you think you're getting. You can remind learners of this when it comes up.

### Lesson 2: Coordinates and units
* Highlight how the WHERE and CONTAINS statements interact - that CONTAINS returns a 1 when something is contained within the defined shape and you're checking WHERE 1=CONTAINS and this is true when something is located in a region
* learners struggle with what we're transforming between and why. Repeat this as many times as possible. The Gaia catalog is in a universal frame (ICRS) but its easier to visualize and create filters specific to GD-1 when the axes of our reference frame are aligned with the GD-1 stream direction. So we will be switching between them often. To the GD-1 frame to build our filter then to the ICRS frame to actually query the database, then back to the GD-1 frame to visualize our results
* Emphasize the basic structure of the SkyCoord object: location along axis 1, location along axis 2, frame. For example in the ICRS frame, axis 1 is `ra` and axis 2 is `dec`. In the GD-1 frame, axis 1 is `phi1` and axis 2 is `phi2`.
* Mention that this is a benefit of using a unified framework like astropy. You can build a custom frame object and then have access to all of the other astropy tools that deal with coordinate transformations
* When you define phi1_min, phi1_max, etc go back to Figure 1 and show learners the region you are defining, connecting the min and max values to the coordinates in the GD-1 frame. This is another place you can mention the benefit of using the GD-1 frame is that we can define a rectangle around the stream.

### Lesson 3: Proper motion 
* This is a really long lesson. But, by the end of it you will have prototyped making the first row of Figure 1. Make sure to keep your eyes on the prize. This is likely the last lesson of Day 1
* Transforming back: the discussion about the motivation to transform back to GD-1 frame is a great time to have student build up some intuition about the proper motion selection and how this relates to the physical picture. Remind them again that GD-1 is a globular cluster that is being pulled into a stream along the phi1 direction. This means that GD-1 stars should have non-zero motion along the phi1 direction and motion close to 0 in the phi2 direction. 
* Don't let learners get too hung up on setting the distance and radial velocity to constants. For distance, the parallax is so close to 0 that they are unreliable and as mentioned, the radial velocity is to avoid errors.
* The plot at the end of the reflex correction is a good example of an intermediate step to determining the best filter for GD-1 stars. You can imagine doing this data exploration yourself and trying first a purely spatial plot, then realizing that there are too many non-GD-1 stars included and that you need another way to filter out foreground stars. Proper motion!
* This lesson takes a slight detour to introduce a few features of Pandas DataFrames. To keep this connected to the story, you can talk about how data exploration is an important part of prototyping your query and making sure you are getting the results you expect in the format you expect them in.
* Starting with selecting the centerline, we do a series of filters on different data frames. Take a minute before you teach this section to make sure you understand what each one represents. We use results_df to build centerline_df. We use centerline_df to determine proper motion limits. We use the proper motion limits determined from centerline_df to select GD-1 stars from results_df. This is selected_df.

* It is also worth noting that we are selecting stars close to the centerline to identify the proper motion cut, so its ok if we exclude some GD-1 stars here. In this case we prefer a more pure sample to a more complete sample. Once we have the proper motion limits from our pure sample (centerline_df), we'll include the full spatial region (results_df) and get all of the GD-1 stars (selected_df).

* Learners may be concerned that we're writing a function for later that doesn't fill the full canvas (when we set `axis('equal')`). You can reassure them that in episode 7 we'll take care of this by learning how to set the figure size, the subplot size, and we'll be using a larger spatial area (that we define in episode 4).

* Notice that the first time we use `DataFrame.to_hdf`, we use the `w` argument to indicate that we want to create a new, empty HDF Store.  For all subsequent uses, we should *not* use the `w` argument, so that we add new Datasets to the existing Store, rather than starting over.
* At the end of Day 1, if a student is lost or struggled with the end of this lesson, point them to the static version of the HDF5 files (TODO: decide where this lives) so that they can read it in with everyone else on Day 2.

* If a learner notices that there is a `/` in front of the keys at the end of the lesson here's what they mean: `/` indicates that the keys  are at the top level of the Dataset hierarchy, and not in a named "group". In future lessons we will add a few more Datasets to this file, but
not so many that we need to organize them into groups.

### Lesson 4: Coordinate transformation and selection
* This is likely the first lesson of Day 2. You should start by having learners start a new notebook. They will need to read in the data they saved yesterday to HDF5 files and the functions they wrote yesterday.

* This lesson is faster than it looks because learners have seen much of the content already

* As in previous lessons, a challenge for learners here is keeping track of what we are doing and why.  Surface periodically to remind them where we are. The big focus of this lesson is that now that we have determined how to filter based on proper motion and eliminated most of the contaminating stars from our prototype, we can repeat the query, this time on the full spatial extent of GD-1 rather than the small portion we used for our prototype.

* Another difficulty is that we re-use several functions from previous notebooks.  If learners are working in the same notebook through multiple lessons, they might already have them.  Otherwise it is probably best to instruct them to [copy and paste from here](https://allendowney.github.io/AstronomicalData/04_select.html) rather than retyping. TODO: Decide how we want to handle this - import?

* This lesson uses the following idiom several times

```
x = ...
y = ...
plt.plot(x, y)
```

This idiom violates the recommendation not to repeat variables names, but since they are defined and used immediately, it might be ok.  If you don't like it, you can inline the expressions.

* The use of ConvexHull in this lesson is a bit of a hack, but some learners will find the computation geometry functions in SciPy useful.

* Make sure to review the query from Lesson 2 (this will happen naturally when you have to retype it) and remind learners what each filter does and why we're using it

* Notice that the definitions of `phi1_min`, `phi1_max`, etc.  are different in this lesson.  Because we are adding more filters, we can select a bigger region without exceeding resource limits.  If learners don't get as many "candidates" as expected, they might be using the old values of these bounds.

* Defining the new region is a good opportunity to go back to the figure and connect these coordinates to the physical picture of GD-1 (e.g. before we were only selecting from this region now we're expanding it to this whole region)

* When you get to `make_dataframe`, you might want to copy and paste it from the notes, rather than retyping.

* Learners may ask why we are using `dict(key1=value1, key2=value2)` rather than `{key1:value1, key2:value2}`. We are using the `dict` syntax so that the key value pairs look like keyword arguments. This may simplify the understanding for learners who are less familiar with dictionaries.

* Recall that the first time we use `DataFrame.to_hdf`, we use the `w` argument to indicate that we want to create a new, empty HDF Store.  For all subsequent uses, we should *not* use the `w` argument, so that we add new Datasets to the existing Store, rather than starting over.

* If learners struggle with the query exercise - remind them that we are defining a polygon in the same way we did before (with the same syntax) but looking at proper motion instead of space.

* If learners get a different number of candidates, its likely they mistyped something.
### Lesson 5: Joining tables

* The early part of this lesson brings back a lot of best practices that learners learned in earlier lessons (e.g. exploring tables, returning the top 5 rows) this is a great opportunity to high-light these and remind learners that they've seen this before and why we're doing it.

* Although we don't immediately use the column definitions in panstarrs1_best_neighbour, we will come back to them in the data exploration part of this lesson

* Throughout this lesson continue to come back to the central theme of starting with a simple query and building up layers of complexity, testing each layer as you go

* The Pan-STARRS join exercise is likely to feel scary to a lot of students - they have only seen one example. Emphasize the building blocks of a join e.g. FROM table1 JOIN table2 ON table1.common_column=table2.common_column. It is also worth pointing out that common_columns will not have the same names in this case

* In case a student asks about the extra `unnamed` column in the CSV section here's the explanation. You may notice that all Pandas `DataFrame`s have an index column which was not part of the original table definition. 
This essentially numbers each row. When we write a `DataFrame` in any other format, the index gets treated like a bonafide column.
For this reason when we write a CSV file and then read it back into a `DataFrame` the index column gets written as an `unnamed` column and then when it is read back in, another index column is created leading to two extraneous columns. 


### Lesson 6: Photometry
* It is easy in this lesson to lose track of the main point: that we want to define a polygon around the main sequence of GD-1 so we can further hone our sample of candidate GD-1 stars. As we spend time on the isochrone, creating the polygon, etc make sure to come back to this big picture often.

* The key take away from the CMD presentation is that GD-1 is a globular cluster which means all of the stars formed at the same time. Therefore we expect the stars in GD-1 to follow a single, tight isochrone, the main sequence of which we can easily identify.

* If you (or someone in the class) are interested in how to calculate the Isochrone from the MIST models, the code can be found TODO: finish this sentence

* In the original paper, they use an idiosyncratic function to define the boundaries of the isochrone filter. The intention is to define a polygon that gets wider as g increases, to reflect increasing uncertainty. For this exercise we will be using a simplified version which is a constant offset from the isochrone.

* Turning the selection we've defined into a polygon (making a loop) can be challenging to explain verbally but fairly simple visually. This is a nice time to return to the slides or draw a picture

* In the manipulation we do to create the polygon, it can get lost that left_color, right_color, and color_loop are all x values for the polygon (and that g and mag_loop are the y values). It is worth coming back to this over and over as you introduce each variable.

* Learners might express concern that the polygon we use to select candidate stars is not closed.  That's ok; the `contains_points` method treats the polygon as if it is closed (by connecting the last point to the first).

* Learners may notice that `contains_points` is not actually a Polygon property. Polygon inherits this method from Patch.

* We mention usetex in the discussion of TeX markup in matplotlib. This option does not work for all operating system/backend combinations (in our first alpha pilot some notebooks required restarting after trying to use this option). But, one error that can be fixed is `LaTeX Error: File 'type1cm.sty' not found.` Learners with this error might have to install a package that contains the fonts LaTeX needs.  On some systems, the packages `texlive-latex-extra` or `cm-super` might be what they need.  [See here for more help with
 this](https://stackoverflow.com/questions/11354149/python-unable-to-render-tex-in-matplotlib).


* Learners may ask why we're initializing an empty array and then creating the columns on the fly. DataFrame initialize with arrays of rows rather than columns, so this is the easiest way without having to do some array manipulation. See https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html for more details.

### Lesson 7: Visualization

* This lesson is largely about showing learners selected capabilities that will make their lives easier. Matplotlib is a huge package with infinite flexibility - this is in no way complete but hopefully gives them some barebones tools to work with and inspired them to explore further. 

* In addition to the links provided in the lesson it is nice to show learners the [list of plotting commands](https://matplotlib.org/stable/api/pyplot_summary.html) and the [examples gallery](https://matplotlib.org/stable/gallery/index.html). 

{% include links.md %}
