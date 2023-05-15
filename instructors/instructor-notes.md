---
title: Instructor Notes
---

## Instructor notes

### Overview

This lesson guides learners through analyzing data from a large database. Scientifically, we are identifying stars in GD-1, a stellar stream in the Milky Way (creating Figure 1 in "[Off the beaten path: Gaia reveals GD-1 stars outside of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian Price-Whelan and Ana Bonaca). The first part of this lesson (episodes 1-6) shows learners how to prototype a query, starting by querying a subset of the sky we ultimately want and then building up stronger and stronger filters locally. With our filters in place, episode 7 performs the full query remote on the remote database, giving us a dataset to visualize in episode 8. Episode 8 demonstrates best practices and tips and tricks to efficiently and effectively visualize data.

Because this episode follows a single dataset throughout, its easy for learners (and instructors) to lose sight of the bigger picture and focus instead on the scientific goals or individual commands. At the beginning of each episode it is recommended that the instructor discuss both the scientific goal of the episode (with frequent references to Figure 1) and highlight the big picture skills that we hope each student takes away from the episode, beyond the specific science case. At the end of the episode the instructors should recap the same information, highlighting the best practices covered.

We have an [onboarding video](https://www.youtube.com/watch?v=gfaNFaKIOrY) and accompanying
[slides](https://docs.google.com/presentation/d/1YosDXx1gBGpBxf6fCEaazFQwZ2dYTWgtYSdPEeD09yo/edit#slide=id.p) available to prepare Instructors to
teach this lesson. After watching this video, please contact [team@carpentries.org](mailto:team@carpentries.org) so that we can record your status as an onboarded Instructor. Instructors who have completed onboarding will be given priority status for teaching at Centrally-Organised
Data Carpentry Foundations of Astronomical Data Science workshops.

### Timing

Unfortunately these episodes are incredibly cumulative and there is not much that can be cut along the way.
If you are running short on time, we recommend eliminating or condensing these sections:

- cutting or stream-lining parts of episode 8 about plotting. This is the most significant way you can make up time.
- stream-line some of the CSV section at the end of episode 6, for instance skip the read back section which demonstrates how extra columns can slip into CSV files
- cut the "Checking the Match" section of episode 6
- cut the part where we check the size of the HDF5 file in episode 4
- cut the brief discussion of context managers in episode 4

### Decisions Made

- We explicitly recommend using Jupyter notebooks rather than Jupyter lab as we have found them more reliable across different platforms and they require less up front explanation.
- As explained in the episodes, we think that it is valuable for learners to know a little bit about both Astropy `Tables` and Pandas `DataFrames` as they both have strengths and weaknesses as both an astronomy software standard and an industry standard.
- There is a lot of typing in this curriculum. To best keep with Carpentries live coding practices - we recommend that you copy and paste from previous cells as needed and not provide pre-filled notebooks. Provide learners with links to the on line lessons and reassure them that they can jump back in at the beginning of any episode if they are too overwhelmed an need a break.

### Introduction

- Astronomical datasets are getting larger and larger. As a result, astronomical analysis is shifting away from working with raw imaging and towards analyzing already reduced data (and simulations) that are stored in large databases. If you know how to ask these databases the right questions, they are an incredibly powerful tool. The ability to interact with these large databases is fundamental to the future of astronomy as well as a skill that is highly valued outside of the astronomical community.

- Motivate the skills taught in this curriculum with your science whether it is transients, exoplanets, or weak lensing (or something else).

- Introduce the tools we will be teaching and why you think these are important for astronomers to know:
  
  - SQL: general purpose language to interact with databases - used widely in astronomy and industry
  - Python: general purpose programming language - used widely in astronomy and industry
  - Astropy: package with common astronomical functions and objects to facilitate analysis of astronomical data
  - Astroquery: python interface to SQL-flavored astronomical databases (this might be a good place to show all of the catalogs that astroquery supports)
  - Pandas: python package to deal with arrays and databases - used widely by astronomy and industry
  - Matplotlib: python visualization package - used widely by astronomy and industry

- Note that this is a balance of astronomy specific tools which allow us to efficiently interact with our specific data and broadly used tools which can prepare learners for a career outside of academia.

- Introduce the skills (and best practices surrounding those skills) we will be teaching:
  
  - Developing and testing SQL and Python
  - Working with remote databases
  - Storing data locally
  - Validating data and analysis pipelines
  - Visualizing results

- Introduce the dataset:
  
  - Science
    - GD-1 is a stream of stars in the Milky Way that is believed to have been a globular cluster that was stretched out by the gravitational potential of the Milky Way
    - Identify stars in GD-1 using Gaia and Pan-STARRs data which is stored in the ESA Gaia database
    - Give a one sentence description of the Gaia and Pan-STARRs surveys
    - Reproduce Figure 1 of [Price-Whelan \& Bonaca, 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...863L..20P/abstract)
  - Why this dataset?
    - Uses data from two broadly used surveys which are accessible through astroquery (which can be used to access a large number of astronomical surveys)
    - Two large catalogs make this too big to be done locally
    - Compelling science case: this analysis extends the known structure of GD-1, discovers new substructure, and identifies a candidate progenitor
    - Compelling use case: the time from data release to publication on the arXiv was less than a week. This models how efficiently this workflow can be and also how one could set up the tools before a data release to quickly turn around a paper

- Show Figure 1 from Price-Whelan \& Bonaca and describe the workflow that we are going to go through
  
  - First we are going to select stars purely based on their location in the sky near where we know GD-1 exists. Until we figure out our full query and have a more robust filter, we are going to limit ourselves to a very small spatial subset of these stars so that we are downloading a manageable dataset. We will also do a rough filter on color (to avoid nearby faint foreground stars which are likely to be red) and on parallax (again to avoid foreground stars)
  - This is a good time to introduce the GD-1 reference frame (coordinates phi1 and phi2) where phi1 is parallel to the stream direction and phi2 is perpendicular to the stream direction. This is the easiest frame to visualize the stream in and often to set query constraints.
  - The first filter we apply will be a filter in proper motion. Because GD-1 is a stretched out globular cluster, the stars should all be moving along the stream (parallel to the x-axis). If we plot the proper motion of all of the stars in the spatial region we selected, we should see a clump of stars with no motion in phi2 and a non-zero motion in phi1. These stars are GD-1. The first filter we apply will be to select for stars with a proper motion that falls in this region of proper motion space (top right plot). The stars selected are shown on the top left plot. We will first identify the range of proper motion we are interested in locally, then we will develop this into a query that can be performed on the Gaia catalog
  - This allows us to see GD-1, but there is still a lot of contamination from non-GD-1 sources. To further refine our filter, we are going to use the fact that we believe GD-1 was a globular cluster and therefore all of the stars formed at the same time. This means we can build a g-i color magnitude diagram (CMD) and select stars that fall on a single isochrone corresponding to the age of GD-1. To build this color magnitude diagram we will need the g and i magnitudes from the Pan-STARRS catalog. The CMD of all of the stars selected based on proper motion is shows in the bottom right panel of Price-Whelan \& Bonaca, 1998 Figure 1, along with the best-fit isochrone. We will then select stars only around the main sequence of the isochrone. The results of this selection can be seen in the bottom left panel of Figure 1. You can see here that GD-1 stands out much more prominently.
  - Once we have prototyped this query on a small subset of the spatial extent of GD-1, we will re-run the query on almost the full extent and recreate Figure 1.

- For each episode start with an overview - often its good to show Figure 1 and talk about what you just did and what you are going to do. Be sure to highlight the learning objectives and discuss how these tools and skills are applicable beyond this dataset and science case.

- End each episode with a summary of what you just did and going over the best-practices. Take this as an opportunity to connect the skills used for this specific science case to the broader skill set that can be applied to any science case.

- Share links to the episode with learners. There is a lot of typing and if a student gets stuck debugging a typo, having the episode can allow them to catch up. Encourage them to type along and only use the episode as a back up.

- It is recommended that at a minimum each day be a new notebook to avoid a single notebook that is so large that it become slow to interact with. While it would be great to demonstrate the best practice of breaking a complex analysis into multiple notebooks by starting each new episode in a new notebook, this will add overhead to your teaching time both to launch a new notebook and get it set up with the previously defined variables, so use your best judgement based on the level of your audience. If you decide to start each episode in a new notebook you can use the "Starting from this Episode" call out at the beginning of each episode to import and define everything you need from previous episodes. If you are using files from the `student_download` directory be sure to add the correct path `student_download/backup-data/<filename>`.

- It is recommended that you do not repeat variable names and that you do not modify queries in place. If a student misses a step then they do not have an easy way to determine that they do not have the same code that you have.

- Although learners should have seen functions in the SWC python episode, its always good to reiterate the process of prototyping what you want to do. Then if you find yourself wanting to do it over and over, make it into a function (a good rule of thumb is if you are copy and pasting, then think about whether you should be using a function).

- Highlight throughout the episode what commands are ADQL specific (`CONTAINS`, `POINTS`, `POLYGON`) and what are general SQL commands (everything else).

- Errors will happen throughout the workshop. In an online setting it can be especially useful for learners to share their screen with the class while you troubleshoot together. This gives other learners the opportunity to see how you troubleshoot and to see common errors. However, be aware of the time and do not derail the entire episode if someone is having an edge case issue.

- Think about how to communicate exercises to learners, especially in an online setting. In our first beta workshop both instructors and learners were really happy with [hackmd](https://hackmd.io).

### Episode 1: Queries

- At the beginning of this episode, possible points of confusion are (1) there are two functions with similar named, `load_tables` and `load_table`, and neither of them actually downloads a table; rather, they download metadata about the tables.

- Note that even though we do not use the parallax column for our analysis, we include it in the results for use in the data exploration section

- For `query2` (and most subsequent queries) if a learner get exactly 2000 results, it is probably because they used `launch_job` rather than `launch_job_async`. This is a good opportunity to remind learners: if you get exactly 2000 results, check that you did an async query and not a synchronous query

- One of the challenges of debugging queries is that `astroquery` provides basically no debugging information other than a generic error message.  So it is important to emphasize careful checking of queries and incremental development: starting with a known-good query and making small changes.

- Throughout the episodes we check the type of our output. This is a good best practice to make sure that you are getting what you think you are getting. You can remind learners of this when it comes up.

### Episode 2: Coordinates and units

- Highlight how the `WHERE` and `CONTAINS` statements interact - that `CONTAINS` returns a `1` when something is contained within the defined shape and you are checking where `CONTAINS` equals `1`. So the statement `WHERE 1=CONTAINS` is true when something is located in a region.

- Learners struggle with what reference frames we are transforming between and why. Repeat this as many times as possible. The Gaia catalog is in a universal frame (ICRS) but it is easier to visualize and create filters specific to GD-1 when the axes of our reference frame are aligned with the GD-1 stream direction. So we will be switching between the ICRS frame and the GD-1 frame often. To the GD-1 frame to build our filter then to the ICRS frame to actually query the database, then back to the GD-1 frame to visualize our results.

- Emphasize the basic structure of the SkyCoord object: location along axis 1, location along axis 2, frame. For example in the ICRS frame, axis 1 is `ra` and axis 2 is `dec`. In the GD-1 frame, axis 1 is `phi1` and axis 2 is `phi2`.

- Mention that we can transform between these different coordinate systems so easily because we are using a unified framework: astropy. This allows learners to build a custom frame object and then have access to all of the other astropy tools that deal with coordinate transformations.

- When you define phi1\_min, phi1\_max, etc go back to Figure 1 and show learners the region you are defining, connecting the min and max values to the coordinates in the GD-1 frame. This is another place you can mention the benefit of using the GD-1 frame is that we can define a rectangle around the stream.

### Episode 3: Tabular Data and Transformations

- Transforming back: the discussion about the motivation to transform back to GD-1 frame is a great time to have student build up some intuition about the proper motion selection and how this relates to the physical picture. Remind them again that GD-1 is a globular cluster that is being pulled into a stream along the phi1 direction. This means that GD-1 stars should have non-zero motion along the phi1 direction and motion close to 0 in the phi2 direction.

- Do not let learners get too hung up on setting the distance and radial velocity to constants. For distance, the parallax is so close to 0 that they are unreliable and as mentioned, the radial velocity is to avoid errors.

- The plot at the end of the reflex correction is a good example of an intermediate step to determining the best filter for GD-1 stars. You can imagine doing this data exploration yourself and trying first a purely spatial plot, then realizing that there are too many non-GD-1 stars included and that you need another way to filter out foreground stars. Proper motion!

- In this episode we are going to define a few functions that we will use in other episodes. This process of prototyping inline and then writing a function is very common. Especially with `make_dataframe` this is a good opportunity to show how the different steps you have performed over the last episode can be elegantly strung together. While defining these functions now takes a little longer, its good to do it when the code is still fresh in learn's minds and will be really useful later.

### Episode 4: Proper motion

- This is likely the final episode of day 1 in a two day workshop.

- This episode takes a slight detour to introduce a few features of Pandas DataFrames. To keep this connected to the story, you can talk about how data exploration is an important part of prototyping your query and making sure you are getting the results you expect in the format you expect them in.

- Starting with selecting the centerline, we do a series of filters on different data frames. Take a minute before you teach this section to make sure you understand what each one represents. We use `results_df` to build `centerline_df`. We use `centerline_df` to determine proper motion limits. We use the proper motion limits determined from `centerline_df` to select GD-1 stars from `results_df`. This is `selected_df`.

- Learners may ask why we did not select this close to the center line from the beginning. It is also worth noting that we are selecting stars close to the center line to identify the proper motion cut, so it is ok if we exclude some GD-1 stars here. In this case we prefer a more pure sample to a more complete sample. Once we have the proper motion limits from our pure sample (`centerline_df`), we will include the full spatial region (`results_df`) and get all of the GD-1 stars (`selected_df`).

- This episode uses the following idiom several times

```
x = ...
y = ...
plt.plot(x, y)
```

This idiom violates the recommendation not to repeat variables names, but since they are defined and used immediately, it should be ok. This syntax simplifies the final `plot` expression making it easier to read.

- Learners may be concerned that we are writing a function for later that does not fill the full canvas (when we set `axis('equal')` in the proper motion selected stars figure). You can reassure them that in episode 7 we will take care of this by learning how to set the figure size, the subplot size, and we will be using a larger spatial area (that we define in episode 4).

- Notice that the previous episode when we used `DataFrame.to_hdf` included the `mode="w"` argument to indicate that we want to create a new, empty HDF Store.  For this (and subsequent) episode(s), we should *not* use the `mode="w"` argument, so that we add new Datasets to the existing Store, rather than starting over.

- At the end of Day 1, if a student is lost or struggled with the end of this episode, point them to their local copy of the HDF5 files so that they can read it in with everyone else on Day 2. They will have downloaded this file as part of the set up in the `student_download/data directory/`.

- If a learner notices that there is a `/` in front of the keys at the end of the episode here is what they mean: `/` indicates that the keys are at the top level of the Dataset hierarchy, and not in a named "group". In future episodes we will add a few more Datasets to this file, but not so many that we need to organize them into groups.

### Episode 5: Coordinate transformation and selection

- This is likely the first episode of Day 2. You should start by having learners start a new notebook. They will need to read in the data they saved yesterday to HDF5 files and the functions they wrote yesterday. You can follow the directions under the collapsable section at the beginning of this episode called "Starting from this episode". If you are using files from the `student_download` directory be sure to add the correct path `student_download/backup-data/<filename>`.

- This episode is faster than it looks because learners have seen much of the content already

- As in previous episodes, a challenge for learners here is keeping track of what we are doing and why. Surface periodically to remind them where we are. The big focus of this episode is that now that we have determined how to filter based on proper motion and eliminated most of the contaminating stars from our prototype, we can repeat the query, this time on the full spatial extent of GD-1 rather than the small portion we used for our prototype.

- The use of ConvexHull in this episode is a bit of a hack, but some learners will find the computation geometry functions in SciPy useful.

- Make sure to review the query `candidate_coord_query_base` from Episode 2 (this will happen naturally when you have to retype it) and remind learners what each filter does and why we are using it

- Notice that the definitions of `phi1_min`, `phi1_max`, etc.  are different in this episode.  Because we are adding more filters, we can select a bigger region without exceeding resource limits.  If learners do not get as many "candidates" as expected, they might be using the old values of these bounds.

- Defining the new region is a good opportunity to go back to the figure and connect these coordinates to the physical picture of GD-1 (e.g. before we were only selecting from this region now we are expanding it to this whole region)

- Learners may ask why we are using `dict(key1=value1, key2=value2)` rather than `{key1:value1, key2:value2}`. We are using the `dict` syntax so that the key value pairs look like keyword arguments. This may simplify the understanding for learners who are less familiar with dictionaries.

- Again, in this episode we should *not* use the `mode="w"` argument, so that we add new Datasets to the existing Store, rather than starting over.

- when writing the `point_series` object to an HDF5 file, learners may see the warning message like the following:
  
  ```/Users/bostroem/opt/anaconda3/envs/AstronomicalData/lib/python3.8/site-packages/pandas/core/generic.py:2434: PerformanceWarning:
  your performance may suffer as PyTables will pickle object types that it cannot
  map directly to c-types [inferred_type->mixed,key->values] [items->None]
  
  pytables.to_hdf(
  ```

This warning is saying that pickle will be slow for large objects. We do not need to worry about this since the
objects we are saving are small.

- If learners struggle with the query exercise - remind them that we are defining a polygon in the same way we did before (with the same syntax) but looking at proper motion instead of space.

- Learners should get the same number of candidates. If they get a different number it is likely they mistyped something.

### Episode 6: Joining tables

- The early part of this episode brings back a lot of best practices that learners learned in previous episodes (e.g. exploring tables, returning the top 5 rows) this is a great opportunity to highlight these and remind learners that they have seen this before and why we are doing it.

- Although we do not immediately use the column definitions in `panstarrs1_best_neighbour`, we will come back to them in the data exploration part of this episode.

- Throughout this episode continue to come back to the central theme of starting with a simple query and building up layers of complexity, testing each layer as you go.

- The Pan-STARRS join exercise is likely to feel scary to a lot of learners - they have only seen one example. Emphasize the building blocks of a `JOIN` e.g. `FROM table1 JOIN table2 ON table1.common_column=table2.common_column`. It is also worth pointing out that common\_column will not have the same names in this case.

- While CSV files are the most basic file format, if you are running short on time, this section can be skipped or abbreviated as learners have often likely encountered CSV files and writing them is not essential to the rest of the curriculum.

- In case a learner asks about the extra `unnamed` column in the CSV section here is the explanation. You may notice that all Pandas `DataFrame`s have an index column which was not part of the original table definition.
  This essentially numbers each row. When we write a `DataFrame` in any other format, the index gets treated like a bonafide column.
  For this reason when we write a CSV file and then read it back into a `DataFrame` the index column gets written as an `unnamed` column and then when it is read back in, another index column is created leading to two extraneous columns.

### Episode 7: Photometry

- It is easy in this episode to lose track of the main point: that we want to define a polygon around the main sequence of GD-1 so we can further hone our sample of candidate GD-1 stars. As we spend time on the isochrone, creating the polygon, etc make sure to come back to this big picture often.

- The key take away from the CMD presentation is that GD-1 is a globular cluster which means all of the stars formed at the same time. Therefore we expect the stars in GD-1 to follow a single, tight isochrone, the main sequence of which we can easily identify.

- We did have to manipulate the MIST isochrone to get it from what is available for download to what we read in. This is too much detail for the curriculum, but the process is detailed in [Making the Isochrone Dataframe](calculating_MIST_isochrone.md) if you or a learner is interested.

- In the original paper, they use an idiosyncratic function to define the boundaries of the isochrone filter. The intention is to define a polygon that gets wider as g increases, to reflect increasing uncertainty. For this exercise we will be using a simplified version which is a constant offset from the isochrone.

- Turning the selection we have defined into a polygon (making a loop) can be challenging to explain verbally but fairly simple visually. This is a nice time to return to the slides or draw a picture

- In the manipulation we do to create the polygon, it can get lost that `left_color`, `right_color`, and `color_loop` are all x values for the polygon (and that `g` and `mag_loop` are the y values). It is worth coming back to this over and over as you introduce each variable.

- Learners might express concern that the polygon we use to select candidate stars is not closed.  That is ok; the `contains_points` method treats the polygon as if it is closed (by connecting the last point to the first).

- Learners may notice that `contains_points` is not actually a `Polygon` property. `Polygon` inherits this method from `Patch`.

- We mention usetex in the discussion of TeX markup in matplotlib. This option does not work for all operating system/backend combinations (in our first alpha pilot some notebooks required restarting after trying to use this option). But, one error that can be fixed is `LaTeX Error: File 'type1cm.sty' not found.` Learners with this error might have to install a package that contains the fonts LaTeX needs.  On some systems, the packages `texlive-latex-extra` or `cm-super` might be what they need.  [See here for more help with
  this](https://stackoverflow.com/questions/11354149/python-unable-to-render-tex-in-matplotlib).

- Learners may ask why we are initializing an empty array and then creating the columns on the fly. DataFrame initializes with arrays of rows rather than columns, so this is the easiest way without having to do some array manipulation. See [https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) for more details.

### Episode 8: Visualization

- This episode is largely about showing learners selected capabilities that will make their lives easier. Matplotlib is a huge package with infinite flexibility - this is in no way complete but hopefully gives them some barebones tools to work with and inspired them to explore further.

- In addition to the links provided in the episode it is nice to show learners the [list of plotting commands](https://matplotlib.org/stable/api/pyplot_summary.html) and the [examples gallery](https://matplotlib.org/stable/gallery/index.html).

- This episode is also a great place to show them how useful the functions we wrote in previous episodes are.

- If you are running short on time, the discussion of customization can be shortened or eliminated and exercises such as adding annotation could be done together as a class.




