---
title: "Basic queries"
teaching: 3000
exercises: 0
questions:
- "How can we select and download the data we want from the Gaia server?"

objectives:
- "Compose a basic query in ADQL/SQL."
- "Use queries to explore a database and its tables."
- "Use queries to download data."
- "Develop, test, and debug a query incrementally."
keypoints:
- "If you can't download an entire dataset (or it's not practical) use queries to select the data you need."

- "Read the metadata and the documentation to make sure you understand the tables, their columns, and what they mean."

- "Develop queries incrementally: start with something simple, test it, and add a little bit at a time."

- "Use ADQL features like `TOP` and `COUNT` to test before you run a query that might return a lot of data."

- "If you know your query will return fewer than 3000 rows, you can 
run it synchronously, which might complete faster (but it doesn't seem to make much difference).  If it might return more than 3000 rows, you should run it asynchronously."

- "ADQL and SQL are not case-sensitive, so you don't have to 
capitalize the keywords, but you should."

- "ADQL and SQL don't require you to break a query into multiple 
lines, but you should."

---

{% include links.md %}

# 1. Queries

This is the first in a series of lessons about working with astronomical data.

As a running example, we will replicate parts of the analysis in a
recent paper, "[Off the beaten path: Gaia reveals GD-1 stars outside
of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian
Price-Whelan and Ana Bonaca.

## Outline

This lesson demonstrates the steps for selecting and downloading data
from the Gaia Database:

1. First we'll make a connection to the Gaia server,

2. We will explore information about the database and the tables it contains,

3. We will write a query and send it to the server, and finally

4. We will download the response from the server.


## Query Language

In order to select data from a database, you have to compose a query,
which is a program written in a "query language".
The query language we'll use is ADQL, which stands for "Astronomical
Data Query Language".

ADQL is a dialect of [SQL](https://en.wikipedia.org/wiki/SQL)
(Structured Query Language), which is by far the most commonly used
query language.  Almost everything you will learn about ADQL also
works in SQL.

[The reference manual for ADQL is
here](http://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html).
But you might find it easier to learn from [this ADQL
Cookbook](https://www.gaia.ac.uk/data/gaia-data-release-1/adql-cookbook).

## Using Jupyter

If you have not worked with Jupyter notebooks before, you might start
with [the tutorial on from Jupyter.org called "Try Classic
Notebook"](https://jupyter.org/try), or [this tutorial from
DataQuest](https://www.dataquest.io/blog/jupyter-notebook-tutorial/).

There are two environments you can use to write and run notebooks: 

* "Jupyter Notebook" is the original, and

* "Jupyter Lab" is a newer environment with more features.

For these lessons, you can use either one.

If you are too impatient for the tutorials, here are the most
important things to know:

1. Notebooks are made up of code cells and text cells (and a few other
less common kinds).  Code cells contain code; text cells, like this
one, contain explanatory text written in
[Markdown](https://www.markdownguide.org/).

2. To run a code cell, click the cell to select it and press
Shift-Enter.  The output of the code should appear below the cell.

3. In general, notebooks only run correctly if you run every code cell
in order from top to bottom.  If you run cells out of order, you are
likely to get errors.

4. You can modify existing cells, but then you have to run them again
to see the effect.

5. You can add new cells, but again, you might have to be careful
about the order you run them in.

6. If you have added or modified cells and the behavior of the
notebook seems strange, you can restart the "kernel", which clears all
of the variables and functions you have defined, and run the cells
again from the beginning.

* If you are using Jupyter notebook, open the `Kernel` menu and select
"Restart and Run All".

* In Jupyter Lab, open the `Kernel` menu and select "Restart Kernel
and Run All Cells"

* In Colab, open the `Runtime` menu and select "Restart and run all"

Before you go on, you might want to explore the other menus and the
toolbar to see what else you can do.

## Connecting to Gaia

The library we'll use to get Gaia data is
[Astroquery](https://astroquery.readthedocs.io/en/latest/).
Astroquery provides `Gaia`, which is an [object that represents a
connection to the Gaia
database](https://astroquery.readthedocs.io/en/latest/gaia/gaia.html).

We can connect to the Gaia database like this:



~~~
from astroquery.gaia import Gaia
~~~
{: .language-python}

~~~
Created TAP+ (v1.2.1) - Connection:
	Host: gea.esac.esa.int
	Use HTTPS: True
	Port: 443
	SSL Port: 443
Created TAP+ (v1.2.1) - Connection:
	Host: geadata.esac.esa.int
	Use HTTPS: True
	Port: 443
	SSL Port: 443

~~~
{: .output}


    

Running this import statement has the effect of creating a
[TAP+](http://www.ivoa.net/documents/TAP/) connection; TAP stands for
"Table Access Protocol", which is a network protocol for sending
queries to the database and getting back the results.

## Databases and Tables

What is a database, anyway?  Most generally, it can be any collection
of data, but when we are talking about ADQL or SQL:

* A database is a collection of one or more named tables.

* Each table is a 2-D array with one or more named columns of data.

We can use `Gaia.load_tables` to get the names of the tables in the
Gaia database.  With the option `only_names=True`, it loads
information about the tables, called "metadata", not the data itself.



~~~
tables = Gaia.load_tables(only_names=True)
~~~
{: .language-python}

~~~
INFO: Retrieving tables... [astroquery.utils.tap.core]
INFO: Parsing tables... [astroquery.utils.tap.core]
INFO: Done. [astroquery.utils.tap.core]

~~~
{: .output}


    



~~~
for table in tables:
    print(table.name)
~~~
{: .language-python}

~~~
external.apassdr9
external.gaiadr2_geometric_distance
external.galex_ais
external.ravedr5_com
external.ravedr5_dr5
external.ravedr5_gra
external.ravedr5_on
external.sdssdr13_photoprimary
external.skymapperdr1_master
external.skymapperdr2_master
external.tmass_xsc
[Output truncated]
~~~
{: .output}


    

So that's a lot of tables.  The ones we'll use are:

* `gaiadr2.gaia_source`, which contains Gaia data from [data release
2](https://www.cosmos.esa.int/web/gaia/data-release-2),

* `gaiadr2.panstarrs1_original_valid`, which contains the photometry
data we'll use from PanSTARRS, and

* `gaiadr2.panstarrs1_best_neighbour`, which we'll use to cross-match
each star observed by Gaia with the same star observed by PanSTARRS.

We can use `load_table` (not `load_tables`) to get the metadata for a
single table.  The name of this function is misleading, because it
only downloads metadata, not the contents of the table.



~~~
meta = Gaia.load_table('gaiadr2.gaia_source')
meta
~~~
{: .language-python}

~~~
Retrieving table 'gaiadr2.gaia_source'
Parsing table 'gaiadr2.gaia_source'...
Done.

<astroquery.utils.tap.model.taptable.TapTableMeta at 0x7f208b57a880>
~~~
{: .output}


    




    



Jupyter shows that the result is an object of type `TapTableMeta`, but
it does not display the contents.

To see the metadata, we have to print the object.



~~~
print(meta)
~~~
{: .language-python}

~~~
TAP Table name: gaiadr2.gaiadr2.gaia_source
Description: This table has an entry for every Gaia observed source as listed in the
Main Database accumulating catalogue version from which the catalogue
release has been generated. It contains the basic source parameters,
that is only final data (no epoch data) and no spectra (neither final
nor epoch).
Num. columns: 96

~~~
{: .output}


    

## Columns

The following loop prints the names of the columns in the table.



~~~
for column in meta.columns:
    print(column.name)
~~~
{: .language-python}

~~~
solution_id
designation
source_id
random_index
ref_epoch
ra
ra_error
dec
dec_error
parallax
parallax_error
[Output truncated]
~~~
{: .output}


    

You can probably guess what many of these columns are by looking at
the names, but you should resist the temptation to guess.
To find out what the columns mean, [read the
documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).

If you want to know what can go wrong when you don't read the
documentation, [you might like this
article](https://www.vox.com/future-perfect/2019/6/4/18650969/married-women-miserable-fake-paul-dolan-happiness).

> ## Exercise
> 
> One of the other tables we'll use is
> `gaiadr2.panstarrs1_original_valid`.  Use `load_table` to get the
> metadata for this table.  How many columns are there and what are
> their names?


>
> > ## Solution
> > 
> > ~~~
> > 
> > meta2 = Gaia.load_table('gaiadr2.panstarrs1_original_valid')
> > print(meta2)
> > 
> > for column in meta2.columns:
> >     print(column.name)
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}



    Retrieving table 'gaiadr2.panstarrs1_original_valid'
    Parsing table 'gaiadr2.panstarrs1_original_valid'...
    Done.
    TAP Table name: gaiadr2.gaiadr2.panstarrs1_original_valid
    Description: The Panoramic Survey Telescope and Rapid Response System (Pan-STARRS) is
    a system for wide-field astronomical imaging developed and operated by
    the Institute for Astronomy at the University of Hawaii. Pan-STARRS1
    (PS1) is the first part of Pan-STARRS to be completed and is the basis
    for Data Release 1 (DR1). The PS1 survey used a 1.8 meter telescope and
    its 1.4 Gigapixel camera to image the sky in five broadband filters (g,
    r, i, z, y).
    
    The current table contains a filtered subsample of the 10 723 304 629
    entries listed in the original ObjectThin table.
    We used only ObjectThin and MeanObject tables to extract
    panstarrs1OriginalValid table, this means that objects detected only in
    stack images are not included here. The main reason for us to avoid the
    use of objects detected in stack images is that their astrometry is not
    as good as the mean objects astrometry: “The stack positions (raStack,
    decStack) have considerably larger systematic astrometric errors than
    the mean epoch positions (raMean, decMean).” The astrometry for the
    MeanObject positions uses Gaia DR1 as a reference catalog, while the
    stack positions use 2MASS as a reference catalog.
    
    In details, we filtered out all objects where:
    
    -   nDetections = 1
    
    -   no good quality data in Pan-STARRS, objInfoFlag 33554432 not set
    
    -   mean astrometry could not be measured, objInfoFlag 524288 set
    
    -   stack position used for mean astrometry, objInfoFlag 1048576 set
    
    -   error on all magnitudes equal to 0 or to -999;
    
    -   all magnitudes set to -999;
    
    -   error on RA or DEC greater than 1 arcsec.
    
    The number of objects in panstarrs1OriginalValid is 2 264 263 282.
    
    The panstarrs1OriginalValid table contains only a subset of the columns
    available in the combined ObjectThin and MeanObject tables. A
    description of the original ObjectThin and MeanObjects tables can be
    found at:
    https://outerspace.stsci.edu/display/PANSTARRS/PS1+Database+object+and+detection+tables
    
    Download:
    http://mastweb.stsci.edu/ps1casjobs/home.aspx
    Documentation:
    https://outerspace.stsci.edu/display/PANSTARRS
    http://pswww.ifa.hawaii.edu/pswww/
    References:
    The Pan-STARRS1 Surveys, Chambers, K.C., et al. 2016, arXiv:1612.05560
    Pan-STARRS Data Processing System, Magnier, E. A., et al. 2016,
    arXiv:1612.05240
    Pan-STARRS Pixel Processing: Detrending, Warping, Stacking, Waters, C.
    Z., et al. 2016, arXiv:1612.05245
    Pan-STARRS Pixel Analysis: Source Detection and Characterization,
    Magnier, E. A., et al. 2016, arXiv:1612.05244
    Pan-STARRS Photometric and Astrometric Calibration, Magnier, E. A., et
    al. 2016, arXiv:1612.05242
    The Pan-STARRS1 Database and Data Products, Flewelling, H. A., et al.
    2016, arXiv:1612.05243
    
    Catalogue curator:
    SSDC - ASI Space Science Data Center
    https://www.ssdc.asi.it/
    Num. columns: 26
    obj_name
    obj_id
    ra
    dec
    ra_error
    dec_error
    epoch_mean
    g_mean_psf_mag
    g_mean_psf_mag_error
    g_flags
    r_mean_psf_mag
    r_mean_psf_mag_error
    r_flags
    i_mean_psf_mag
    i_mean_psf_mag_error
    i_flags
    z_mean_psf_mag
    z_mean_psf_mag_error
    z_flags
    y_mean_psf_mag
    y_mean_psf_mag_error
    y_flags
    n_detections
    zone_id
    obj_info_flag
    quality_flag


## Writing queries

By now you might be wondering how we download the actual data.  With
tables this big, you generally don't.  Instead, you use queries to
select only the data you want.

A query is a string written in a query language like SQL; for the Gaia
database, the query language is a dialect of SQL called ADQL.

Here's an example of an ADQL query.



~~~
query1 = """SELECT 
TOP 10
source_id, ra, dec, parallax 
FROM gaiadr2.gaia_source
"""
~~~
{: .language-python}

**Python note:** We use a [triple-quoted
string](https://docs.python.org/3/tutorial/introduction.html#strings)
here so we can include line breaks in the query, which makes it easier
to read.

The words in uppercase are ADQL keywords:

* `SELECT` indicates that we are selecting data (as opposed to adding
or modifying data).

* `TOP` indicates that we only want the first 10 rows of the table,
which is useful for testing a query before asking for all of the data.

* `FROM` specifies which table we want data from.

The third line is a list of column names, indicating which columns we want.  

In this example, the keywords are capitalized and the column names are
lowercase.  This is a common style, but it is not required.  ADQL and
SQL are not case-sensitive.

Also, the query is broken into multiple lines to make it more
readable.  This is a common style, but not required.  Line breaks
don't affect the behavior of the query.

To run this query, we use the `Gaia` object, which represents our
connection to the Gaia database, and invoke `launch_job`:



~~~
job = Gaia.launch_job(query1)
job
~~~
{: .language-python}

~~~
<astroquery.utils.tap.model.job.Job at 0x7f208b5a49a0>
~~~
{: .output}





    



The result is an object that represents the job running on a Gaia server.

If you print it, it displays metadata for the forthcoming table.



~~~
print(job)
~~~
{: .language-python}

~~~
<Table length=10>
   name    dtype  unit                            description                             n_bad
--------- ------- ---- ------------------------------------------------------------------ -----
source_id   int64      Unique source identifier (unique within a particular Data Release)     0
       ra float64  deg                                                    Right ascension     0
      dec float64  deg                                                        Declination     0
 parallax float64  mas                                                           Parallax     3
Jobid: None
Phase: COMPLETED
Owner: None
Output file: sync_20210113094733.xml.gz
[Output truncated]
~~~
{: .output}


    

Don't worry about `Results: None`.  That does not actually mean there
are no results.

However, `Phase: COMPLETED` indicates that the job is complete, so we
can get the results like this:



~~~
results = job.get_results()
type(results)
~~~
{: .language-python}

~~~
astropy.table.table.Table
~~~
{: .output}





    



**Optional detail:**  Why is `table` repeated three times?  The first
is the name of the module, the second is the name of the submodule,
and the third is the name of the class.  Most of the time we only care
about the last one.  It's like the Linnean name for gorilla, which is
*Gorilla gorilla gorilla*.

The result is an [Astropy
Table](https://docs.astropy.org/en/stable/table/), which is similar to
a table in an SQL database except:

* SQL databases are stored on disk drives, so they are persistent;
that is, they "survive" even if you turn off the computer.  An Astropy
`Table` is stored in memory; it disappears when you turn off the
computer (or shut down this Jupyter notebook).

* SQL databases are designed to process queries.  An Astropy `Table`
can perform some query-like operations, like selecting columns and
rows.  But these operations use Python syntax, not SQL.

Jupyter knows how to display the contents of a `Table`.



~~~
results1
~~~
{: .language-python}

~~~
<Table length=10>
     source_id              ra         ...        parallax       
                           deg         ...          mas          
       int64             float64       ...        float64        
------------------- ------------------ ... ----------------------
4036442223689093632  268.1632016769232 ...     1.1699115050505098
4036432465518023424  268.4099737354832 ...    0.09728878869604544
4036431056767171072  268.4987929943058 ...                     --
4036427586432136576  268.1550746929837 ...    -1.6517420537894711
4036422295031001856   268.295159195894 ...                     --
4036435454861834880   268.402863428808 ...    0.47166978762943135
[Output truncated]
~~~
{: .output}





<i>Table length=10</i>
<table id="table139777753632192" class="table-striped table-bordered table-condensed">
<thead><tr><th>source_id</th><th>ra</th><th>dec</th><th>parallax</th></tr></thead>
<thead><tr><th></th><th>deg</th><th>deg</th><th>mas</th></tr></thead>
<thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>4036442223689093632</td><td>268.1632016769232</td><td>-39.45184716473943</td><td>1.1699115050505098</td></tr>
<tr><td>4036432465518023424</td><td>268.4099737354832</td><td>-39.48255132974542</td><td>0.09728878869604544</td></tr>
<tr><td>4036431056767171072</td><td>268.4987929943058</td><td>-39.546378864610695</td><td>--</td></tr>
<tr><td>4036427586432136576</td><td>268.1550746929837</td><td>-39.640672063878434</td><td>-1.6517420537894711</td></tr>
<tr><td>4036422295031001856</td><td>268.295159195894</td><td>-39.75110784405333</td><td>--</td></tr>
<tr><td>4036435454861834880</td><td>268.402863428808</td><td>-39.4490745170462</td><td>0.47166978762943135</td></tr>
<tr><td>4036423497622557184</td><td>268.22982540766145</td><td>-39.695491626875175</td><td>-0.0018594810744405113</td></tr>
<tr><td>4036432706082489344</td><td>268.53765532143115</td><td>-39.49836075294007</td><td>0.294780085311733</td></tr>
<tr><td>4036431984481125248</td><td>268.4176240506027</td><td>-39.52535476978414</td><td>-0.23707477297436713</td></tr>
<tr><td>4036439784146944768</td><td>268.118736733789</td><td>-39.55860223906007</td><td>--</td></tr>
</table>



Each column has a name, units, and a data type.

For example, the units of `ra` and `dec` are degrees, and their data
type is `float64`, which is a 64-bit floating-point number, used to
store measurements with a fraction part.

This information comes from the Gaia database, and has been stored in
the Astropy `Table` by Astroquery.

> ## Exercise
> 
> Read [the documentation of this
> table](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html)
> and choose a column that looks interesting to you.  Add the column
> name to the query and run it again.  What are the units of the column
> you selected?  What is its data type?


>
> > ## Solution
> > 
> > ~~~
> > 
> > # Let's add
> > #
> > # radial_velocity : Radial velocity (double, Velocity[km/s] )
> > #
> > # Spectroscopic radial velocity in the solar barycentric 
> > # reference frame.
> > #
> > # The radial velocity provided is the median value of the 
> > # radial velocity measurements at all epochs.
> > 
> > query = """SELECT 
> > TOP 10
> > source_id, ra, dec, parallax, radial_velocity
> > FROM gaiadr2.gaia_source
> > """
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}



## Asynchronous queries

`launch_job` asks the server to run the job "synchronously", which
normally means it runs immediately.  But synchronous jobs are limited
to 2000 rows.  For queries that return more rows, you should run
"asynchronously", which mean they might take longer to get started.

If you are not sure how many rows a query will return, you can use the
SQL command `COUNT` to find out how many rows are in the result
without actually returning them.  We'll see an example in the next
lesson.

The results of an asynchronous query are stored in a file on the
server, so you can start a query and come back later to get the
results.
For anonymous users, files are kept for three days.

As an example, let's try a query that's similar to `query1`, with these changes:

* It selects the first 3000 rows, so it is bigger than we should run
synchronously.

* It selects two additional columns, `pmra` and `pmdec`, which are
proper motions along the axes of `ra` and `dec`.

* It uses a new keyword, `WHERE`.



~~~
query2 = """SELECT 
TOP 3000
source_id, ra, dec, pmra, pmdec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1
"""
~~~
{: .language-python}

A `WHERE` clause indicates which rows we want; in this case, the query
selects only rows "where" `parallax` is less than 1.  This has the
effect of selecting stars with relatively low parallax, which are
farther away.
We use this clause to exclude nearby stars that are unlikely to be part of GD-1.

`WHERE` is one of the most common clauses in ADQL/SQL, and one of the
most useful, because it allows us to download only the rows we need
from the database.

We use `launch_job_async` to submit an asynchronous query.



~~~
job = Gaia.launch_job_async(query)
job
~~~
{: .language-python}

~~~
INFO: Query finished. [astroquery.utils.tap.core]

<astroquery.utils.tap.model.job.Job at 0x7f208c75e7c0>
~~~
{: .output}


    




    



And here are the results.



~~~
results = job.get_results()
results
~~~
{: .language-python}

~~~
<Table length=3000>
     source_id              ra         ...        parallax       
                           deg         ...          mas          
       int64             float64       ...        float64        
------------------- ------------------ ... ----------------------
4036432465518023424  268.4099737354832 ...    0.09728878869604544
4036427586432136576  268.1550746929837 ...    -1.6517420537894711
4036435454861834880   268.402863428808 ...    0.47166978762943135
4036423497622557184 268.22982540766145 ... -0.0018594810744405113
4036432706082489344 268.53765532143115 ...      0.294780085311733
4036431984481125248  268.4176240506027 ...   -0.23707477297436713
[Output truncated]
~~~
{: .output}





<i>Table length=3000</i>
<table id="table139777750656528" class="table-striped table-bordered table-condensed">
<thead><tr><th>source_id</th><th>ra</th><th>dec</th><th>pmra</th><th>pmdec</th><th>parallax</th></tr></thead>
<thead><tr><th></th><th>deg</th><th>deg</th><th>mas / yr</th><th>mas / yr</th><th>mas</th></tr></thead>
<thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>4036432465518023424</td><td>268.4099737354832</td><td>-39.48255132974542</td><td>-6.18591186141027</td><td>-9.901482693050884</td><td>0.09728878869604544</td></tr>
<tr><td>4036427586432136576</td><td>268.1550746929837</td><td>-39.640672063878434</td><td>2.902531307436574</td><td>-2.6372995461983555</td><td>-1.6517420537894711</td></tr>
<tr><td>4036435454861834880</td><td>268.402863428808</td><td>-39.4490745170462</td><td>-2.9204452670058174</td><td>-0.6600330794265028</td><td>0.47166978762943135</td></tr>
<tr><td>4036423497622557184</td><td>268.22982540766145</td><td>-39.695491626875175</td><td>-3.3502947606249847</td><td>-3.9994189219114302</td><td>-0.0018594810744405113</td></tr>
<tr><td>4036432706082489344</td><td>268.53765532143115</td><td>-39.49836075294007</td><td>2.2117535185619137</td><td>0.26989064159601095</td><td>0.294780085311733</td></tr>
<tr><td>4036431984481125248</td><td>268.4176240506027</td><td>-39.52535476978414</td><td>-4.657776446841011</td><td>-1.5298512820469425</td><td>-0.23707477297436713</td></tr>
<tr><td>4036433942986772864</td><td>268.4446322062673</td><td>-39.45434881922685</td><td>-1.268163073772277</td><td>-3.143326774839096</td><td>0.4383827477396896</td></tr>
<tr><td>4036434767622765568</td><td>268.3397669855326</td><td>-39.533206030665674</td><td>-2.570383074581705</td><td>-7.8309562837625695</td><td>0.9802160615520338</td></tr>
<tr><td>4036426658716007680</td><td>268.2325720330626</td><td>-39.636792360242836</td><td>-3.483222377252905</td><td>-11.518985840909323</td><td>0.3402866531323366</td></tr>
<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
<tr><td>4036429643676361088</td><td>268.1940481392251</td><td>-39.53484810691843</td><td>-0.5399180238167743</td><td>-0.43950577474250674</td><td>0.09673690064167999</td></tr>
<tr><td>4036434389663288448</td><td>268.51391606591034</td><td>-39.40538040808745</td><td>-2.6667782256956745</td><td>-3.3321195693502235</td><td>0.614141340458613</td></tr>
<tr><td>4036423291461004800</td><td>268.2446471322681</td><td>-39.698336039409675</td><td>-3.525450125069836</td><td>-4.404888528602299</td><td>0.6173594385857027</td></tr>
<tr><td>4036440196464371328</td><td>268.1093614056419</td><td>-39.511003642540814</td><td>-4.149726381202896</td><td>-2.4847916296292567</td><td>0.477742037030516</td></tr>
<tr><td>4036428720300382720</td><td>268.29821418217523</td><td>-39.57587285989207</td><td>-5.998673555253831</td><td>-9.23082959079349</td><td>-0.47420031990483885</td></tr>
<tr><td>4036442601646439168</td><td>268.14514023769215</td><td>-39.44636747310117</td><td>-0.6481888801608459</td><td>-4.502398437965107</td><td>0.4456397029007401</td></tr>
<tr><td>4036433393230070016</td><td>268.54224417983687</td><td>-39.45140942552161</td><td>1.3749953846556249</td><td>-9.205922550553257</td><td>0.22018046862051377</td></tr>
<tr><td>4036429407495950464</td><td>268.22710083295215</td><td>-39.54661989184268</td><td>-3.802526532062553</td><td>-3.227275729455894</td><td>-0.14612958706070997</td></tr>
<tr><td>4036439200034820608</td><td>268.05859875239844</td><td>-39.61886564077019</td><td>-2.1172181651309834</td><td>-3.7754096391993404</td><td>-0.5167970688281044</td></tr>
<tr><td>4036441772690773632</td><td>268.16153516094687</td><td>-39.48158028904507</td><td>-4.015920758803185</td><td>-3.4190173584838024</td><td>0.45270857335252485</td></tr>
</table>



You might notice that some values of `parallax` are negative.  As
[this FAQ
explains](https://www.cosmos.esa.int/web/gaia/archive-tips#negative%20parallax),
"Negative parallaxes are caused by errors in the observations."
Negative parallaxes have "no physical meaning," but they can be a
"useful diagnostic on the quality of the astrometric solution."

> ## Exercise
> 
> The clauses in a query have to be in the right order.  Go back and
> change the order of the clauses in `query2` and run it again.
> 
> The modified query should fail, but notice that you don't get much
> useful debugging information.
> 
> For this reason, developing and debugging ADQL queries can be really
> hard.  A few suggestions that might help:
> 
> * Whenever possible, start with a working query, either an example you
> find online or a query you have used in the past.
> 
> * Make small changes and test each change before you continue.
> 
> * While you are debugging, use `TOP` to limit the number of rows in
> the result.  That will make each test run faster, which reduces your
> development time.
> 
> * Launching test queries synchronously might make them start faster, too.


>
> > ## Solution
> > 
> > ~~~
> > 
> > # In this example, the WHERE clause is in the wrong place
> > 
> > query = """SELECT 
> > TOP 3000
> > WHERE parallax < 1
> > source_id, ref_epoch, ra, dec, parallax
> > FROM gaiadr2.gaia_source
> > """
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}



## Operators

In a `WHERE` clause, you can use any of the [SQL comparison
operators](https://www.w3schools.com/sql/sql_operators.asp); here are
the most common ones:

| Symbol | Operation
|--------| :---
| `>` | greater than
| `<` | less than
| `>=` | greater than or equal
| `<=` | less than or equal
| `=` | equal
| `!=` or `<>` | not equal

Most of these are the same as Python, but some are not.  In
particular, notice that the equality operator is `=`, not `==`.
Be careful to keep your Python out of your ADQL!

You can combine comparisons using the logical operators:

* AND: true if both comparisons are true
* OR: true if either or both comparisons are true

Finally, you can use `NOT` to invert the result of a comparison. 

> ## Exercise
> 
> [Read about SQL operators
> here](https://www.w3schools.com/sql/sql_operators.asp) and then modify
> the previous query to select rows where `bp_rp` is between `-0.75` and
> `2`.


>
> > ## Solution
> > 
> > ~~~
> > 
> > # Here's a solution using > and < operators
> > 
> > query = """SELECT 
> > TOP 10
> > source_id, ref_epoch, ra, dec, parallax
> > FROM gaiadr2.gaia_source
> > WHERE parallax < 1 
> >   AND bp_rp > -0.75 AND bp_rp < 2
> > """
> > 
> > # And here's a solution using the BETWEEN operator
> > 
> > query = """SELECT 
> > TOP 10
> > source_id, ref_epoch, ra, dec, parallax
> > FROM gaiadr2.gaia_source
> > WHERE parallax < 1 
> >   AND bp_rp BETWEEN -0.75 AND 2
> > """
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}



`bp_rp` contains BP-RP color, which is the difference between two
other columns, `phot_bp_mean_mag` and `phot_rp_mean_mag`.
You can [read about this variable
here](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).

This [Hertzsprung-Russell
diagram](https://sci.esa.int/web/gaia/-/60198-gaia-hertzsprung-russell-diagram)
shows the BP-RP color and luminosity of stars in the Gaia catalog
(Copyright: ESA/Gaia/DPAC, CC BY-SA 3.0 IGO).

<img width="300"
src="https://github.com/AllenDowney/AstronomicalData/raw/main/images/1567214809100-ESA_Gaia_DR2_HRD_Gaia_625.jpg">

Selecting stars with `bp-rp` less than 2 excludes many [class M dwarf
stars](https://xkcd.com/2360/), which are low temperature, low
luminosity.  A star like that at GD-1's distance would be hard to
detect, so if it is detected, it it more likely to be in the
foreground.

## Formatting queries

The queries we have written so far are string "literals", meaning that
the entire string is part of the program.
But writing queries yourself can be slow, repetitive, and error-prone.

It is often a good idea to write Python code that assembles a query
for you.  One useful tool for that is the [string `format`
method](https://www.w3schools.com/python/ref_string_format.asp).

As an example, we'll divide the previous query into two parts; a list
of column names and a "base" for the query that contains everything
except the column names.

Here's the list of columns we'll select.  



~~~
columns = 'source_id, ra, dec, pmra, pmdec, parallax'
~~~
{: .language-python}

And here's the base; it's a string that contains at least one format
specifier in curly brackets (braces).



~~~
query3_base = """SELECT 
TOP 10 
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2
"""
~~~
{: .language-python}

This base query contains one format specifier, `{columns}`, which is a
placeholder for the list of column names we will provide.

To assemble the query, we invoke `format` on the base string and
provide a keyword argument that assigns a value to `columns`.



~~~
query3 = query3_base.format(columns=columns)
~~~
{: .language-python}

The result is a string with line breaks.  If you display it, the line
breaks appear as `\n`.



~~~
query3
~~~
{: .language-python}

~~~
'SELECT \nTOP 10 \nsource_id, ra, dec, pmra, pmdec, parallax\nFROM gaiadr2.gaia_source\nWHERE parallax < 1\n  AND bp_rp BETWEEN -0.75 AND 2\n'
~~~
{: .output}





    



But if you print it, the line breaks appear as... line breaks.



~~~
print(query3)
~~~
{: .language-python}

~~~
SELECT 
TOP 10 
source_id, ra, dec, pmra, pmdec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2


~~~
{: .output}


    

Notice that the format specifier has been replaced with the value of `columns`.

Let's run it and see if it works:



~~~
job = Gaia.launch_job(query3)
print(job)
~~~
{: .language-python}

~~~
<Table length=10>
   name    dtype    unit                              description                            
--------- ------- -------- ------------------------------------------------------------------
source_id   int64          Unique source identifier (unique within a particular Data Release)
       ra float64      deg                                                    Right ascension
      dec float64      deg                                                        Declination
     pmra float64 mas / yr                         Proper motion in right ascension direction
    pmdec float64 mas / yr                             Proper motion in declination direction
 parallax float64      mas                                                           Parallax
Jobid: None
Phase: COMPLETED
[Output truncated]
~~~
{: .output}


    



~~~
results = job.get_results()
results
~~~
{: .language-python}

~~~
<Table length=10>
     source_id              ra         ...       parallax      
                           deg         ...         mas         
       int64             float64       ...       float64       
------------------- ------------------ ... --------------------
4036432495535155072 268.41729223407117 ... -0.11176468032974321
4036441536493921920   268.163280168269 ... -0.06466314453980221
4036442601646438400  268.1387976914471 ...  -1.1448803059362858
4036429201337340672  268.2935748074385 ...   0.4772565168467984
4036435209981391360  268.3235008543401 ...  -0.2260844472181426
4036430781889616640  268.4323413972199 ... 0.021914470221676332
[Output truncated]
~~~
{: .output}





<i>Table length=10</i>
<table id="table139777753221728" class="table-striped table-bordered table-condensed">
<thead><tr><th>source_id</th><th>ra</th><th>dec</th><th>pmra</th><th>pmdec</th><th>parallax</th></tr></thead>
<thead><tr><th></th><th>deg</th><th>deg</th><th>mas / yr</th><th>mas / yr</th><th>mas</th></tr></thead>
<thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>4036432495535155072</td><td>268.41729223407117</td><td>-39.46918929346299</td><td>-1.4657595289399132</td><td>-5.748127423097124</td><td>-0.11176468032974321</td></tr>
<tr><td>4036441536493921920</td><td>268.163280168269</td><td>-39.49035741409665</td><td>-2.133016188076998</td><td>-0.117511099698601</td><td>-0.06466314453980221</td></tr>
<tr><td>4036442601646438400</td><td>268.1387976914471</td><td>-39.451546875392324</td><td>-0.9017024454865927</td><td>-6.79073640236943</td><td>-1.1448803059362858</td></tr>
<tr><td>4036429201337340672</td><td>268.2935748074385</td><td>-39.511082764583186</td><td>0.23803118725105005</td><td>-2.619620943565854</td><td>0.4772565168467984</td></tr>
<tr><td>4036435209981391360</td><td>268.3235008543401</td><td>-39.47076322705205</td><td>-5.269315172794177</td><td>-6.311972787283363</td><td>-0.2260844472181426</td></tr>
<tr><td>4036430781889616640</td><td>268.4323413972199</td><td>-39.568388828520796</td><td>0.2101300515915823</td><td>0.5020913360695279</td><td>0.021914470221676332</td></tr>
<tr><td>4036434595821886080</td><td>268.51291280069273</td><td>-39.39090803101428</td><td>-2.0281461390898623</td><td>-1.09405424552038</td><td>0.42083898443584056</td></tr>
<tr><td>4036426933645905664</td><td>268.2033082657327</td><td>-39.62426797607084</td><td>0.40080910282280036</td><td>-4.320453719127832</td><td>0.6494251580592608</td></tr>
<tr><td>4036425108189475968</td><td>268.26633179055267</td><td>-39.648450617498646</td><td>-2.0326696514976903</td><td>-5.868948034915663</td><td>0.23435564036004816</td></tr>
<tr><td>4036440368305153280</td><td>268.0417228205479</td><td>-39.55154954991695</td><td>-7.807636960289898</td><td>-7.1469993883954475</td><td>-0.06901432561194154</td></tr>
</table>



Good so far.

> ## Exercise
> 
> This query always selects sources with `parallax` less than 1.  But
> suppose you want to take that upper bound as an input.
> 
> Modify `query3_base` to replace `1` with a format specifier like
> `{max_parallax}`.  Now, when you call `format`, add a keyword argument
> that assigns a value to `max_parallax`, and confirm that the format
> specifier gets replaced with the value you provide.


>
> > ## Solution
> > 
> > ~~~
> > 
> > query_base = """SELECT 
> > TOP 10
> > {columns}
> > FROM gaiadr2.gaia_source
> > WHERE parallax < {max_parallax} AND 
> > bp_rp BETWEEN -0.75 AND 2
> > """
> > 
> > query = query_base.format(columns=columns,
> >                             max_parallax=0.5)
> > print(query)
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}



    SELECT 
    TOP 10
    source_id, ra, dec, pmra, pmdec, parallax
    FROM gaiadr2.gaia_source
    WHERE parallax < 0.5 AND 
    bp_rp BETWEEN -0.75 AND 2
    


## Summary

This notebook demonstrates the following steps:

1. Making a connection to the Gaia server,

2. Exploring information about the database and the tables it contains,

3. Writing a query and sending it to the server, and finally

4. Downloading the response from the server as an Astropy `Table`.

In the next lesson we will extend these queries to select a particular
region of the sky.

## Best practices

* If you can't download an entire dataset (or it's not practical) use
queries to select the data you need.

* Read the metadata and the documentation to make sure you understand
the tables, their columns, and what they mean.

* Develop queries incrementally: start with something simple, test it,
and add a little bit at a time.

* Use ADQL features like `TOP` and `COUNT` to test before you run a
query that might return a lot of data.

* If you know your query will return fewer than 3000 rows, you can run
it synchronously, which might complete faster.  If it might return
more than 3000 rows, you should run it asynchronously.

* ADQL and SQL are not case-sensitive, so you don't have to capitalize
the keywords, but you should.

* ADQL and SQL don't require you to break a query into multiple lines,
but you should.

Jupyter notebooks can be good for developing and testing code, but
they have some drawbacks.  In particular, if you run the cells out of
order, you might find that variables don't have the values you expect.

There are a few things you can do to mitigate these problems:

* Make each section of the notebook self-contained.  Try not to use
the same variable name in more than one section.

* Keep notebooks short.  Look for places where you can break your
analysis into phases with one notebook per phase.



~~~

~~~
{: .language-python}
