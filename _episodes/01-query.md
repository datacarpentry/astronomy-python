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
FIXME

{% include links.md %}

# Queries

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
"Table Access Protocol".  It is a network protocol for sending queries
to the database and getting back the results.

It looks like it connects to two servers, `gea.esac.esa.int` and
`geadata.esac.esa.int`; we don't know why, but possibly one of them is
used for metadata and the other for data.

## Databases and Tables

What is a database, anyway?  Most generally, it can be any collection
of data, but when we are talking about ADQL or SQL:

* A database is a collection of one or more named tables.

* Each table is a 2-D array with one or more named columns of data.

We can use `Gaia.load_tables` to get the names of the tables in the
Gaia database.  With the option `only_names=True`, it loads
information about the tables, called the "metadata", not the data
itself.

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
external.tmass_xsc
public.hipparcos
public.hipparcos_newreduction
public.hubble_sc
public.igsl_source
public.igsl_source_catalog_ids
public.tycho2
public.dual
tap_config.coord_sys
tap_config.properties
tap_schema.columns
tap_schema.key_columns
tap_schema.keys
tap_schema.schemas
tap_schema.tables
gaiadr1.aux_qso_icrf2_match
gaiadr1.ext_phot_zero_point
gaiadr1.allwise_best_neighbour
gaiadr1.allwise_neighbourhood
gaiadr1.gsc23_best_neighbour
gaiadr1.gsc23_neighbourhood
gaiadr1.ppmxl_best_neighbour
gaiadr1.ppmxl_neighbourhood
gaiadr1.sdss_dr9_best_neighbour
gaiadr1.sdss_dr9_neighbourhood
gaiadr1.tmass_best_neighbour
gaiadr1.tmass_neighbourhood
gaiadr1.ucac4_best_neighbour
gaiadr1.ucac4_neighbourhood
gaiadr1.urat1_best_neighbour
gaiadr1.urat1_neighbourhood
gaiadr1.cepheid
gaiadr1.phot_variable_time_series_gfov
gaiadr1.phot_variable_time_series_gfov_statistical_parameters
gaiadr1.rrlyrae
gaiadr1.variable_summary
gaiadr1.allwise_original_valid
gaiadr1.gsc23_original_valid
gaiadr1.ppmxl_original_valid
gaiadr1.sdssdr9_original_valid
gaiadr1.tmass_original_valid
gaiadr1.ucac4_original_valid
gaiadr1.urat1_original_valid
gaiadr1.gaia_source
gaiadr1.tgas_source
gaiadr2.aux_allwise_agn_gdr2_cross_id
gaiadr2.aux_iers_gdr2_cross_id
gaiadr2.aux_sso_orbit_residuals
gaiadr2.aux_sso_orbits
gaiadr2.dr1_neighbourhood
gaiadr2.allwise_best_neighbour
gaiadr2.allwise_neighbourhood
gaiadr2.apassdr9_best_neighbour
gaiadr2.apassdr9_neighbourhood
gaiadr2.gsc23_best_neighbour
gaiadr2.gsc23_neighbourhood
gaiadr2.hipparcos2_best_neighbour
gaiadr2.hipparcos2_neighbourhood
gaiadr2.panstarrs1_best_neighbour
gaiadr2.panstarrs1_neighbourhood
gaiadr2.ppmxl_best_neighbour
gaiadr2.ppmxl_neighbourhood
gaiadr2.ravedr5_best_neighbour
gaiadr2.ravedr5_neighbourhood
gaiadr2.sdssdr9_best_neighbour
gaiadr2.sdssdr9_neighbourhood
gaiadr2.tmass_best_neighbour
gaiadr2.tmass_neighbourhood
gaiadr2.tycho2_best_neighbour
gaiadr2.tycho2_neighbourhood
gaiadr2.urat1_best_neighbour
gaiadr2.urat1_neighbourhood
gaiadr2.sso_observation
gaiadr2.sso_source
gaiadr2.vari_cepheid
gaiadr2.vari_classifier_class_definition
gaiadr2.vari_classifier_definition
gaiadr2.vari_classifier_result
gaiadr2.vari_long_period_variable
gaiadr2.vari_rotation_modulation
gaiadr2.vari_rrlyrae
gaiadr2.vari_short_timescale
gaiadr2.vari_time_series_statistics
gaiadr2.panstarrs1_original_valid
gaiadr2.gaia_source
gaiadr2.ruwe

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
only downloads metadata.

~~~
meta = Gaia.load_table('gaiadr2.gaia_source')
meta
~~~
{: .language-python}

~~~
Retrieving table 'gaiadr2.gaia_source'
Parsing table 'gaiadr2.gaia_source'...
Done.

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

In `meta`, the name of the table appears as
`gaiadr2.gaiadr2.gaia_source`, which is the "qualified name", but when
we load the metadata, we refer to it as `gaiadr2.gaia_source`.

~~~
# Solution

# The error message, last time we tried, was

# Retrieving table 'gaiadr2.gaiadr2.gaia_source'
# 500 Error 500:
# esavo.tap.TAPException: esavo.tap.TAPException: Schema cannot be null

# Which is not remotely helpful.

# The point of this exercise is to alert the participants to the difficulty
# of debugging queries with VERY limited feedback.  So developing and testing
# incrementally is very important.
~~~
{: .language-python}
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
parallax_over_error
pmra
pmra_error
pmdec
pmdec_error
ra_dec_corr
ra_parallax_corr
ra_pmra_corr
ra_pmdec_corr
dec_parallax_corr
dec_pmra_corr
dec_pmdec_corr
parallax_pmra_corr
parallax_pmdec_corr
pmra_pmdec_corr
astrometric_n_obs_al
astrometric_n_obs_ac
astrometric_n_good_obs_al
astrometric_n_bad_obs_al
astrometric_gof_al
astrometric_chi2_al
astrometric_excess_noise
astrometric_excess_noise_sig
astrometric_params_solved
astrometric_primary_flag
astrometric_weight_al
astrometric_pseudo_colour
astrometric_pseudo_colour_error
mean_varpi_factor_al
astrometric_matched_observations
visibility_periods_used
astrometric_sigma5d_max
frame_rotator_object_type
matched_observations
duplicated_source
phot_g_n_obs
phot_g_mean_flux
phot_g_mean_flux_error
phot_g_mean_flux_over_error
phot_g_mean_mag
phot_bp_n_obs
phot_bp_mean_flux
phot_bp_mean_flux_error
phot_bp_mean_flux_over_error
phot_bp_mean_mag
phot_rp_n_obs
phot_rp_mean_flux
phot_rp_mean_flux_error
phot_rp_mean_flux_over_error
phot_rp_mean_mag
phot_bp_rp_excess_factor
phot_proc_mode
bp_rp
bp_g
g_rp
radial_velocity
radial_velocity_error
rv_nb_transits
rv_template_teff
rv_template_logg
rv_template_fe_h
phot_variable_flag
l
b
ecl_lon
ecl_lat
priam_flags
teff_val
teff_percentile_lower
teff_percentile_upper
a_g_val
a_g_percentile_lower
a_g_percentile_upper
e_bp_min_rp_val
e_bp_min_rp_percentile_lower
e_bp_min_rp_percentile_upper
flame_flags
radius_val
radius_percentile_lower
radius_percentile_upper
lum_val
lum_percentile_lower
lum_percentile_upper
datalink_url
epoch_photometry_url

~~~
{: .output}

You can probably guess what many of these columns are by looking at
the names, but you should resist the temptation to guess.
To find out what the columns mean, [read the
documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).

If you want to know what can go wrong when you don't read the
documentation, [you might like this
article](https://www.vox.com/future-perfect/2019/6/4/18650969/married-women-miserable-fake-paul-dolan-happiness).

### Exercise

One of the other tables we'll use is
`gaiadr2.panstarrs1_original_valid`.  Use `load_table` to get the
metadata for this table.  How many columns are there and what are
their names?

~~~
# Solution

meta2 = Gaia.load_table('gaiadr2.panstarrs1_original_valid')
print(meta2)
~~~
{: .language-python}

~~~
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

~~~
{: .output}

~~~
# Solution

for column in meta2.columns:
    print(column.name)
~~~
{: .language-python}

~~~
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

~~~
{: .output}

## Writing queries

By now you might be wondering how we actually download the data.  With
tables this big, you generally don't.  Instead, you use queries to
select only the data you want.

A query is a string written in a query language like SQL; for the Gaia
database, the query language is a dialect of SQL called ADQL.

Here's an example of an ADQL query.

~~~
query1 = """SELECT 
TOP 10
source_id, ref_epoch, ra, dec, parallax 
FROM gaiadr2.gaia_source"""
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

To run this query, we use the `Gaia` object, which represents our
connection to the Gaia database, and invoke `launch_job`:

~~~
job1 = Gaia.launch_job(query1)
job1
~~~
{: .language-python}

~~~
~~~
{: .output}

The result is an object that represents the job running on a Gaia server.

If you print it, it displays metadata for the forthcoming table.

~~~
print(job1)
~~~
{: .language-python}

~~~
<Table length=10>
   name    dtype  unit                            description                             n_bad
--------- ------- ---- ------------------------------------------------------------------ -----
source_id   int64      Unique source identifier (unique within a particular Data Release)     0
ref_epoch float64   yr                                                    Reference epoch     0
       ra float64  deg                                                    Right ascension     0
      dec float64  deg                                                        Declination     0
 parallax float64  mas                                                           Parallax     3
Jobid: None
Phase: COMPLETED
Owner: None
Output file: sync_20201117154748.xml.gz
Results: None

~~~
{: .output}

Don't worry about `Results: None`.  That does not actually mean there
are no results.

However, `Phase: COMPLETED` indicates that the job is complete, so we
can get the results like this:

~~~
results1 = job1.get_results()
type(results1)
~~~
{: .language-python}

~~~
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
~~~
{: .output}

Each column has a name, units, and a data type.

For example, the units of `ra` and `dec` are degrees, and their data
type is `float64`, which is a 64-bit floating-point number, used to
store measurements with a fraction part.

This information comes from the Gaia database, and has been stored in
the Astropy `Table` by Astroquery.

### Exercise

Read [the documentation of this
table](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html)
and choose a column that looks interesting to you.  Add the column
name to the query and run it again.  What are the units of the column
you selected?  What is its data type?

## Asynchronous queries

`launch_job` asks the server to run the job "synchronously", which
normally means it runs immediately.  But synchronous jobs are limited
to 2000 rows.  For queries that return more rows, you should run
"asynchronously", which mean they might take longer to get started.

If you are not sure how many rows a query will return, you can use the
SQL command `COUNT` to find out how many rows are in the result
without actually returning them.  We'll see an example of this later.

The results of an asynchronous query are stored in a file on the
server, so you can start a query and come back later to get the
results.

For anonymous users, files are kept for three days.

As an example, let's try a query that's similar to `query1`, with two changes:

* It selects the first 3000 rows, so it is bigger than we should run
synchronously.

* It uses a new keyword, `WHERE`.

~~~
query2 = """SELECT TOP 3000
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1
"""
~~~
{: .language-python}
A `WHERE` clause indicates which rows we want; in this case, the query
selects only rows "where" `parallax` is less than 1.  This has the
effect of selecting stars with relatively low parallax, which are
farther away.  We'll use this clause to exclude nearby stars that are
unlikely to be part of GD-1.

`WHERE` is one of the most common clauses in ADQL/SQL, and one of the
most useful, because it allows us to select only the rows we need from
the database.

We use `launch_job_async` to submit an asynchronous query.

~~~
job2 = Gaia.launch_job_async(query2)
print(job2)
~~~
{: .language-python}

~~~
INFO: Query finished. [astroquery.utils.tap.core]
<Table length=3000>
   name    dtype  unit                            description                            
--------- ------- ---- ------------------------------------------------------------------
source_id   int64      Unique source identifier (unique within a particular Data Release)
ref_epoch float64   yr                                                    Reference epoch
       ra float64  deg                                                    Right ascension
      dec float64  deg                                                        Declination
 parallax float64  mas                                                           Parallax
Jobid: 1605646069281O
Phase: COMPLETED
Owner: None
Output file: async_20201117154749.vot
Results: None

~~~
{: .output}

And here are the results.

~~~
results2 = job2.get_results()
results2
~~~
{: .language-python}

~~~
~~~
{: .output}

You might notice that some values of `parallax` are negative.  As
[this FAQ
explains](https://www.cosmos.esa.int/web/gaia/archive-tips#negative%20parallax),
"Negative parallaxes are caused by errors in the observations."
Negative parallaxes have "no physical meaning," but they can be a
"useful diagnostic on the quality of the astrometric solution."

Later we will see an example where we use `parallax` and
`parallax_error` to identify stars where the distance estimate is
likely to be inaccurate.

### Exercise

The clauses in a query have to be in the right order.  Go back and
change the order of the clauses in `query2` and run it again.

The query should fail, but notice that you don't get much useful
debugging information.

For this reason, developing and debugging ADQL queries can be really
hard.  A few suggestions that might help:

* Whenever possible, start with a working query, either an example you
find online or a query you have used in the past.

* Make small changes and test each change before you continue.

* While you are debugging, use `TOP` to limit the number of rows in
the result.  That will make each attempt run faster, which reduces
your testing time.

* Launching test queries synchronously might make them start faster, too.

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

### Exercise

[Read about SQL operators
here](https://www.w3schools.com/sql/sql_operators.asp) and then modify
the previous query to select rows where `bp_rp` is between `-0.75` and
`2`.

You can [read about this variable
here](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).

~~~
# Solution

# This is what most people will probably do

query = """SELECT TOP 10
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp > -0.75 AND bp_rp < 2
"""
~~~
{: .language-python}
~~~
# Solution

# But if someone notices the BETWEEN operator, 
# they might do this

query = """SELECT TOP 10
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp BETWEEN -0.75 AND 2
"""
~~~
{: .language-python}
This [Hertzsprung-Russell
diagram](https://sci.esa.int/web/gaia/-/60198-gaia-hertzsprung-russell-diagram)
shows the BP-RP color and luminosity of stars in the Gaia catalog.

Selecting stars with `bp-rp` less than 2 excludes many [class M dwarf
stars](https://xkcd.com/2360/), which are low temperature, low
luminosity.  A star like that at GD-1's distance would be hard to
detect, so if it is detected, it it more likely to be in the
foreground.

## Cleaning up

Asynchronous jobs have a `jobid`.

~~~
job1.jobid, job2.jobid
~~~
{: .language-python}

~~~
~~~
{: .output}

Which you can use to remove the job from the server.

~~~
Gaia.remove_jobs([job2.jobid])
~~~
{: .language-python}

~~~
Removed jobs: '['1605646069281O']'.

~~~
{: .output}

If you don't remove it job from the server, it will be removed
eventually, so don't feel too bad if you don't clean up after
yourself.

## Formatting queries

So far the queries have been string "literals", meaning that the
entire string is part of the program.
But writing queries yourself can be slow, repetitive, and error-prone.

It is often a good idea to write Python code that assembles a query
for you.  One useful tool for that is the [string `format`
method](https://www.w3schools.com/python/ref_string_format.asp).

As an example, we'll divide the previous query into two parts; a list
of column names and a "base" for the query that contains everything
except the column names.

Here's the list of columns we'll select.  

~~~
columns = 'source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity'
~~~
{: .language-python}
And here's the base; it's a string that contains at least one format
specifier in curly brackets (braces).

~~~
query3_base = """SELECT TOP 10 
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
~~~
{: .output}

But if you print it, the line breaks appear as... line breaks.

~~~
print(query3)
~~~
{: .language-python}

~~~
SELECT TOP 10 
source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2


~~~
{: .output}

Notice that the format specifier has been replaced with the value of `columns`.

Let's run it and see if it works:

~~~
job3 = Gaia.launch_job(query3)
print(job3)
~~~
{: .language-python}

~~~
<Table length=10>
      name       dtype    unit                              description                             n_bad
--------------- ------- -------- ------------------------------------------------------------------ -----
      source_id   int64          Unique source identifier (unique within a particular Data Release)     0
             ra float64      deg                                                    Right ascension     0
            dec float64      deg                                                        Declination     0
           pmra float64 mas / yr                         Proper motion in right ascension direction     0
          pmdec float64 mas / yr                             Proper motion in declination direction     0
       parallax float64      mas                                                           Parallax     0
 parallax_error float64      mas                                         Standard error of parallax     0
radial_velocity float64   km / s                                                    Radial velocity    10
Jobid: None
Phase: COMPLETED
Owner: None
Output file: sync_20201117154752.xml.gz
Results: None

~~~
{: .output}

~~~
results3 = job3.get_results()
results3
~~~
{: .language-python}

~~~
~~~
{: .output}

Good so far.

### Exercise

This query always selects sources with `parallax` less than 1.  But
suppose you want to take that upper bound as an input.

Modify `query3_base` to replace `1` with a format specifier like
`{max_parallax}`.  Now, when you call `format`, add a keyword argument
that assigns a value to `max_parallax`, and confirm that the format
specifier gets replaced with the value you provide.

~~~
# Solution

query4_base = """SELECT TOP 10
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < {max_parallax} AND 
bp_rp BETWEEN -0.75 AND 2
"""
~~~
{: .language-python}
~~~
# Solution

query4 = query4_base.format(columns=columns,
                          max_parallax=0.5)
print(query)
~~~
{: .language-python}

~~~
SELECT TOP 10
source_id, ref_epoch, ra, dec, parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp BETWEEN -0.75 AND 2


~~~
{: .output}

**Style note:**  You might notice that the variable names in this
notebook are numbered, like `query1`, `query2`, etc.

The advantage of this style is that it isolates each section of the
notebook from the others, so if you go back and run the cells out of
order, it's less likely that you will get unexpected interactions.

A drawback of this style is that it can be a nuisance to update the
notebook if you add, remove, or reorder a section.

What do you think of this choice?  Are there alternatives you prefer?

## Summary

This notebook demonstrates the following steps:

1. Making a connection to the Gaia server,

2. Exploring information about the database and the tables it contains,

3. Writing a query and sending it to the server, and finally

4. Downloading the response from the server as an Astropy `Table`.

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
it synchronously, which might complete faster (but it doesn't seem to
make much difference).  If it might return more than 3000 rows, you
should run it asynchronously.

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
