---
title: "Coordinate Transformations"
teaching: 3000
exercises: 0
questions:

- "How do we transform celestial coordinates from one frame to another and save results in files?"

objectives:

- "Use Python string formatting to compose more complex ADQL queries."

- "Work with coordinates and other quantities that have units."

- "Download the results of a query and store them in a file."

keypoints:

- "For measurements with units, use `Quantity` objects that represent units explicitly and check for errors."

- "Use the `format` function to compose queries; it is often faster and less error-prone."

- "Develop queries incrementally: start with something simple, test it, and add a little bit at a time."

- "Once you have a query working, save the data in a local file.  If you shut down the notebook and come back to it later, you can reload the file; you don't have to run the query again."

---

{% include links.md %}

# Coordinates and units

This is the second in a series of notebooks related to astronomy data.

As a running example, we are replicating parts of the analysis in a
recent paper, "[Off the beaten path: Gaia reveals GD-1 stars outside
of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian M.
Price-Whelan and Ana Bonaca.

In the first notebook, we wrote ADQL queries and used them to select
and download data from the Gaia server.

In this notebook, we'll pick up where we left off and write a query to
select stars from the region of the sky where we expect GD-1 to be.

## Outline

We'll start with an example that does a "cone search"; that is, it
selects stars that appear in a circular region of the sky.

Then, to select stars in the vicinity of GD-1, we'll:

* Use `Quantity` objects to represent measurements with units.

* Use Astropy to convert coordinates from one frame to another.

* Use the ADQL keywords `POLYGON`, `CONTAINS`, and `POINT` to select
stars that fall within a polygonal region.

* Submit a query and download the results.

* Store the results in a FITS file.

After completing this lesson, you should be able to

* Use Python string formatting to compose more complex ADQL queries.

* Work with coordinates and other quantities that have units.

* Download the results of a query and store them in a file.

## Selecting a region

One of the most common ways to restrict a query is to select stars in
a particular region of the sky.

For example, here's a query from the [Gaia archive
documentation](https://gea.esac.esa.int/archive-help/adql/examples/index.html)
that selects "all the objects ... in a circular region centered at
(266.41683, -29.00781) with a search radius of 5 arcmin (0.08333
deg)."

~~~
query = """
SELECT 
TOP 10 source_id
FROM gaiadr2.gaia_source
WHERE 1=CONTAINS(
  POINT(ra, dec),
  CIRCLE(266.41683, -29.00781, 0.08333333))
"""
~~~
{: .language-python}
This query uses three keywords that are specific to ADQL (not SQL):

* `POINT`: a location in [ICRS
coordinates](https://en.wikipedia.org/wiki/International_Celestial_Reference_System),
specified in degrees of right ascension and declination.

* `CIRCLE`: a circle where the first two values are the coordinates of
the center and the third is the radius in degrees.

* `CONTAINS`: a function that returns `1` if a `POINT` is contained in
a shape and `0` otherwise.

Here is the [documentation of
`CONTAINS`](http://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html#tth_sEc4.2.12).

A query like this is called a cone search because it selects stars in a cone.

Here's how we run it.

~~~
from astroquery.gaia import Gaia

job = Gaia.launch_job(query)
job
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

<astroquery.utils.tap.model.job.Job at 0x7f59bd93e490>
~~~
{: .output}

~~~
result = job.get_results()
result
~~~
{: .language-python}

~~~
<Table length=10>
     source_id     
       int64       
-------------------
4057468321929794432
4057468287575835392
4057482027171038976
4057470349160630656
4057470039924301696
4057469868125641984
4057468351995073024
[Output truncated]
~~~
{: .output}

> ## Exercise
> 
> When you are debugging queries like this, you can use `TOP` to limit
> the size of the results, but then you still don't know how big the
> results will be.
> 
> An alternative is to use `COUNT`, which asks for the number of rows
> that would be selected, but it does not return them.
> 
> In the previous query, replace `TOP 10 source_id` with
> `COUNT(source_id)` and run the query again.  How many stars has Gaia
> identified in the cone we searched?
> > 
> > ~~~
> > 
> > query = """
> > SELECT 
> > COUNT(source_id)
> > FROM gaiadr2.gaia_source
> > WHERE 1=CONTAINS(
> >   POINT(ra, dec),
> >   CIRCLE(266.41683, -29.00781, 0.08333333))
> > """
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}


## Getting GD-1 Data

From the Price-Whelan and Bonaca paper, we will try to reproduce
Figure 1, which includes this representation of stars likely to belong
to GD-1:

<img
src="https://github.com/datacarpentry/astronomy-python/raw/gh-pages/fig/gd1-4.png">

The axes of this figure are defined so the x-axis is aligned with the
stars in GD-1, and the y-axis is perpendicular.

* Along the x-axis ($\phi_1$) the figure extends from -100 to 20 degrees.

* Along the y-axis ($\phi_2$) the figure extends from about -8 to 4 degrees.

Ideally, we would select all stars from this rectangle, but there are
more than 10 million of them, so

* That would be difficult to work with,

* As anonymous Gaia users, we are limited to 3 million rows in a
single query, and

* While we are developing and testing code, it will be faster to work
with a smaller dataset.

So we'll start by selecting stars in a smaller rectangle near the
center of GD-1, from -55 to -45 degrees $\phi_1$ and -8 to 4 degrees
$\phi_2$.

But first we let's see how to represent quantities with units like degrees.

## Working with coordinates

Coordinates are physical quantities, which means that they have two
parts, a value and a unit.

For example, the coordinate $30^{\circ}$ has value 30 and its units are degrees.

Until recently, most scientific computation was done with values only;
units were left out of the program altogether, [often with
catastrophic
results](https://en.wikipedia.org/wiki/Mars_Climate_Orbiter#Cause_of_failure).

Astropy provides tools for including units explicitly in computations,
which makes it possible to detect errors before they cause disasters.

To use Astropy units, we import them like this:

~~~
import astropy.units as u
~~~
{: .language-python}
`u` is an object that contains most common units and all SI units.

You can use `dir` to list them, but you should also [read the
documentation](https://docs.astropy.org/en/stable/units/).

~~~
dir(u)
~~~
{: .language-python}

~~~
['A',
 'AA',
 'AB',
 'ABflux',
 'ABmag',
 'AU',
 'Angstrom',
 'B',
 'Ba',
 'Barye',
 'Bi',
[Output truncated]
~~~
{: .output}

To create a quantity, we multiply a value by a unit.

~~~
quantity = 30 * u.degree
type(quantity)
~~~
{: .language-python}

~~~
astropy.units.quantity.Quantity
~~~
{: .output}

The result is a `Quantity` object.

Jupyter knows how to display `Quantities` like this:

~~~
quantity
~~~
{: .language-python}

~~~
<Quantity 30. deg>
~~~
{: .output}

## Transforming coordinates

Astropy provides a `SkyCoord` object that represents sky coordinates
relative to a specified frame.

The following example creates a `SkyCoord` object that represents the
approximate coordinates of
[Betelgeuse](http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=Betelgeuse)
(alf Ori) in the ICRS frame.

~~~
from astropy.coordinates import SkyCoord

ra = 88.8 * u.degree
dec = 7.4 * u.degree
coord_icrs = SkyCoord(ra=ra, dec=dec, frame='icrs')

coord_icrs
~~~
{: .language-python}

~~~
<SkyCoord (ICRS): (ra, dec) in deg
    (88.8, 7.4)>
~~~
{: .output}

`SkyCoord` provides a function that transforms to other frames.
For example, we can transform `coords_icrs` to Galactic coordinates like this:

~~~
coord_galactic = coord_icrs.transform_to('galactic')
coord_galactic
~~~
{: .language-python}

~~~
<SkyCoord (Galactic): (l, b) in deg
    (199.79693102, -8.95591653)>
~~~
{: .output}

To transform to and from GD-1 coordinates, we'll use a frame defined
by [Gala](https://gala-astro.readthedocs.io/en/latest/), which is an
Astropy-affiliated library that provides tools for galactic dynamics.

Gala provides `GD1Koposov10`, which is "[a Heliocentric spherical
coordinate system defined by the orbit of the GD-1
stream](https://gala-astro.readthedocs.io/en/latest/_modules/gala/coordinates/gd1.html)"

~~~
from gala.coordinates import GD1Koposov10

gd1_frame = GD1Koposov10()
gd1_frame
~~~
{: .language-python}

~~~
<GD1Koposov10 Frame>
~~~
{: .output}

We can use it to find the coordinates of Betelgeuse in the GD-1 frame,
like this:

~~~
coord_gd1 = coord_icrs.transform_to(gd1_frame)
coord_gd1
~~~
{: .language-python}

~~~
<SkyCoord (GD1Koposov10): (phi1, phi2) in deg
    (-94.97222038, 34.5813813)>
~~~
{: .output}

> ## Exercise
> 
> Let's find the location of GD-1 in ICRS coordinates.
> 
> 1. Create a `SkyCoord` object at 0°, 0° in the GD-1 frame.
> 
> 2. Transform it to the ICRS frame.
> 
> Hint: Because ICRS is built into Astropy, you can specify it by name,
> `icrs` (as we did with `galactic`).
> > 
> > ~~~
> > 
> > coord_gd1 = SkyCoord(0*u.degree, 0*u.degree, frame=gd1_frame)
> > 
> > # Note: because ICRS is built into Astropy, 
> > # we can identify it by name
> > coord_gd1.transform_to('icrs')
> > 
> > # More formally, we could instantiate it
> > from astropy.coordinates import ICRS
> > icrs_frame = ICRS()
> > coord_gd1.transform_to(icrs_frame)
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}


## Selecting a rectangle

Now we'll use these coordinate transformations to define a rectangle
in the GD-1 frame and transform it to ICRS.

The following variables define the boundaries of the rectangle in
$\phi_1$ and $\phi_2$.

~~~
phi1_min = -55 * u.degree
phi1_max = -45 * u.degree
phi2_min = -8 * u.degree
phi2_max = 4 * u.degree
~~~
{: .language-python}
To represent a rectangle, we'll use two lists of coordinates and
multiply by their units.

~~~
def make_rectangle(x1, x2, y1, y2):
    """Return the corners of a rectangle."""
    xs = [x1, x1, x2, x2, x1]
    ys = [y1, y2, y2, y1, y1]
    return xs, ys
~~~
{: .language-python}
~~~
phi1_rect, phi2_rect = make_rectangle(
    phi1_min, phi1_max, phi2_min, phi2_max)
~~~
{: .language-python}
`phi1_rect` and `phi2_rect` represent the coordinates of the corners
of a rectangle in the GD-1 frame.

In order to use them in a Gaia query, we have to convert them to ICRS.

~~~
import gala.coordinates as gc

corners = SkyCoord(phi1=phi1_rect, phi2=phi2_rect, frame=gd1_frame)
corners
~~~
{: .language-python}

~~~
<SkyCoord (GD1Koposov10): (phi1, phi2) in deg
    [(-55., -8.), (-55.,  4.), (-45.,  4.), (-45., -8.), (-55., -8.)]>
~~~
{: .output}

Now we can use `transform_to` to convert to ICRS coordinates.

~~~
import astropy.coordinates as coord

corners_icrs = corners.transform_to('icrs')
corners_icrs
~~~
{: .language-python}

~~~
<SkyCoord (ICRS): (ra, dec) in deg
    [(146.27533314, 19.26190982), (135.42163944, 25.87738723),
     (141.60264825, 34.3048303 ), (152.81671045, 27.13611254),
     (146.27533314, 19.26190982)]>
~~~
{: .output}

Notice that a rectangle in one coordinate system is not necessarily a
rectangle in another.  In this example, the result is a polygon.

## Selecting a polygon

In order to use this polygon as part of an ADQL query, we have to
convert it to a string with a comma-separated list of coordinates, as
in this example:

```
"""
POLYGON(143.65, 20.98, 
        134.46, 26.39, 
        140.58, 34.85, 
        150.16, 29.01)
"""
```

`corners_icrs` behaves like a list, so we can use a `for` loop to
iterate through the points.

~~~
for point in corners_icrs:
    print(point)
~~~
{: .language-python}

~~~
<SkyCoord (ICRS): (ra, dec) in deg
    (146.27533314, 19.26190982)>
<SkyCoord (ICRS): (ra, dec) in deg
    (135.42163944, 25.87738723)>
<SkyCoord (ICRS): (ra, dec) in deg
    (141.60264825, 34.3048303)>
<SkyCoord (ICRS): (ra, dec) in deg
    (152.81671045, 27.13611254)>
<SkyCoord (ICRS): (ra, dec) in deg
    (146.27533314, 19.26190982)>

~~~
{: .output}

From that, we can select the coordinates `ra` and `dec`:

~~~
for point in corners_icrs:
    print(point.ra, point.dec)
~~~
{: .language-python}

~~~
146d16m31.1993s 19d15m42.8754s
135d25m17.902s 25d52m38.594s
141d36m09.5337s 34d18m17.3891s
152d49m00.1576s 27d08m10.0051s
146d16m31.1993s 19d15m42.8754s

~~~
{: .output}

The results are quantities with units, but if we select the `value`
part, we get a dimensionless floating-point number.

~~~
for point in corners_icrs:
    print(point.ra.value, point.dec.value)
~~~
{: .language-python}

~~~
146.27533313607782 19.261909820533692
135.42163944306296 25.87738722767213
141.60264825107333 34.304830296257144
152.81671044675923 27.136112541397996
146.27533313607782 19.261909820533692

~~~
{: .output}

We can use string `format` to convert these numbers to strings.

~~~
point_base = "{point.ra.value}, {point.dec.value}"

t = [point_base.format(point=point)
     for point in corners_icrs]
t
~~~
{: .language-python}

~~~
['146.27533313607782, 19.261909820533692',
 '135.42163944306296, 25.87738722767213',
 '141.60264825107333, 34.304830296257144',
 '152.81671044675923, 27.136112541397996',
 '146.27533313607782, 19.261909820533692']
~~~
{: .output}

The result is a list of strings, which we can join into a single
string using `join`.

~~~
point_list = ', '.join(t)
point_list
~~~
{: .language-python}

~~~
'146.27533313607782, 19.261909820533692, 135.42163944306296, 25.87738722767213, 141.60264825107333, 34.304830296257144, 152.81671044675923, 27.136112541397996, 146.27533313607782, 19.261909820533692'
~~~
{: .output}

Notice that we invoke `join` on a string and pass the list as an argument.

Before we can assemble the query, we need `columns` again (as we saw
in the previous notebook).

~~~
columns = 'source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity'
~~~
{: .language-python}
Here's the base for the query, with format specifiers for `columns`
and `point_list`.

~~~
query_base = """SELECT {columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({point_list}))
"""
~~~
{: .language-python}
And here's the result:

~~~
query = query_base.format(columns=columns, 
                          point_list=point_list)
print(query)
~~~
{: .language-python}

~~~
SELECT source_id, ra, dec, pmra, pmdec, parallax, parallax_error, radial_velocity
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON(146.27533313607782, 19.261909820533692, 135.42163944306296, 25.87738722767213, 141.60264825107333, 34.304830296257144, 152.81671044675923, 27.136112541397996, 146.27533313607782, 19.261909820533692))


~~~
{: .output}

As always, we should take a minute to proof-read the query before we launch it.

The result will be bigger than our previous queries, so it will take a
little longer.

~~~
job = Gaia.launch_job_async(query)
print(job)
~~~
{: .language-python}

~~~
INFO: Query finished. [astroquery.utils.tap.core]
<Table length=140340>
      name       dtype    unit                              description                             n_bad 
--------------- ------- -------- ------------------------------------------------------------------ ------
      source_id   int64          Unique source identifier (unique within a particular Data Release)      0
             ra float64      deg                                                    Right ascension      0
            dec float64      deg                                                        Declination      0
           pmra float64 mas / yr                         Proper motion in right ascension direction      0
          pmdec float64 mas / yr                             Proper motion in declination direction      0
       parallax float64      mas                                                           Parallax      0
 parallax_error float64      mas                                         Standard error of parallax      0
[Output truncated]
~~~
{: .output}

Here are the results.

~~~
results = job.get_results()
len(results)
~~~
{: .language-python}

~~~
140340
~~~
{: .output}

There are more than 100,000 stars in this polygon, but that's a
manageable size to work with.

## Saving results

This is the set of stars we'll work with in the next step.  But since
we have a substantial dataset now, this is a good time to save it.

Storing the data in a file means we can shut down this notebook and
pick up where we left off without running the previous query again.

Astropy `Table` objects provide `write`, which writes the table to disk.

~~~
filename = 'gd1_results.fits'
results.write(filename, overwrite=True)
~~~
{: .language-python}
Because the filename ends with `fits`, the table is written in the
[FITS format](https://en.wikipedia.org/wiki/FITS), which preserves the
metadata associated with the table.

If the file already exists, the `overwrite` argument causes it to be
overwritten.

To see how big the file is, we can use `ls` with the `-lh` option,
which prints information about the file including its size in
human-readable form.

~~~
!ls -lh gd1_results.fits
~~~
{: .language-python}

~~~
-rw-rw-r-- 1 downey downey 8.6M Dec 29 11:47 gd1_results.fits

~~~
{: .output}

The file is about 8.6 MB.  If you are using Windows, `ls` might not
work; in that case, try:

```
!dir gd1_results.fits
```

## Summary

In this notebook, we composed more complex queries to select stars
within a polygonal region of the sky.  Then we downloaded the results
and saved them in a FITS file.

In the next notebook, we'll reload the data from this file and
replicate the next step in the analysis, using proper motion to
identify stars likely to be in GD-1.

## Best practices

* For measurements with units, use `Quantity` objects that represent
units explicitly and check for errors.

* Use the `format` function to compose queries; it is often faster and
less error-prone.

* Develop queries incrementally: start with something simple, test it,
and add a little bit at a time.

* Once you have a query working, save the data in a local file.  If
you shut down the notebook and come back to it later, you can reload
the file; you don't have to run the query again.
