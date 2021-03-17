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

# 2. Coordinates and Units

In the previous lesson, we wrote ADQL queries and used them to select
and download data from the Gaia server.

In this lesson, we'll pick up where we left off and write a query to
select stars from a particular region of the sky.

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

## Working with Units

The measurements we will work with are physical quantities, which
means that they have two parts, a value and a unit.
For example, the coordinate $30^{\circ}$ has value 30 and its units are degrees.

Until recently, most scientific computation was done with values only;
units were left out of the program altogether, [often with
catastrophic
results](https://en.wikipedia.org/wiki/Mars_Climate_Orbiter#Cause_of_failure).

Astropy provides tools for including units explicitly in computations,
which makes it possible to detect errors before they cause disasters.

To use Astropy units, we import them like this:


```python

~~~
import astropy.units as u
~~~
{: .language-python}
```

`u` is an object that contains most common units and all SI units.

You can use `dir` to list them, but you should also [read the
documentation](https://docs.astropy.org/en/stable/units/).


```python

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

```




    



To create a quantity, we multiply a value by a unit.


```python

~~~
angle = 10 * u.degree
type(angle)
~~~
{: .language-python}

~~~
astropy.units.quantity.Quantity
~~~
{: .output}

```




    



The result is a `Quantity` object.
Jupyter knows how to display `Quantities` like this:


```python

~~~
angle
~~~
{: .language-python}

~~~
<Quantity 10. deg>
~~~
{: .output}

```




$10 \; \mathrm{{}^{\circ}}$



Quantities provide a method called `to` that converts to other units.
For example, we can compute the number of arcminutes in `angle`:


```python

~~~
angle_arcmin = angle.to(u.arcmin)
angle_arcmin
~~~
{: .language-python}

~~~
<Quantity 600. arcmin>
~~~
{: .output}

```




$600 \; \mathrm{{}^{\prime}}$



If you add quantities, Astropy converts them to compatible units, if possible:


```python

~~~
angle + 30 * u.arcmin
~~~
{: .language-python}

~~~
<Quantity 10.5 deg>
~~~
{: .output}

```




$10.5 \; \mathrm{{}^{\circ}}$



If the units are not compatible, you get an error.
For example:

```
angle + 5 * u.second
```

causes a `UnitConversionError`.

> ## Exercise
> 
> Create a quantity that represents 5
> [arcminutes](https://en.wikipedia.org/wiki/Minute_and_second_of_arc)
> and assign it to a variable called `radius`.
> 
> Then convert it to degrees.


```python

~~~
## Solution

radius = 5 * u.arcmin
print(radius)

radius.to(u.degree)
~~~
{: .language-python}

~~~
5.0 arcmin

<Quantity 0.08333333 deg>
~~~
{: .output}

```

    




$0.083333333 \; \mathrm{{}^{\circ}}$



## Selecting a Region

One of the most common ways to restrict a query is to select stars in
a particular region of the sky.
For example, here's a query from the [Gaia archive
documentation](https://gea.esac.esa.int/archive-help/adql/examples/index.html)
that selects objects in a circular region centered at (88.8, 7.4) with
a search radius of 5 arcmin (0.08333 deg).


```python

~~~
query_cone = """SELECT 
TOP 10 
source_id
FROM gaiadr2.gaia_source
WHERE 1=CONTAINS(
  POINT(ra, dec),
  CIRCLE(88.8, 7.4, 0.08333333))
"""
~~~
{: .language-python}
```

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


```python

~~~
from astroquery.gaia import Gaia

job = Gaia.launch_job(query_cone)
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

<astroquery.utils.tap.model.job.Job at 0x7f277785fa30>
~~~
{: .output}

```

    




    




```python

~~~
results = job.get_results()
results
~~~
{: .language-python}

~~~
<Table length=10>
     source_id     
       int64       
-------------------
3322773965056065536
3322773758899157120
3322774068134271104
3322773930696320512
3322774377374425728
3322773724537891456
3322773724537891328
[Output truncated]
~~~
{: .output}

```




<i>Table length=10</i>
<table id="table139807485721280" class="table-striped table-bordered table-condensed">
<thead><tr><th>source_id</th></tr></thead>
<thead><tr><th>int64</th></tr></thead>
<tr><td>3322773965056065536</td></tr>
<tr><td>3322773758899157120</td></tr>
<tr><td>3322774068134271104</td></tr>
<tr><td>3322773930696320512</td></tr>
<tr><td>3322774377374425728</td></tr>
<tr><td>3322773724537891456</td></tr>
<tr><td>3322773724537891328</td></tr>
<tr><td>3322773930696321792</td></tr>
<tr><td>3322773724537890944</td></tr>
<tr><td>3322773930696322176</td></tr>
</table>



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


```python
>
> > ## Solution
> > 
> > ~~~
> > 
> > query = """SELECT 
> > COUNT(source_id)
> > FROM gaiadr2.gaia_source
> > WHERE 1=CONTAINS(
> >   POINT(ra, dec),
> >   CIRCLE(88.8, 7.4, 0.08333333))
> > """
> > 
> > job = Gaia.launch_job(query)
> > results = job.get_results()
> > results
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}
```

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

But first we let's see how to represent these coordinates with Astropy.

## Transforming coordinates

Astropy provides a `SkyCoord` object that represents sky coordinates
relative to a specified frame.

The following example creates a `SkyCoord` object that represents the
approximate coordinates of
[Betelgeuse](http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=Betelgeuse)
(alf Ori) in the ICRS frame.

[ICRS](https://www.iers.org/IERS/EN/Science/ICRS/ICRS.html) is the
"International Celestial Reference System", adopted in 1997 by the
International Astronomical Union.


```python

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

```




    



`SkyCoord` provides a function that transforms to other frames.
For example, we can transform `coords_icrs` to Galactic coordinates like this:


```python

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

```




    



Notice that in the Galactic frame, the coordinates are called `l` and
`b`, not `ra` and `dec`.

To transform to and from GD-1 coordinates, we'll use a frame defined
by [Gala](https://gala-astro.readthedocs.io/en/latest/), which is an
Astropy-affiliated library that provides tools for galactic dynamics.

Gala provides
[`GD1Koposov10`](https://gala-astro.readthedocs.io/en/latest/_modules/gala/coordinates/gd1.html),
which is "a Heliocentric spherical coordinate system defined by the
orbit of the GD-1 stream".


```python

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

```




    



We can use it to find the coordinates of Betelgeuse in the GD-1 frame,
like this:


```python

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

```




    



Notice that the coordinates are called `phi1` and `phi2`.
These are the coordinates shown in the figure from the paper, above.

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


```python
>
> > ## Solution
> > 
> > ~~~
> > 
> > origin_gd1 = SkyCoord(0*u.degree, 0*u.degree, frame=gd1_frame)
> > 
> > # OR
> > 
> > origin_gd1 = SkyCoord(phi1=0*u.degree, 
> >                       phi2=0*u.degree, 
> >                       frame=gd1_frame)
> > 
> > # Note: because ICRS is built into Astropy, 
> > # we can identify it by string name
> > origin_gd1.transform_to('icrs')
> > 
> > # More formally, we could instantiate it
> > from astropy.coordinates import ICRS
> > icrs_frame = ICRS()
> > origin_gd1.transform_to(icrs_frame)
> > ~~~
> > {: .language-python}
> {: .solution}
{: .challenge}
```

Notice that the origin of the GD-1 frame maps to `ra=200`, exactly, in
ICRS.  That's by design.

## Selecting a rectangle

Now we'll use these coordinate transformations to define a rectangle
in the GD-1 frame and transform it to ICRS.

The following variables define the boundaries of the rectangle in
$\phi_1$ and $\phi_2$.


```python

~~~
phi1_min = -55 * u.degree
phi1_max = -45 * u.degree
phi2_min = -8 * u.degree
phi2_max = 4 * u.degree
~~~
{: .language-python}
```

To create a rectangle, we'll use the following function, which takes
the lower and upper bounds as parameters.


```python

~~~
def make_rectangle(x1, x2, y1, y2):
    """Return the corners of a rectangle."""
    xs = [x1, x1, x2, x2, x1]
    ys = [y1, y2, y2, y1, y1]
    return xs, ys
~~~
{: .language-python}
```

The return value is a tuple containing a list of coordinates in `phi1`
followed by a list of coordinates in `phi2`.


```python

~~~
phi1_rect, phi2_rect = make_rectangle(
    phi1_min, phi1_max, phi2_min, phi2_max)
~~~
{: .language-python}
```

`phi1_rect` and `phi2_rect` contains the coordinates of the corners of
a rectangle in the GD-1 frame.

In order to use them in a Gaia query, we have to convert them to ICRS.
First we'll put them into a `SkyCoord` object.


```python

~~~
corners = SkyCoord(phi1=phi1_rect, phi2=phi2_rect, frame=gd1_frame)
corners
~~~
{: .language-python}

~~~
<SkyCoord (GD1Koposov10): (phi1, phi2) in deg
    [(-55., -8.), (-55.,  4.), (-45.,  4.), (-45., -8.), (-55., -8.)]>
~~~
{: .output}

```




    



Now we can use `transform_to` to convert to ICRS coordinates.


```python

~~~
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

```




    



Notice that a rectangle in one coordinate system is not necessarily a
rectangle in another.  In this example, the result is a
(non-rectangular) polygon.

## Defining a polygon

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

`SkyCoord` provides `to_string`, which produces a list of strings.


```python

~~~
t = corners_icrs.to_string()
t
~~~
{: .language-python}

~~~
['146.275 19.2619',
 '135.422 25.8774',
 '141.603 34.3048',
 '152.817 27.1361',
 '146.275 19.2619']
~~~
{: .output}

```




    



We can use the Python string function `join` to join `t` into a single
string (with spaces between the pairs):


```python

~~~
s = ' '.join(t)
s
~~~
{: .language-python}

~~~
'146.275 19.2619 135.422 25.8774 141.603 34.3048 152.817 27.1361 146.275 19.2619'
~~~
{: .output}

```




    



That's almost what we need, but we have to replace the spaces with commas.


```python

~~~
s.replace(' ', ', ')
~~~
{: .language-python}

~~~
'146.275, 19.2619, 135.422, 25.8774, 141.603, 34.3048, 152.817, 27.1361, 146.275, 19.2619'
~~~
{: .output}

```




    



The following function combines these steps.


```python

~~~
def skycoord_to_string(skycoord):
    """Convert SkyCoord to string."""
    t = skycoord.to_string()
    s = ' '.join(t)
    return s.replace(' ', ', ')
~~~
{: .language-python}
```

Here's how we use it.


```python

~~~
point_list = skycoord_to_string(corners_icrs)
point_list
~~~
{: .language-python}

~~~
'146.275, 19.2619, 135.422, 25.8774, 141.603, 34.3048, 152.817, 27.1361, 146.275, 19.2619'
~~~
{: .output}

```




    



## Assembling the query

Now we're ready to assemble the query. 
We need `columns` again (as we saw in the previous lesson).


```python

~~~
columns = 'source_id, ra, dec, pmra, pmdec, parallax'
~~~
{: .language-python}
```

And here's the query base we used in the previous lesson:


```python

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
```

Now we'll add a `WHERE` clause to select stars in the polygon we defined.


```python

~~~
query4_base = """SELECT
TOP 10
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({point_list}))
"""
~~~
{: .language-python}
```

The query base contains format specifiers for `columns` and `point_list`.

We'll use `format` to fill in these values.


```python

~~~
query4 = query4_base.format(columns=columns, 
                          point_list=point_list)
print(query4)
~~~
{: .language-python}

~~~
SELECT
TOP 10
source_id, ra, dec, pmra, pmdec
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON(146.275, 19.2619, 135.422, 25.8774, 141.603, 34.3048, 152.817, 27.1361, 146.275, 19.2619))


~~~
{: .output}

```

    

As always, we should take a minute to proof-read the query before we launch it.


```python

~~~
job = Gaia.launch_job_async(query4)
print(job)
~~~
{: .language-python}

~~~
INFO: Query finished. [astroquery.utils.tap.core]
<Table length=10>
   name    dtype    unit                              description                            
--------- ------- -------- ------------------------------------------------------------------
source_id   int64          Unique source identifier (unique within a particular Data Release)
       ra float64      deg                                                    Right ascension
      dec float64      deg                                                        Declination
     pmra float64 mas / yr                         Proper motion in right ascension direction
    pmdec float64 mas / yr                             Proper motion in declination direction
Jobid: 1615815873808O
Phase: COMPLETED
[Output truncated]
~~~
{: .output}

```

    

Here are the results.


```python

~~~
results = job.get_results()
results
~~~
{: .language-python}

~~~
<Table length=10>
    source_id              ra         ...        pmdec       
                          deg         ...       mas / yr     
      int64             float64       ...       float64      
------------------ ------------------ ... -------------------
637987125186749568 142.48301935991023 ...   2.941813096629439
638285195917112960 142.25452941346344 ... -12.165984395577347
638073505568978688 142.64528557468074 ...  -7.950659620550862
638086386175786752 142.57739430926034 ...  -2.584105480335548
638049655615392384 142.58913564478618 ...  -4.941079187010136
638267565075964032 141.81762228999614 ...  1.8838892877285924
[Output truncated]
~~~
{: .output}

```




<i>Table length=10</i>
<table id="table139807251332016" class="table-striped table-bordered table-condensed">
<thead><tr><th>source_id</th><th>ra</th><th>dec</th><th>pmra</th><th>pmdec</th></tr></thead>
<thead><tr><th></th><th>deg</th><th>deg</th><th>mas / yr</th><th>mas / yr</th></tr></thead>
<thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>637987125186749568</td><td>142.48301935991023</td><td>21.75771616932985</td><td>-2.5168384683875766</td><td>2.941813096629439</td></tr>
<tr><td>638285195917112960</td><td>142.25452941346344</td><td>22.476168171141378</td><td>2.6627020143457996</td><td>-12.165984395577347</td></tr>
<tr><td>638073505568978688</td><td>142.64528557468074</td><td>22.16693224953078</td><td>18.30674739454163</td><td>-7.950659620550862</td></tr>
<tr><td>638086386175786752</td><td>142.57739430926034</td><td>22.22791951401365</td><td>0.9877856720147953</td><td>-2.584105480335548</td></tr>
<tr><td>638049655615392384</td><td>142.58913564478618</td><td>22.110783166677418</td><td>0.24443878227817095</td><td>-4.941079187010136</td></tr>
<tr><td>638267565075964032</td><td>141.81762228999614</td><td>22.375696125322275</td><td>-3.413174589660796</td><td>1.8838892877285924</td></tr>
<tr><td>638028902333511168</td><td>143.18339801317677</td><td>22.2512465812369</td><td>7.848511762712128</td><td>-21.391145547787154</td></tr>
<tr><td>638085767700610432</td><td>142.9347319464589</td><td>22.46244080823965</td><td>-3.6585960944321476</td><td>-12.486419770278376</td></tr>
<tr><td>638299863230178304</td><td>142.26769745823267</td><td>22.640183776884836</td><td>-1.8168370892218297</td><td>1.0537342990941316</td></tr>
<tr><td>637973067758974208</td><td>142.89551292869012</td><td>21.612824100339875</td><td>-8.645166256559063</td><td>-44.41164172204947</td></tr>
</table>



Finally, we can remove `TOP 10` run the query again.

The result is bigger than our previous queries, so it will take a little longer.


```python

~~~
query5_base = """SELECT
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({point_list}))
"""
~~~
{: .language-python}
```


```python

~~~
query5 = query5_base.format(columns=columns, 
                          point_list=point_list)
print(query5)
~~~
{: .language-python}

~~~
SELECT
source_id, ra, dec, pmra, pmdec
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON(146.275, 19.2619, 135.422, 25.8774, 141.603, 34.3048, 152.817, 27.1361, 146.275, 19.2619))


~~~
{: .output}

```

    


```python

~~~
job = Gaia.launch_job_async(query5)
print(job)
~~~
{: .language-python}

~~~
INFO: Query finished. [astroquery.utils.tap.core]
<Table length=140339>
   name    dtype    unit                              description                            
--------- ------- -------- ------------------------------------------------------------------
source_id   int64          Unique source identifier (unique within a particular Data Release)
       ra float64      deg                                                    Right ascension
      dec float64      deg                                                        Declination
     pmra float64 mas / yr                         Proper motion in right ascension direction
    pmdec float64 mas / yr                             Proper motion in declination direction
Jobid: 1615815886707O
Phase: COMPLETED
[Output truncated]
~~~
{: .output}

```

    


```python

~~~
results = job.get_results()
len(results)
~~~
{: .language-python}

~~~
140339
~~~
{: .output}

```




    



There are more than 100,000 stars in this polygon, but that's a
manageable size to work with.

## Saving results

This is the set of stars we'll work with in the next step.  But since
we have a substantial dataset now, this is a good time to save it.

Storing the data in a file means we can shut down this notebook and
pick up where we left off without running the previous query again.

Astropy `Table` objects provide `write`, which writes the table to disk.


```python

~~~
filename = 'gd1_results.fits'
results.write(filename, overwrite=True)
~~~
{: .language-python}
```

Because the filename ends with `fits`, the table is written in the
[FITS format](https://en.wikipedia.org/wiki/FITS), which preserves the
metadata associated with the table.

If the file already exists, the `overwrite` argument causes it to be
overwritten.

We can use `getsize` to confirm that the file exists and check the size:


```python

~~~
from os.path import getsize

MB = 1024 * 1024
getsize(filename) / MB
~~~
{: .language-python}

~~~
5.36407470703125
~~~
{: .output}

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

* Use the `format` function to compose queries; code written this way
is easier to read and less error-prone.

* Develop queries incrementally: start with something simple, test it,
and add a little bit at a time.

* Once you have a query working, save the data in a local file.  If
you shut down the notebook and come back to it later, you can reload
the file; you don't have to run the query again.
