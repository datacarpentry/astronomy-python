---
title: "Photometry"
teaching: 3000
exercises: 0
questions:

- "How do we use Matplotlib to define a polygon and select points that fall inside it?"

objectives:

- "Use isochrone data to specify a polygon and determine which points fall inside it."

- "Use Matplotlib features to customize the appearance of figures."

keypoints:

- "Matplotlib provides operations for working with points, polygons, and other geometric entities, so it's not just for making figures."

- "Use Matplotlib options to control the size and aspect ratio of figures to make them easier to interpret."

- "Record every element of the data analysis pipeline that would be needed to replicate the results."

---

{% include links.md %}

# 6. Photometry

This is the sixth in a series of notebooks related to astronomy data.

As a continuing example, we will replicate part of the analysis in a
recent paper, "[Off the beaten path: Gaia reveals GD-1 stars outside
of the main stream](https://arxiv.org/abs/1805.00425)" by Adrian M.
Price-Whelan and Ana Bonaca.

In the previous lesson we downloaded photometry data from Pan-STARRS,
which is available from the same server we've been using to get Gaia
data.

The next step in the analysis is to select candidate stars based on
the photometry data.
The following figure from the paper is a color-magnitude diagram
showing the stars we previously selected based on proper motion:

<img width="300"
src="https://github.com/datacarpentry/astronomy-python/raw/gh-pages/fig/gd1-3.png">

In red is a theoretical isochrone, showing where we expect the stars
in GD-1 to fall based on the metallicity and age of their original
globular cluster.

By selecting stars in the shaded area, we can further distinguish the
main sequence of GD-1 from mostly younger background stars.

## Outline

Here are the steps in this notebook:

1. We'll reload the data from the previous notebook and make a
color-magnitude diagram.

2. We'll use an isochrone computed by MIST to specify a polygonal
region in the color-magnitude diagram and select the stars inside it.

After completing this lesson, you should be able to

* Use Matplotlib to specify a `Polygon` and determine which points
fall inside it.

* Use Pandas to merge data from multiple `DataFrames`, much like a
database `JOIN` operation.

## Reload the data

The following cell downloads the photometry data we created in the
previous notebook.


```python

~~~
import os
from wget import download

filename = 'gd1_data.hdf'
filepath = 'https://github.com/AllenDowney/AstronomicalData/raw/main/data/'

if not os.path.exists(filename):
    print(download(filepath+filename))
~~~
{: .language-python}
```


```python

~~~
import pandas as pd

candidate_df = pd.read_hdf(filename, 'candidate_df')
~~~
{: .language-python}
```

## Plotting photometry data

Now that we have photometry data from Pan-STARRS, we can replicate the
[color-magnitude
diagram](https://en.wikipedia.org/wiki/Galaxy_color%E2%80%93magnitude_diagram)
from the original paper:

<img width="300"
src="https://github.com/datacarpentry/astronomy-python/raw/gh-pages/fig/gd1-3.png">

The y-axis shows the apparent magnitude of each source with the [g
filter](https://en.wikipedia.org/wiki/Photometric_system).

The x-axis shows the difference in apparent magnitude between the g
and i filters, which indicates color.

Stars with lower values of (g-i) are brighter in g-band than in
i-band, compared to other stars, which means they are bluer.

Stars in the lower-left quadrant of this diagram are less bright than
the others, and have lower metallicity, which means they are [likely
to be
older](http://spiff.rit.edu/classes/ladder/lectures/ordinary_stars/ordinary.html).

Since we expect the stars in GD-1 to be older than the background
stars, the stars in the lower-left are more likely to be in GD-1.

The following function takes a table containing photometry data and
draws a color-magnitude diagram.
The input can be an Astropy `Table` or Pandas `DataFrame`, as long as
it has columns named `g_mean_psf_mag` and `i_mean_psf_mag`.




```python

~~~
import matplotlib.pyplot as plt

def plot_cmd(table):
    """Plot a color magnitude diagram.
    
    table: Table or DataFrame with photometry data
    """
    y = table['g_mean_psf_mag']
    x = table['g_mean_psf_mag'] - table['i_mean_psf_mag']

    plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)

    plt.xlim([0, 1.5])
    plt.ylim([14, 22])
    plt.gca().invert_yaxis()

    plt.ylabel('$Magnitude (g)$')
    plt.xlabel('$Color (g-i)$')
~~~
{: .language-python}
```

`plot_cmd` uses a new function, `invert_yaxis`, to invert the `y`
axis, which is conventional when plotting magnitudes, since lower
magnitude indicates higher brightness.

`invert_yaxis` is a little different from the other functions we've
used.  You can't call it like this:

```
plt.invert_yaxis()          # doesn't work
```

You have to call it like this:

```
plt.gca().invert_yaxis()          # works
```

`gca` stands for "get current axis".  It returns an object that
represents the axes of the current figure, and that object provides
`invert_yaxis`.

**In case anyone asks:** The most likely reason for this inconsistency
in the interface is that `invert_yaxis` is a lesser-used function, so
it's not made available at the top level of the interface.

Here's what the results look like.


```python

~~~
plot_cmd(candidate_df)
~~~
{: .language-python}

~~~
<Figure size 432x288 with 1 Axes>
~~~
{: .output}

```


    
![png](06-photo_files/06-photo_11_0.png)
    


Our figure does not look exactly like the one in the paper because we
are working with a smaller region of the sky, so we don't have as many
stars.  But we can see an overdense region in the lower left that
contains stars with the photometry we expect for GD-1.

In the next section we'll use an isochrone to specify a polygon that
contains this overdense regioin.

## Isochrone

Based on our best estimates for the ages of the stars in GD-1 and
their metallicity, we can compute a [stellar
isochrone](https://en.wikipedia.org/wiki/Stellar_isochrone) that
predicts the relationship between their magnitude and color.

In fact, we can use [MESA Isochrones & Stellar
Tracks](http://waps.cfa.harvard.edu/MIST/) (MIST) to compute it for
us.

Using the [MIST Version 1.2 web
interface](http://waps.cfa.harvard.edu/MIST/interp_isos.html), we
computed an isochrone with the following parameters:
    
* Rotation initial v/v_crit = 0.4

* Single age, linear scale = 12e9

* Composition [Fe/H] = -1.35

* Synthetic Photometry, PanStarrs

* Extinction av = 0

The following cell downloads the results:


```python

~~~
import os
from wget import download

filename = 'MIST_iso_5fd2532653c27.iso.cmd'
filepath = 'https://github.com/AllenDowney/AstronomicalData/raw/main/data/'

if not os.path.exists(filename):
    print(download(filepath+filename))
~~~
{: .language-python}
```

To read this file we'll download a Python module [from this
repository](https://github.com/jieunchoi/MIST_codes).


```python

~~~
import os
from wget import download

filename = 'read_mist_models.py'
filepath = 'https://github.com/jieunchoi/MIST_codes/raw/master/scripts/'

if not os.path.exists(filename):
    print(download(filepath+filename))
~~~
{: .language-python}
```

Now we can read the file:


```python

~~~
import read_mist_models

filename = 'MIST_iso_5fd2532653c27.iso.cmd'
iso = read_mist_models.ISOCMD(filename)
~~~
{: .language-python}

~~~
Reading in: MIST_iso_5fd2532653c27.iso.cmd

~~~
{: .output}

```

    

The result is an `ISOCMD` object.


```python

~~~
type(iso)
~~~
{: .language-python}

~~~
read_mist_models.ISOCMD
~~~
{: .output}

```




    



It contains a list of arrays, one for each isochrone.


```python

~~~
type(iso.isocmds)
~~~
{: .language-python}

~~~
list
~~~
{: .output}

```




    



We only got one isochrone.


```python

~~~
len(iso.isocmds)
~~~
{: .language-python}

~~~
1
~~~
{: .output}

```




    



So we can select it like this:


```python

~~~
iso_array = iso.isocmds[0]
~~~
{: .language-python}
```

It's a NumPy array:


```python

~~~
type(iso_array)
~~~
{: .language-python}

~~~
numpy.ndarray
~~~
{: .output}

```




    



But it's an unusual NumPy array, because it contains names for the columns.


```python

~~~
iso_array.dtype
~~~
{: .language-python}

~~~
dtype([('EEP', '<i4'), ('isochrone_age_yr', '<f8'), ('initial_mass', '<f8'), ('star_mass', '<f8'), ('log_Teff', '<f8'), ('log_g', '<f8'), ('log_L', '<f8'), ('[Fe/H]_init', '<f8'), ('[Fe/H]', '<f8'), ('PS_g', '<f8'), ('PS_r', '<f8'), ('PS_i', '<f8'), ('PS_z', '<f8'), ('PS_y', '<f8'), ('PS_w', '<f8'), ('PS_open', '<f8'), ('phase', '<f8')])
~~~
{: .output}

```




    



Which means we can select columns using the bracket operator:


```python

~~~
iso_array['phase']
~~~
{: .language-python}

~~~
array([0., 0., 0., ..., 6., 6., 6.])
~~~
{: .output}

```




    



We can use `phase` to select the part of the isochrone for stars in
the main sequence and red giant phases.


```python

~~~
phase_mask = (iso_array['phase'] >= 0) & (iso_array['phase'] < 3)
phase_mask.sum()
~~~
{: .language-python}

~~~
354
~~~
{: .output}

```




    




```python

~~~
main_sequence = iso_array[phase_mask]
len(main_sequence)
~~~
{: .language-python}

~~~
354
~~~
{: .output}

```




    



The other two columns we'll use are `PS_g` and `PS_i`, which contain
simulated photometry data for stars with the given age and
metallicity, based on a model of the Pan-STARRS sensors.

We'll use these columns to superimpose the isochrone on the
color-magnitude diagram, but first we have to use a [distance
modulus](https://en.wikipedia.org/wiki/Distance_modulus) to scale the
isochrone based on the estimated distance of GD-1.

We can use the `Distance` object from Astropy to compute the distance modulus.


```python

~~~
import astropy.coordinates as coord
import astropy.units as u

distance = 7.8 * u.kpc
distmod = coord.Distance(distance).distmod.value
distmod
~~~
{: .language-python}

~~~
14.4604730134524
~~~
{: .output}

```




    



Now we can compute the scaled magnitude and color of the isochrone.


```python

~~~
mag_g = main_sequence['PS_g'] + distmod
color_g_i = main_sequence['PS_g'] - main_sequence['PS_i']
~~~
{: .language-python}
```

Now we can plot it on the color-magnitude diagram like this.


```python

~~~
plot_cmd(candidate_df)
plt.plot(color_g_i, mag_g);
~~~
{: .language-python}

~~~
<Figure size 432x288 with 1 Axes>
~~~
{: .output}

```


    
![png](06-photo_files/06-photo_41_0.png)
    


The theoretical isochrone passes through the overdense region where we
expect to find stars in GD-1.

Let's save this result so we can reload it later without repeating the
steps in this section.

So we can save the data in an HDF5 file, we'll put it in a Pandas
`DataFrame` first:


```python

~~~
import pandas as pd

iso_df = pd.DataFrame()
iso_df['mag_g'] = mag_g
iso_df['color_g_i'] = color_g_i

iso_df.head()
~~~
{: .language-python}

~~~
       mag_g  color_g_i
0  28.294743   2.195021
1  28.189718   2.166076
2  28.051761   2.129312
3  27.916194   2.093721
4  27.780024   2.058585
~~~
{: .output}

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>mag_g</th>
      <th>color_g_i</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>28.294743</td>
      <td>2.195021</td>
    </tr>
    <tr>
      <th>1</th>
      <td>28.189718</td>
      <td>2.166076</td>
    </tr>
    <tr>
      <th>2</th>
      <td>28.051761</td>
      <td>2.129312</td>
    </tr>
    <tr>
      <th>3</th>
      <td>27.916194</td>
      <td>2.093721</td>
    </tr>
    <tr>
      <th>4</th>
      <td>27.780024</td>
      <td>2.058585</td>
    </tr>
  </tbody>
</table>
</div>



And then save it.


```python

~~~
filename = 'gd1_isochrone.hdf5'
iso_df.to_hdf(filename, 'iso_df')
~~~
{: .language-python}
```

## Making a polygon

The following cell downloads the isochrone we made in the previous
section, if necessary.


```python

~~~
import os
from wget import download

filename = 'gd1_isochrone.hdf5'
filepath = 'https://github.com/AllenDowney/AstronomicalData/raw/main/data/'

if not os.path.exists(filename):
    print(download(filepath+filename))
~~~
{: .language-python}
```

Now we can read it back in.


```python

~~~
iso_df = pd.read_hdf(filename, 'iso_df')
iso_df.head()
~~~
{: .language-python}

~~~
       mag_g  color_g_i
0  28.294743   2.195021
1  28.189718   2.166076
2  28.051761   2.129312
3  27.916194   2.093721
4  27.780024   2.058585
~~~
{: .output}

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>mag_g</th>
      <th>color_g_i</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>28.294743</td>
      <td>2.195021</td>
    </tr>
    <tr>
      <th>1</th>
      <td>28.189718</td>
      <td>2.166076</td>
    </tr>
    <tr>
      <th>2</th>
      <td>28.051761</td>
      <td>2.129312</td>
    </tr>
    <tr>
      <th>3</th>
      <td>27.916194</td>
      <td>2.093721</td>
    </tr>
    <tr>
      <th>4</th>
      <td>27.780024</td>
      <td>2.058585</td>
    </tr>
  </tbody>
</table>
</div>



Here's what the isochrone looks like on the color-magnitude diagram.


```python

~~~
plot_cmd(candidate_df)
plt.plot(iso_df['color_g_i'], iso_df['mag_g']);
~~~
{: .language-python}

~~~
<Figure size 432x288 with 1 Axes>
~~~
{: .output}

```


    
![png](06-photo_files/06-photo_51_0.png)
    


In the bottom half of the figure, the isochrone passes through the
overdense region where the stars are likely to belong to GD-1.

In the top half, the isochrone passes through other regions where the
stars have higher magnitude and metallicity than we expect for stars
in GD-1.

So we'll select the part of the isochrone that lies in the overdense region.

`g_mask` is a Boolean Series that is `True` where `g` is between 18.0 and 21.5.


```python

~~~
g = iso_df['mag_g']

g_mask = (g > 18.0) & (g < 21.5)
g_mask.sum()
~~~
{: .language-python}

~~~
117
~~~
{: .output}

```




    



We can use it to select the corresponding rows in `iso_df`:


```python

~~~
iso_masked = iso_df[g_mask]
iso_masked.head()
~~~
{: .language-python}

~~~
        mag_g  color_g_i
94  21.411746   0.692171
95  21.322466   0.670238
96  21.233380   0.648449
97  21.144427   0.626924
98  21.054549   0.605461
~~~
{: .output}

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>mag_g</th>
      <th>color_g_i</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>94</th>
      <td>21.411746</td>
      <td>0.692171</td>
    </tr>
    <tr>
      <th>95</th>
      <td>21.322466</td>
      <td>0.670238</td>
    </tr>
    <tr>
      <th>96</th>
      <td>21.233380</td>
      <td>0.648449</td>
    </tr>
    <tr>
      <th>97</th>
      <td>21.144427</td>
      <td>0.626924</td>
    </tr>
    <tr>
      <th>98</th>
      <td>21.054549</td>
      <td>0.605461</td>
    </tr>
  </tbody>
</table>
</div>



Now, to select the stars in the overdense region, we have to define a
polygon that includes stars near the isochrone.

The original paper uses the following formulas to define the left and
right boundaries.


```python

~~~
g = iso_masked['mag_g']
left_color = iso_masked['color_g_i'] - 0.4 * (g/28)**5
right_color = iso_masked['color_g_i'] + 0.8 * (g/28)**5
~~~
{: .language-python}
```

The intention is to define a polygon that gets wider as `g` increases,
to reflect increasing uncertainty.

But we can do about as well with a simpler formula:


```python

~~~
g = iso_masked['mag_g']
left_color = iso_masked['color_g_i'] - 0.06
right_color = iso_masked['color_g_i'] + 0.12
~~~
{: .language-python}
```

Here's what these boundaries look like:


```python

~~~
plot_cmd(candidate_df)

plt.plot(left_color, g, label='left color')
plt.plot(right_color, g, label='right color')

plt.legend();
~~~
{: .language-python}

~~~
<Figure size 432x288 with 1 Axes>
~~~
{: .output}

```


    
![png](06-photo_files/06-photo_61_0.png)
    


## Which points are in the polygon?

Matplotlib provides a `Polygon` object that we can use to check which
points fall in the polygon we just constructed.

To make a `Polygon`, we need to assemble `g`, `left_color`, and
`right_color` into a loop, so the points in `left_color` are connected
to the points of `right_color` in reverse order.

We'll use the following function, which takes two arrays and joins
them front-to-back:


```python

~~~
import numpy as np

def front_to_back(first, second):
    """Join two arrays front to back."""
    return np.append(first, second[::-1])
~~~
{: .language-python}
```

`front_to_back` uses a "slice index" to reverse the elements of `second`.

As explained in the [NumPy
documentation](https://numpy.org/doc/stable/reference/arrays.indexing.html),
a slice index has three parts separated by colons:

* `start`: The index of the element where the slice starts.

* `stop`: The index of the element where the slice ends.

* `step`: The step size between elements.

In this example, `start` and `stop` are omitted, which means all
elements are selected.

And `step` is `-1`, which means the elements are in reverse order.

We can use `front_to_back` to make a loop that includes the elements
of `left_color` and `right_color`:


```python

~~~
color_loop = front_to_back(left_color, right_color)
color_loop.shape
~~~
{: .language-python}

~~~
(234,)
~~~
{: .output}

```




    



And a corresponding loop with the elements of `g` in forward and reverse order.


```python

~~~
mag_loop = front_to_back(g, g)
mag_loop.shape
~~~
{: .language-python}

~~~
(234,)
~~~
{: .output}

```




    



Here's what the loop looks like.


```python

~~~
plot_cmd(candidate_df)
plt.plot(color_loop, mag_loop);
~~~
{: .language-python}

~~~
<Figure size 432x288 with 1 Axes>
~~~
{: .output}

```


    
![png](06-photo_files/06-photo_69_0.png)
    


To make a `Polygon`, it will be convenient to put `color_loop` and
`mag_loop` into a `DataFrame`:


```python

~~~
loop_df = pd.DataFrame()
loop_df['color_loop'] = color_loop
loop_df['mag_loop'] = mag_loop
loop_df.head()
~~~
{: .language-python}

~~~
   color_loop   mag_loop
0    0.632171  21.411746
1    0.610238  21.322466
2    0.588449  21.233380
3    0.566924  21.144427
4    0.545461  21.054549
~~~
{: .output}

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>color_loop</th>
      <th>mag_loop</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.632171</td>
      <td>21.411746</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.610238</td>
      <td>21.322466</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.588449</td>
      <td>21.233380</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.566924</td>
      <td>21.144427</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.545461</td>
      <td>21.054549</td>
    </tr>
  </tbody>
</table>
</div>



Now we can pass `loop_df` to `Polygon`:


```python

~~~
from matplotlib.patches import Polygon

polygon = Polygon(loop_df)
polygon
~~~
{: .language-python}

~~~
<matplotlib.patches.Polygon at 0x7fd8a96084f0>
~~~
{: .output}

```




    



The result is a `Polygon` object , which provides `contains_points`,
which figures out which points are inside the polygon.

To test it, we'll create a list with two points, one inside the
polygon and one outside.


```python

~~~
points = [(0.4, 20), 
          (0.4, 16)]
~~~
{: .language-python}
```

Now we can make sure `contains_points` does what we expect.


```python

~~~
inside = polygon.contains_points(points)
inside
~~~
{: .language-python}

~~~
array([ True, False])
~~~
{: .output}

```




    



The result is an array of Boolean values.

We are almost ready to select stars whose photometry data falls in
this polygon.  But first we need to do some data cleaning.

## Save the polygon

[Reproducibile
research](https://en.wikipedia.org/wiki/Reproducibility#Reproducible_research)
is "the idea that ... the full computational environment used to
produce the results in the paper such as the code, data, etc. can be
used to reproduce the results and create new work based on the
research."

This Jupyter notebook is an example of reproducible research because
it contains all of the code needed to reproduce the results, including
the database queries that download the data and and analysis.

In this lesson we used an isochrone to derive a polygon, which we used
to select stars based on photometry.
So it is important to record the polygon as part of the data analysis pipeline.

Here's how we can save it in an HDF file.


```python

~~~
filename = 'gd1_data.hdf'
loop_df.to_hdf(filename, 'loop_df')
~~~
{: .language-python}
```

## Selecting based on photometry

Now let's see how many of the candidate stars are inside the polygon we chose.
We'll put color and magnitude data from `candidate_df` into a new `DataFrame`:


```python

~~~
points = pd.DataFrame()

points['color'] = candidate_df['g_mean_psf_mag'] - candidate_df['i_mean_psf_mag']
points['mag'] = candidate_df['g_mean_psf_mag']

points.head()
~~~
{: .language-python}

~~~
    color      mag
0  0.3804  17.8978
1  1.6092  19.2873
2  0.4457  16.9238
3  1.5902  19.9242
4  1.4853  16.1516
~~~
{: .output}

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>color</th>
      <th>mag</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.3804</td>
      <td>17.8978</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1.6092</td>
      <td>19.2873</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.4457</td>
      <td>16.9238</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1.5902</td>
      <td>19.9242</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1.4853</td>
      <td>16.1516</td>
    </tr>
  </tbody>
</table>
</div>



Which we can pass to `contains_points`:


```python

~~~
inside = polygon.contains_points(points)
inside
~~~
{: .language-python}

~~~
array([False, False, False, ..., False, False, False])
~~~
{: .output}

```




    



The result is a Boolean array.  We can use `sum` to see how many stars
fall in the polygon.


```python

~~~
inside.sum()
~~~
{: .language-python}

~~~
454
~~~
{: .output}

```




    



Now we can use `inside` as a mask to select stars that fall inside the polygon.


```python

~~~
winner_df = candidate_df[inside]
~~~
{: .language-python}
```

Let's make a color-magnitude plot one more time, highlighting the
selected stars with green markers.


```python

~~~
plot_cmd(candidate_df)
plt.plot(color_g_i, mag_g)
plt.plot(color_loop, mag_loop)

x = winner_df['g_mean_psf_mag'] - winner_df['i_mean_psf_mag']
y = winner_df['g_mean_psf_mag']
plt.plot(x, y, 'go', markersize=0.5, alpha=0.5);
~~~
{: .language-python}

~~~
<Figure size 432x288 with 1 Axes>
~~~
{: .output}

```


    
![png](06-photo_files/06-photo_90_0.png)
    


It looks like the selected stars are, in fact, inside the polygon,
which means they have photometry data consistent with GD-1.

Finally, we can plot the coordinates of the selected stars:


```python

~~~
plt.figure(figsize=(10,2.5))

x = winner_df['phi1']
y = winner_df['phi2']
plt.plot(x, y, 'ko', markersize=0.7, alpha=0.9)

plt.xlabel('ra (degree GD1)')
plt.ylabel('dec (degree GD1)')

plt.axis('equal');
~~~
{: .language-python}

~~~
<Figure size 720x180 with 1 Axes>
~~~
{: .output}

```


    
![png](06-photo_files/06-photo_92_0.png)
    


This example includes two new Matplotlib commands:

* `figure` creates the figure.  In previous examples, we didn't have
to use this function; the figure was created automatically.  But when
we call it explicitly, we can provide arguments like `figsize`, which
sets the size of the figure.

* `axis` with the parameter `equal` sets up the axes so a unit is the
same size along the `x` and `y` axes.

In an example like this, where `x` and `y` represent coordinates in
space, equal axes ensures that the distance between points is
represented accurately.

## Write the data

Finally, let's write the selected stars to a file.


```python

~~~
filename = 'gd1_data.hdf'
winner_df.to_hdf(filename, 'winner_df')
~~~
{: .language-python}
```


```python

~~~
from os.path import getsize

MB = 1024 * 1024
getsize(filename) / MB
~~~
{: .language-python}

~~~
2.512819290161133
~~~
{: .output}

```




    



## Summary

In this lesson, we used photometry data from Pan-STARRS to draw a
color-magnitude diagram.
We used an isochrone to define a polygon and select stars we think are
likely to be in GD-1.  Plotting the results, we have a clearer picture
of GD-1, similar to Figure 1 in the original paper.

## Best practices

* Matplotlib provides operations for working with points, polygons,
and other geometric entities, so it's not just for making figures.

* Use Matplotlib options to control the size and aspect ratio of
figures to make them easier to interpret.  In this example, we scaled
the axes so the size of a degree is equal along both axes.

* Record every element of the data analysis pipeline that would be
needed to replicate the results.


```python

~~~

~~~
{: .language-python}
```
