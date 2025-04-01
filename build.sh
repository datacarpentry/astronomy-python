# copy notebooks from jb

cp ~/AstronomicalData/01_query.ipynb 01-query.ipynb
cp ~/AstronomicalData/02_coords.ipynb 02-coords.ipynb

# convert tags to comments
python convert.py

# convert ipynb to md
jupyter nbconvert --to markdown 01-query.ipynb
jupyter nbconvert --to markdown 02-coords.ipynb


# push it to GutHub
#git add 01-query.md
#git add 02-coords.md
#git commit -m "Automated push from AstronomicalData"
#git push
