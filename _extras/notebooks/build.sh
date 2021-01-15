# copy notebooks from source repo

cp ~/AstronomicalData/01_query.ipynb 01-query.ipynb
cp ~/AstronomicalData/02_coords.ipynb 02-coords.ipynb
cp ~/AstronomicalData/03_motion.ipynb 03-motion.ipynb
cp ~/AstronomicalData/04_select.ipynb 04-select.ipynb
cp ~/AstronomicalData/05_join.ipynb 05-join.ipynb
cp ~/AstronomicalData/06_photo.ipynb 06-photo.ipynb
cp ~/AstronomicalData/07_plot.ipynb 07-plot.ipynb
git add [0-7]*.ipynb
git commit -m "Updating notebooks"

# add Carpentries markup
python convert.py

# convert ipynb to md
jupyter nbconvert --to markdown 01-query.ipynb
jupyter nbconvert --to markdown 02-coords.ipynb
jupyter nbconvert --to markdown 03-motion.ipynb
jupyter nbconvert --to markdown 04-select.ipynb
jupyter nbconvert --to markdown 05-join.ipynb
jupyter nbconvert --to markdown 06-photo.ipynb
jupyter nbconvert --to markdown 07-plot.ipynb
mv [0-7]*.md ../../_episodes
rsync -a [0-7]*files ../../_episodes

# push it to GutHub
cd ../../_episodes
git add 0*.md
git add 0*files
git commit -m "Automated translation to md"
git push
