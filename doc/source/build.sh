#!/bin/bash

set -e

STATIC_NB_DIR=../build/html/notebooks/

if [[ -n $READTHEDOCS ]]
then
    echo "Settings for Readthedocs"
    # This line makes ghostscript read the Fontmap.GS present in this
    # directory.  Without it, weblogo fails on Readthedocs because it cannot
    # find the ArialMT font.
    export GS_LIB=$PWD
    STATIC_NB_DIR=./_build/html/notebooks/
fi

mkdir -p $STATIC_NB_DIR

pandoc --version

pip install MDAnalysis

# Prepare the notebooks
echo "Jupyter-nbconvert version: `jupyter-nbconvert --version`"
for notebook in ./notebooks/*.ipynb
do
    name=${notebook%.ipynb}
    echo $notebook $name
    # Convert the notebook to rst readable by sphinx
    # Jupyter can convert directly to rst, yet there is an issue with the
    # titles <https://github.com/ipython/ipython/issues/8674>.
    # Therefore, we do the conversion in two steps.
    jupyter-nbconvert --to markdown --execute $notebook --output ${name}.md \
        --ExecutePreprocessor.timeout=180 \
        --ExecutePreprocessor.allow_errors=True
    # Because of the 'implicit_figures' extension, pandoc uses the alt text
    # of images as figure caption. Yet, IPython sets the alt text as 'png'
    # for all figures. Therefore, we disable the extension. The filter
    # changes the alt text to a different one for each figure because, without
    # the 'implicit_figures', pandoc writes the images with the alt text as
    # replacement text, and replacement text cannot repeat themselves.
    # The filter also fix the image path because the origin directory for the
    # relative path is not the same between where pandoc runs and where
    # ipython runs.
    pandoc --from markdown-implicit_figures \
           --filter ./pandoc_fix_img.py \
           -i ${name}.md -o ${name}_.rst

    # Clean and adapt the rst file to have better rendering by sphinx
    cat >> header.txt << EOF

.. contents:: Table of Contents
   :local:
   :backlinks: none

.. note::

   This page is initialy a jupyter notebook. You can see a \`notebook HTML
   render <$(basename $name)_.html>\`__ of it or download the \`notebook
   itself <$(basename $name).ipynb>\`__.

EOF
    sed -Ee 's/:([^:]+):``([^`]+)``/:\1:`\2`/g' ${name}_.rst | \
        sed -Ee '/^===+/ r header.txt' | \
        sed -Ee 's@\.\. figure:: \./notebooks/@.. figure:: @g' > ${name}.rst

    # Clean unnecessary files
    rm ${name}.md ${name}_.rst header.txt

    # Move the notebook images to the build directory
    if [[ -d ${name}_files ]]
    then
        rm -rf $STATIC_NB_DIR/${name}_files
        cp -r ${name}_files $STATIC_NB_DIR/
    fi

    # Make available a direct HTML render of the notebook,
    # and the notebook itself
    jupyter nbconvert --to html $notebook --output $STATIC_NB_DIR/$(basename $name)_.html
    cp $notebook $STATIC_NB_DIR
done

# Build the doc
#make html
