#!/bin/bash

STATIC_NB_DIR=../build/html/notebooks/
mkdir -p $STATIC_NB_DIR

pip install MDAnalysis

# Prepare the notebooks
for notebook in ./notebooks/*.ipynb
do
    name=${notebook%.ipynb}
    echo $notebook $name
    # Convert the notebook to rst readable by sphinx
    # Jupyter can convert directly to rst, yet there is an issue with the
    # titles <https://github.com/ipython/ipython/issues/8674>.
    # Therefore, we do the conversion in two steps.
    jupyter-nbconvert --to markdown --execute $notebook --output ${name}.md
    pandoc -i ${name}.md -o ${name}_.rst

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
        sed -Ee '/^===+/ r header.txt' > ${name}.rst

    # Clean unnecessary files
    rm ${name}.md ${name}_.rst header.txt

    # Make available a direct HTML render of the notebook,
    # and the notebook itself
    jupyter nbconvert --to html $notebook --output $STATIC_NB_DIR/$(basename $name)_.html
    cp $notebook $STATIC_NB_DIR
done

# Build the doc
#make html
