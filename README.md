This branch maintains the *development* notes for the project.
These rendered notes can be viewed at [http://jeetsukumaran.github.io/inphest](http://jeetsukumaran.github.io/inphest).

To edit the documentation, clone this repository:

    $ git clone https://github.com/jeetsukumaran/inphest.git

and checkout out the "``gh-pages``" branch:

    $ git checkout -t -b gh-pages origin/gh-pages

The documentation file is "``src/inphest-development-journal.tex``".
After editing the "``src/inphest-development-journal.tex``", run "``make``" to update the "``index.html``" file, and then commit the changes:

    $ vim src/inphest-development-journal.tex
    $ make
    $ git commit -a

The build product, "``index.html``", is a tracked, version-controlled file.
This is necessary to have it properly [served by GitHub](http://jeetsukumaran.github.io/inphest).
