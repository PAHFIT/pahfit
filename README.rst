PAHFIT
======

``PAHFIT`` is a decomposition model and tool for astronomical infrared spectra, focusing on dust and gas emission features from the interstellar medium (see `Smith, J.D.T., Draine B.T., et al., 2007, ApJ, 656, 770 <http://tir.astro.utoledo.edu/jdsmith/research/pahfit.php>`_).

This package provides an updated python implementation of ``PAHFIT``.  While the original versions of ``PAHFIT`` (``v1.x``) were written in IDL and focused mainly on Spitzer/IRS spectroscopic observations, the newer python-based versions (``>=v2.0``) will expand instrument coverage to other existing (e.g., AKARI) and planned (e.g., JWST) facilities, and will offer a more flexible modeling framework suitable for modeling a wider range of astrophysical sources.

Based on and inspired by the IDL code ``PAHFIT`` by JD Smith and Bruce Draine.

Build checks/status
-------------------

.. image:: http://readthedocs.org/projects/pahfit/badge/?version=latest
   :target: http://dust-extinction.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://github.com/PAHFIT/pahfit/workflows/Python%20Tests/badge.svg
   :target: https://github.com/PAHFIT/pahfit/actions/
   :alt: Test Status

.. image:: https://codecov.io/gh/PAHFIT/pahfit/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/PAHFIT/pahfit
   :alt: Test Coverage Status

.. image:: https://img.shields.io/lgtm/grade/python/g/PAHFIT/pahfit.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/PAHFIT/pahfit/context:python
   :alt: LGTM Status

.. image:: https://app.codacy.com/project/badge/Grade/01a75df3279e45609906c1f28a4ca867
   :target: https://www.codacy.com/gh/PAHFIT/pahfit?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=PAHFIT/pahfit&amp;utm_campaign=Badge_Grade
   :alt: Codacy Status

Packaging
---------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

Documentation
-------------

Hosted by readthedocs: <http://pahfit.readthedocs.io/en/latest/>

In Development!
---------------

This code is currently in active development.
Contributions welcome (see below).

Contributors
------------
* Thomas Lai
* Alexandros Maragkoudakis
* J.D. Smith
* Karl Gordon
* Henrik Spoon
* Els Peeters
* Karin Sandstrom

License
-------

This project is Copyright (c) PAHFit Developers and licensed under
the terms of the GNU GPL v3+ license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.

Contributing
------------

Please open a new issue or new pull request for bugs, feedback, or new features
you would like to see.   If there is an issue you would like to work on, please
leave a comment and we will be happy to assist.   New contributions and
contributors are very welcome!

New to github or open source projects?  If you are unsure about where to start
or haven't used github before, please feel free to contact `@karllark`.
Want more information about how to make a contribution?  Take a look at
the astropy `contributing`_ and `developer`_ documentation.

Feedback and feature requests?   Is there something missing you would like
to see?  Please open an issue or send an email to  `@karllark`.
PAHFIT follows the `Astropy Code of Conduct`_ and strives to provide a
welcoming community to all of our users and contributors.

We love contributions! pahfit is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

*This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
pahfit based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.*

.. _AstroPy: https://www.astropy.org/
.. _contributing: https://docs.astropy.org/en/stable/index.html#contributing
.. _developer: https://docs.astropy.org/en/stable/index.html#developer-documentation
.. _Astropy Code of Conduct:  https://www.astropy.org/about.html#codeofconduct
