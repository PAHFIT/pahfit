##############
How to Install
##############

From source
===========

``pahfit`` can be installed from the source code in the normal
python fashion after downloading it from the git repo::

    pip install -e .

Using pip
=========

``pahfit`` can also be installed using pip::

    # from the master trunk on the repository, considered developmental code
    pip install git+https://github.com/PAHFIT/pahfit.git

If you get an error that is includes `SSLError(SSLCertVerificationError`, the
following command may work to remove this error::

    pip install --trusted-host=pypi.org --trusted-host=files.pythonhosted.org git+https://github.com/PAHFIT/pahfit.git
