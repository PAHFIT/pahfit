"""
pahfit.features

Manage PAHFIT features, their parameters, and parameter attributes.
The PAHFIT model contains a flexible combination of features, each of
which have associated parameters and parameter attributes.  These
parameters are specified by a combination of a named science pack, and
instrument-specific information such as wavelength-dependent line
resolution.

The main class, Features, inherits from astropy.table.Table.  All the
grouping, sorting, selection, and indexing operations from astropy
tables are therefore also available for pahfit.features.Features.

  Usage:

  TBD
"""

import os
import numpy as np
from astropy.table import vstack, Table, TableAttribute
from astropy.io.misc.yaml import yaml
from importlib import resources
from pahfit.errors import PAHFITFeatureError
from pahfit.features.features_format import BoundedParTableFormatter
import pahfit.units

# Feature kinds and associated parameters
KIND_PARAMS = {'starlight': {'temperature', 'tau'},
               'dust_continuum': {'temperature', 'tau'},
               'line': {'wavelength', 'power'},  # 'fwhm', Instrument Pack detail!
               'dust_feature': {'wavelength', 'fwhm', 'power'},
               'attenuation': {'model', 'tau', 'geometry'},
               'absorption': {'wavelength', 'fwhm', 'tau', 'geometry'}}

# Parameter default units: flux density/intensity/power (other units determined on fit)
# Note: power is actually amplitude for now. Change this to
# intensity_power when the power features are implemented.

PARAM_UNITS = {'temperature': pahfit.units.temperature,
               'wavelength': pahfit.units.wavelength,
               'fwhm': pahfit.units.wavelength,
               'power': pahfit.units.intensity_power}


class UniqueKeyLoader(yaml.SafeLoader):
    def construct_mapping(self, node, deep=False):
        mapping = set()
        for key_node, _ in node.value:
            key = self.construct_object(key_node, deep=deep)
            if key in mapping:
                raise PAHFITFeatureError(f"Duplicate {key!r} key found in YAML.")
            mapping.add(key)
        return super().construct_mapping(node, deep)


def value_bounds(val, bounds):
    """Compute bounds for a bounded value.

    Arguments:
    ----------
      val: The value to bound.

      bounds: Either None for no relevant bounds (i.e. fixed), or a
        two element iterable specifying (min, max) bounds.  Each of
        min and max can be a numerical value, None (for infinite
        bounds, either negative or positive, as appropriate), or a
        string ending in:

          #: an absolute offset from the value
          %: a percentage offset from the value

        Offsets are necessarily negative for min bound, positive for
        max bounds.

    Examples:
    ---------

      A bound of ('-1.5%', '0%') would indicate a minimum bound
        1.5% below the value, and a max bound at the value itself.

      A bound of ('-0.1#', None) would indicate a minimum bound 0.1 below
        the value, and infinite maximum bound.

    Returns:
    -------

      A 3 element tuple (value, min, max).
        Any missing bound is replaced with the numpy.nan value.

    Raises:
    -------

      ValueError: if bounds are specified and the value does not fall
        between them.

    """
    if val is None:
        val = np.ma.masked
    if not bounds:
        return (val,) + 2 * (np.nan,)  # (val,nan,nan) indicates fixed
    ret = [val]
    for i, b in enumerate(bounds):
        if isinstance(b, str):
            if b.endswith('%'):
                b = val * (1. + float(b[:-1]) / 100.)
            elif b.endswith('#'):
                b = val + float(b[:-1])
            else:
                raise PAHFITFeatureError(f"Incorrectly formatted bound {b}")
        elif b is None:
            b = np.inf if i else -np.inf
        ret.append(b)

    if (val < ret[1] or val > ret[2]):
        raise PAHFITFeatureError(f"Value <{ret[0]}> is not between bounds: {ret[1:]}")
    return tuple(ret)


class Features(Table):
    """A class for holding a table of PAHFIT features and associated
    parameter information.

    Note that each parameter has an associated `kind', and that each
    kind has an associated set of allowable parameters (see
    `KIND_PARAMS`).

    See Also
    --------
    `~astropy.table.Table`: The parent table class.
    """

    TableFormatter = BoundedParTableFormatter

    param_covar = TableAttribute(default=[])
    _group_attrs = set(('bounds', 'features', 'kind'))  # group-level attributes
    _param_attrs = set(('value', 'bounds'))  # Each parameter can have these attributes
    _no_bounds = set(('name', 'group', 'kind', 'geometry', 'model'))  # str attributes (no bounds)
    _bounds_dtype = np.dtype([("val", float), ("min", float), ("max", float)])  # bounded param type

    @classmethod
    def read(cls, file, *args, **kwargs):
        """Read a table from file.

        If reading a YAML file, read it in as a science pack and
        return the new table. Otherwise, use astropy's normal Table
        reader.
        """
        if file.endswith(".yaml") or file.endswith(".yml"):
            return cls._read_scipack(file)
        else:
            table = super().read(file, *args, **kwargs)
            table.add_index('name')
            return table

    @classmethod
    def _read_scipack(cls, file):
        """Read a science pack specification from YAML file.

        Arguments:
        ----------

          file: the name of the file, either a full valid path, or
            named file in the PAHFIT science_packs directory.!

        Returns:
        --------

          table: A filled pahfit.features.Features table.
        """

        feat_tables = dict()

        if not os.path.isfile(file):
            file = resources.files("pahfit") / "packs/science" / file
        try:
            with open(file) as fd:
                scipack = yaml.load(fd, Loader=UniqueKeyLoader)
        except IOError as e:
            raise PAHFITFeatureError("Error reading science pack file\n"
                                     f"\t{file}\n\t{repr(e)}")
        for (name, elem) in scipack.items():
            try:
                keys = elem.keys()
            except AttributeError:
                raise PAHFITFeatureError("Invalid science pack"
                                         f" format at {name}\n\t{file}")

            try:
                kind = elem.pop('kind')
            except KeyError:
                raise PAHFITFeatureError(f"No kind found for {name}\n\t{file}")

            try:
                valid_params = KIND_PARAMS[kind]
            except KeyError:
                raise PAHFITFeatureError(f"Unknown kind {kind} for {name}\n\t{file}")
            unknown_params = [x for x in keys
                              if not (x in valid_params or x in cls._group_attrs)]
            if unknown_params:
                raise PAHFITFeatureError(f"Unknown {kind} parameters:"
                                         f" {', '.join(unknown_params)}\n\t{file}")

            hasFeatures = 'features' in elem
            hasLists = any(k not in cls._group_attrs
                           and (isinstance(v, (tuple, list))
                                or (isinstance(v, dict)
                                    and cls._param_attrs.isdisjoint(v.keys())))
                           for (k, v) in elem.items())
            if hasFeatures and hasLists:
                raise PAHFITFeatureError("A single group cannot contain both 'features'"
                                         f" and parameter list(s): {name}\n\t{file}")
            isGroup = (hasFeatures or hasLists)
            bounds = None
            if isGroup:  # A named group of features
                if 'bounds' in elem:
                    if not isinstance(elem['bounds'], dict):
                        for p in cls._no_bounds:
                            if p in elem:
                                raise PAHFITFeatureError(f"Parameter {p} cannot have "
                                                         f"bounds: {name}\n\t{file}")
                        if sum(p in elem for p in valid_params) > 1:
                            raise PAHFITFeatureError("Groups with simple bounds "
                                                     "can only specify a single "
                                                     f"parameter: {name}\n\t{file}")
                        if hasFeatures:
                            raise PAHFITFeatureError("Groups with simple bounds "
                                                     "cannot specify "
                                                     f"'features': {name}\n\t{file}")
                    bounds = elem.pop('bounds')
                if hasFeatures:  # our group uses a features dict
                    for n, v in elem['features'].items():
                        if bounds and 'bounds' not in v:  # inherit bounds
                            v['bounds'] = bounds
                        cls._add_feature(kind, feat_tables, n, group=name, **v)
                elif hasLists:  # a "shortcut" feature group, using lists
                    llen = []
                    for k, v in elem.items():
                        if k in cls._group_attrs:
                            continue
                        if not isinstance(v, (tuple, list, dict)):
                            raise PAHFITFeatureError(f"All non-group parameters in {name} "
                                                     f"must be lists or dicts:\n\t{file}")
                        llen.append(len(v))

                    if not all(x == llen[0] for x in llen):
                        raise PAHFITFeatureError(f"All parameter lists in group {name} "
                                                 f"must be the same length:\n\t{file}")
                    ngroup = llen[0]
                    feat_names = None
                    for k, v in elem.items():
                        if isinstance(elem[k], dict):
                            if not feat_names:  # First names win
                                feat_names = list(elem[k].keys())
                            elem[k] = list(elem[k].values())  # turn back into a value list
                    if not feat_names:  # no names: construct one for each group feature
                        feat_names = [f"{name}{x:02}" for x in range(ngroup)]
                    for i in range(ngroup):  # Iterate over list(s) adding feature
                        v = {k: elem[k][i] for k in valid_params if k in elem}
                        cls._add_feature(kind, feat_tables, feat_names[i],
                                         group=name, bounds=bounds, **v)
                else:
                    raise PAHFITFeatureError(f"Group {name} needs either 'features' or"
                                             f"parameter list(s):\n\t{file}")
            else:  # Just one standalone feature
                cls._add_feature(kind, feat_tables, name, **elem)
        return cls._construct_table(feat_tables)

    @classmethod
    def _add_feature(cls, kind: str, t: dict, name: str, *,
                     bounds=None, group='_none_', **pars):
        """Adds an individual feature to the passed dictionary t."""
        if kind not in t:
            t[kind] = {}  # group by kind
        if name not in t[kind]:
            t[kind][name] = {}
        t[kind][name]['group'] = group
        t[kind][name]['kind'] = kind
        for (param, val) in pars.items():
            if param not in KIND_PARAMS[kind]:
                continue
            if isinstance(val, dict):  # A param attribute dictionary
                unknown_attrs = [x for x in val.keys() if x not in cls._param_attrs]
                if unknown_attrs:
                    raise PAHFITFeatureError("Unknown parameter attributes for"
                                             f" {name} ({kind}, {group}): "
                                             f"{', '.join(unknown_attrs)}")
                if 'value' not in val:
                    raise PAHFITFeatureError("Missing 'value' attribute for "
                                             f"{name} ({kind}, {group})")
                value = val['value']
                if 'bounds' in val:  # individual param bounds
                    if param in cls._no_bounds:
                        raise PAHFITFeatureError("Parameter {param} cannot have bounds: "
                                                 f"{name} ({kind}, {group})")
                    bounds = val['bounds']
            else:
                value = val  # a bare value
            if isinstance(bounds, dict):
                b = bounds.get(param)
                if b and param in cls._no_bounds:
                    raise PAHFITFeatureError("Parameter {param} cannot have bounds: "
                                             f"{name} ({kind}, {group})")
            else:  # Simple bounds
                b = bounds
            try:
                t[kind][name][param] = (value if param in cls._no_bounds
                                        else value_bounds(value, b))
            except ValueError as e:
                raise PAHFITFeatureError("Error initializing value and bounds for"
                                         f" {name} ({kind}, {group}):\n\t{e}")

    @classmethod
    def _construct_table(cls, inp: dict):
        """Construct a masked table with units from input dictionary
        INP.  INP is a dictionary with feature names as the key, and a
        dictionary of feature parameters as value.  Each value in the
        feature parameter dictionary is either a value or tuple of 3
        values for bounds.
        """
        tables = []
        for (kind, features) in inp.items():
            if kind == "_ratios":
                continue
            kp = KIND_PARAMS[kind]  # All params for this kind
            rows = []
            for (name, params) in features.items():
                for missing in kp - params.keys():
                    if missing in cls._no_bounds:
                        params[missing] = 0.0
                    else:
                        params[missing] = value_bounds(0.0, bounds=(0.0, None))
                rows.append(dict(name=name, **params))
            param_names = rows[0].keys()
            dtypes = [str if x in cls._no_bounds else cls._bounds_dtype for x in param_names]
            t = cls(rows, names=param_names, dtype=dtypes)
            tables.append(t)
        tables = vstack(tables)
        for cn, col in tables.columns.items():
            if cn in PARAM_UNITS:
                col.unit = PARAM_UNITS[cn]
        cls._index_table(tables)

        if '_ratios' in inp:
            tables.meta['_ratios'] = inp['_ratios']
        return tables

    @staticmethod
    def _index_table(tbl):
        for indx in ('name', 'group'):
            tbl.add_index(indx)

    def mask_feature(self, name, mask_value=True):
        """Mask all the parameters of a feature.

        The masks of the bounds are left intact, so we don't lose the
        fixed / not fixed distinction.

        This is used to indicate that the parameter values of this
        feature were not fit. This mask should not affect the model
        constructor. It is purely a way to indicate to the user that the
        parameter values are meaningless.

        mask_value : bool
            Set this to False to undo the mask

        """
        row = self.loc[name]
        relevant_params = KIND_PARAMS[row['kind']]
        for col_name in relevant_params:
            if col_name in self._no_bounds:
                # these are all strings, so can't mask
                pass
            else:
                # mask only the value, not the bounds
                row[col_name].mask['val'] = mask_value

    def unmask_feature(self, name):
        """Remove the mask for all parameters of a feature."""
        self.mask_feature(name, mask_value=False)

    def _base_repr_(self, *args, **kwargs):
        """Omit dtype on self-print."""
        return super()._base_repr_(*args, ** kwargs | dict(show_dtype=False))
