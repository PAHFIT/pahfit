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
import astropy.units as u
from pkg_resources import resource_filename
from pahfit.errors import PAHFITFeatureError
from pahfit.features.features_format import BoundedMaskedColumn, BoundedParTableFormatter
from pahfit.units import UNITS

# Feature kinds and associated parameters
KIND_PARAMS = {'starlight': {'temperature', 'tau'},
               'dust_continuum': {'temperature', 'tau'},
               'line': {'wavelength', 'power'},  # 'fwhm', Instrument Pack detail!
               'dust_feature': {'wavelength', 'fwhm', 'power'},
               'attenuation': {'model', 'tau', 'geometry'},
               'absorption': {'wavelength', 'fwhm', 'tau', 'geometry'}}

# Parameter default units: flux density/intensity/power (other units determined on fit)
PARAM_UNITS = {'temperature': UNITS.temperature.value,
               'wavelength': UNITS.wavelength.value,
               'fwhm': UNITS.wavelength.value}


class UniqueKeyLoader(yaml.SafeLoader):
    def construct_mapping(self, node, deep=False):
        mapping = set()
        for key_node, _ in node.value:
            key = self.construct_object(key_node, deep=deep)
            if key in mapping:
                raise PAHFITFeatureError(f"Duplicate {key!r} key found in YAML.")
            mapping.add(key)
        return super().construct_mapping(node, deep)


def value_bounds(val, bounds, no_masked=False):
    """Compute bounds for a bounded value.

    Parameters
    ----------
      val : float
          The value to bound.

      bounds : float or str iterable, or None
          Either None for no relevant bounds (i.e. fixed), or a two
          element iterable specifying (min, max) bounds.  Each of min
          and max can be a numerical value, None (for infinite bounds,
          either negative or positive, as appropriate), or a string
          ending in:

            #: an absolute offset from the value
            %: a percentage offset from the value

          Offsets are necessarily negative for min bound, positive for
          max bounds.

    Keywords
    --------

      no_masked: bool, default: False
          If set, the value ``None`` will be use instead of `np.ma.masked`.

    Returns:
    -------

      The value and bounds as a 3 element tuple (value, min, max).
          Any missing bound is replaced with the numpy `masked' value.

    Notes
    --------

      A bound of ('-1.5%', '0%') would indicate a minimum bound
        1.5% below the value, and a max bound at the value itself.

      A bound of ('-0.1#', None) would indicate a minimum bound 0.1 below
        the value, and infinite maximum bound.

    Raises:
    -------

      ValueError: if bounds are specified and the value does not fall
          between them.
    """

    ma_val = None if no_masked else np.ma.masked
    if val is None or not bounds:
        return (val,) + 2 * (ma_val,)  # Fixed
    ret = [val]
    for i, b in enumerate(bounds):
        if isinstance(b, str) and val:
            if b.endswith('%'):
                b = val * (1. + float(b[:-1]) / 100.)
            elif b.endswith('#'):
                b = val + float(b[:-1])
            else:
                raise PAHFITFeatureError(f"Incorrectly formatted bound: {b}")
        elif b is None:
            b = np.inf if i else -np.inf  # lower/upper bound missing
        ret.append(b)

    if (val < ret[1] or val > ret[2]):
        raise ValueError(f"Value <{ret[0]}> is not between bounds: {ret[1:]}")
    return tuple(ret)


class Features(Table):
    """A class for holding a table of PAHFIT features and associated
    parameter information.

    Note that each parameter has an associated `kind', and that each
    kind has an associated set of allowable parameters (see
    `KIND_PARAMS`.

    See Also
    --------
    `~astropy.table.Table`: The parent table class.
    """

    TableFormatter = BoundedParTableFormatter
    MaskedColumn = BoundedMaskedColumn

    param_covar = TableAttribute(default=[])
    _param_attrs = set(('value', 'bounds', 'tied'))  # params can have these attributes
    _group_attrs = set(('bounds', 'features', 'kind', 'tied'))  # group-level attributes
    _no_bounds = set(('name', 'group', 'geometry', 'model'))  # String attributes (no bounds)

    def __repr__(self):
        repr = super().__repr__()
        if '_ratios' in self.meta:
            good = [k for k in self.meta['_ratios'].keys()
                    if k in self['name'] or k in self['group']]
            if good:
                st = tuple(x + " ("
                           + ", ".join(self.meta['_ratios'][x].keys()) + ")"
                           for x in good)
                repr += "\n\n" + "Tied:\t" + "\n\t".join(st)
        return repr

    @classmethod
    def read(cls, file, *args, **kwargs):
        """Read a table from file.

        Parameters
        ----------

          file : str
              The name of the file to read, either a full valid path,
              or named file in the PAHFIT science_packs directory.

        Returns
        -------
          table : Features
              A filled `~pahfit.features.Features` table.

        Notes
        -----
        If reading a YAML file, reads it in as a science pack and
        return the new table. Otherwise, uses astropy's normal Table
        reader.
        """
        if file.endswith(".yaml") or file.endswith(".yml"):
            return cls._read_scipack(file)
        else:
            table = super().read(file, *args, **kwargs)
            cls._index_table(table)
            return table

    @classmethod
    def _read_scipack(cls, file):
        """Read a science pack specification from YAML file.

        For parameters and return, see `read`.
        """

        feat_tables = dict()

        if not os.path.isfile(file):
            pack_path = resource_filename("pahfit", "packs/science")
            file = os.path.join(pack_path, file)
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
                                or (isinstance(v, dict)  # names: values dict
                                    and cls._param_attrs.isdisjoint(v.keys())))
                           for (k, v) in elem.items())
            if hasFeatures and hasLists:
                raise PAHFITFeatureError("A single group cannot contain both 'features'"
                                         f" and parameter list(s): {name}\n\t{file}")
            isGroup = (hasFeatures or hasLists)
            bounds = None
            if isGroup:  # A named group of features
                if 'bounds' in elem:  # group-level bounds
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
                if 'tied' in elem:  # group-level ties
                    for denom, tie in elem.pop('tied').items():
                        if not ('param' in tie):
                            raise PAHFITFeatureError(f"Group-level ties require param: {name}")
                        try:
                            cls._add_tie_ratio(feat_tables, name, tie['param'], denom, tie)
                        except ValueError as e:
                            raise PAHFITFeatureError(f"{e}: {name}, {tie['param']}")
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
    def _add_tie_ratio(cls, t: dict, num: str, param: str, denom: str, tie):
        """Add tie ratios to dictionary t.

        Creates a ``_ratios`` entry for ratio (with bounds) between
        numerator ``num`` to denominator ``denom`` (both feature or
        group names), for parameter ``param``.  Pass the ``tie`` dict
        or string.  If ``tie`` is a string, it is ignored.
        """
        if '_ratios' not in t:
            t['_ratios'] = {}
        if num not in t['_ratios']:
            t['_ratios'][num] = {}
        if isinstance(tie, dict) and 'ratio' in tie:
            # Table doesn't like masked in .meta
            ratio = value_bounds(*cls._parse_value(tie['ratio']), no_masked=True)
        else:
            ratio = None
        t['_ratios'][num][param] = dict(denom=denom, ratio=ratio)

    @classmethod
    def _parse_value(cls, val):
        """Parse a value for param with optional bounds.

        Parameters
        ----------

        val : dict or float
            The value to parse

        Returns
        -------

        (value, bounds) : (float, None or tuple)
        """
        bounds = None
        if isinstance(val, dict):  # A param attribute dictionary
            unknown_attrs = [x for x in val.keys() if x not in cls._param_attrs]
            if unknown_attrs:
                raise ValueError(f"Unknown parameter attributes {', '.join(unknown_attrs)}")
            if 'value' not in val:
                raise ValueError("Missing 'value' attribute")
            value = val['value']
            if 'bounds' in val:  # individual param bounds
                bounds = val['bounds']
        else:
            value = val  # a bare value
        return (value, bounds)

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
            if isinstance(val, dict) and 'tied' in val:
                tie = val['tied']
                if isinstance(tie, dict):
                    try:
                        feat = tie['feature']
                    except KeyError:
                        raise PAHFITFeatureError(f"Tie requires feature: {name}, {param}")
                else:
                    feat = tie  # Just a name of a feature or group: no ratio!
                try:
                    cls._add_tie_ratio(t, name, param, feat, tie)
                except ValueError as e:
                    raise PAHFITFeatureError(f"{e}: {name}, {param}")
            else:
                try:
                    (value, b) = cls._parse_value(val)
                except ValueError as e:
                    raise PAHFITFeatureError(f"{name} ({kind}, {group}):\n\t{e}")
                if isinstance(bounds, dict):  # group-level bounds
                    if param in bounds:
                        if b:
                            raise PAHFITFeatureError("Cannot set both group and feature bounds"
                                                     f" for {param}: {name} ({kind=}, {group=})")
                        b = bounds[param]
                else:
                    b = b or bounds  # maybe none
                if b and param in cls._no_bounds:
                    raise PAHFITFeatureError(f"Parameter {param} cannot have bounds: "
                                             f"{name} ({kind=}, {group=})")
                try:
                    t[kind][name][param] = (value if param in cls._no_bounds
                                            else value_bounds(value, b))
                except ValueError as e:
                    raise PAHFITFeatureError("Error initializing value and bounds for"
                                             f" {name} ({kind}, {group}):\n\t{e}")

    @classmethod
    def _construct_table(cls, inp: dict):
        """Construct a masked table from input dictionary INP.
        INP is a dictionary with feature names as the keys, and a
        dictionary of feature parameters as value.  Each value in the
        feature parameter dictionary is either a value or tuple of 3
        values, expressing bounds.
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
                    else:  # Any missing required parameters: [0,] bounds
                        params[missing] = value_bounds(0.0, bounds=(0.0, None))
                rows.append(dict(name=name, **params))
            table_columns = rows[0].keys()
            t = cls(rows, names=table_columns)
            for p in KIND_PARAMS[kind]:
                if p not in cls._no_bounds:
                    t[p].info.format = "0.4g"  # Nice format (customized by Formatter)
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
                row[col_name].mask[0] = mask_value

    def unmask_feature(self, name):
        """Remove the mask for all parameters of a feature."""
        self.mask_feature(name, mask_value=False)
