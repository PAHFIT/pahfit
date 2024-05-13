import numpy as np
from astropy.table.pprint import TableFormatter


# * Special table formatting for bounded (val, min, max) values
def fmt_func(fmt: str):
    """Format bounded variables specially."""
    if fmt.startswith('%'):
        fmt = fmt[1:]

    def _fmt(x):
        ret = f"{x['val']:{fmt}}"
        if np.isnan(x['min']) and np.isnan(x['max']):
            return ret + " (fixed)"
        else:
            mn = ("-∞" if np.isnan(x['min']) or x['min'] == -np.inf
                  else f"{x['min']:{fmt}}")
            mx = ("∞" if np.isnan(x['max']) or x['max'] == np.inf
                  else f"{x['max']:{fmt}}")
            return f"{ret} ({mn}, {mx})"
    return _fmt


class BoundedParTableFormatter(TableFormatter):
    """Format bounded parameters.
    Bounded parameters are 3-field structured arrays, with fields
    'val', 'min', and 'max'.  To be set as Table's `TableFormatter'.
    """
    def _pformat_table(self, table, *args, **kwargs):
        bpcols = []
        tlfmt = table.meta.get('pahfit_format')
        try:
            for col in table.columns.values():
                if len(col.dtype) == 3:  # bounded!
                    bpcols.append((col, col.info.format))
                    fmt = col.meta.get('pahfit_format') or tlfmt or "g"
                    col.info.format = fmt_func(fmt)
            return super()._pformat_table(table, *args, **kwargs)
        finally:
            for col, fmt in bpcols:
                col.info.format = fmt

    def _name_and_structure(self, name, *args):
        "Simplified column name: no val, min, max needed."
        return name
